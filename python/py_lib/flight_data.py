"""
==========================================================================
FLIGHT DATA - Betaflight Controller Analysis DATA IMPORT CLASS
==========================================================================
Betaflight Controller Tuning Analysis Script
Purpose: Read Data for further calculations

Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
Supervisor: [Michael Peter]
Date: [28.02.2026]
==========================================================================
"""

import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
import time

from .pidtuninglib import header_info


@dataclass
class ColumnIndices:
    """Storage for column indices from flight log."""
    time: int = 0
    gyroUnfilt: np.ndarray = field(default_factory=lambda: np.array([]))
    gyroADC: np.ndarray = field(default_factory=lambda: np.array([]))
    setpoint: np.ndarray = field(default_factory=lambda: np.array([]))
    rcCommand: np.ndarray = field(default_factory=lambda: np.array([]))
    axisP: np.ndarray = field(default_factory=lambda: np.array([]))
    axisI: np.ndarray = field(default_factory=lambda: np.array([]))
    axisD: np.ndarray = field(default_factory=lambda: np.array([]))
    axisSum: np.ndarray = field(default_factory=lambda: np.array([]))
    motor: np.ndarray = field(default_factory=lambda: np.array([]))
    debug: np.ndarray = field(default_factory=lambda: np.array([]))
    axisSumPI: np.ndarray = field(default_factory=lambda: np.array([]))
    sinarg: int = 0


class FlightData:
    """Flight data loader and preprocessor for Betaflight logs."""
    
    def __init__(self, file_path: str):
        """
        Initialize flight data object.
        
        Parameters
        ----------
        file_path : str
            Path to the flight log file (.csv)
        """
        self.file_path = Path(file_path)
        self.linewidth = 1.2
        
        # Data storage
        self.time: Optional[np.ndarray] = None
        self.data: Optional[np.ndarray] = None
        self.ind: Optional[ColumnIndices] = None
        self.para: Optional[header_info] = None
        
        # Calculated data
        self.Ts_log: Optional[float] = None
        self.Ts_cntr: Optional[float] = None
        
    def get_data(self) -> 'FlightData':
        """
        Load and process flight log data.
        
        Returns
        -------
        FlightData
            Self with loaded data
        """
        # Extract header information
        self.para, nheader, self.ind, ind_cntr = self._extract_header_information()
        
        # Load data from CSV or cached pickle file
        print("Loading flight data...")
        start_time = time.time()
        
        pickle_path = self.file_path.with_suffix('.pkl')
        
        try:
            # Try to load from pickle cache
            with open(pickle_path, 'rb') as f:
                self.data = pickle.load(f)
            print(f"Loaded from cache in {time.time() - start_time:.3f}s")
        except (FileNotFoundError, pickle.PickleError):
            # Load from CSV
            self.data = pd.read_csv(
                self.file_path,
                skiprows=nheader,
                header=None
            ).values
            
            # Save to pickle for faster loading next time
            with open(pickle_path, 'wb') as f:
                pickle.dump(self.data, f)
            print(f"Loaded from CSV and cached in {time.time() - start_time:.3f}s")
        
        # Data preprocessing
        self._preprocess_data(ind_cntr)
        
        return self
    
    def _extract_header_information(self):
        """Extract header parameters and column indices from log file."""
        # Get header parameters
        para = header_info.get_header(self.file_path)
        
        # Find data header line and count header rows
        nheader = 0
        ind = ColumnIndices()
        ind_cntr = 0
        
        with open(self.file_path, 'r') as f:
            for line_num, line in enumerate(f):
                nheader = line_num + 1
                if 'loopIteration' in line:
                    # Parse column names
                    ind, ind_cntr = self._parse_column_names(line)
                    break
        
        return para, nheader, ind, ind_cntr
    
    def _parse_column_names(self, header_line: str):
        """Parse column names from data header line."""
        ind = ColumnIndices()
        
        # Split by comma and remove quotes
        columns = [col.strip('"') for col in header_line.strip().split(',')]
        
        # Lists to store multi-column indices
        gyro_unfilt = []
        gyro_adc = []
        setpoint = []
        rc_command = []
        axis_p = []
        axis_i = []
        axis_d = []
        axis_sum = []
        motor = []
        debug = []
        
        for idx, col in enumerate(columns):
            if col == 'time (us)':
                ind.time = idx
            elif 'gyroUnfilt' in col:
                gyro_unfilt.append(idx)
            elif 'gyroADC' in col:
                gyro_adc.append(idx)
            elif col.startswith('setpoint['):
                setpoint.append(idx)
            elif col.startswith('rcCommand['):
                rc_command.append(idx)
            elif col.startswith('axisP['):
                axis_p.append(idx)
            elif col.startswith('axisI['):
                axis_i.append(idx)
            elif col.startswith('axisD['):
                axis_d.append(idx)
            elif col.startswith('axisSum['):
                axis_sum.append(idx)
            elif col.startswith('motor['):
                motor.append(idx)
            elif col.startswith('debug['):
                debug.append(idx)
        
        # Convert to numpy arrays
        ind.gyroUnfilt = np.array(gyro_unfilt)
        ind.gyroADC = np.array(gyro_adc)
        ind.setpoint = np.array(setpoint)
        ind.rcCommand = np.array(rc_command)
        ind.axisP = np.array(axis_p)
        ind.axisI = np.array(axis_i)
        ind.axisD = np.array(axis_d)
        ind.axisSum = np.array(axis_sum)
        ind.motor = np.array(motor)
        ind.debug = np.array(debug)
        
        ind_cntr = len(columns)
        
        return ind, ind_cntr
    
    def _preprocess_data(self, ind_cntr: int):
        """Preprocess loaded data."""
        # Expand indices for additional data columns
        self.ind.axisSumPI = np.arange(ind_cntr, ind_cntr + 3)
        self.ind.sinarg = self.ind.debug[0]
        
        # Convert microseconds to seconds for time vector
        self.time = (self.data[:, self.ind.time] - self.data[0, self.ind.time]) * 1.0e-6
        
        # Unscale high resolution gain if enabled
        if self.para.get_int('blackbox_high_resolution'):
            blackbox_high_resolution_scale = 10.0
            ind_bb_high_res = np.concatenate([
                self.ind.gyroADC,
                self.ind.gyroUnfilt,
                self.ind.rcCommand,
                self.ind.setpoint[:3]
            ])
            self.data[:, ind_bb_high_res] = (
                self.data[:, ind_bb_high_res] / blackbox_high_resolution_scale
            )
        
        # Unscale and remap sinarg
        sinarg_scale = 5.0e3
        self.data[:, self.ind.sinarg] = self.data[:, self.ind.sinarg] / sinarg_scale
        
        # Create additional entry for PI sum
        pi_sum = self.data[:, self.ind.axisP] + self.data[:, self.ind.axisI]
        self.data = np.column_stack([self.data, pi_sum])
        
        # Calculate sampling times
        Ts = self.para.get_int('looptime') * 1.0e-6  # Gyro loop
        self.Ts_cntr = self.para.get_int('pid_process_denom') * Ts  # Control loop
        self.Ts_log = self.para.get_int('frameIntervalPDenom') * self.Ts_cntr  # Logging loop