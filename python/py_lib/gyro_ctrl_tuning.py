"""
==========================================================================
GYRO CTRL TUNING - Betaflight Controller Analysis GYRO TUNING CLASS
==========================================================================
Betaflight Controller Tuning Analysis Script
Purpose: Calculation for Gyro Tuning

Author: [Janick Dort, Yuri Bianchi, Dario Jurietti]
Supervisor: [Michael Peter]
Date: [28.02.2026]
==========================================================================
"""

import numpy as np
import control as ct
from scipy.signal import hann
from typing import Optional, List
import warnings
import time

from .pidtuninglib import (
    header_info,
    ClosedLoop,
    apply_rotfiltfilt,
    estimate_frequency_response,
    calculate_transfer_functions,
    calculate_closed_loop,
    downsample_frd,
    calculate_step_response_from_frd,
    get_pid_scale,
)


def get_ind_eval(sinarg: np.ndarray, gyro_data: np.ndarray, threshold: float = 500.0) -> np.ndarray:
    """
    Select evaluation intervals based on phase window and signal variance.
    
    Parameters
    ----------
    sinarg : np.ndarray
        Phase signal
    gyro_data : np.ndarray
        Gyro measurement data
    threshold : float
        Variance threshold for valid intervals
    
    Returns
    -------
    np.ndarray
        Boolean mask of valid evaluation indices
    """
    ndata = len(sinarg)
    ind_eval_candidate = sinarg > 0
    
    signal = np.zeros(ndata)
    signal[ind_eval_candidate] = 1
    dsignal = np.diff(signal, prepend=0)
    
    ind_eval_start = np.where(dsignal > 0.9)[0]
    ind_eval_end = np.where(dsignal < -0.9)[0] - 1
    
    ind_eval = np.zeros(ndata, dtype=bool)
    
    for start, end in zip(ind_eval_start, ind_eval_end):
        ind_verify = slice(start, end + 1)
        if np.var(gyro_data[ind_verify]) > threshold:
            ind_eval[ind_verify] = True
    
    return ind_eval


class GyroCtrlTuning:
    """Gyro controller tuning and analysis."""
    
    def __init__(
        self,
        data: np.ndarray,
        ind,
        Ts_log: float,
        para: header_info,
        Ts_cntr: float
    ):
        """
        Initialize gyro controller tuning object.
        
        Parameters
        ----------
        data : np.ndarray
            Flight data
        ind : ColumnIndices
            Column indices
        Ts_log : float
            Logging sample time
        para : header_info
            Parameter structure
        Ts_cntr : float
            Control loop sample time
        """
        self.data = data
        self.ind = ind
        self.Ts_log = Ts_log
        self.para = para
        self.Ts_cntr = Ts_cntr
        self.linewidth = 1.2
        
        # Storage for results
        self.T: List[ct.FRD] = []
        self.Guw: List[ct.FRD] = []
        self.Gvw: List[ct.FRD] = []
        self.P_gef: List[ct.FRD] = []
        self.P: List[ct.FRD] = []
        self.Cpi: List[ct.FRD] = []
        self.Cd: List[ct.FRD] = []
        self.Gf_ana: List[ct.FRD] = []
        self.Cpi_ana: List[ct.FRD] = []
        self.Cd_ana: List[ct.FRD] = []
        self.PID: List[np.ndarray] = []
        self.para_used: List[header_info] = []
        self.omega_bode: Optional[np.ndarray] = None
        self.C_T: List[ct.FRD] = []
        self.C_Guw: List[ct.FRD] = []
        self.Coh: List[ct.FRD] = []
        self.throttle_avg: Optional[float] = None
        self.Nest: Optional[int] = None
        self.ind_ax: Optional[int] = None
        
        # New controller results
        self.Cpi_ana_new: Optional[ct.FRD] = None
        self.Gf_ana_new: Optional[ct.FRD] = None
        self.Cd_ana_new: Optional[ct.FRD] = None
        self.CloLoAan: Optional[ClosedLoop] = None
        self.CloLoAanNew: Optional[ClosedLoop] = None
        self.step_time: Optional[np.ndarray] = None
        self.step_resp_tra: Optional[np.ndarray] = None
        self.step_resp_com: Optional[np.ndarray] = None
    
    def calculate_transfer_func(
        self,
        Nestfatra: float,
        koverlaptra: float
    ) -> 'GyroCtrlTuning':
        """
        Calculate transfer functions from flight data.
        
        Parameters
        ----------
        Nestfatra : float
            Window length in seconds
        koverlaptra : float
            Overlap factor (0-1)
        
        Returns
        -------
        GyroCtrlTuning
            Self with calculated transfer functions
        """
        # Analysis window parameters
        self.Nest = int(np.round(Nestfatra / self.Ts_log))
        Noverlap = int(np.floor(koverlaptra * self.Nest))
        window = hann(self.Nest, sym=False)
        
        # Design linear filter for zero phase excitation filter
        Dlp = np.sqrt(3) / 2
        wlp = 2 * np.pi * 10
        Glp_cont = ct.tf([wlp**2], [1, 2*Dlp*wlp, wlp**2])
        Glp = ct.c2d(Glp_cont, self.Ts_log, method='tustin')
        
        # Process each axis
        n_axes = 3
        self.T = [None] * n_axes
        self.Guw = [None] * n_axes
        self.Gvw = [None] * n_axes
        self.P = [None] * n_axes
        self.Cpi = [None] * n_axes
        self.Cd = [None] * n_axes
        self.C_T = [None] * n_axes
        self.C_Guw = [None] * n_axes
        self.P_gef = [None] * n_axes
        self.Coh = [None] * n_axes
        self.Gf_ana = [None] * n_axes
        self.Cpi_ana = [None] * n_axes
        self.Cd_ana = [None] * n_axes
        self.PID = [None] * n_axes
        self.para_used = [None] * n_axes
        
        for ind_axis in range(n_axes):
            # Get evaluation indices
            ind_eval = get_ind_eval(
                self.data[:, self.ind.sinarg],
                self.data[:, self.ind.gyroADC[ind_axis]]
            )
            
            sinarg_ax = self.data[:, self.ind.sinarg].copy()
            sinarg_ax[~ind_eval] = 0
            
            # Calculate average throttle
            self.throttle_avg = np.median(self.data[ind_eval, self.ind.setpoint[3]]) / 1.0e3
            
            # Input signal: w (filtered setpoint)
            w = self.data[:, self.ind.setpoint[ind_axis]]
            inp = apply_rotfiltfilt(Glp, sinarg_ax, w[:, np.newaxis])
            
            # Output y: filtered gyro
            y = self.data[:, self.ind.gyroADC[ind_axis]]
            out_y = apply_rotfiltfilt(Glp, sinarg_ax, y[:, np.newaxis])
            
            # Calculate complementary sensitivity T
            T_ax, C_T_ax, freq_hz, _ = estimate_frequency_response(
                inp[ind_eval].ravel(),
                out_y[ind_eval].ravel(),
                window,
                Noverlap,
                self.Nest,
                self.Ts_log
            )
            
            # Calculate control sensitivity Guw
            u = self.data[:, self.ind.axisSum[ind_axis]]
            out_u = apply_rotfiltfilt(Glp, sinarg_ax, u[:, np.newaxis])
            Guw_ax, C_Guw_ax, _, _ = estimate_frequency_response(
                inp[ind_eval].ravel(),
                out_u[ind_eval].ravel(),
                window,
                Noverlap,
                self.Nest,
                self.Ts_log
            )
            
            # Calculate PI controller response Gvw
            v = self.data[:, self.ind.axisSumPI[ind_axis]]
            out_v = apply_rotfiltfilt(Glp, sinarg_ax, v[:, np.newaxis])
            Gvw_ax, _, _, _ = estimate_frequency_response(
                inp[ind_eval].ravel(),
                out_v[ind_eval].ravel(),
                window,
                Noverlap,
                self.Nest,
                self.Ts_log
            )
     
            # Calculate plant response (indirect method)
            P_gef_ax = T_ax / Guw_ax
            
            # Calculate controller frequency responses
            Cpi_ax = Gvw_ax / (1 - T_ax)
            Cd_ax = Guw_ax * Gvw_ax / T_ax * (1 / Guw_ax - 1 / Gvw_ax)
            
            # Frequency vector for Bode plots
            omega_bode_ax = 2 * np.pi * T_ax.omega / (2 * np.pi)
            
            # Store results
            self.T[ind_axis] = T_ax
            self.Guw[ind_axis] = Guw_ax
            self.Gvw[ind_axis] = Gvw_ax
            self.P_gef[ind_axis] = P_gef_ax
            self.Cpi[ind_axis] = Cpi_ax
            self.Cd[ind_axis] = Cd_ax
            self.omega_bode = omega_bode_ax
            self.C_T[ind_axis] = C_T_ax
            self.C_Guw[ind_axis] = C_Guw_ax
            
            # Analytical controller transfer functions
            Cpi_ana_ax, Cd_ana_ax, Gf_ana_ax, PID_ax, para_used_ax = \
                calculate_transfer_functions(
                    self.para,
                    ind_axis,
                    self.throttle_avg,
                    self.Ts_cntr
                )
            
            # Downsample to logging grid if necessary
            if Gf_ana_ax.dt < self.Ts_log:
                freq_target = T_ax.omega / (2 * np.pi)
                Gf_ana_ax = downsample_frd(Gf_ana_ax, self.Ts_log, freq_target)
                Cpi_ana_ax = downsample_frd(Cpi_ana_ax, self.Ts_log, freq_target)
                Cd_ana_ax = downsample_frd(Cd_ana_ax, self.Ts_log, freq_target)
            
            # Store analytical results
            self.Gf_ana[ind_axis] = Gf_ana_ax
            self.Cpi_ana[ind_axis] = Cpi_ana_ax
            self.Cd_ana[ind_axis] = Cd_ana_ax
            self.PID[ind_axis] = PID_ax
            self.para_used[ind_axis] = para_used_ax
            self.P[ind_axis] = P_gef_ax / Gf_ana_ax
            self.Coh[ind_axis] = C_T_ax * C_Guw_ax
        
        return self
    
    def calculate_new_controller(
        self,
        ind_ax: int,
        P_new: float,
        I_new: float,
        D_new: float,
        default_parameters: bool,
        para_new: header_info
    ) -> 'GyroCtrlTuning':
        """
        Calculate new controller parameters.
        
        Parameters
        ----------
        ind_ax : int
            Axis index (0=roll, 1=pitch, 2=yaw)
        P_new : float
            New P gain
        I_new : float
            New I gain
        D_new : float
            New D gain
        default_parameters : bool
            Use default parameters
        para_new : header_info
            New parameter structure
        
        Returns
        -------
        GyroCtrlTuning
            Self with new controller
        """
        start_time = time.time()
        
        pid_axis = ['rollPID', 'pitchPID', 'yawPID']
        
        # Handle default parameter case
        if default_parameters:
            para_new = self.para.make_copy()
            pid_vec = self.para.get_list(pid_axis[ind_ax])
            P_new = pid_vec[0]
            I_new = pid_vec[1]
            D_new = pid_vec[2]
        
        # Display used PID configuration
        print('   used PID parameters are:')
        print(f'      {pid_axis[ind_ax]}: {self.para.get_list(pid_axis[ind_ax])[:3]}')
        
        # Display used parameters
        pu = self.para_used[ind_ax]
        print('   used parameters are:')
        for key in pu.data.keys():
            val = pu.data[key]
            try:
                val_num = float(val.replace(',', '').split()[0])
                print(f'      {key}: {int(val_num)}')
            except:
                print(f'      {key}: {val}')
        
        # Get axis-specific scaling
        pid_scale = get_pid_scale(ind_ax) + [1]
        
        # Calculate new gains
        PID_new = np.zeros(4)
        PID_new[0] = P_new * pid_scale[0]
        fI = self.PID[ind_ax][1] / (2 * np.pi * self.PID[ind_ax][0])
        PID_new[1] = I_new * pid_scale[1]
        fI_new = PID_new[1] / (2 * np.pi * PID_new[0])
        PID_new[2] = D_new * pid_scale[2]
        PID_new[3] = 0  # No feedforward
        
        print(f'   used fI is: {fI:.2f} Hz\n')
        
        # Update parameters
        print('   new PID parameters are:')
        para_new_pid = np.round(PID_new / pid_scale).astype(int)
        para_new_list = [para_new_pid[0], para_new_pid[1], para_new_pid[2], 
                        para_new_pid[2], para_new_pid[3]]
        para_new.data[pid_axis[ind_ax]] = ','.join(map(str, para_new_list))
        print(f'      {pid_axis[ind_ax]}: {para_new_pid[:3]}')
        
        # Generate new transfer functions
        self.Cpi_ana_new, self.Cd_ana_new, self.Gf_ana_new, PID_new, para_used_new = \
            calculate_transfer_functions(para_new, ind_ax, self.throttle_avg, self.Ts_cntr)
        
        # Display new parameters
        print('   new parameters are:')
        for key in para_used_new.data.keys():
            val = para_used_new.data[key]
            try:
                val_num = float(val.replace(',', '').split()[0])
                print(f'      {key}: {int(val_num)}')
            except:
                print(f'      {key}: {val}')
        
        print(f'   new used fI is: {fI_new:.2f} Hz\n')
        
        # Downsample if necessary
        if self.Gf_ana_new.dt < self.Ts_log:
            freq_target = self.P_gef[ind_ax].omega / (2 * np.pi)
            self.Gf_ana_new = downsample_frd(self.Gf_ana_new, self.Ts_log, freq_target)
            self.Cpi_ana_new = downsample_frd(self.Cpi_ana_new, self.Ts_log, freq_target)
            self.Cd_ana_new = downsample_frd(self.Cd_ana_new, self.Ts_log, freq_target)
        
        self.ind_ax = ind_ax
        
        print(f'Calculation time: {time.time() - start_time:.3f}s')
        
        return self
    
    def get_tuning_data(self, do_compensate_iterm: bool) -> 'GyroCtrlTuning':
        """
        Calculate closed loop responses and step responses.
        
        Parameters
        ----------
        do_compensate_iterm : bool
            Apply I-term compensation
        
        Returns
        -------
        GyroCtrlTuning
            Self with tuning data
        """
        start_time = time.time()
        
        # Calculate closed loop responses
        CL_ana = calculate_closed_loop(
            self.Cpi_ana[self.ind_ax],
            ct.tf([1], [1], self.Ts_log),
            self.P[self.ind_ax],
            self.Gf_ana[self.ind_ax],
            self.Cd_ana[self.ind_ax]
        )
        
        CL_ana_new = calculate_closed_loop(
            self.Cpi_ana_new,
            ct.tf([1], [1], self.Ts_log),
            self.P[self.ind_ax],
            self.Gf_ana_new,
            self.Cd_ana_new
        )
        
        # Apply I-term compensation if enabled
        if do_compensate_iterm:
            Cpi_com = self.Cpi[self.ind_ax] / self.Cpi_ana[self.ind_ax]
            
            CL_ana_ = calculate_closed_loop(
                self.Cpi_ana[self.ind_ax] * Cpi_com,
                ct.tf([1], [1], self.Ts_log),
                self.P[self.ind_ax],
                self.Gf_ana[self.ind_ax],
                self.Cd_ana[self.ind_ax]
            )
            
            CL_ana_new_ = calculate_closed_loop(
                self.Cpi_ana_new * Cpi_com,
                ct.tf([1], [1], self.Ts_log),
                self.P[self.ind_ax],
                self.Gf_ana_new,
                self.Cd_ana_new
            )

            CL_ana.T = CL_ana_.T
            CL_ana_new.T = CL_ana_new_.T
        
        self.CloLoAan = CL_ana
        self.CloLoAanNew = CL_ana_new
        
        # Step response analysis
        f_max = min(
            self.para.get_int('dyn_notch_min_hz'),
            self.para.get_int('gyro_rpm_notch_min')
        )
        
        T_mean = 0.1 * np.array([-1, 1]) + (self.Nest * self.Ts_log) / 2
        self.step_time = np.arange(self.Nest) * self.Ts_log
        
        # Tracking performance
        step_resp = np.column_stack([
            calculate_step_response_from_frd(CL_ana.T, f_max),
            calculate_step_response_from_frd(CL_ana_new.T, f_max),
            calculate_step_response_from_frd(self.T[self.ind_ax], f_max)
        ])
        
        # Normalize
        mask = (self.step_time > T_mean[0]) & (self.step_time < T_mean[1])
        step_resp_mean = np.mean(step_resp[mask, :], axis=0)
        step_resp = step_resp / step_resp_mean
        self.step_resp_tra = step_resp
        
        # Disturbance rejection
        step_resp = np.column_stack([
            calculate_step_response_from_frd(CL_ana.SP, f_max),
            calculate_step_response_from_frd(CL_ana_new.SP, f_max)
        ])
        
        step_resp_mean = np.mean(step_resp[mask, :], axis=0)
        step_resp = step_resp - step_resp_mean
        self.step_resp_com = step_resp
        
        print(f'Tuning data calculation time: {time.time() - start_time:.3f}s')
        
        return self