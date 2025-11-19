"""
Betaflight Controller Tuning - Step Response Analysis

This script analyzes Betaflight blackbox log data to extract and visualize
frequency and step responses from chirp excitation signals.
"""
import pickle
import time
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal.windows import hann

import mylib2

# =========================================================================
# CONSTANTS
# =========================================================================
MICROSECONDS_TO_SECONDS = 1.0e-6
SAMPLE_TIME_BASE = 125e-6  # Base sampling time in seconds
OVERLAP_RATIO = 0.9
MAX_FREQUENCY_HZ = 200

# =========================================================================
# FILE CONFIGURATION
# =========================================================================
DATA_FILE = Path("./20250908/20250908_flipmini_00.bbl.csv")


# =========================================================================
# HELPER FUNCTIONS
# =========================================================================
def find_header_line(file_path: Path) -> int:
    """
    Find the line number where the CSV header containing 'loopIteration' is located.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        Zero-indexed line number of the header
    """
    with open(file_path, "r") as f:
        for line_count, line in enumerate(f):
            if "loopIteration" in line:
                return line_count
    raise ValueError("Could not find 'loopIteration' in CSV file")


def load_or_preprocess_data(csv_path: Path) -> Tuple[np.ndarray, list]:
    """
    Load preprocessed pickle data or read and preprocess CSV file.
    
    Args:
        csv_path: Path to the CSV file
        
    Returns:
        Tuple of (data array, column names list)
    """
    pickle_path = csv_path.with_suffix(".pkl")
    
    # Try loading preprocessed data
    if pickle_path.exists():
        with open(pickle_path, "rb") as f:
            data_dict = pickle.load(f)
        print(f"Loaded preprocessed data from {pickle_path}")
        return data_dict["data"], data_dict["columns"]
    
    # Read and preprocess CSV
    print("Reading CSV file...")
    header_line = find_header_line(csv_path)
    df = pd.read_csv(csv_path, skiprows=header_line)
    
    column_names = df.columns.tolist()
    data = df.values
    
    # Save for future use
    with open(pickle_path, "wb") as f:
        pickle.dump({"data": data, "columns": column_names}, f)
    print(f"Saved preprocessed data to {pickle_path}")
    
    return data, column_names


def get_column_index(column_names: list, col_name: str) -> int:
    """
    Get column index by name, raising error if not found.
    
    Args:
        column_names: List of column names
        col_name: Name of the column to find
        
    Returns:
        Index of the column
        
    Raises:
        ValueError: If column is not found
    """
    try:
        return column_names.index(col_name)
    except ValueError:
        raise ValueError(f"Column '{col_name}' not found in data")


def detect_chirp_boundaries(chirp_signal: np.ndarray) -> Tuple[int, int]:
    """
    Detect the start and end indices of the chirp signal.
    
    Args:
        chirp_signal: Array containing the chirp signal values
        
    Returns:
        Tuple of (start_index, end_index)
    """
    chirp_active = chirp_signal > 0
    # Find rising edge (start) and falling edge (end)
    edges = np.diff(chirp_active.astype(int), prepend=0)
    start_indices = np.where(edges == 1)[0]
    end_indices = np.where(edges == -1)[0]
    
    if len(start_indices) == 0 or len(end_indices) == 0:
        raise ValueError("Could not detect chirp boundaries")
    
    return start_indices[0], end_indices[0]


def calculate_frequency_response_params(ts_log: float) -> dict:
    """
    Calculate parameters for frequency response estimation.
    
    Args:
        ts_log: Logging sample time in seconds
        
    Returns:
        Dictionary containing calculation parameters
    """
    n_estimate = round(2.0 / ts_log)
    n_overlap = round(OVERLAP_RATIO * n_estimate)
    window = hann(n_estimate, sym=False)
    
    return {
        "n_estimate": n_estimate,
        "n_overlap": n_overlap,
        "window": window,
        "fs": 1 / (SAMPLE_TIME_BASE),
    }


def plot_raw_data(time: np.ndarray, setpoint: np.ndarray, gyro: np.ndarray) -> None:
    """Plot raw setpoint and gyro data."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
    
    ax1.plot(time, setpoint)
    ax1.set_title("Setpoint")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Roll Rate (deg/s)")
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(time, gyro)
    ax2.set_title("Gyro Unfiltered")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Roll Rate (deg/s)")
    ax2.grid(True, alpha=0.3)
    
    fig.suptitle("Raw Data", fontsize=14, fontweight='bold')
    fig.tight_layout()


def plot_chirp_data(
    time: np.ndarray, setpoint: np.ndarray, gyro: np.ndarray
) -> None:
    """Plot data over chirp excitation period."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
    
    ax1.plot(time, setpoint)
    ax1.set_title("Setpoint")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Roll Rate (deg/s)")
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(time, gyro)
    ax2.set_title("Gyro Unfiltered")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Roll Rate (deg/s)")
    ax2.grid(True, alpha=0.3)
    
    fig.suptitle("Data Over Chirp", fontsize=14, fontweight='bold')
    fig.tight_layout()


def plot_frequency_response(freq: np.ndarray, G_fresp: np.ndarray, f_max: float) -> None:
    """Plot frequency response magnitude."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    G_magnitude_db = 20 * np.log10(np.abs(np.squeeze(G_fresp)))
    ax.plot(freq, G_magnitude_db)
    ax.set_xlim(0, f_max)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.set_title("Frequency Response")
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()


def plot_step_response(step_response: np.ndarray, f_max: float) -> None:
    """Plot step response."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(step_response)
    ax.set_xlim(0, f_max)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Amplitude")
    ax.set_title("Step Response")
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()


# =========================================================================
# MAIN SCRIPT
# =========================================================================
def main() -> None:
    """Main analysis workflow."""
    print("=" * 70)
    print("Betaflight Controller Tuning - Step Response Analysis")
    print("=" * 70)
    
    # Load data
    start_time = time.time()
    data, column_names = load_or_preprocess_data(DATA_FILE)
    elapsed_time = time.time() - start_time
    print(f"Data loaded: {data.shape[0]} rows in {elapsed_time:.2f} seconds\n")
    
    # Get column indices
    col_time = get_column_index(column_names, "time")
    col_setpoint_roll = get_column_index(column_names, "setpoint[0]")
    col_gyro_unfilt_roll = get_column_index(column_names, "gyroUnfilt[0]")
    col_debug_chirp = get_column_index(column_names, "debug[0]")
    
    # Extract time vector and convert to seconds
    time_vec = (data[:, col_time] - data[0, col_time]) * MICROSECONDS_TO_SECONDS
    print(f"Extracted {len(time_vec)} samples ({time_vec[-1]:.2f} seconds of data)\n")
    
    # Detect chirp boundaries
    chirp_signal = data[:, col_debug_chirp]
    chirp_start, chirp_end = detect_chirp_boundaries(chirp_signal)
    print(f"Chirp detected: start={chirp_start}, end={chirp_end}\n")
    
    # Extract roll data
    setpoint_roll = data[:, col_setpoint_roll]
    gyro_unfilt_roll = data[:, col_gyro_unfilt_roll]
    
    # Extract data over chirp period
    time_over_chirp = time_vec[chirp_start:chirp_end]
    setpoint_over_chirp = setpoint_roll[chirp_start:chirp_end]
    gyro_over_chirp = gyro_unfilt_roll[chirp_start:chirp_end]
    
    # Calculate sampling times
    ts_base = SAMPLE_TIME_BASE
    ts_controller = 2 * ts_base
    ts_log = 2 * ts_controller
    
    # Calculate frequency response
    params = calculate_frequency_response_params(ts_log)
    print(f"Frequency response parameters:")
    print(f"  N_estimate: {params['n_estimate']}")
    print(f"  N_overlap: {params['n_overlap']}")
    print(f"  Window length: {len(params['window'])}\n")
    
    G, C, freq, _ = mylib2.estimate_frequency_response(
        setpoint_over_chirp,
        gyro_over_chirp,
        params["window"],
        params["n_overlap"],
        params["n_estimate"],
        ts_log,
    )
    
    # Calculate step response
    step_resp = mylib2.calculate_step_response_from_frd(G, MAX_FREQUENCY_HZ)
    
    # Visualization
    print("Generating plots...")
    plot_raw_data(time_vec, setpoint_roll, gyro_unfilt_roll)
    plot_chirp_data(time_over_chirp, setpoint_over_chirp, gyro_over_chirp)
    plot_frequency_response(freq, G.fresp, MAX_FREQUENCY_HZ)
    plot_step_response(step_resp, MAX_FREQUENCY_HZ)
    
    print("\nAnalysis complete! Displaying plots...")
    plt.show()


if __name__ == "__main__":
    main()
