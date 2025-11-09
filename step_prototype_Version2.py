"""
Simple CSV Reader - Extract Setpoint and GyroUnfilt Data
Improved with SciPy for signal processing
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as time_module
import pickle
from scipy import signal

# =========================================================================
# 1. DEFINE FILE PATH
# =========================================================================

file_path = "./20250908/20250908_flipmini_00.bbl.csv"

# =========================================================================
# 2. READ THE ENTIRE FILE
# =========================================================================

print("Reading CSV file...")
start_time = time_module.time()

# Read CSV file - pandas will automatically handle the header
# First, we need to find where the actual data starts
with open(file_path, "r") as f:
    line_count = 0
    for line in f:
        line_count += 1
        # Look for loopIteration with or without quotes
        if "loopIteration" in line:
            header_line = line_count - 1  # 0-indexed for pandas
            break

# Load or read data
mat_file_path = file_path.replace(".bbl.csv", ".pkl")

try:
    # Try to load preprocessed data
    with open(mat_file_path, "rb") as f:
        data_dict = pickle.load(f)
        data = data_dict["data"]
        column_names = data_dict["columns"]
    print(f"Loaded preprocessed data from {mat_file_path}")
except FileNotFoundError:
    # Read the CSV file
    df = pd.read_csv(file_path, skiprows=header_line)

    # Store column names and convert to numpy array
    column_names = df.columns.tolist()
    data = df.values

    # Save for future use
    with open(mat_file_path, "wb") as f:
        pickle.dump({"data": data, "columns": column_names}, f)
    print(f"Saved preprocessed data to {mat_file_path}")

elapsed_time = time_module.time() - start_time
print(f"Data loaded: {data.shape[0]} rows in {elapsed_time:.2f} seconds")


# =========================================================================
# 3. FIND COLUMN INDICES
# =========================================================================
def find_column(col_name):
    """Find column index by name"""
    try:
        return column_names.index(col_name)
    except ValueError:
        return None


col_time = find_column("time")
col_setpoint_roll = find_column("setpoint[0]")
col_setpoint_pitch = find_column("setpoint[1]")
col_setpoint_yaw = find_column("setpoint[2]")
col_gyro_unfilt_roll = find_column("gyroUnfilt[0]")
col_gyro_unfilt_pitch = find_column("gyroUnfilt[1]")
col_gyro_unfilt_yaw = find_column("gyroUnfilt[2]")
col_debug_chirp = find_column("debug[0]")


print("\nColumn indices:")
print(f"  time: {col_time}")
print(f"  debug[0] (chirp): {col_debug_chirp}")
print(f"  setpoint: [{col_setpoint_roll}, {col_setpoint_pitch}, {col_setpoint_yaw}]")
print(
    f"  gyroUnfilt: [{col_gyro_unfilt_roll}, {col_gyro_unfilt_pitch}, {col_gyro_unfilt_yaw}]"
)

# =========================================================================
# 4. EXTRACT THE DATA
# =========================================================================

# Time vector (convert from microseconds to seconds)
time = (data[:, col_time] - data[0, col_time]) * 1.0e-6

# Extract chirp signal
chirp = data[:, col_debug_chirp]

# Detect chirp start and end (transition from 0 to non-zero)
chirp_binary = (chirp > 0).astype(int)
chirp_diff = np.diff(np.concatenate([[0], chirp_binary]))  # prepend 0 to align indices
chirp_start_idx = np.where(chirp_diff == 1)[0]  # rising edge
chirp_end_idx = np.where(chirp_diff == -1)[0]  # falling edge

# Handle edge case where chirp is still active at end of recording
if len(chirp_start_idx) > 0 and (
    len(chirp_end_idx) == 0 or chirp_end_idx[-1] < chirp_start_idx[-1]
):
    chirp_end_idx = np.append(chirp_end_idx, len(time) - 1)

print(f"\nFound {len(chirp_start_idx)} chirp events:")
for i in range(len(chirp_start_idx)):
    print(
        f"  Chirp {i + 1}: Start at {time[chirp_start_idx[i]]:.3f} s (idx {chirp_start_idx[i]}), "
        f"End at {time[chirp_end_idx[i]]:.3f} s (idx {chirp_end_idx[i]}), "
        f"Duration: {time[chirp_end_idx[i]] - time[chirp_start_idx[i]]:.3f} s"
    )

# Setpoint data (deg/sec)
setpoint_roll = data[:, col_setpoint_roll]
setpoint_pitch = data[:, col_setpoint_pitch]
setpoint_yaw = data[:, col_setpoint_yaw]

# Gyro unfiltered data (deg/sec)
gyro_unfilt_roll = data[:, col_gyro_unfilt_roll]
gyro_unfilt_pitch = data[:, col_gyro_unfilt_pitch]
gyro_unfilt_yaw = data[:, col_gyro_unfilt_yaw]

print(f"\nExtracted {len(time)} samples ({time[-1]:.2f} seconds of data)")

# =========================================================================
# 5. PROCESS EACH CHIRP EVENT
# =========================================================================

# Calculate sampling time
dt = np.median(np.diff(time))  # Use median for robustness
fs = 1 / dt
print(f"\nSampling rate: {fs:.2f} Hz (dt = {dt:.6f} s)")


def estimate_frequency_response(x, y, fs, nperseg):
    """
    Estimates the frequency response using SciPy's Welch method.
    Returns the frequency, transfer function (H), and coherence (C).
    
    Parameters:
    -----------
    x : array_like
        Input signal
    y : array_like
        Output signal
    fs : float
        Sampling frequency
    nperseg : int
        Length of each segment for Welch's method
        
    Returns:
    --------
    freq : ndarray
        Frequency vector
    H : ndarray (complex)
        Transfer function estimate
    C : ndarray (real)
        Coherence function
    """
    # Use hanning window, 50% overlap is standard
    noverlap = nperseg // 2
    window = signal.windows.hann(nperseg)

    # Cross-spectral density (input to output)
    freq, Pxy = signal.csd(x, y, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap)
    
    # Power-spectral density of the input
    _, Pxx = signal.welch(x, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap)
    
    # Power-spectral density of the output (for coherence)
    _, Pyy = signal.welch(y, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap)

    # Transfer function H(f) = Pxy(f) / Pxx(f)
    H = Pxy / (Pxx + 1e-10)
    
    # Coherence C(f) = |Pxy(f)|^2 / (Pxx(f) * Pyy(f))
    # Values close to 1 indicate high correlation at that frequency
    C = np.abs(Pxy)**2 / (Pxx * Pyy + 1e-10)

    return freq, H, C


# Process each chirp event
for chirp_num in range(len(chirp_start_idx)):
    print(f"\n--- Processing Chirp Event {chirp_num + 1} ---")

    # Extract data for this chirp period
    idx_start = chirp_start_idx[chirp_num]
    idx_end = chirp_end_idx[chirp_num]

    # Extract time-aligned data
    time_chirp = time[idx_start : idx_end + 1] - time[idx_start]  # Zero-referenced time

    u_roll = setpoint_roll[idx_start : idx_end + 1].copy()
    u_pitch = setpoint_pitch[idx_start : idx_end + 1].copy()
    u_yaw = setpoint_yaw[idx_start : idx_end + 1].copy()

    y_roll = gyro_unfilt_roll[idx_start : idx_end + 1].copy()
    y_pitch = gyro_unfilt_pitch[idx_start : idx_end + 1].copy()
    y_yaw = gyro_unfilt_yaw[idx_start : idx_end + 1].copy()

    # Remove DC offset
    u_roll -= np.mean(u_roll)
    u_pitch -= np.mean(u_pitch)
    u_yaw -= np.mean(u_yaw)

    y_roll -= np.mean(y_roll)
    y_pitch -= np.mean(y_pitch)
    y_yaw -= np.mean(y_yaw)

    print(f"Chirp duration: {time_chirp[-1]:.3f} s ({len(time_chirp)} samples)")

    # =====================================================================
    # 6. PLOT RAW INPUT/OUTPUT SIGNALS
    # =====================================================================

    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    fig.suptitle(f"Chirp {chirp_num + 1} - Raw Signals", fontsize=14, fontweight="bold")

    # Roll
    axes[0].plot(time_chirp, u_roll, "b", linewidth=1.5, label="Setpoint (Input)")
    axes[0].plot(time_chirp, y_roll, "r", linewidth=1.5, label="Gyro Unfilt (Output)")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Rate (deg/s)")
    axes[0].set_title(f"Roll - Chirp Event {chirp_num + 1}")
    axes[0].legend()
    axes[0].grid(True)

    # Pitch
    axes[1].plot(time_chirp, u_pitch, "b", linewidth=1.5, label="Setpoint (Input)")
    axes[1].plot(time_chirp, y_pitch, "r", linewidth=1.5, label="Gyro Unfilt (Output)")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Rate (deg/s)")
    axes[1].set_title(f"Pitch - Chirp Event {chirp_num + 1}")
    axes[1].legend()
    axes[1].grid(True)

    # Yaw
    axes[2].plot(time_chirp, u_yaw, "b", linewidth=1.5, label="Setpoint (Input)")
    axes[2].plot(time_chirp, y_yaw, "r", linewidth=1.5, label="Gyro Unfilt (Output)")
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Rate (deg/s)")
    axes[2].set_title(f"Yaw - Chirp Event {chirp_num + 1}")
    axes[2].legend()
    axes[2].grid(True)

    plt.tight_layout()

    # =====================================================================
    # 7. ESTIMATE TRANSFER FUNCTION USING SCIPY
    # =====================================================================

    # Window settings for spectral estimation
    window_length = min(512, len(u_roll) // 4)

    # Estimate transfer functions and coherence using SciPy
    f, H_roll, C_roll = estimate_frequency_response(u_roll, y_roll, fs, nperseg=window_length)
    _, H_pitch, C_pitch = estimate_frequency_response(u_pitch, y_pitch, fs, nperseg=window_length)
    _, H_yaw, C_yaw = estimate_frequency_response(u_yaw, y_yaw, fs, nperseg=window_length)

    # =====================================================================
    # 8. PLOT FREQUENCY RESPONSE (BODE-STYLE)
    # =====================================================================

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    fig.suptitle(
        f"Chirp {chirp_num + 1} - Frequency Response", fontsize=14, fontweight="bold"
    )

    # Magnitude
    axes[0].semilogx(
        f, 20 * np.log10(np.abs(H_roll) + 1e-10), "b", linewidth=1.5, label="Roll"
    )
    axes[0].semilogx(
        f, 20 * np.log10(np.abs(H_pitch) + 1e-10), "r", linewidth=1.5, label="Pitch"
    )
    axes[0].semilogx(
        f, 20 * np.log10(np.abs(H_yaw) + 1e-10), "g", linewidth=1.5, label="Yaw"
    )
    axes[0].set_xlabel("Frequency (Hz)")
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].set_title(f"Transfer Function Magnitude - Chirp Event {chirp_num + 1}")
    axes[0].legend()
    axes[0].grid(True, which="both", alpha=0.3)

    # Phase
    axes[1].semilogx(
        f, np.unwrap(np.angle(H_roll)) * 180 / np.pi, "b", linewidth=1.5, label="Roll"
    )
    axes[1].semilogx(
        f, np.unwrap(np.angle(H_pitch)) * 180 / np.pi, "r", linewidth=1.5, label="Pitch"
    )
    axes[1].semilogx(
        f, np.unwrap(np.angle(H_yaw)) * 180 / np.pi, "g", linewidth=1.5, label="Yaw"
    )
    axes[1].set_xlabel("Frequency (Hz)")
    axes[1].set_ylabel("Phase (deg)")
    axes[1].set_title("Transfer Function Phase")
    axes[1].legend()
    axes[1].grid(True, which="both", alpha=0.3)

    plt.tight_layout()

    # =====================================================================
    # 8b. PLOT COHERENCE
    # =====================================================================

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    fig.suptitle(
        f"Chirp {chirp_num + 1} - Coherence", fontsize=14, fontweight="bold"
    )

    ax.semilogx(f, C_roll, "b", linewidth=1.5, label="Roll")
    ax.semilogx(f, C_pitch, "r", linewidth=1.5, label="Pitch")
    ax.semilogx(f, C_yaw, "g", linewidth=1.5, label="Yaw")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Coherence")
    ax.set_title(f"Coherence (1 = perfect correlation) - Chirp Event {chirp_num + 1}")
    ax.set_ylim([0, 1.05])
    ax.axhline(y=0.5, color='k', linestyle='--', linewidth=0.8, alpha=0.5, label='Threshold (0.5)')
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)

    plt.tight_layout()

    # =====================================================================
    # 9. COMPUTE IMPULSE RESPONSE AND STEP RESPONSE
    # =====================================================================

    # Convert frequency response to impulse response
    # Create full spectrum (negative frequencies are complex conjugate)
    H_roll_full = np.concatenate([H_roll, np.conj(H_roll[-2:0:-1])])
    H_pitch_full = np.concatenate([H_pitch, np.conj(H_pitch[-2:0:-1])])
    H_yaw_full = np.concatenate([H_yaw, np.conj(H_yaw[-2:0:-1])])

    h_roll = np.real(np.fft.ifft(H_roll_full))
    h_pitch = np.real(np.fft.ifft(H_pitch_full))
    h_yaw = np.real(np.fft.ifft(H_yaw_full))

    # Truncate to reasonable length
    n_impulse = min(500, len(h_roll) // 2)
    h_roll = h_roll[:n_impulse]
    h_pitch = h_pitch[:n_impulse]
    h_yaw = h_yaw[:n_impulse]

    # Compute step response (cumulative sum of impulse response)
    step_roll = np.cumsum(h_roll) * dt
    step_pitch = np.cumsum(h_pitch) * dt
    step_yaw = np.cumsum(h_yaw) * dt

    time_step = np.arange(n_impulse) * dt

    # =====================================================================
    # 10. PLOT STEP RESPONSES
    # =====================================================================

    fig, axes = plt.subplots(3, 1, figsize=(12, 8))
    fig.suptitle(
        f"Chirp {chirp_num + 1} - Step Response", fontsize=14, fontweight="bold"
    )

    axes[0].plot(time_step, step_roll, "b", linewidth=1.5)
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Amplitude")
    axes[0].set_title(f"Roll Step Response - Chirp Event {chirp_num + 1}")
    axes[0].grid(True)

    axes[1].plot(time_step, step_pitch, "r", linewidth=1.5)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Amplitude")
    axes[1].set_title(f"Pitch Step Response - Chirp Event {chirp_num + 1}")
    axes[1].grid(True)

    axes[2].plot(time_step, step_yaw, "g", linewidth=1.5)
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Amplitude")
    axes[2].set_title(f"Yaw Step Response - Chirp Event {chirp_num + 1}")
    axes[2].grid(True)

    plt.tight_layout()

print("\n=== Processing Complete ===")

# Show all plots
plt.show()