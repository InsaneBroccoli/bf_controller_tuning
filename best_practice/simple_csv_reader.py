"""
Simple CSV Reader - Extract Setpoint and GyroUnfilt Data
No external functions required (only pandas, numpy, matplotlib)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# from pathlib import Path
import time as time_module
import pickle

# =========================================================================
# 1. DEFINE FILE PATH
# =========================================================================

file_path = "../20250908/20250908_flipmini_00.bbl.csv"

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


def find_column(col_names_to_try):
    """Find column index by trying multiple possible names"""
    for col_name in col_names_to_try:
        try:
            return column_names.index(col_name)
        except ValueError:
            continue
    return None


# Print all column names for debugging
print(f"\nAvailable columns ({len(column_names)} total):")
print(column_names[:20])  # Print first 20 columns
print("...")

# Try both quoted and unquoted versions
col_time = find_column(["time", '"time"'])
col_setpoint_roll = find_column(["setpoint[0]", '"setpoint[0]"'])
col_setpoint_pitch = find_column(["setpoint[1]", '"setpoint[1]"'])
col_setpoint_yaw = find_column(["setpoint[2]", '"setpoint[2]"'])
col_gyro_unfilt_roll = find_column(["gyroUnfilt[0]", '"gyroUnfilt[0]"'])
col_gyro_unfilt_pitch = find_column(["gyroUnfilt[1]", '"gyroUnfilt[1]"'])
col_gyro_unfilt_yaw = find_column(["gyroUnfilt[2]", '"gyroUnfilt[2]"'])
col_debug_chirp = find_column(["debug[0]", '"debug[0]"'])

print("\nColumn indices:")
print(f"  time: {col_time}")
print(f"  debug[0] (chirp): {col_debug_chirp}")
print(f"  setpoint: [{col_setpoint_roll}, {col_setpoint_pitch}, {col_setpoint_yaw}]")
print(
    f"  gyroUnfilt: [{col_gyro_unfilt_roll}, {col_gyro_unfilt_pitch}, {col_gyro_unfilt_yaw}]"
)

# Check if we found all required columns
required_cols = {
    "time": col_time,
    "debug[0]": col_debug_chirp,
    "setpoint[0]": col_setpoint_roll,
    "setpoint[1]": col_setpoint_pitch,
    "setpoint[2]": col_setpoint_yaw,
    "gyroUnfilt[0]": col_gyro_unfilt_roll,
    "gyroUnfilt[1]": col_gyro_unfilt_pitch,
    "gyroUnfilt[2]": col_gyro_unfilt_yaw,
}

missing_cols = [name for name, idx in required_cols.items() if idx is None]
if missing_cols:
    print(f"\nERROR: Could not find the following columns: {missing_cols}")
    print("\nPlease check the column names in your CSV file.")
    exit(1)

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


def tfestimate(x, y, window, noverlap, nfft, fs):
    """
    Estimate transfer function using Welch's method
    Similar to MATLAB's tfestimate
    """
    # Number of segments
    step = len(window) - noverlap
    n_segments = (len(x) - noverlap) // step

    # Initialize accumulators (use complex dtype)
    Pxy = np.zeros(nfft // 2 + 1, dtype=complex)
    Pxx = np.zeros(nfft // 2 + 1, dtype=complex)

    # Process each segment
    for i in range(n_segments):
        start = i * step
        end = start + len(window)

        if end > len(x):
            break

        # Extract and window the segments
        x_seg = x[start:end] * window
        y_seg = y[start:end] * window

        # Compute FFT
        X = np.fft.fft(x_seg, nfft)
        Y = np.fft.fft(y_seg, nfft)

        # Accumulate cross and auto power spectral densities
        # Only keep positive frequencies
        Pxy += X[: nfft // 2 + 1] * np.conj(Y[: nfft // 2 + 1])
        Pxx += X[: nfft // 2 + 1] * np.conj(X[: nfft // 2 + 1])

    # Average
    Pxy /= n_segments
    Pxx /= n_segments

    # Transfer function estimate
    # Pxx should be real (it's an auto-correlation), so we take the real part
    H = Pxy / (np.real(Pxx) + 1e-10)  # Add small value to avoid division by zero

    # Frequency vector
    f = np.fft.fftfreq(nfft, 1 / fs)[: nfft // 2 + 1]

    return H, f


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
    # 7. ESTIMATE TRANSFER FUNCTION
    # =====================================================================

    # Window settings for spectral estimation
    window_length = min(512, len(u_roll) // 4)
    window = np.hanning(window_length)
    noverlap = window_length // 2
    nfft = max(1024, 2 ** int(np.ceil(np.log2(len(u_roll)))))

    # Estimate transfer functions
    H_roll, f = tfestimate(u_roll, y_roll, window, noverlap, nfft, fs)
    H_pitch, _ = tfestimate(u_pitch, y_pitch, window, noverlap, nfft, fs)
    H_yaw, _ = tfestimate(u_yaw, y_yaw, window, noverlap, nfft, fs)

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

# Show all plots
plt.show(block=False)

print("\n=== Processing Complete ===")

time_module.sleep(60)
plt.close("all")
