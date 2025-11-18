import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time as time_module
import pickle
import mylib
from scipy.signal.windows import hann
from control import frequency_response
from typing import Optional

# =========================================================================
# 1. DEFINE FILE PATH
# =========================================================================

file_path = "./20250908/20250908_flipmini_00.bbl.csv"

# =========================================================================
# 2. READ THE ENTIRE FILE
# =========================================================================
header_line: int = 0  # Line number where data starts
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


def find_column(col_name: str) -> Optional[int]:
    """Find column index by name"""
    try:
        return column_names.index(col_name)
    except ValueError:
        return None


col_time = find_column("time")
col_setpoint_roll = find_column("setpoint[0]")
col_gyro_unfilt_roll = find_column("gyroUnfilt[0]")
col_debug_chirp = find_column("debug[0]")

# =========================================================================
# 3. EXTRACT THE DATA
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

chirp_start = chirp_start_idx[0]
chirp_end = chirp_end_idx[0]


print(chirp_start, chirp_end)

# extract raw data
setpoint_roll = data[:, col_setpoint_roll]
gyro_unfilt_roll = data[:, col_gyro_unfilt_roll]

# Extract data over first chirp
time_over_chirp = time[chirp_start:chirp_end]
setpoint_over_chirp = setpoint_roll[chirp_start:chirp_end]
gyro_over_chirp = gyro_unfilt_roll[chirp_start:chirp_end]

print(f"\nExtracted {len(time)} samples ({time[-1]:.2f} seconds of data)")

# =========================================================================
# 4. PROCESS DATA
# =========================================================================
# Calculate sampling time
Ts = 125 * 1.0e-6
Ts_cntr = 2 * Ts
Ts_log = 2 * Ts_cntr

# Parameters
fs = 1 / Ts
f_max_hz = 650
Nest = round(2.0 / Ts_log)
koverlap = 0.9
Noverlap = round(koverlap * Nest)
window = hann(Nest)

G, C, freq, _ = mylib.estimate_frequency_response(
    setpoint_over_chirp, gyro_over_chirp, window, Noverlap, Nest, Ts_log
)
step_resp = mylib.calculate_step_response_from_frd(G, f_max_hz)

# =========================================================================
# 5. PLOT DATA
# =========================================================================
# PLOT 1
fig1, (ax1, ax2) = plt.subplots(2, 1)

ax1.plot(time, setpoint_roll)
ax1.set_title("setpoint")

ax2.plot(time, gyro_unfilt_roll)
ax2.set_title("gyro unfiltered")

fig1.suptitle("raw data")

# PLOT 2
fig2, (ax1, ax2) = plt.subplots(2, 1)

ax1.plot(time_over_chirp, setpoint_over_chirp)
ax1.set_title("setpoint")

ax2.plot(time_over_chirp, gyro_over_chirp)
ax2.set_title("gyro unfiltered")

fig2.suptitle("data over chirp")

# PLOT 3
G_squeezed = np.squeeze(G.fresp)
plt.figure()
plt.plot(freq, 20 * np.log10(G_squeezed))
plt.xlim(0, f_max_hz)
# plt.ylim(-10.6, 1.2)
plt.suptitle("Frequency Response")

# PLOT 4
plt.figure()
plt.plot(step_resp)
plt.xlim(0, f_max_hz)
plt.suptitle("Frequency Response")

plt.tight_layout()
plt.show()
