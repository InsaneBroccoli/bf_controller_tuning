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


def find_column(col_name: str): -> int
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

# Handle edge case where chirp is still active at end of recording
if len(chirp_start_idx) > 0 and (
    len(chirp_end_idx) == 0 or chirp_end_idx[-1] < chirp_start_idx[-1]
):
    chirp_end_idx = np.append(chirp_end_idx, len(time) - 1)

# Setpoint data (deg/sec)
setpoint_roll = data[:, col_setpoint_roll]

# Gyro unfiltered data (deg/sec)
gyro_unfilt_roll = data[:, col_gyro_unfilt_roll]

print(f"\nExtracted {len(time)} samples ({time[-1]:.2f} seconds of data)")

# =========================================================================
# 4. PROCESS DATA
# =========================================================================

# =========================================================================
# 5. PROCESS EACH CHIRP EVENT
# =========================================================================

# Calculate sampling time
