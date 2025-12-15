import pickle
import time

import numpy as np
import pandas as pd
import control as ct
import scipy.signal as signal
import py_lib.pidtuninglib as tlib

from matplotlib import pyplot as plt
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple
from scipy.io import loadmat

bbl_data_path = Path("../20250908/20250908_flipmini_00.bbl.csv")


@dataclass
class FlightData:
    """Container for flight log data with column-based access."""

    values: np.ndarray
    columns: List[str]

    def __getitem__(self, col_name: str) -> np.ndarray:
        """Access column data by name."""
        if col_name not in self.columns:
            raise ValueError(f"Column {col_name} not found")
        idx = self.columns.index(col_name)
        return self.values[:, idx]

    @property
    def n_rows(self) -> int:
        return self.values.shape[0]


def get_header_index(file_path: Path) -> int:
    """Find the line number where CSV header starts (contains 'loopIteration')."""
    with file_path.open("r", encoding="utf-8") as f:
        for line_num, line in enumerate(f):
            if "loopIteration" in line:
                return line_num
    return 0


def load_flight_data(file_path: Path) -> FlightData:
    """Load flight data from CSV, using pickle cache if available."""
    cache_path = Path(str(file_path).replace(".bbl.csv", ".pkl"))

    # Load from cache if exists
    if cache_path.exists():
        print(f"loading cached data from {cache_path}")
        cache = pickle.load(cache_path.open("rb"))
        return FlightData(values=cache["data"], columns=cache["columns"])

    # Parse CSV and create cache
    print(f"Reading CSV file from {file_path}")
    time_start = time.time()

    idx_header = get_header_index(file_path)
    df = pd.read_csv(file_path, skiprows=idx_header)

    data_obj = FlightData(values=df.values, columns=df.columns.tolist())

    with cache_path.open("wb") as f:
        pickle.dump({"data": data_obj.values, "columns": data_obj.columns}, f)

    print(f"Cached data to {cache_path}")
    print(
        f"Data loaded: {data_obj.n_rows} rows in {time.time() - time_start:.2f} seconds"
    )

    return data_obj


def get_chirp_indices(chirp_signal: np.ndarray) -> Tuple[int, int]:
    """Detect start and end indices of chirp signal (rising and falling edges)."""
    chirp_binary = (chirp_signal > 0).astype(int)
    edge_diff = np.diff(chirp_binary, prepend=0)

    start_idx = np.where(edge_diff == 1)[0]  # Rising edge
    end_idx = np.where(edge_diff == -1)[0]  # Falling edge

    return start_idx[0], end_idx[0]


def test_output(
    output_data: Path, time_step: np.ndarray, step_resp: np.ndarray
) -> None:
    """Validate Python results against MATLAB reference output."""
    data_out = loadmat(output_data)

    np.testing.assert_allclose(
        time_step, data_out["t_step"].squeeze(), rtol=1e-9, atol=1e-9
    )
    np.testing.assert_allclose(
        step_resp, data_out["step_resp"].squeeze(), rtol=1e-9, atol=1e-9
    )

    print("Python output matches MATLAB output")


def main():
    start_time = time.time()

    # Load flight log data
    flight_data = load_flight_data(bbl_data_path)

    print("\n--- Available Columns ---")
    for col in flight_data.columns:
        print(f"- {col}")
    print("-------------------------\n")

    # Extract time series (convert from microseconds to seconds)
    t = (flight_data["time"] - flight_data["time"][0]) * 1.0e-6
    setpoint_roll = flight_data["setpoint[0]"]
    gyro_unfilt_roll = flight_data["gyroUnfilt[0]"]

    # Plot full dataset
    fig1, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(t, setpoint_roll)
    ax1.set_title("setpoint")
    ax1.grid(True)
    ax2.plot(t, gyro_unfilt_roll)
    ax2.set_title("gyro unfiltered")
    ax2.grid(True)
    fig1.suptitle("RAW DATA")

    # Extract chirp interval for frequency response analysis
    chirp_signal = flight_data["debug[0]"]
    chirp_start, chirp_end = get_chirp_indices(chirp_signal)

    time_chirp = t[chirp_start:chirp_end]
    setpoint_roll_chirp = setpoint_roll[chirp_start:chirp_end]
    gyro_unfilt_roll_chirp = gyro_unfilt_roll[chirp_start:chirp_end]

    # Plot chirp interval data
    fig2, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(time_chirp, setpoint_roll_chirp)
    ax1.set_title("setpoint")
    ax1.grid(True)
    ax2.plot(time_chirp, gyro_unfilt_roll_chirp)
    ax2.set_title("gyro unfiltered")
    ax2.grid(True)
    fig2.suptitle("DATA OVER CHIRP")

    # Define sampling periods for control loop and logging
    Ts = 125 * 1.0e-6  # Base sampling period (125 μs)
    Ts_cntr = 2 * Ts  # Control loop period
    Ts_log = 2 * Ts_cntr  # Logging period

    # Welch method parameters for frequency response estimation
    f_max_hz = 200
    Nest = round(2.0 / Ts_log)  # FFT window size
    koverlap = 0.9
    Noverlap = int(koverlap * Nest)  # Window overlap
    window = signal.windows.hann(Nest, sym=False)

    # Estimate frequency response from input/output data
    G, _, _, _ = tlib.estimate_frequency_response(
        setpoint_roll_chirp, gyro_unfilt_roll_chirp, window, Noverlap, Nest, Ts_log
    )

    # Plot Bode diagram
    plt.figure()
    ct.bode_plot(
        G,
        title="FREQUENCY RESPONSE",
        dB=True,
        Hz=True,
    )

    # Calculate step response from frequency response data
    step_resp = tlib.calculate_step_response_from_frd(G, f_max_hz)
    time_step = np.arange(len(step_resp)) * Ts_log

    # Plot step response
    plt.figure()
    plt.plot(time_step, step_resp)
    plt.xlim(0, 0.2)
    plt.grid(True)
    plt.title("STEP RESPONSE")

    time_elapsed = time.time() - start_time
    print(f"\nTotal execution time: {time_elapsed:.2f} seconds")

    # Verify results against MATLAB reference
    output_data = Path("./testing/output.mat")
    test_output(output_data, time_step, step_resp)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
