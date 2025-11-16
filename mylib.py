import numpy as np
from scipy import signal
from control import frd
from typing import Tuple


def estimate_frequency_response(
    inp: np.ndarray,
    out: np.ndarray,
    fs: float,
    window: np.ndarray,
    nperseg: int,
    noverlap: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Estimate cross spectral density and power spectral densities
    freq, Pxy = signal.csd(
        inp, out, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap
    )
    _, Pxx = signal.csd(
        inp, inp, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap
    )
    _, Pyy = signal.csd(
        out, out, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap
    )

    # Calculate frequency response function
    g = Pxy / Pxx

    # Calculate coherence
    c = np.abs(Pxy) ** 2 / (Pxx * Pyy)

    # Truncate DC (freq=0) to avoid divide-by-zero issues
    return freq[1:], g[1:], c[1:]


def calculate_step_response_from_frd(G_frd: frd, f_max_hz: float) -> np.ndarray:
    # Extract complex frequency response
    g = G_frd.magnitude.flatten() * np.exp(1j * G_frd.phase.flatten())

    # Reconstruct DC (simulate symmetry at zero freq)
    g_dc = g[0]  # Use g[0] again as a placeholder for DC
    g = np.insert(g, 0, g_dc)  # Prepend DC component

    # Extend frequency vector accordingly
    freq = G_frd.frequency / (2 * np.pi)
    freq = np.insert(freq, 0, 0.0)

    # Zero out above f_max_hz
    g[freq > f_max_hz] = 0

    # Construct full symmetric spectrum
    g_full = np.concatenate([g, np.conj(g[-2:0:-1])])

    # Step response is cumulative sum of real part of IFFT
    step_resp = np.cumsum(np.real(np.fft.ifft(g_full)))

    return step_resp
