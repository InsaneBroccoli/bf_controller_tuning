import numpy as np
from scipy import signal
from control import frd, FRD
from typing import Tuple


def estimate_frequency_response(
    inp: np.ndarray,
    out: np.ndarray,
    window: np.ndarray,
    Noverlap: int,
    Nest: int,
    Ts: float,
    delta: float = 0.0,
) -> Tuple[FRD, FRD, np.ndarray, np.ndarray]:
    # Todo: check if it is usefull to remove mean here
    inp = inp - np.mean(inp)
    out = out - np.mean(out)

    Ndata = inp.shape[0]

    # factor 2 so that the magnitude corresponds to a single sided spectrum
    # 2.3*sin(2*pi*f0*time) <=> sqrt(puu(f0)) = 2.3
    W = np.sum(window) / Nest / 2

    Pavg = np.zeros((Nest, 3), dtype=np.complex128)

    Navg = 0
    ind_start = 0
    ind_end = Nest
    Ndelta = Nest - Noverlap
    while ind_end <= Ndata:
        ind = np.arange(ind_start, ind_end)

        inp_act = inp[ind]
        out_act = out[ind]

        # Todo: check if it is usefull to remove mean here
        inp_act = inp_act - np.mean(inp_act)
        out_act = out_act - np.mean(out_act)

        # needed squeeze to make element wise multiplication work correctly
        # inp_act = np.squeeze(inp_act)
        # out_act = np.squeeze(out_act)

        inp_act = window * inp_act  # type: ignore
        out_act = window * out_act  # type: ignore

        U = np.fft.fft(inp_act) / (Nest * W)
        Y = np.fft.fft(out_act) / (Nest * W)

        Pavg += np.vstack([U * np.conj(U), Y * np.conj(U), Y * np.conj(Y)]).T
        Navg += 1

        ind_start += Ndelta
        ind_end += Ndelta

    Pavg /= Navg

    g, c = _calc_freqresp_and_cohere(Pavg, delta)
    df = 1 / (Nest * Ts)
    freq = np.arange(0, 1 / Ts, df)

    # Create the FRD (Frequency Response Data) objects
    G = frd(g, freq)
    C = frd(c, freq)

    return G, C, freq, Pavg


def _calc_freqresp_and_cohere(
    P: np.ndarray, delta: float
) -> Tuple[np.ndarray, np.ndarray]:
    P[:, 0] = P[:, 0] + delta
    g = P[:, 1] / P[:, 0]
    c = np.abs(P[:, 1]) ** 2 / (P[:, 0] * P[:, 2])

    return g, c


def estimate_frf_and_coherence(
    inp: np.ndarray,
    out: np.ndarray,
    fs: float,
    window: np.ndarray,
    nperseg: int,
    noverlap: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Estimate cross spectral density and power spectral densities
    freq, Pxy = signal.csd(
        inp,
        out,
        fs=fs,
        window=window,  # type: ignore
        nperseg=nperseg,
        noverlap=noverlap,
    )
    _, Pxx = signal.csd(
        inp,
        inp,
        fs=fs,
        window=window,  # type: ignore
        nperseg=nperseg,
        noverlap=noverlap,
    )
    _, Pyy = signal.csd(
        out,
        out,
        fs=fs,
        window=window,  # type: ignore
        nperseg=nperseg,
        noverlap=noverlap,
    )

    # Calculate frequency response function
    g = Pxy / Pxx

    # Calculate coherence
    c = np.abs(Pxy) ** 2 / (Pxx * Pyy)

    # Truncate DC (freq=0) to avoid divide-by-zero issues
    return freq[1:], g[1:], c[1:]


def calculate_step_response_from_frd(G_frd: FRD, f_max_hz: float) -> np.ndarray:
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
