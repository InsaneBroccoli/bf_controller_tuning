import numpy as np
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
    """
    Welch single-sided X/Y spectra with amplitude-calibrated scaling.
    Identical implementation to the MATLAB version.
    """
    # Assumptions
    assert Nest % 2 == 0, "This implementation assumes even Nest."
    assert len(window) == Nest, "window length must equal Nest."

    # Global mean removal
    inp = inp - np.mean(inp)
    out = out - np.mean(out)

    Ndata = len(inp)

    # Frequency axis (0 .. Fs/2), length Nfreq = Nest/2+1
    fs = 1 / Ts
    freq = np.arange(Nest // 2 + 1) * (fs / Nest)
    Nfreq = len(freq)

    # Normalization (bakes single-sided doubling into all bins)
    # W = sum(window)/Nest/2  => dividing by Nest*W == dividing by (sum(window)/2)
    W = np.sum(window) / Nest / 2

    # Columns: [Suu, Syu, Syy]
    Pavg = np.zeros((Nfreq, 3), dtype=complex)

    Navg = 0
    ind_start = 0
    ind_end = Nest
    Ndelta = Nest - Noverlap

    while ind_end <= Ndata:
        ind = slice(ind_start, ind_end)

        inp_act = inp[ind]
        out_act = out[ind]

        # Optional per-segment mean removal (matches MATLAB script)
        inp_act = inp_act - np.mean(inp_act)
        out_act = out_act - np.mean(out_act)

        inp_act = window * inp_act
        out_act = window * out_act

        U = np.fft.fft(inp_act) / (Nest * W)
        Y = np.fft.fft(out_act) / (Nest * W)

        # Two-sided spectra
        # Suu = U*conj(U), Syu = Y*conj(U), Syy = Y*conj(Y)
        # Note: In Python, we stack them.
        # U * np.conj(U) results in real values, but we keep complex dtype for consistency
        Pact = np.column_stack([U * np.conj(U), Y * np.conj(U), Y * np.conj(Y)])

        # One-sided conversion with DC & Nyquist fix (power / 4)
        Pseg = Pact[:Nfreq, :]
        Pseg[0, :] = Pseg[0, :] / 4  # DC
        Pseg[-1, :] = Pseg[-1, :] / 4  # Nyquist (exists since Nest is even)

        Pavg += Pseg
        Navg += 1

        # Next segment
        ind_start += Ndelta
        ind_end += Ndelta

    if Navg > 0:
        Pavg /= Navg

    g, c = _calc_freqresp_and_cohere(Pavg, delta)

    # Return frd with frequency in Hz
    # Note: control.frd expects omega (rad/s) usually, but we can store Hz if we are consistent.
    # However, standard control library usage prefers rad/s.
    # The MATLAB script returns an FRD object with 'Units', 'Hz'.
    # The Python control library's FRD object stores omega.
    # To match MATLAB's output structure (freq in Hz), we return freq separately.
    # We will store omega in the FRD object to be compatible with other control lib functions.
    omega = freq * 2 * np.pi

    # FIX: Pass Ts as a keyword argument 'dt'
    G = frd(g, omega)  # Python FRD uses rad/s for the frequency vector
    C = frd(c, omega)

    return G, C, freq, Pavg


def _calc_freqresp_and_cohere(
    P: np.ndarray, delta: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Helper for G and coherence.
    """
    # P is [Nfreq x 3] (Suu, Syu, Syy)
    # We need to be careful not to modify P in place if it's used elsewhere,
    # but here it's a local copy/result.

    Suu = P[:, 0] + delta
    Syu = P[:, 1]
    Syy = P[:, 2]

    g = Syu / Suu
    c = np.abs(Syu) ** 2 / (Suu * Syy)

    # c should be real, but division of complex numbers might leave tiny imag part
    return g, np.real(c)


def calculate_step_response_from_frd(G: FRD, f_max_hz: float) -> np.ndarray:
    """
    Step response from frequency response data.
    Identical implementation to the MATLAB version.

    Args:
        G: FRD object containing complex frequency response samples
        f_max_hz: scalar cutoff frequency [Hz]; frequencies above are zeroed

    Returns:
        step_resp: real-valued step response approximation in time domain
    """
    # 1. Extract frequency response data
    # G.fresp is usually (n_outputs, n_inputs, n_freq). Squeeze to 1D.
    g = np.squeeze(G.fresp)

    # Extract frequency vector in Hz
    # In estimate_frequency_response, we stored omega = freq * 2*pi
    freq = G.omega / (2 * np.pi)

    # 2. Replace missing DC (NaN) with the second frequency point if needed
    # MATLAB: if isnan(abs(g(1))) g(1) = g(2); end
    if np.isnan(np.abs(g[0])):
        g[0] = g[1]

    # 3. Zero out response above f_max
    # MATLAB: g(freq > f_max) = 0;
    g[freq > f_max_hz] = 0

    # 4. Construct full symmetric spectrum
    # MATLAB: g_full = [g; conj(g(end-1:-1:2))];
    # Python: g is 0..Nyquist.
    # We want indices: Nyquist-1 down to 1 (skipping DC).
    # Python indices: -2 down to 1. Slice: [-2 : 0 : -1]
    g_mirror = np.conj(g[-2:0:-1])
    g_full = np.concatenate((g, g_mirror))

    # 5. Step response is cumulative sum of real part of IFFT
    # MATLAB: step_resp = cumsum(real(ifft(g_full)));
    step_resp = np.cumsum(np.real(np.fft.ifft(g_full)))

    return step_resp
