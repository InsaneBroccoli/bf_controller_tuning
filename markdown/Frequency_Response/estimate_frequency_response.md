<p align="right">
  <img src="./Images/zhaw_logo.jpg"
     alt="ZHAW Logo"
     width="100"
     style="float:right; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>

# Function `estimate_frequency_response`

The `estimate_frequency_response` function estimates the **frequency response** and **coherence** between a measured input and output signal using the **Welch method** for spectral averaging [1][2].  
It returns an amplitude-calibrated, single-sided frequency response with correct phase and energy scaling.

As an example, consider a system excited by a random or chirp signal, where both input and output are measured to identify its linear dynamics.
<p align="center">
  <img src="./Images/inp_outp_welch.jpg"
     alt="Original noisy signals"
     width="800"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>

---

## 1) Spectral Estimation and Theoretical Basis

### Derivation of the Frequency Response Formula

In the frequency domain, the relation between input and output of a linear time-invariant (LTI) system can be expressed as

\[
G(s) = \frac{Y(s)}{U(s)},
\]

where \( U(s) \) and \( Y(s) \) are the Fourier transforms of the input \( u(t) \) and output \( y(t) \).  
This ratio describes how each frequency component of the input is scaled and phase-shifted by the system.

In practice, however, direct division of \( Y(s) \) and \( U(s) \) is sensitive to noise. 
To obtain a more reliable estimate, we use averaged spectral quantities:

\[
S_{yu}(s) = E\{Y(s)U^*(s)\}, \qquad S_{uu}(s) = E\{|U(s)|^2\},
\]

which leads to the practical and statistically robust estimator

\[
G(f) = \frac{S_{yu}(s)}{S_{uu}(s)}.
\]

This formulation provides a consistent estimate of the system’s amplitude and phase response across all frequencies [1][2].

<p align="center">
  <img src="./Images/bode_Syu_Suu.jpg"
     alt="Original noisy signals"
     width="800"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>

---

## 2) The Welch Method for Averaging

The function implements the **Welch averaging method** [1, pp. 356–364], which provides a smoother and statistically robust estimate by dividing the data into overlapping, windowed segments.

1. **Segmentation:**  
   The signals are divided into short overlapping sections of length \( N_\text{est} \) with overlap.

2. **Windowing:**  
   Each segment is multiplied by a window (e.g. Hann) to reduce spectral leakage.

3. **FFT and normalization:**  
   Each windowed segment is transformed using the FFT and normalized by  
   \( W = \sum w(n)/2 \), ensuring correct amplitude calibration [1, p. 358].

4. **Power spectra formation:**  
   For each segment:
   \[
   S_{uu,k}(s) = U_k(s)U_k^*(s), \quad
   S_{yu,k}(s) = Y_k(s)U_k^*(s), \quad
   S_{yy,k}(s) = Y_k(s)Y_k^*(s).
   \]
   The spectra are then converted to **one-sided** form by doubling all interior bins and dividing the DC and Nyquist bins by 4 [1, p. 360].

5. **Averaging:**  
   The spectra from all segments are averaged to yield \( \overline{S}_{uu}, \overline{S}_{yu}, \overline{S}_{yy} \).

<p align="center">
  <img src="./Images/bode_welch.jpg"
     alt="Welch spectral estimation illustration"
     width="750"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>

---

## 3) Regularization and Robustness

To avoid division by near-zero values of \( S_{uu}(f) \), a small positive constant `delta` is added:
\[
G(s) = \frac{S_{yu}(s)}{S_{uu}(s) + \delta}.
\]
This **regularized spectral inversion** improves numerical stability [2, pp. 208–210].

---

## References

[1] P. M. Djuric, S. M. Kay, *Statistical Digital Signal Processing and Modeling*, Prentice Hall, 1993, **pp. 356–375.**  
[2] R. Pintelon, J. Schoukens, *System Identification: A Frequency Domain Approach*, 2nd ed., IEEE Press, 2012, **pp. 57–63, 205–212.**
