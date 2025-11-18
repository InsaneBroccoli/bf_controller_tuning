# Spectral Analysis of Gyro and Control Signals

To better understand the dynamic behavior of a system, it is often useful to analyze how it responds to inputs at various frequencies. In this implementation, the `estimate_spectra` function is used to perform this frequency-domain analysis [1, pp. 45–48].

## Overview of the Algorithm

The function processes each signal column independently and performs the following steps:

### 1. Mean Removal

A global mean is removed from every signal. Additionally, each segment undergoes per-segment mean removal, which suppresses low-frequency drift and numerical bias.(DC components are therefore intentionally attenuated.)

### 2) Windowing of Each Segment

Before the FFT is computed, each data segment is multiplied by the chosen analysis window (e.g., a periodic Hann window). The window reduces unwanted frequency spreading by making the signal smoothly fade in and out at the edges of each segment. Without this tapering, sudden jumps at the segment boundaries would create artificial frequency components that do not actually exist in the original signal.

\[
x_{\mathrm{win}}(n) = x(n) \cdot w(n), \qquad n = 0, \ldots, \text{Nest}-1
\]

A periodic Hann window is commonly used because it provides:

- smooth tapering with minimal discontinuities at the segment boundaries  
- reduced leakage compared to rectangular windows  
- good balance of main-lobe width and side-lobe suppression  

After windowing, the segment is ready for amplitude-correct spectral estimation.

### 3. Segmentation with Overlap

The function divides the input signal into multiple analysis windows (segments) of length  \( \text{Nest} \). These segments can overlap to a configurable degree.  By default, a **90% overlap** is used to ensure smooth and consistent spectral estimates. 
The number of overlapping samples is defined as \( \text{N}_{\text{Overlap}} \).  Accordingly, the shift between two consecutive segments is:

\[
\text{N}_{\text{Shift}} = \text{Nest} - \text{N}_{\text{Overlap}}.
\]

Overlaps in the range of **50–90%** yield smoother and statistically more stable spectra,  which aligns with standard recommendations for Welch averaging [2, pp. 455–457].

### 3. Windowing and Amplitude-Correct FFT Scaling

The function uses a custom normalization:

\[
W = \frac{\sum w}{2}
\]

The FFT output is then scaled by:

\[
\frac{1}{Nest \cdot W}
\]

This scaling intentionally pre-doubles all bins so that interior bins already contain the correct single-sided amplitude representation.  
However, the **DC** and **Nyquist** bins must not be doubled and are therefore corrected:

\[
P_{\mathrm{DC}} = \frac{P_{\mathrm{DC}}}{4}, 
\qquad
P_{\mathrm{Nyq}} = \frac{P_{\mathrm{Nyq}}}{4}
\]

This correction ensures:

- correct amplitude level for sinusoidal bins  
- correct energy conservation  
- proper handling of real-valued signals with even FFT length
- 
These characteristics are well documented in frequency-domain signal analysis literature [2, pp. 50–53].

---

### 4) One-Sided Power Spectrum Construction

From the two-sided FFT power spectrum, only the positive half is used:

\[
P_{\text{1-sided}}(k) = |U(k)|^{2}
\]

The power is averaged across all segments:

\[
P_{\text{avg}}(f) =
\frac{1}{N_{\text{segments}}}
\sum P_{\text{seg}}(f)
\]

The output is:

- **Pavg**: one-sided power spectrum  
- **freq**: matching frequency vector  

To obtain single-sided amplitude spectra:

\[
A(f) = \sqrt{P_{\text{avg}}(f)}
\]

This form is used throughout the analysis for gyro data, filtered gyro signals, and control-loop sums [1, pp. 45–50].

---

## References

[1] Pintelon, R., & Schoukens, J. (2012). *System Identification: A Frequency Domain Approach*  
    (2nd ed.). Wiley-IEEE Press, pp. 45–60, 118–121.

[2] Hayes, M. H. (1996). *Statistical Digital Signal Processing and Modeling*.  
    Wiley, pp. 455–457.