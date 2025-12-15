# Spectrogram Estimation

The function `estimate_spectrogram` extends the spectral estimation as described in [Spectra Analysis](./Spectra_Analysis.md) by computing a **single-sided, amplitude-correct spectrogram** over an additional coordinate $y$.  
This allows the frequency content of a signal to be analyzed not only over time, but also along another dimension such as height, position, angle, or system state.

Instead of forming a single averaged spectrum, the algorithm groups all FFT segments into $N_{\text{res}}$ bins along the thrust-axis and computes one spectrum per bin.  
The amplitude calibration, FFT scaling, and DC/Nyquist correction follow the same method as in `estimate_spectra`.

---

## Function Overview

### 1. Binning Along the Secondary Coordinate

The input vector $y$ is divided into $N_{\text{res}}$ equally spaced bins:

$$y_{\text{axis}} = \left[\, y_{\min},\, y_{\min} + \Delta y, \ldots, y_{\max} - \Delta y \,\right]$$
$$\Delta y = \frac{y_{\max} - y_{\min}}{N_{\text{res}}}$$

Each FFT segment is assigned to one or more of these bins depending on its corresponding \( y \)-values.

---

### 2. Segmentation, Windowing, and FFT

Segmentation, windowing with a tapered analysis window (e.g., Hann), and amplitude-correct  FFT scaling follow the same steps described in the [Spectra Analysis](./Spectra_Analysis.md) chapter. [2, pp. 455–457]

Each segment produces a one-sided power vector:

$$P_{\text{seg}}(f) = |U(f)|^2$$

with DC and Nyquist bins corrected as:

$$P_{\mathrm{DC}} = \frac{P_{\mathrm{DC}}}{4}$$

$$P_{\mathrm{Nyq}} = \frac{P_{\mathrm{Nyq}}}{4}$$

---

### 3. Accumulation Into $y$-Bins

For each segment, all samples of $y$ belonging to that segment are mapped to their corresponding $y$-bins.  
The resulting segment spectrum is added to those bins, weighted by the number of matching samples:

$$P_{\text{avg}}(y_i, f) = \frac{1}{N_i} \sum_{k \in \mathcal{S}_i} P_{\text{seg},k}(f)$$

where $N_i$ is the number of contributing segment-samples for bin $i$.

This creates a 2-D matrix:

$$P_{\text{avg}} \in \mathbb{R}^{N_{\text{res}} \times N_{\text{freq}}}$$

representing the power spectrum as a function of both frequency and \( y \) [1, pp. 45–48][2, pp. 455–457].

---

### 4. 2-D Smoothing of the Spectrum

A small weighted 3×3 convolution kernel is applied to reduce noise and produce a visually smooth spectrogram.  
The kernel preserves energy by normalizing the sum of all weights.

$$K =
\begin{bmatrix}
1 & 3 & 1 \\
3 & 5 & 3 \\
1 & 3 & 1
\end{bmatrix}$$

$$K_\text{norm} = \frac{K}{\sum K}$$

The smoothed spectrogram is obtained as:

$$P_{\mathrm{smooth}} = \frac{P_{\mathrm{avg}} * K_\text{norm}}{\mathbf{1} * K_\text{norm}}$$


## Outputs

- **$P_{avg}$** — the averaged single-sided power spectrogram  
- **$freq$** — corresponding frequency axis ($0 … \frac{f_s}{2}$)  
- **$y_{axis}$** — center coordinates of the $y$-bins  

Amplitude spectrograms can be obtained via:

$$A(y,f) = \sqrt{P_{\text{avg}}(y,f)}$$

---

## References

[1] R. Pintelon, J. Schoukens, *System Identification: A Frequency Domain Approach*,  
2nd ed., IEEE Press, 2012, **pp. 45–60.**

[2] M. H. Hayes, *Statistical Digital Signal Processing and Modeling*,  
John Wiley & Sons, 1996, **pp. 51–52, 455–457.**

# ========================================

# Spectrogram Estimation

The function `estimate_spectrogram` extends the spectral estimation as described in [Spectra Analysis](./Spectra_Analysis. md) by computing a **single-sided, amplitude-correct spectrogram** over an additional dimension such as throttle, height, or system state.

A spectrogram is a time-frequency representation that shows how the frequency content of a signal evolves over time or along another variable. Unlike a standard spectrum that averages across the entire signal, a spectrogram reveals dynamic changes in vibration patterns, making it particularly useful for identifying throttle-dependent resonances or transient issues in FPV drones.

This technique is widely used in both offline analysis tools like this MATLAB implementation and in real-time viewers such as the Betaflight Blackbox Explorer, where it helps pilots visualize and diagnose mechanical issues or tune filters.

---

## How Spectrograms Work

### Overview of the Process

1. **Divide the signal into overlapping time segments**  
   The input data is split into short windows, each analyzed independently.

2. **Apply a windowing function**  
   A tapered window (e.g., Hann) reduces spectral leakage at segment boundaries. 

3. **Compute the FFT for each segment**  
   Each windowed segment is transformed into the frequency domain. 

4. **Calculate magnitude**  
   The magnitude or power at each frequency is computed from the FFT output.

5. **Bin the results by a secondary variable**  
   Segments are grouped by throttle, time, or another parameter to create a 2D representation.

6. **Average and smooth**  
   Multiple measurements in each bin are averaged, and optional smoothing is applied for better visualization.

---

## Detailed Algorithm

### 1. Binning Along the Secondary Coordinate

The input vector $y$ (e.g., throttle percentage) is divided into $N_{\text{res}}$ equally spaced bins:

$$y_{\text{axis}} = \left[\, y_{\min},\, y_{\min} + \Delta y, \ldots, y_{\max} - \Delta y \,\right]$$

$$\Delta y = \frac{y_{\max} - y_{\min}}{N_{\text{res}}}$$

Each FFT segment is assigned to one or more bins based on the average value of $y$ during that segment.  This allows the spectrogram to show how frequency content varies with throttle or other flight parameters [1, pp.  45-48].

**Implementation Note:**  
Betaflight Blackbox Explorer typically uses 100 throttle bins spanning 0-100% throttle.  Each FFT chunk's average throttle determines which bin receives its frequency data.

---

### 2. Segmentation and Windowing

The signal is divided into overlapping segments: 

- **Segment length:** Determined by the desired frequency resolution and sampling rate.  Betaflight typically uses chunks of a few hundred milliseconds. 
  
- **Overlap:** A moving window approach with overlap (e.g., 50-90%) increases the number of data points and improves time resolution.

- **FFT size:** Rounded to the nearest power of 2 for computational efficiency. 

A **Hann window** or similar taper is applied to each segment to minimize spectral leakage: 

$$w_{\text{Hann}}(n) = 0.5 \left(1 - \cos\left(\frac{2\pi n}{N-1}\right)\right)$$

The windowed signal becomes: 

$$x_{\text{windowed}}(n) = x(n) \cdot w(n)$$

This step is critical for accurate frequency estimation, as sharp transitions at segment edges would otherwise spread energy across multiple frequency bins [2, pp. 51-52].

---

### 3. FFT Computation and Magnitude Calculation

Each windowed segment is transformed via FFT:

$$X(k) = \sum_{n=0}^{N-1} x_{\text{windowed}}(n) \cdot e^{-j2\pi kn/N}$$

Since the gyro data is real-valued, only the **first half** of the FFT output (one-sided spectrum) is used:

$$N_{\text{freq}} = \lfloor N_{\text{FFT}}/2 \rfloor$$

The **magnitude** at each frequency bin is computed as: 

$$|X(k)| = \sqrt{\text{Re}(X(k))^2 + \text{Im}(X(k))^2}$$

or equivalently, as implemented in Betaflight: 

```javascript
magnitudes[i] = Math.hypot(re, im);
```

The **power** is the squared magnitude: 

$$P_{\text{seg}}(f) = |X(f)|^2$$

**DC and Nyquist Correction:**  
The DC and Nyquist frequency bins appear only once in the one-sided spectrum, so their power is divided by 4 for amplitude-correct scaling:

$$P_{\mathrm{DC}} = \frac{P_{\mathrm{DC}}}{4}, \quad P_{\mathrm{Nyq}} = \frac{P_{\mathrm{Nyq}}}{4}$$

This correction ensures the amplitude spectrum matches the original signal amplitude [1, pp. 45-48]. 

---

### 4. Accumulation Into Bins

For each FFT segment, the corresponding $y$-values (e.g., throttle samples) are mapped to their bins.  The segment's power spectrum is added to those bins and averaged:

$$P_{\text{avg}}(y_i, f) = \frac{1}{N_i} \sum_{k \in \mathcal{S}_i} P_{\text{seg},k}(f)$$

where:
- $N_i$ is the number of FFT segments contributing to bin $i$
- $\mathcal{S}_i$ is the set of segments assigned to bin $i$

This creates a 2D matrix:

$$P_{\text{avg}} \in \mathbb{R}^{N_{\text{res}} \times N_{\text{freq}}}$$

representing power as a function of both **frequency** (horizontal) and **throttle or time** (vertical) [1, pp. 45-48][2, pp. 455-457]. 

---

### 5. 2D Smoothing (Optional)

To reduce noise and improve visual clarity, a small weighted convolution kernel is applied:

$$K = \begin{bmatrix}
1 & 3 & 1 \\
3 & 5 & 3 \\
1 & 3 & 1
\end{bmatrix}$$

The kernel is normalized to preserve energy:

$$K_{\text{norm}} = \frac{K}{\sum K}$$

The smoothed spectrogram is:

$$P_{\text{smooth}} = \frac{P_{\text{avg}} * K_{\text{norm}}}{\mathbf{1} * K_{\text{norm}}}$$

Betaflight applies a similar blur filter to smooth imperfections in the heatmap visualization.

---

## Visualization:  Heat Map

The 2D power matrix is displayed as a **heat map** or **spectrogram plot**:

- **Horizontal axis:** Frequency (0 to $f_s/2$ Hz)
- **Vertical axis:** Secondary variable (e.g., throttle 0-100%)
- **Color:** Amplitude or power level

### Color Mapping

In Betaflight Blackbox Explorer, colors are mapped using HSL color space: 

```javascript
const fftColorScale = 100 / (zoomY * SCALE_HEATMAP);
valuePlot = Math.round(Math.min(valuePlot * fftColorScale, 100));
canvasCtx.fillStyle = `hsl(360, 100%, ${valuePlot}%)`;
```

- **Brighter colors** (white/yellow) indicate higher vibration amplitude
- **Darker colors** (red/black) indicate lower amplitude
- Scale factor of 1.1 adjusts the dynamic range for better visibility

---

## Frequency Resolution

The frequency resolution is determined by: 

$$\Delta f = \frac{f_s}{2 \cdot N_{\text{freq}}}$$

or equivalently:

$$\text{Frequency step} = \frac{0.5 \cdot f_s}{N_{\text{FFT}}}$$

where $f_s$ is the sampling rate (e.g., blackbox logging rate).

Higher frequency resolution requires longer time segments, but this reduces time resolution.  This is a fundamental tradeoff in spectral analysis [2, pp. 51-52].

---

## Outputs

- **Pavg** — Averaged single-sided power spectrogram matrix
- **freq** — Frequency axis vector (0 to $f_s/2$ Hz)
- **y_axis** — Center coordinates of the secondary variable bins

**Amplitude spectrograms** can be obtained via:

$$A(y, f) = \sqrt{P_{\text{avg}}(y, f)}$$

---

## Practical Use in Drone Tuning

Spectrograms are particularly useful for:

1. **Identifying throttle-dependent resonances**  
   Vertical streaks indicate frequencies that appear at specific throttle levels. 

2. **Detecting motor issues**  
   Uneven patterns across motors suggest mechanical problems such as bent shafts or unbalanced props.

3. **Filter tuning**  
   Shows which frequency ranges need notch filtering or dynamic filtering.

4. **Validating fixes**  
   Before and after comparison shows whether mechanical or software changes reduced vibrations.

5. **Diagnosing prop wash**  
   Characteristic patterns can indicate aerodynamic disturbances.

---

## References

[1] R. Pintelon, J. Schoukens, *System Identification:  A Frequency Domain Approach*,  
2nd ed., IEEE Press, 2012, **pp.  45-60.**

[2] M.  H. Hayes, *Statistical Digital Signal Processing and Modeling*,  
John Wiley & Sons, 1996, **pp. 51-52, 455-457.**

---

## Additional Resources

- [Betaflight Blackbox Explorer Source Code](https://github.com/betaflight/blackbox-log-viewer)
- [Betaflight Blackbox Logging Guide](https://www.betaflight.com/docs/wiki/guides/current/Black-Box-logging-and-usage)
- [Spectra Analysis Documentation](./Spectra_Analysis.md)