# Determination of the Transfer Functions
For the calculation of the transfer functions, corresponding input and output signals were selected from the measurement data.  
To reduce the influence of noise and disturbances, the functions [`apply_rotfiltfilt`](./ApplyRotRotfilt.md) and [`estimate_frequency_response`](./estimate_frequency_response.md) were used. (More information in [Function apply_rotfilfil](./ApplyRotfiltfilt.md) and [Function estimate_frequency_response](./estimate_frequency_response.md))
The first function shifts the signal to the baseband and applies zero-phase filtering to suppress unwanted components, while the second performs the frequency response estimation using the Welch method with amplitude-correct scaling [1][2].  

Additionally, care was taken to use as few signals as possible from the closed loop, since noise in those signals is fed back and can therefore be amplified.  
With these measures, a robust and noise-resistant estimation of the transfer functions could be achieved.

---

## Measured Transfer Paths

In principle, the transfer function $P(j\omega)$ can be determined by measuring the input and output signals and calculating the spectrum using the **Fast Fourier Transform (FFT)**:
$G_{yv}(j\omega) = P(j\omega) = \frac{Y(j\omega)}{U(j\omega)}$

However, the resulting measurement is affected by disturbances that are present throughout the closed loop. For this reason, measurements are usually taken **outside the closed loop** to minimize the influence of noise on the signals.

The individual transfer functions are related through the closed-loop relationships:

$$G_{yu}(j\omega)= \frac{P(j\omega)}{1 + P(j\omega)C_{D}(j\omega)G_{fD}}$$
$$G_{yw}(j\omega) = T(j\omega) = \frac{G_{yu}(j\omega)\,C_{PI}(j\omega)}{1 + G_{yu}(j\omega)\,C_{PI}(j\omega)}$$
$$G_{uw}(j\omega) = \frac{C_{PI}(j\omega)}{1 + G_{yu}(j\omega)\,C_{PI}(j\omega)}$$

$$\Rightarrow P(j\omega) = \frac{T(j\omega)}{G_{wu}(j\omega)}$$

These relationships make it possible to **reconstruct the plant transfer function \( P(j\omega) \) from closed-loop measurement data** without opening the control loop during the experiment.
<p align="center">
  <img src="./Images/Regelstrecke.png"
     alt="ZHAW Logo"
     width="800"
     style="float:center; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>
---

## 2) Compensation of Gyro Path Filters

Since the measured signal $y(t)$ includes the effect of all filters in the gyro signal path (e.g., DLPF, notch, or software low-pass), the resulting $P(j\omega)$ represents the **filtered plant**:

$P(j\omega) = P_\text{real}(j\omega) \cdot G_f(j\omega)$ where $G_f(j\omega)$ denotes the gyro filter dynamics.

To obtain the **true plant**, the known filter model $G_{f}(j\omega)$ is divided out:
$P_\text{real}(j\omega) = \frac{P(j\omega)}{G_{f}(j\omega)}.$

This step restores the unfiltered system dynamics of the motor and frame, while maintaining consistency with the measured closed-loop responses.

---

## 3) Coherence and Measurement Quality

The coherence is derived from the **auto** and **cross power spectra** of the input $U(f)$ and output $Y(f)$.  
They are defined as

$$S_{UU}(j\omega) = U(j\omega) \cdot \overline{U(j\omega)}$$
$$S_{YU}(j\omega) = Y(j\omega) \cdot \overline{U(j\omega)}$$
$$S_{YY}(j\omega) = Y(j\omega) \cdot \overline{Y(j\omega)}$$

Using these spectra, the **magnitude-squared coherence** is obtained as $\gamma^2(f) = \frac{|S_{YU}(f)|^2}{S_{YY}(f)\,S_{UU}(f)}$.

A high coherence $\gamma^2(f) \approx 1$ indicates a strong linear relationship between the input and output signals,  
confirming that the frequency response estimation is reliable [1, pp. 128–131].  
If the coherence is close to zero, it indicates a very weak or non-existent relationship between the signals.  
Therefore, it serves as an indicator of an **inaccurate or noise-dominated frequency response estimate** at that frequency.

## References

[1] R. Pintelon, J. Schoukens, *System Identification: A Frequency Domain Approach*, 2nd ed., IEEE Press, 2012, **pp. 45–60**.  
[2] P. M. Djuric, S. M. Kay, *Statistical Digital Signal Processing and Modeling*, Prentice Hall, 1993, **pp. 25–35**.
