<p align="right">
  <img src="./Images/zhaw_logo.jpg"
     alt="ZHAW Logo"
     width="100"
     style="float:right; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>

# Determination of the Transfer Functions
For the calculation of the transfer functions, corresponding input and output signals were selected from the measurement data.  
To reduce the influence of noise and disturbances, the functions [`apply_rotfiltfilt`](./ApplyRotRotfilt.md) and [`estimate_frequency_response`](./estimate_frequency_response.md) were used.  
The first function shifts the signal to the baseband and applies zero-phase filtering to suppress unwanted components, while the second performs the frequency response estimation using the Welch method with amplitude-correct scaling [1][2].  

Additionally, care was taken to use as few signals as possible from the closed loop, since noise in those signals is fed back and can therefore be amplified.  
With these measures, a robust and noise-resistant estimation of the transfer functions could be achieved.

---

## 1) Measured Transfer Paths

In principle, the transfer function \( P(s) \) can be determined by measuring the input and output signals and calculating the spectrum using the **Fast Fourier Transform (FFT)**:
\[
G_{yv}(s) = P(s) = \frac{Y(s)}{U(s)}
\]

However, the resulting measurement is affected by disturbances that are present throughout the closed loop. For this reason, measurements are usually taken **outside the closed loop** to minimize the influence of noise on the signals.

The individual transfer functions are related through the closed-loop relationships:
\[
G_{yu}(s)= \frac{P(s)}{1 + P(s)C_{D}(s)}, \qquad
G_{yw}(s) = T(s) = \frac{G_{yu}(s)\,C_{PI}(s)}{1 + G_{yu}(s)\,C_{PI}(s)}, \qquad
G_{uw}(s) = \frac{C_{PI}(s)}{1 + G_{yu}(s)\,C_{PI}(s)}
\]
\[
\Rightarrow P(s) = \frac{T(s)}{G_{wu}(s)}
\]

These relationships make it possible to **reconstruct the plant transfer function \( P(s) \) from closed-loop measurement data** without opening the control loop during the experiment.



---

## 2) Compensation of Gyro Path Filters

Since the measured signal \( y(t) \) includes the effect of all filters in the gyro signal path (e.g., DLPF, notch, or software low-pass), the resulting \( P(s) \) represents the **filtered plant**:
\[
P(s) = P_\text{real}(s) \cdot G_f(s),
\]
where \( G_f(s) \) denotes the gyro filter dynamics.

To obtain the **true plant**, the known filter model \( G_{f}(s) \) is divided out:
\[
P_\text{real}(s) = \frac{P(s)}{G_{f,\text{ana}}(s)}.
\]

This step restores the unfiltered system dynamics of the motor and frame, while maintaining consistency with the measured closed-loop responses.

---

## References

[1] R. Pintelon, J. Schoukens, *System Identification: A Frequency Domain Approach*, 2nd ed., IEEE Press, 2012, **pp. 45–60**.  
[2] P. M. Djuric, S. M. Kay, *Statistical Digital Signal Processing and Modeling*, Prentice Hall, 1993, **pp. 25–35**.
