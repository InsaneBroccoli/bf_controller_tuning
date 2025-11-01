# Function `ApplyRotRotfilt`

<p align="right">
  <img src="../Images/zhaw_logo.jpg"
       alt="ZHAW Logo"
       width="250"
       style="float:right; margin-left:20px; border-radius:10px;">
</p>

The apply_rotfiltfilt function mathematically shifts the constantly changing frequency signal (the logarithmic chirp) so that it appears as a signal with a constant frequency. In this baseband, interfering unwanted signal components can be effectively removed with a low-pass filter without any time distortion. The filtered signal is then ‘shifted back’ to its original frequency curve. This results in a smoothed, phase-correct signal that contains the essential information of the original chirp but is free of noise and high-frequency interference components.

## 1) Time-Dependent Phase Rotation (`sinarg`)

The function **`apply_rotfiltfilt`** uses a **time-dependent phasor**
\[
p(t) = e^{i\,\text{sinarg}(t)},
\]
where `sinarg` represents the **phase progression** of the carrier signal – for example, the phase of a **chirp** (a signal whose frequency changes continuously over time).  
`sinarg` is given in **radians** and can increase either **linearly** or **nonlinearly**, depending on the frequency behavior of the input signal.

By multiplying the input signal \(x(t)\) with this phasor \(p(t)\) or its complex conjugate \(p^*(t)\),
\[
y_R(t) = x(t)\,p(t), \qquad y_Q(t) = x(t)\,p^*(t),
\]
the signal is mathematically **shifted in frequency**.  
You can imagine this as "rotating" the current carrier frequency of the signal down to **0 Hz (baseband)**.  
Once the signal is centered around 0 Hz, it becomes much easier and more precise to filter.

In mathematical terms, this corresponds to a **frequency shift**, as described by the well-known **Fourier relation**:
\[
e^{i2\pi\xi_0 t}\,f(t) \;\xleftrightarrow{\mathcal{F}}\; \hat{f}(\xi - \xi_0).
\]
The difference here is that \(\xi_0\) is **not constant** — it changes **over time**.  
This means that the instantaneous frequency of the chirp is dynamically shifted to the baseband at every moment.
