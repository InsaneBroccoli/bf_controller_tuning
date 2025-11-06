<p align="right">
  <img src="./Images/zhaw_logo.jpg"
     alt="ZHAW Logo"
     width="100"
     style="float:right; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>

# Filter and Controller Structure

The functions `get_filter` and `calculate_transfer_functions` together form the basis for constructing the controller and filter structure in the closed control loop.  
Both follow the filter and controller architecture used in **Betaflight**.

---

## Function `get_filter`

This function generates a **discrete filter** section of a specified type and returns it as both a state-space model and a transfer function. The filters are implemented according to the Betaflight filter architecture and frequency characteristics.

According to **Betaflight**, the user can choose which of these different filter types to apply depending on the axis and signal path.
- **PT1, PT2, PT3:** Low-pass filters with Butterworth (3dB point is the same for all PT's ) correction.  
  They attenuate high-frequency noise while maintaining a smooth phase response.  
  Higher-order variants (pt2, pt3) provide steeper roll-off at the cutoff frequency.

- **Biquad:** 2nd-order low-pass filter with \( Q = 1/\sqrt{2} \) (Butterworth characteristic).  
  Offers a clean –12 dB/octave slope and minimal overshoot in the time domain.

- **Notch:** Band-stop (notch) filter with quality factor \( Q \) computed by `get_notch_Q`.  
  Used to suppress narrow-band resonances such as frame vibrations or motor noise.

- **Lead / Lag :** Phase compensation or lead–lag network.  
  These filters modify the phase response (lead or lag) to improve control stability or compensate for sensor delay.

All filters are discretized with the sample time \( T_s \).  The function `get_filter` is called internally by `calculate_transfer_functions` to build the gyro, D-term, and P-term filter paths.

---

## Function `calculate_transfer_functions`

This function constructs the **filter chains and controllers** for a selected axis.  
It creates three filter paths:
1. **Gyro path \( G_f \):** Low-pass filters, dynamic low-pass filters, notch filters, optional phase compensation  
2. **D-term path \( G_d \):** Filter chain for the derivative path, cascaded with \( C_D \)  
3. **P-term path \( G_{f_p} \):** Optional path with phase compensation or yaw low-pass

From the axis parameters, the discrete **PI(+P)** and **D** controllers are built as:
\[
C_{PI}(z) = K_p\,G_{f_p} + K_i\,T_s\frac{z}{z-1}, \qquad
C_D(z) = \frac{K_d}{T_s}\frac{1 - z^{-1}}{z^{-1}} \cdot G_d
\]

The function returns the state-space models \( C_{PI} \), \( C_D \), \( G_f \), and the effective PID gain vector.  
