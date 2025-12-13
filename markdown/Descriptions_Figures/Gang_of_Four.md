# Gang of four

The Gang of Four is a standard set of four transfer functions that together fully characterize a feedback control system. In this case your drone. The plots will help you to evaluate stability, tracking performance, disturbance rejection and controller workload. Each plot highlights a different aspect of the system, so it is strongly recommended to work with all of them.

<p align="center">
  <img src="./Images/Gang_of_Four.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>

## Bode Diagram (T) Tracking / Closed Loop Response

This plot shows how well the controller follows the setpoint across frequencies. Good tracking at low frequencies means the drone holds its angle accurately. A roll-off at higher frequencies is necessary so the controller does not react to motor noise.

**How it helps tuning:**

- If the closed-loop gain is too low → increase P or I for sharper response.
- If tracking causes overshoot or oscillation → reduce P or increase D.
- If the curve falls off too early → controller is too weak → raise P.

**Goal**

- Stable at low frequencies (**0 dB**)
- Bandwidth between **20-60 Hz**
- Strong noise suppression at high frequencies (**< -20 dB**)


## Sensitivity S – Disturbance Sensitivity & Stability Margin

The Sensitivity indicates how strongly a system reacts to external disturbances like wind. Low sensitivity means the controller suppresses these disturbances well, high sensitivity indicates a more fragile system. At high frequencies, it is totally normal and expected that the Sensitivity approaches 0 dB. In this region, the controller cannot react anymore, so disturbances are neither amplified nor actively suppressed.
A small overshoot in the middle frequency range is also normal. As long as the peak is not higher than 6 dB, the system should still be stable.

**How it helps tuning**.

- A strong peak → system near instability → reduce P or increase D
- A smooth, low sensitivity curve indicates stable and robust tuning
- This plot shows how aggressively the system can be tuned before losing stability

**Goal**

- Low sensitivity at low frequencies (strong disturbance rejection)
- Peak in the mid-frequency range stays below **6 dB**
- Sensitivity returns smoothly towards **0 dB** at high frequencies (normal behavior)
  
## Controller Effort SC

This plot shows how much control activity the PI and D loops generate. Large values indicate a heavy workload or that the controller is responding strongly to noise, especially in the high-frequency range.

**How it helps tuning:**

- High controller effort at high frequencies → too much noise → increase D-term filtering
- Excessive effort at mid frequencies → D too high, resulting in hot motors
- A balanced controller effort curve shows efficient controller behavior

**Goal**

- Controller effort stays within a moderate level **(±0–3 dB)** in the useful frequency range
- Mid-frequency rise remains below **~6 dB** to prevent hot motors and excessive D-term activity
- High-frequency controller effort drops to below **–20 dB** for effective noise suppression

## Compliance SP – Disturbance Transmission / Flexibility

Compliance indicates how much external disturbance is transmitted to the system output. Low values are preferred, as they show effective rejection of vibrations or external forces. High values suggest that disturbances are not being adequately suppressed.

**How it helps tuning**

- A peak in compliance reveals a frequency where disturbances pass through strongly → potential frame resonance or insufficient D-term.
- Broad elevation across frequencies → controller too soft → increase P or D.
- Filter effects can also be seen here; good filtering reduces compliance at high frequencies.

**Goal**

- Low compliance at low frequencies **(< –20 dB)** for strong disturbance rejection.
- Any mid-frequency peaks stay below **~6 dB** (avoids resonance amplification).
- Compliance continues to decrease toward high frequencies **(≤ –20 dB)** to prevent vibrations from passing through the system.
