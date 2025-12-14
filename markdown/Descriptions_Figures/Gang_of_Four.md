# Gang of Four

The Gang of Four is a standard set of four transfer functions that together fully characterize a feedback control system. In this case, your drone. The plots will help you to evaluate stability, tracking performance, disturbance rejection and controller workload. Each plot highlights a different aspect of the system, so it is strongly recommended to work with all of them.

<p align="center">
  <img src="./Images/Gang_of_Four.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>

The key frequency regions for your drone are:

- **~ 0 – 100 Hz:** Actual flight movements and response to pilot commands occur here.
- **~ 100 – 600 Hz:** Dominated by vibrations from motors and propellers.
- **~ >600 Hz:** Primarily high-frequency sensor and electronic noise.

## Bode Diagram (T) – Tracking / Closed Loop Response

This plot shows how well the controller follows the setpoint at different frequencies. Good tracking at lower frequencies (0 – 100 Hz, particularly 0 – 20 Hz) means the drone holds its angle accurately and responds to control inputs. A roll-off at higher frequencies (above 100 Hz) is necessary so the controller does not react to vibrations or noise.

**How it helps tuning:**

- If the closed-loop gain is too low in the movement range (0 – 100 Hz) → increase P or I for sharper response.
- If tracking causes overshoot or oscillation in the movement range → reduce P or increase D.
- If the gain curve falls off already below 20 Hz → controller is too weak → raise P.

**Goal**

- Gain flat and close to **0 dB** from 0 to ~20 Hz (precise tracking and stable control)
- Gradual roll-off from ~20 Hz to 100 Hz (responsive and agile, but not over-driven)
- Strong attenuation above **100 Hz** (gain **< –20 dB**) to block motor/prop vibrations and high-frequency noise

## Sensitivity S – Disturbance Sensitivity & Stability Margin

Sensitivity indicates how strongly the system reacts to external disturbances, such as wind. Low sensitivity in the **0 – 100 Hz** region means the controller effectively suppresses disturbances. At frequencies above **100 Hz**, it's normal and expected for sensitivity to approach 0 dB, since the controller cannot react to such fast disturbances. A small overshoot (peak) in the **100 – 600 Hz** range is common and usually acceptable if it stays below 6 dB.

**How it helps tuning:**

- A strong peak in the movement or vibration range (0 – 600 Hz) → system near instability → reduce P or increase D.
- A smooth, low sensitivity curve in the 0 – 100 Hz range indicates stable and robust tuning.
- Overshoot in the 100 – 600 Hz range should be minimized, but minor peaks below 6 dB are acceptable.

**Goal**

- Sensitivity low (**<0 dB**) from 0 – 100 Hz (strong disturbance rejection during flight)
- Peak in 100 – 600 Hz region stays below **6 dB** (prevents instability or motor/prop oscillation)
- Sensitivity returns smoothly toward **0 dB** above 600 Hz (controller disengages from high-frequency noise)

## Controller Effort SC

This plot shows how much control activity (from P, I, and D terms) is generated in each frequency band. Large values, especially in the **100 – 600 Hz** or **>600 Hz** ranges, can indicate the controller is overreacting to vibrations or noise.

**How it helps tuning:**

- High controller effort above 100 Hz → too much noise passed to motors, try increasing D-term filtering.
- Excessive effort in 100 – 600 Hz → D term too high, leading to hot motors.
- Balanced controller effort curve: active in the 0 – 100 Hz movement range, drops off above 100 Hz.

**Goal**

- Controller effort is moderate (**±0–3 dB**) from 0 – 100 Hz (active control where it matters)
- Any rise in 100 – 600 Hz stays below **~6 dB** (prevents excessive workload and motor heat)
- Above 600 Hz, controller effort drops to **<–20 dB** (filters out all high-frequency noise)

## Compliance SP – Disturbance Transmission / Flexibility

Compliance indicates how much disturbance is transmitted to the system output. Low values in the **0 – 100 Hz** region mean effective rejection of flight disturbances. Peaks in the **100 – 600 Hz** region point to resonances or weaknesses in vibration handling.

**How it helps tuning:**

- A peak in compliance in the 100 – 600 Hz range → disturbances pass through, indicating possible resonance or not enough D term.
- Broad elevation from 0 – 100 Hz → controller is too “soft” → try increasing P or D.
- Low compliance above 600 Hz verifies good high-frequency noise filtering.

**Goal**

- Compliance kept low (**< –20 dB**) from 0 – 100 Hz for strong flight disturbance rejection
- Any peak in 100 – 600 Hz remains below **~6 dB** (avoids amplifying vibrational resonance)
- Compliance continues to decrease and stays **≤ –20 dB** above 600 Hz (suppresses sensor and electronic noise)
Ask anything

## Summary

In genaral, the sensitivity, controller effort, and compliance should follow a smooth, rounded shape. This indicates a stable and robust system that effectively rejects disturbances at low frequencies (0–100 Hz), filters out high-frequency noise (above 100 Hz), and ensures the controller works efficiently without overreacting or generating sharp peaks that might cause instability, resonance, or excessive workload.