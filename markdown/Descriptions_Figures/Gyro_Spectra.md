# Gyro Spectra

## Filtersettings

Filters in flight controllers are essential for reducing sensor noise and providing the PID controller with a clean signal. Different filters serve different purposes: some remove unwanted frequencies, while others specifically target vibration-related noise during flight. Using the right combination of filters is crucial for achieving stable flight and efficient motor behavior.
One important filter is the low-pass filter, which should act as an anti-aliasing filter. An anti-aliasing filter ensures that very high frequencies do not fold back into the measurable signal during logging. However, it does not remove motor noise seen by the PID controller during flight, because this noise lies well below the anti-aliasing cutoff. Therefore, additional filters are necessary to clean up the gyro signal for flight control.
Recommended filters:

- RPM Filter: Should always be enabled. It uses exact motor RPM values to remove motor harmonics and their multiples very effectively, making the filtering process highly precise.
- Dynamic Notch Filter: Strongly recommended. It automatically adjusts to vibration peaks during flight, which is useful for handling changing vibration patterns.
- Fixed Gyro Notch Filter: Optional but essential for a perfect tune. It is only useful when all three gyro axes show a strong peak at the same frequency, usually indicating a frame or prop resonance.
- D-Term Low-Pass Filter: Should always be active. Since the D-term amplifies high-frequency noise, a cutoff around 80–120 Hz helps keep motors cool and prevents oscillations.

## Information of Spectra

The gyro spectra are very important because they show how much noise and vibration the drone produces at different frequencies. This helps you understand which parts of the signal are useful for control and which parts need to be filtered out. By looking at the spectra, you can clearly see the motor harmonics, frame vibrations, and other noise peaks that could disturb the PID controller.
When analysing the spectra, you should pay attention to:

<p align="center">
  <img src="./Images/Gyro_Spectra.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>

- **Resonance peaks:** If all three gyro axes show a peak at the same frequency, it often means a frame or propeller resonance. In this case, a fixed notch filter can be very helpful.
- **Noise floor:** A high noise floor means a lot of general vibration or bad filtering. A clean noise floor indicates that the filters are working well.
- **Axis differences:** If one axis is much noisier than the others, it may indicate mechanical issues like a bent motor shaft, unbalanced prop, or loose screws.

By combining the spectra with the filter settings, you can choose filters that keep the noise low while still allowing fast and responsive control behaviour.

## Axis Sum Spectra

The Axis Sum spectra display, for each frequency, the **summed PID controller outputs for each axis (Roll, Pitch, Yaw)**. In Blackbox log analysis, the Axis Sum is **not** a direct recording from the flight controller, but a computed field created by adding all relevant control terms associated with each axis. Specifically, for each axis, the following terms are summed:

- Proportional (P)
- Integral (I)
- Derivative (D)
- Feedforward (F)
- Setpoint weighting or smoothing (S, if present)

The Axis Sum represents the total control effort the flight controller applies to each axis (Roll, Pitch, Yaw) at any moment. Looking at the spectra, these traces show how much the controller needs to correct for vibrations, oscillations, or aggressive maneuvers. Low values mean the drone is running smoothly with little intervention, while high values at certain frequencies highlight issues where the controller is fighting disturbances or instability.

The formula used by Blackbox log viewers is:
```
axisSum[axis] = axisP[axis] + axisI[axis] + axisD[axis] + axisF[axis] + axisS[axis]
```
For a well-tuned drone, Axis Sum amplitudes should be low across all frequencies except near the motor bands, where more control output is naturally needed. Spikes or an elevated noise floor in the Axis Sum indicate that the flight controller is working unusually hard—this can point to tuning problems, mechanical issues, or poor filtering.

## Source 
The calculation and field logic for Axis Sum are implemented in the Betaflight Blackbox Log Viewer project, specifically in the [`injectComputedFields()`](https://github.com/betaflight/blackbox-log-viewer/blob/7de89d7f92904b3aab8d20a47f432c30e2871ea9/src/flightlog.js#L859-L889) function:
> axisSum[axis] = axisP[axis] + axisI[axis] + axisD[axis] + axisF[axis] + axisS[axis]  
> (with a PID sum limit if defined; computed per axis for every log entry)