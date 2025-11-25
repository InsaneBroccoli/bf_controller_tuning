# Gyro Spectra

## Filtersettings

A Lowpassfilter is used as an anti-alising Filter. An anti-aliasing filter is mainly important for data recording and analysis. It makes sure that very high frequencies do not fold back into the measurable signal during logging. However, it does not remove the motor noise that the PID controller sees during flight, because this noise lies well below the anti-aliasing cutoff.
The RPM filter should always be turned on. It uses the exact motor RPM values to remove motor harmonics and their multiples very effectively. This makes the whole filtering process more precise.
The dynamic notch filter is also recommendet to use. It automatically adjusts itself to shifting vibration peaks during flight, keeping noise low without filtering too much.
A fixed gyro notch filter is optional, but for a perfect tune essential. It is only useful when all three gyro axes show a strong peak at the same frequency, which usually means a frame or prop resonance.
The D-term low-pass filter should also be active. Since the D-term amplifies high-frequency noise, a cutoff around 80–120 Hz helps keep the motors cool and prevents oscillations. After flying:
- If motors are warm or hot, lower the D-term LPF frequency.
- If motors stay cool, you can increase it slightly for better response.

## Information of Spectra
The gyro spectra are very important because they show how much noise the drone produces at different frequencies. This helps you understand which parts of the signal are useful for control and which parts need to be filtered out. By looking at the spectra, you can clearly see the motor harmonics, frame vibrations, and other noise peaks that could disturb the PID controller.
When analysing the spectra, you should pay attention to:

- Motor harmonics: These usually appear between 200–600 Hz and often have strong peaks. They help you decide where the LPF cutoffs should be placed.
- Resonance peaks: If all three gyro axes show a peak at the same frequency, it often means a frame or propeller resonance. In this case, a fixed notch filter can be very helpful.
- Noise floor: A high noise floor means a lot of general vibration or bad filtering. A clean noise floor indicates that the filters are working well.
- Chirp response: During chirp excitation, you can see how the system reacts across many frequencies. This helps you understand how the filters and the PID loop behave in real flight conditions.
- Axis differences: If one axis is much noisier than the others, it may indicate mechanical issues like a bent motor shaft, unbalanced prop, or loose screws.

By combining the spectra with the filter settings, you can choose filters that keep the noise low while still allowing fast and responsive control behaviour.

## Axis Sum Spectra

A well-filtered axis-sum should drop below 1 shortly after the motor frequency range and stay there. This shows that the drone is running smoothly and that vibrations and noise are well under control. You should also look for any strong peaks or a raised noise floor, as these can indicate resonances, unbalanced motors, or general issues with vibrations. A clean axis-sum means the system is stable and the filters are working correctly.

<p align="center">
  <img src="./Images/Gyro_Spectra.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>