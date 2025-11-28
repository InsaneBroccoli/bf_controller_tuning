# Gyro Spectogram

The spectrograms show how the gyro noise is applyed over frequency and throttle for all three axes, once without and once with filtering. This makes the plots very useful to understand how the drone behaves over the whole throttle range and how effective the filters are.
Spectrograms are useful because they show how noise develops across throttle and frequency. In the unfiltered data, the characteristic noise patterns become visible, including the areas where motor harmonics and other disturbances dominate. These regions are important to identify because they represent the frequencies that can interfere with the PID controller and therefore need to be handled by the filters.
After filtering, spectrograms make it easy to see how effectively the filter configuration suppresses these noise regions. They help determine whether the filters are strong enough to reduce high-frequency disturbances while still preserving the lower-frequency components that contain the meaningful control information. This makes spectrograms an important tool for checking filter settings and ensuring a good balance between noise reduction and control responsiveness.

## What You Should Look For

You should pay attention to:

- **Where the noise band appear and how they move with throttle.** This will help you to identifiy the frequency ranges that filters must convert.
- **At which frequencies the filtering becomes effective.** You can see whether the filters reduce the noice as expected.
- **How strog the noise is across the throttle range.** This indicates if addtional filtering is needed.

## How the Spectograms Guide Filter Seeting

These plot will help you to decide:

- **IWhether an RPM filter is needed:** If you see motor-related harmonics — smooth, upward-moving, parallel lines that shift with throttle, you should enable an RPM filter.
- **How wide the dynamic notch filter should be** since the spectrogram shows how far the noise band moves under load.
- **Whether a fixed gyro notch is useful or not** if a strong noise band stays at the same frequency regardless of throttle.

## Why This Plot is Relevant

Spectograms give the user a full overview of how the drone behaves at different frequency and different thrust. This helps to ensure that filters are set when need, but not to over filtering the system. This lead to a smoother flight, cooler motors and more responsive control.

<p align="center">
  <img src="./Images/Gyro_Spectograms.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>