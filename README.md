# Betaflight Controller Tuning

## Overview
This repository provides an **offline controller tuning framework** for Betaflight.
It enables analysis of **Blackbox (BBL) flight logs** with automatic chirp excitation,
to tune **linear filters** and **PID controllers** offline.

The workflow:

1. Enable and configure the **chirp signal generator** in Betaflight ([see implementation below](#implementation-of-the-chirp-signal-generator-in-betaflight)).
2. Perform flights with chirp excitation and log data (`.bbl` files).
3. Convert `.bbl` files to `.csv` using the [Betaflight Blackbox Explorer](https://blackbox.betaflight.com/)
4. Use the Python script [`main.py`](./main.py) together with the functions in [`py_lib/`](./py_lib) to:
   - Extract log data
   - Compute frequency responses and spectra
   - Display spectra and spectrograms
   - Estimate closed-loop behavior
   - Extract plant (quad) dynamics seen from the pids
   - Extract measured controller dynamics and compare to the theoretical (split in PI and D  2-DOF controller)
   - Enables offline PID and filter tuning
   - Display step responses of setpoint tracking and input disturbance rejection (compare actual tune to a new user specified tune)
   - Display Bode plots of closed-loop system, e.g. Tracking, Sensitivity, Controller Effort, Compliance (compare actual tune to a new user specified tune)

This repository is intended for use with [Betaflight PR #13105](https://github.com/betaflight/betaflight/pull/13105), which adds the chirp generator mode.

The configurator and the blackbox explorer can be found here:

- [Betaflight Configurator](https://app.betaflight.com)
- [Betaflight Blackbox Explorer](https://blackbox.betaflight.com/)

## Requirements
- **Python**
- **activated [conda environment](./envs/tuning.yml)**

## Quickstart

The repo ships the Blackbox binary log at [`logs/20260417/6_inch_drone.TXT`](./logs/20260417/6_inch_drone.TXT). Open it in the [Betaflight Blackbox Explorer](https://blackbox.betaflight.com/) and export to CSV. **Save the exported `6_inch_drone.TXT.csv` in the same folder as the `.TXT` file** (i.e. `logs/20260417/`) — `main.py` looks for it there. Then run:

```bash
conda env create -f envs/tuning.yml   # first time
conda activate tuning
python main.py
```

## Implementation of the CHIRP Signal Generator in betaflight
The chirp signal generator in Betaflight provides an automated excitation input for analyzing the quadcopter’s dynamic response. It produces a sweep signal whose frequency rises exponentially from a defined starting point to a chosen final value over a preset duration.

This generator outputs a signal that is added directly to currentPidSetpoint inside pid.c. Because of this, the pilot is still in full control over the aircraft during the measurement, and the test can be performed in both `ACRO` and `ANGLE` flight modes.

Because a typical rate-controlled closed-loop system exhibits differentiating behavior from `pidSetpoint` to `pidSum` at low frequencies (up to around 30 Hz), the chirp signal is shaped by a Lag Filter before being injected into the loop.

In order to bring the `CHIRP` signal generator onto your drone, you must include it as a custom define while flashing the Betaflight 2025.12.xx Firmware. You can just type CHIRP into the Field under User Definitions (see picture below). After flashing the Firmware, the `CHIRP` mode appears in the modes tab.

It is recommended to assign the `CHIRP` mode to a non-momentary switch. During the first activation, the chirp is applied to the roll axis. After toggling the switch off and on, the signal is applied to the pitch axis and then to the yaw axis. After one full cycle it resets to the roll axis again. This allows repeated analysis of all axes as needed. The mode can be turned off and reactivated at any time.

When `CHIRP` mode is active, the goggles show CHIR as flight-mode label. Once the chirp sequence is complete (after 20 sec), a blinking message CHIRP IS FINISHED is displayed and the flight-mode label changes back to current flight mode.

If you want to test the mode before flight, it is also possible to do this safely on the bench. To prevent damages to your components, lower the PID gains (e.g., P = 10, I = 0, D = 0), and set the values for chirp_amplitude_roll, chirp_amplitude_pitch, and chirp_amplitude_yaw to 10 while the quadcopter is in `ACRO` mode.

## CLI Paramters

| Name                            | Default Value        | Explanation                                               |
| ------------------------------- | -------------------- | --------------------------------------------------------- |
| `chirp_lag_freq_hz`             | 3 Hz                 | leadlag1Filter cutoff/pole to shape the excitation signal |
| `chirp_lead_freq_hz`            | 30 Hz                | leadlag1Filter cutoff/zero                                |
| `chirp_amplitude_roll`          | 230 deg/sec          | amplitude roll in degree/second                           |
| `chirp_amplitude_pitch`         | 230 deg/sec          | amplitude pitch in degree/second                          |
| `chirp_amplitude_yaw`           | 180 deg/sec          | amplitude yaw in degree/second                            |
| `chirp_frequency_start_deci_hz` | 0.2 Hz / 2 deciHz    | start frequency in units of 0.1 hz                        |
| `chirp_frequency_end_deci_hz`   | 600 Hz / 6000 deciHz | end frequency in units of 0.1 hz                          |
| `chirp_time_seconds`            | 20 sec               | excitation time                                           |

## Assumptions for Offline Tuning
- Dynamic Notch filters are tuned (time-variant)
- RPM filters are tuned (time-variant)
- Thrust Linear is tuned (nonlinear)
- Iterm Relax is tuned (nonlinear)
- Feedforward (FF) is disabled (nonlinear)
- Dynamic Damping is disabled (Dmax = D or 0) (nonlinear)
- Debug mode is set to CHIRP (set debug_mode = CHIRP)
- Blackbox high-resolution logging is enabled (set blackbox_high_resolution = ON)
- Let chirp run for the full chirp_time_seconds

## Tuning Recommendations
- Ensure motors/outputs do not saturate during chirp.
- Adjust amplitude and lag filter if saturation occurs.

## Recommended procedure within one log file per flight
1. Perform two to three throttle sweeps.
2. Complete at least one full sequence of chirp signal excitation, covering roll, pitch, and yaw axes. It is preferable to
cycle through all axes twice. Whether you choose ACRO or ANGLE mode does not matter. Fly in an open space and try to maintain altitude. Be prepared to adjust the throttle as the chirp generator runs. Aim for a smooth and steady flight during the chirp excitation. Ideally, the quadcopter should maintain a steady position and orientation (appart from the axes thats excited).
3. Conduct some maneuvers to test propwash handling, including 180-degree and 360-degree flips. Enjoy yourself and have fun!

## Data for evaluation
To evaluate your flight data, the log file from the Blackbox is required. This has to be converted to a .csv file. You can do this with the [Betaflight Blackbox Explorer](https://blackbox.betaflight.com/).

# Tuning

1. Open the main.py file. This is the only file you need to make changes to tune
2. With `do_insert_legends` you can decide if you would like to have a legend on your plots
3. Define a path to the csv file. For example your file is saved under `bf_controller_tuning\logs\20260417` you have to define `log_folder = 'logs';`, `flight_folder = 20260417` and `log_name = '6_inch_drone.TXT.csv'` your filename.
4. With `ind_ax` you can choose which axis (Roll = 0, Pitch = 1, Yaw = 2) you want to tune
5. `do_compensate_iterm` defines if you want to tune your drone with I-Term Relax (recommended) 
6. For the first flight it is recommended to set default_parameters on `true`. In this case the plot uses the same parameters for the new tune as the old.
7. After this you can define the filter types and you can change their frequencies
8. Below that, you can enter the new PID parameters for Roll, Pitch and Yaw. You can either multiply them or just enter the new value. The values matching the values in Betaflight
9. Now press on **run**, that the first calculation can start. After you pressed run, it should open several figures. If you want to know, how to tune with them go to [Descriptions Figures](./markdown/Descriptions_Figures/).
10. Now enter your new parameters in the main.py file, hit run and repeat as many times as desired.

# How the Calculations Work

You might wonder how we get our results an plots. We can totally understand your wish. Therefore, we created a [detailed documentations](./markdown/). With this documentation and little prior knowledge you might be able to understand our calculations and how the code works.

## Example Flight
- [YouTube Example](https://www.youtube.com/watch?v=bU63eY66QX0)

## Related Theory
- [Chirp](https://github.com/InsaneBroccoli/bf_controller_tuning/tree/PA_final/markdown/Sinarg)
- [Data](https://github.com/InsaneBroccoli/bf_controller_tuning/tree/PA_final/markdown/Data)
- [Frequency Response](https://github.com/InsaneBroccoli/bf_controller_tuning/tree/PA_final/markdown/Frequency_Response)
- [Spectrum and Spectrogram](https://github.com/InsaneBroccoli/bf_controller_tuning/tree/PA_final/markdown/Spectrum_Spectogram)

