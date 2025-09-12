# BF Controller Tuning

## Overview

This repository provides an **offline controller tuning framework** for Betaflight.
It enables analysis of **Blackbox (BBL) flight logs** with automatic chirp excitation,
to tune **linear filters** and **PID controllers** offline.

The workflow:

1. Enable and configure the **chirp signal generator** in Betaflight.
2. Perform flights with chirp excitation and log data (`.bbl` files).
3. Convert `.bbl` files to `.csv` using the [Betaflight Blackbox Explorer]
4. Use the MATLAB script [`bf_controller_tuning.m`](./bf_controller_tuning.m) together with the functions in [`lib/`](./lib) to:
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

- [Betaflight Configurator](https://master.app.betaflight.com/)
- [Betaflight Blackbox Explorer](https://master.blackbox.betaflight.com/)

## Requirements

- **MATLAB** (last tested with R2024a)
- **Control System Toolbox**
- **Signal Processing Toolbox**

## Chirp Signal Generator in Betaflight

### Concept

An automatic chirp signal generator is used to excite the quadcopter during flight. The chirp frequency increases exponentially over a defined period, starting from a chosen initial frequency and ending at a specified final frequency.

The chirp signal generator is the module responsible for producing this excitation signal.

The chirp signal is added to the currentPidSetpoint in the `pid.c` file. This means the quadcopter can still be controlled while the measurement is running. The system’s behavior can be tested in both `ACRO` and `ANGLE` modes.

Because a typical rate-controlled closed-loop system exhibits differentiating behavior from `pidSetpoint` to `pidSum` at low frequencies (up to around 30 Hz), the chirp signal is shaped by a Lag Filter before being injected into the loop.

The chirp generator is implemented as a feature. To enable it, include chirp in the custom defines or add the `CHIRP` flag when building locally.

It is recommended to assign the `CHIRP` mode to a non-momentary switch. The first time `CHIRP` mode is activated, the chirp is applied to the roll axis. If the switch is toggled off and on again, it is applied to the pitch axis, then to the yaw axis, and then cycles back to roll. This makes it possible to cycle through all three axes multiple times if needed. The `CHIRP` mode can be disabled and re-enabled at any point.

When `CHIRP` mode is enabled, `CHIR` is displayed as the active flight mode in the goggles. Once the chirp signal finishes, a blinking warning `CHIRP IS FINISHED` is shown, and the `CHIR` label disappears.

The chirp mode can also be tested safely on the bench without propellers. In this case, reduce the PID gains (e.g. P = 10, I = 0, D = 0) and set `chirp_amplitude_roll = chirp_amplitude_pitch = chirp_amplitude_yaw = 10` while in `ACRO` mode.

## CLI Parameters
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
- Debug mode is set to `CHIRP` (`set debug_mode = CHIRP`)
- Blackbox high-resolution logging is enabled (`set blackbox_high_resolution = ON`)
- Let chirp run for the full `chirp_time_seconds`

## Tuning Recommendations

- Ensure motors/outputs do not saturate during chirp.
- Adjust amplitude and lag filter if saturation occurs.

## Typical testing procedure within one BB-log file / flight

1. Perform two to three throttle sweeps while in ``ACRO`` mode.
2. Complete at least one full sequence of chirp signal excitation, covering roll, pitch, and yaw axes. It is preferable to cycle through all axes twice. Whether you choose ``ACRO`` or ``ANGLE`` mode does not matter. Fly in an open space and try to maintain altitude. Be prepared to adjust the throttle as the chirp generator runs. Aim for a smooth and steady flight during the chirp excitation. Ideally, the quadcopter should maintain a steady position and orientation (appart from the axes thats excited).
3. Conduct some maneuvers to test propwash handling, including 180-degree and 360-degree flips. Enjoy yourself and have fun!

## Data Required

- A `.bbl` log (with chirp excitation enabled). This can be converted to `.csv` using the [Betaflight Blackbox Explorer](https://master.blackbox.betaflight.com/).

## Example Flight

- [YouTube Example](https://www.youtube.com/watch?v=bU63eY66QX0)

## Related Theory

- [Chirp (Wikipedia)](https://en.wikipedia.org/wiki/Chirp)
- [MATLAB chirp()](https://ch.mathworks.com/help/signal/ref/chirp.html)

## Repository Structure (as of 21.09.2025)

```tree
.
├── bf_controller_tuning.m           # Main analysis script
├── dev/
│   └── bf_controller_tuning.m       # Old evaluation script (reference only)
├── examples/                        # Example scripts
│   ├── bf_biquad_delay.m
│   ├── bf_biquad_step.m
│   ├── bf_chirp_signals.m
│   ├── bf_filter_compare_sampling_rates.m
│   └── bf_notch_Q.m
├── lib/                             # Core analysis functions
│   ├── apply_rotfiltfilt.m
│   ├── calculate_closed_loop.m
│   ├── calculate_controllers.m
│   ├── calculate_step_response_from_frd.m
│   ├── calculate_transfer_functions.m
│   ├── downsample_frd.m
│   ├── estimate_frequency_response.m
│   ├── estimate_spectra.m
│   ├── estimate_spectrogram.m
│   ├── expand_multiple_figure_nr.m
│   ├── extract_header_information.m
│   ├── get_chirp_signals.m
│   ├── get_fcut_from_D_and_fcenter.m
│   ├── get_fcut_from_exp.m
│   ├── get_filter.m
│   ├── get_ind_eval.m
│   ├── get_my_colors.m
│   ├── get_notch_Q.m
│   ├── get_pid_scale.m
│   └── get_switch_case_text_from_para.m
├── LICENSE
└── README.md
```
