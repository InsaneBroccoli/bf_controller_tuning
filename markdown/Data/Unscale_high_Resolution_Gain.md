# Unscaling High Resolution Gains

Betaflight provides a **High Resolution Mode** for the Blackbox logger. In this mode, various signals are scaled up by a factor of 10. This allows the values to be stored as integers, which saves memory and simultaneously provides higher numerical precision.  
The affected signals are:

- **gyroADC** (filtered gyro data)
- **gyroUnfilt** (unfiltered gyro data)
- **rcCommand** (control inputs)
- **setpoint** (controller target values)

## Rescaling for Analysis

To correctly analyze these data — for example, for FFT analysis or PID tuning — the affected channels must be divided by the same scaling factor. Only then do the values represent real gyro and RC measurements in physical units.  
Without this rescaling, all affected values would be 10 times too large, which would lead to severe errors in the analysis:

- incorrect gyro amplitude → wrong FFT peaks
- incorrect setpoints → distorted tracking analysis
- incorrect RC values → incorrect stick response models
- faulty PID calculations
- overall: any frequency response analysis would be unusable
