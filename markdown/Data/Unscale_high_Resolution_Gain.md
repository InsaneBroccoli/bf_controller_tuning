# Unscaling High Resolution Gains in Betaflight Blackbox Logs

## Overview

Betaflight provides a **High Resolution Mode** for the Blackbox logger that improves data precision for critical flight signals.  Understanding this feature is essential for accurate flight analysis, PID tuning, and frequency domain analysis.

## What is High Resolution Mode?

High Resolution Mode is an optional Blackbox setting that multiplies specific signals by a factor of 10 before storing them. This improves precision while keeping the storage format efficient.

### Why Does This Feature Exist?

The flight controller has limited storage bandwidth and CPU resources.  Blackbox logs store data as integers rather than floating-point numbers because:
- Integer operations are faster on the flight controller
- Integer storage uses less memory bandwidth
- Integer writes to SD cards are more efficient

However, storing values like 45.7°/s as integers would lose the decimal precision.  By multiplying by 10 first, the value becomes 457, preserving one decimal place while still using integer storage.

This is particularly important for:
- Low-amplitude signals where rounding errors would dominate
- Precise tracking of small oscillations
- Accurate noise floor measurements
- Detailed PID performance analysis

## Affected Signals

High Resolution Mode scales these four signals by a factor of **10**:

| Signal Name | Description |
|------------|-------------|
| **gyroADC** | Filtered gyro data (after all digital filters) |
| **gyroUnfilt** | Raw gyro data before filtering |
| **rcCommand** | Processed RC stick commands |
| **setpoint** | Final controller target values | 

## Rescaling for Analysis

To correctly analyze the data, for example FFT analysis or PID tuning, the affected channels must be divided by the same scaling factor. Only then do the values represent real gyro and RC measurements in physical units.  
Without this rescaling, all affected values would be 10 times too large, which would lead to severe errors in the analysis:

- incorrect gyro amplitude → wrong FFT peaks
- incorrect setpoints → distorted tracking analysis
- incorrect RC values → incorrect stick response models
- faulty PID calculations
- overall: any frequency response analysis would be unusable

## Enabling High Resolution Mode

### Via CLI
```
set blackbox_high_resolution = ON
save
```

Check current setting with:
```
get blackbox_high_resolution
```

## Technical Details

### Storage Impact
High resolution mode increases log file size because:
- The stored integer values are larger (more digits)
- The compression algorithms have slightly less redundancy to exploit
- The trade-off is generally worthwhile for the precision gained

### Processing Workflow
The rescaling must happen at a specific point in your data pipeline: 

```
1. Load CSV or MAT file
2. Parse header and extract blackbox_high_resolution flag
3. Build signal index mapping
4. → RESCALE if high_resolution = 1 (divide by 10)
5. Convert time units
6. Calculate sampling intervals
7. Apply other unit conversions
8. Proceed with analysis
```

Rescaling must occur **after loading** but **before any mathematical operations** on the data.

## References and Sources

This documentation is based on the following sources: 

### Primary Sources

**Betaflight Firmware Source Code:**
- Blackbox configuration structure definition:   
  [betaflight/blackbox. h](https://github.com/betaflight/betaflight/blob/7f8e4d5b8af17bee89b482d7f4aa1ad40e772496/src/main/blackbox/blackbox.h#L8-L95)  
  Defines the `blackboxConfig_t` structure including the `high_resolution` field.

- Blackbox implementation and scaling:   
  [betaflight/blackbox.c](https://github.com/betaflight/betaflight/blob/7f8e4d5b8af17bee89b482d7f4aa1ad40e772496/src/main/blackbox/blackbox.c#L412-L500)  
  Shows usage of `blackboxHighResolutionScale` variable.

**Betaflight Blackbox Log Viewer:**
- Log parser implementation:  
  [blackbox-log-viewer/flightlog_parser.js](https://github.com/betaflight/blackbox-log-viewer/blob/7de89d7f92904b3aab8d20a47f432c30e2871ea9/src/flightlog_parser.js#L153-L216)  
  Demonstrates parsing of `blackbox_high_resolution` parameter from log headers.

### Related Project Documentation

- [Data Import and Header Parsing](./Dataimport.md) - Loading Betaflight logs and parsing header information
- [Time Conversion and Sampling](./Convert_and_evolution_Time.md) - Converting timestamps and handling sampling intervals

### Additional Resources

- [Betaflight GitHub Repository](https://github.com/betaflight/betaflight) - Official firmware source code
- [Betaflight Blackbox Wiki](https://github.com/betaflight/betaflight/wiki/Blackbox-logging) - Official documentation on blackbox logging features
