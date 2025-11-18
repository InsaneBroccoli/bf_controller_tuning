# Loading and Header Parsing for Betaflight Blackbox Logs

## Why Header Information Matters

A Betaflight Blackbox log consists of two logically separate sections:  
the **header**, which contains configuration and system information, and the **data block**, which contains the actual recorded sensor values.  
The measurement data alone are not sufficient for correct analysis — the header provides essential context that determines how the data must be interpreted. For this reason, both sections are read and processed systematically.

## Efficient Loading of Large Log Files in MATLAB

The first part of the code focuses on loading the log file as efficiently as possible.  
When the file is processed for the first time, the CSV data are fully read and then saved as a `.mat` file. For all subsequent analyses, this `.mat` version is loaded directly.  
This avoids the time-consuming CSV parsing step and significantly speeds up the workflow, especially when working with large logs.

## Extracting Parameter Data from the Blackbox Header

The function `extract_header_information()` reads the textual header section of the log file.  
This part of the log contains all relevant Betaflight parameters that were active during the flight — such as loop frequencies, filter configurations, logging settings, and other internal system values.

The extraction begins once the line containing `frameIntervalI` is found.  
From this point onward, Betaflight lists its parameter lines sequentially.  
All lines up to the line containing `loopIteration` are interpreted as parameters, parsed, and stored in the structure `para`. This creates a clean and structured collection of all metadata needed for downstream analysis.

## Reading the Names of the Measurement Signals

The line containing the keyword `loopIteration` marks both the end of the header and the beginning of the data section.  
This line contains the complete list of all recorded measurement signals in the exact order they appear in the CSV data.

The function processes this line and builds a structure `ind`, which assigns each signal name to its corresponding column index.  
This makes later access to the data clearer and more robust, because signals can be referenced by name rather than by numeric column position.

## Reading the Measurement Data

Once both the parameters and the signal names have been identified, the actual measurement data can be loaded.  
Because the exact number of header lines (`Nheader`) has been determined, MATLAB can reliably read all measurement values by skipping the header and interpreting the following lines as a numerical matrix.

This matrix, stored in the variable `data`, contains one measurement per row and one signal per column.  
The previously created index structure `ind` determines which column corresponds to which signal.  
This allows for clear and direct access, for example:

- `data(:, ind.gyroADC(1))` – Gyro X  
- `data(:, ind.setpoint(3))` – Yaw setpoint  
- `data(:, ind.motor)` – Motor outputs  

In this way, all recorded sensor data are available in a clean, structured form and can be used reliably in subsequent analysis steps.
