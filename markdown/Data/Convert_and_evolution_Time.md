# Time Conversion and Sampling Intervals in Betaflight Blackbox Logs

In Betaflight Blackbox logs, time is not stored in seconds but as a continuously increasing counter measured in microseconds. For any correct analysis of the log, these raw timestamps must first be converted into a usable time axis. At the same time, the actual sampling intervals between log entries must be determined, as they form the basis for all subsequent signal processing.

The log steps in a Betaflight Blackbox recording are not always the same size because the logger does not operate with a guaranteed constant sampling interval. This is due to several technical factors related to how Betaflight collects sensor data and writes to the SD card. Betaflight is not a hard real-time system — tasks such as gyro sampling, filtering, PID calculations, motor updates, telemetry, and OSD run with different priorities. The Blackbox logger has relatively low priority and may therefore be executed with slight delays. Additionally, SD card write operations introduce timing variations because their internal processes require different amounts of time.

By calculating the time differences between log entries, one can determine how much time passed between consecutive samples. These differences represent the actual sampling intervals. Converting these intervals back into microseconds provides a high-resolution view of the timing structure of the log.

Analyzing these sampling intervals enables several important evaluations:

- **Determine the effective logging rate**  
  This makes it possible to verify whether Betaflight is truly logging at the intended rate — for example, 1 kHz — or if the actual rate is slightly higher or lower.

- **Detect jitter**  
  Variations in sampling time can result from SD card latency, CPU load, or Betaflight's scheduling behavior. These variations become visible when examining the time differences.

- **Identify errors such as dropped frames**  
  Large jumps in the recorded time for example, 10 to 20 milliseconds indicate delays or missing data due to SD card write stalls or system load.

- **Enable numerically correct analysis**  
  Accurate sampling intervals are essential for correct FFT calculations, filter design, frequency analysis, and controller evaluation. These methods rely on the true time base provided by the log.

Together, these steps convert the raw Betaflight timestamp data into a physically correct time structure, enabling reliable and precise analysis of the recorded signals.
