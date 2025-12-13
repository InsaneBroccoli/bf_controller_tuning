# Flight Overview - Gyro Data

This figure provides an overview of the most important signals recorded during the flight.

### Top subplot (Gyro)

The angular rates for roll, pitch, and yaw are shown. You can clearly see the chirp patterns, where each axis is driven by a signal that slowly changes its frequency over time. These chirps help reveal how the system reacts to different frequencies.

### Second subplot (AxisSum)

The AxisSum plot shows the sum of the absolute gyro rates for roll, pitch, and yaw. It provides a single value representing the total rotational activity of the drone. This plot is useful for detecting phases with high overall motion, vibration, or control effort, for example during aggressive inputs or unstable flight conditions. Larger values indicate higher total gyro activity. AxisSum does not provide information about axis coupling or interactions between axes, such effects must be evaluated by comparing the individual roll, pitch, and yaw gyro signals. Overall, AxisSum is an indicator of total rotational intensity, not a measure of axis independence or control quality.

### Third subplot (Motors):

This plot shows the motor speeds throughout the entire flight. Each motor is plotted separately to allow for comparison between them. The y-axis represents motor speed in revolutions per minute (RPM). Differences between motor traces can indicate control corrections, load changes, or possible imbalances during flight.

### Bottom subplot (Throttle)

This plot shows the throttle stick input throughout the flight. It represents the amount of throttle commanded by the pilot. The throttle signal ranges from 0 to 1000, where 0 corresponds to minimum throttle and 1000 corresponds to 100% throttle.

<p align="center">
  <img src="./Images/Flight_Overview.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>