# Step Response

## Tracking T - Interpretation & Tuning Relevance


The tracking response shows how quickly and accurately the drone follows the step command. A good tracking curve rises steeply, reaches the target without excessive overshoot, and then holds stable.

**What matters for tuning**

- **Fast rise time:**  Too slow $\rightarrow$ P or I too low
- **Minimal overshoot:** Overshoot too high $\rightarrow$ Reduce P or increase D

If the measured and the actual calculated values do not match, check if `do_compensate_iterm` is activated in main.

## Compliance SP – Interpretation & Tuning Relevance


The compliance response describes how the system reacts to a disturbance step rather than a command step. It shows how effectively the controller rejects external disturbances.

**What matters for tuning**

- **Low peak amplitudes:** If peak is over 0.75 $\rightarrow$ P too low or D not high enough
- **Slow decrease:** If decrease is slower than 0.3s $\rightarrow$ increase D-Term
- **Unsteady course:** If oscillations appear $\rightarrow$ D too low or P too high

<p align="center">
  <img src="./Images/Step_response.jpg"
     alt="Original noisy signals"
     width="1000"
     style="float:center; margin-left:10px; margin-right:10px;">
</p>