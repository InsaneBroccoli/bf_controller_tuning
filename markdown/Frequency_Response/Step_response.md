<p align="right">
  <img src="./Images/zhaw_logo.jpg"
     alt="ZHAW Logo"
     width="100"
     style="float:right; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>

# Step Response from Frequency Response Data

The function `calculate_step_response_from_frd` approximates the **time-domain step response** of the newly calculated transfer function in [Closed-Loop Analysis](./CalculateClosedLoop.md).  
The **step response** illustrates how the system reacts to a sudden change during a unit step and therefore provides a direct view of the **transient behavior** of the control loop. From this response, key performance indicators such as **rise time**, **overshoot**, **settling time**, and **steady-state error** can be observed.  
A smooth and fast rise indicates a well-tuned controller, whereas overshoot or oscillations reveal underdamped dynamics or insufficient phase margin [1, pp. 118–121][2, pp. 455–457].

---

## Method

Through the calculations in the previous steps, the complete **closed-loop transfer function** of the system has already been obtained.  
Based on this transfer function, it is now possible to determine the **step response**, which describes how the system reacts to a sudden change in the input signal.  

To compute this response, the transfer function is transformed from the **frequency domain** back into the **time domain** using the **inverse Fourier transform**.  
This transformation yields the system’s **impulse response** \( g(t) \), which characterizes how the system reacts to an instantaneous impulse.  
Since the step input \( x(t) \) represents a unit step (Heaviside function), two of the three signals in the convolution relationship are known.  
By convolving the known step input with the impulse response, the output \( y(t) \) — i.e., the step response — can be obtained:  

\[
y(t) = \int x(t)\,g(t)\,dt
\]

In this form, the step response reflects the combined effects of all filters, delays, and feedback mechanisms within the control loop.  
A well-damped and smooth curve indicates a stable and responsive system,  
whereas oscillations or overshoot reveal underdamped or poorly tuned controller dynamics [1, pp. 118–121][2, pp. 455–457].


---

**References**

[1] Pintelon, R., & Schoukens, J. (2012). *System Identification: A Frequency Domain Approach* (2nd ed.). Wiley-IEEE Press, pp. 118–121.  
[2] Hayes, M. H. (1996). *Statistical Digital Signal Processing and Modeling*. Wiley, pp. 455–457.
