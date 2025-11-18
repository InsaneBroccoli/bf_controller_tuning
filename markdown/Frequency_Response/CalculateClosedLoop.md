# Closed-Loop Analysis

The function `calculate_closed_loop` determines the main transfer relationships of the closed-loop control system, including how the controller, plant, and filters — as described in [Filter and Controller](./Filter_and_Controller.md) and [Frequency Response Estimation](./estimate_frequency_response.md) interact to shape the overall system behavior.  
This structure separates the control actions into an **inner loop** (e.g., rate or current control) and an **outer loop** (e.g., attitude or position control), allowing for a flexible distribution of stability and tracking performance — a concept also applied in the Betaflight control architecture.

## Sensitivity

The **sensitivity function** describes how the closed-loop system reacts to disturbances and noise.  
Low sensitivity at low frequencies ensures good disturbance rejection and accurate reference tracking,  
while higher sensitivity at high frequencies helps to avoid noise amplification and maintain robustness [1, pp. 35–39][2, pp. 192–196].

In the Bode plot, \( |S(j\omega)| \) should ideally be **small at low frequencies** (below the controller bandwidth)  
to minimize steady-state errors and suppress disturbances.  
Around the **crossover frequency**, a moderate peak may appear — typically limited to about **6 dB** —  
indicating the trade-off between tracking performance and stability margin.  
At **high frequencies**, \( |S(j\omega)| \) should approach 0 dB, ensuring that the controller does not react to sensor noise or unmodelled dynamics.  

A well-shaped sensitivity curve therefore combines **good disturbance rejection**, **sufficient stability margin**,  
and **robust behavior** against uncertainties in the plant model.


## Loop Calculation

As already mentioned, the control loop is split into an **inner** and an **outer** loop.  
First, we form the inner loop (open-loop gain, sensitivity, and the closed inner path from `u` to `y`).

\[L_{in}(j\omega) = P(j\omega) \cdot G_{f}(j\omega) \qquad S_{in}(j\omega) = \frac{1}{1+L_{in}(j\omega)} \qquad G_{yu}(j\omega) = \frac{L_{in}(j\omega)}{1+L_{in}(j\omega) \cdot G_{fD}(j\omega) \cdot C_{D}(j\omega)}\] 

With this Information we are able to calculate the hole closed loop system.

\[L_{out}(j\omega) = G_{yu}(j\omega) \cdot C_{PI}(j\omega) \qquad S_{out}(j\omega) = \frac{1}{1+L_{out}(j\omega)} \qquad G_{yw}(j\omega) = \frac{L_{out}(j\omega)}{1+L_{out}(j\omega)}\]

<p align="center">
  <img src="./Images/Regelstrecke.png"
     alt="ZHAW Logo"
     width="800"
     style="float:center; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>
---

**References**

[1] Pintelon, R., & Schoukens, J. (2012). *System Identification: A Frequency Domain Approach* (2nd ed.). Wiley-IEEE Press.  
[2] Hayes, M. H. (1996). *Statistical Digital Signal Processing and Modeling*. Wiley.
