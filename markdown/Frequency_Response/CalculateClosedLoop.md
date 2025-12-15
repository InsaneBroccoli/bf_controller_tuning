# Closed-Loop Analysis

The goal of the function `calculate_closed_loop` is to calculate the differnt parts of the complete closed-loop control system, including its transfer functions, sensitivity functions, compliance and controller effort. This allows us to analyze the different parts of the control system and understand how they interact to shape the overall system behavior.  
This structure separates the control actions into an **inner loop** and an **outer loop**

---

## Transfer Functions

The transfer functions describe how a signal is transformed as it passes through the different components of the control system.
To calculate the transferfunction of the system, we have to split the control loop into an **inner** and an **outer** loop.
First wie have to calculate the transferfunction of the inner loop (path from `u` to `y`).

$$L_{in}(\omega) = P(\omega) \cdot G_{f}(\omega)$$ $$\qquad G_{yu}(\omega) = \frac{L_{in}(\omega)}{1+L_{in}(\omega) \cdot G_{fD}(\omega) \cdot C_{D}(\omega)}$$

With this Information we are able to calculate the hole closed loop system.

$$L_{out}(\omega) = G_{yu}(omega) \cdot C_{PI}(\omega)$$ $$\qquad G_{yw}(\omega) = \frac{L_{out}(\omega)}{1+L_{out}(\omega)}$$

<p align="center">
  <img src="./Images/Regelstrecke.png"
     alt="ZHAW Logo"
     width="800"
     style="float:center; margin-left:20px; margin-right:20px;
     margin-top:20px;">
</p>

---

## Sensitivity

The sensitivity describes how the closed-loop system reacts to disturbances and noise. To achieve good disturbance rejection and accurate reference tracking, the sensitivity should be low at low frequencies, while higher sensitivity at high frequencies helps to avoid noise amplification and maintain robustness [1][2].

To calculate the sensitivity functions of the inner and outer loop, we use the following formulas:

$$S_{in}\(omega) = \frac{1}{1+L_{in}(\omega)}$$ $$\qquad$$ $$S_{out}(\omega) = \frac{1}{1+L_{out}(\omega)}$$

---

## Controller Effort

The controller effort indicates how much effort the controller has to exert to maintain the desired system performance. It is important to monitor this value to ensure that the controller does not exceed its physical limits, which could lead to instability or damage to the system [1][2].

The controller effort can be calculated using the following formula:

$$SC(\omega) = C_{PI}(\omega) \cdot S_{out}(\omega)$$

---

## Compliance

The compliance describes how the system responds to external disturbances or changes in the reference input. A high compliance indicates that the system is able to adapt quickly to changes, while a low compliance indicates that the system is more rigid and less responsive [1][2].

The compliance can be calculated using the following formula:

$$SP(\omega) = G_{yu}(\omega) \cdot S_{out}(\omega)$$

---

## References

[1] R. Pintelon, J. Schoukens, *System Identification: A Frequency Domain Approach*, 2nd ed., Wiley-IEEE Press, 2012, **pp. 35–39, 150–180.**  
[2] M. H. Hayes, *Statistical Digital Signal Processing and Modeling*, Wiley, 1996, **pp. 192–196, 250–280.**
