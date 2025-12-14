# PID Controller Construction

Before performing control-theoretic analysis on Betaflight PID loops, the PID gains must be converted from their Betaflight-specific representation into physically meaningful controller gains. Betaflight applies internal scaling factors to its P, I, and D terms, so the raw values shown in the configurator or CLI are not the actual gains used inside the control loop [1].

To obtain correct controller behavior in simulations, frequency-domain models, or closed-loop analysis, these gains have to be rescaled and then used to construct discrete-time PI and D controllers.

---

## 1. Betaflight PID Gain Scaling

Betaflight internally scales its PID gains using fixed numerical factors defined in the firmware implementation [1].  
These scaling factors are applied internally before the PID terms act on the control error.

To convert Betaflight PID values into physically meaningful controller gains, the corresponding scale factors must therefore be applied explicitly when constructing analytical or simulation models.

This scaling is derived directly from Betaflight’s internal controller implementation and ensures that the converted PID values correspond to the effective gains used in the control loop [1].

---

## 2. Scaling of New PID Gains

Given the original PID gain vector as configured in Betaflight,

$PID = [K_p,\; K_i,\; K_d],$

the scaled PID values are obtained by applying the corresponding Betaflight scaling factors [1].

This step converts the user-facing Betaflight PID parameters into controller gains that are suitable for control-theoretic analysis and model-based evaluation.

---

## 3. Construction of PI and D Controllers

Using the scaled gains, the discrete-time controllers are constructed with the sample time $T_s$ according to standard discrete control formulations [2].

### PI Controller

$C_{PI}(z) = K_p + K_i\,T_s \frac{1}{1 - z^{-1}}$  [2]

### D Controller

$C_D(z) = \frac{K_d}{T_s} (1 - z^{-1})$   [2]

The resulting controllers are implemented as:

- **\(C_{PI}\)** — proportional–integral controller  
- **\(C_D\)** — discrete derivative controller  

With these newly computed controllers, it becomes possible to accurately predict the closed-loop behavior of the system under updated Betaflight PID gains using analytical, simulation, or frequency-domain methods.

---

## References

[1] Betaflight Development Team, *Betaflight Flight Controller Firmware – PID Controller Implementation*,  
GitHub Repository, https://github.com/betaflight/betaflight (accessed: 2025-12-14).

[2] H. Lutz, *Taschenbuch der Regelungstechnik*, Springer Vieweg, aktuelle Auflage.
