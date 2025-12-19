# PID Controller Construction

Before performing control-theoretic analysis on Betaflight PID loops, the PID gains must be converted from their Betaflight-specific representation into physically meaningful controller gains. Betaflight applies internal scaling factors to its P, I and D terms, so the raw values shown in the configurator or CLI are not the actual gains used inside the control loop [1].

To obtain correct controller behavior in simulations, frequency-domain models, or closed-loop analysis, these gains have to be rescaled and then used to construct discrete-time PI and D controllers.

---

## Betaflight PID Gain Scaling

Betaflight internally scales its PID gains using fixed numerical factors [1].  
To convert Betaflight PID values into physically meaningful controller gains, the following scale factors are applied:

This scaling comes **directly from Betaflight’s internal controller implementation** and ensures that the converted PID values correspond to real control gains [1].

---

## Scaling of New PID Gains

Given the original PID vector:

$$PID = [K_p, K_i, K_d]$$

the scaled PID values are obtained as:

This step converts user-facing Betaflight PID values into controller gains suitable for analysis.

---

## Construction of PI and D Controllers

The discrete-time controllers (with sample time $T_s$) are constructed as:

### PI Controller

$$C_{PI} = K_p\ \cdot G_f + K_i \cdot T_s \cdot \frac{1}{1-z^{-1}}$$

### D Controller

$$C_D = \frac{K_d}{T_s}(1 - z^{-1})$$

The output of the implementation is:

- **$C_{PI}$** — proportional + integral controller  
- **$C_D$** — discrete derivative controller

With these newly computed controllers, it becomes possible to predict how the system will respond when operated under the updated control gains. [2]

---

## References

[1] Betaflight Development Team, *Betaflight PID Controller Implementation*, GitHub Repository.  
https://github.com/betaflight/betaflight  

[2] H. Lutz and W. Wendt, Taschenbuch der Regelungstechnik, 6th ed., Springer Vieweg, Berlin, Heidelberg, pp. 183–186.
