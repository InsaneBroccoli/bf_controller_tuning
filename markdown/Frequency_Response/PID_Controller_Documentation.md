# PID Controller Construction

Before performing control-theoretic analysis on Betaflight PID loops, the PID gains must be converted from their Betaflight-specific representation into physically meaningful controller gains. Betaflight applies internal scaling factors to its P, I, and D terms, so the raw values shown in the configurator or CLI are not the actual gains used inside the control loop.

To obtain correct controller behavior in simulations, frequency-domain models, or closed-loop analysis, these gains have to be rescaled and then used to construct discrete-time PI and D controllers. 

---

## 1. Betaflight PID Gain Scaling

Betaflight internally scales its PID gains using fixed numerical factors.  
To convert Betaflight PID values into physically meaningful controller gains, the following scale factors are applied:

```
PTERM_SCALE = 0.032029  
ITERM_SCALE = 0.244381  
DTERM_SCALE = 0.000529
```

This scaling comes **directly from Betaflight’s internal controller implementation** and ensures that the converted PID values correspond to real control gains.

---

## 2. Scaling of New PID Gains

Given the original PID vector:

$$PID = [K_p,\; K_i,\; K_d]$$

the scaled PID values are obtained as:

```
PID_new(1) = P_new * pid_scale(1)
PID_new(2) = I_new * pid_scale(2)
PID_new(3) = D_new * pid_scale(3)
PID_new(4) = 0   % no feedforward
```

This step converts user-facing Betaflight PID values into controller gains suitable for analysis.

---

## 3. Construction of PI and D Controllers

The discrete-time controllers (with sample time $T_s$) are constructed as:

### PI Controller
$C_{PI} = K_p\,G_f + K_i\,T_s \cdot \frac{1}{1 - z^{-1}}$

### D Controller
$C_D = \frac{K_d}{T_s}(1 - z^{-1})$

The output of the implementation is:

- **Cpi** — proportional + integral controller 
- **Cd** — discrete derivative controller

With these newly computed controllers, it becomes possible to predict how the system will respond when operated under the updated control gains.
