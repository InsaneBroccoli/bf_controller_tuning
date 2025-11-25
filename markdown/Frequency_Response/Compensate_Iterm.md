# Compensate Iterm

## Iterm Relax in Betaflight

The I-term in a PID controller is responsible for compensating long-term deviations. For example, if the drone slightly drifts away from its target angle or if an axis does not reach the desired position due to external influences, the I-term accumulates these errors and ensures that the drone eventually moves into this position.
Iterm Relax, on the other hand, ensures that the I-term is weakened during rapid changes. During fast stick movements, large deviations occur between the stick and the gyro. However, these are completely normal and often intended by the pilot. In such cases, the I-term integrates this intentional error and overloads itself. This causes the drone to overshoot the target value.
Iterm Relax prevents this by suppressing the I-term during fast changes, but allowing it again during slow changes. This way, the drone can perform quick angle changes, while still correcting static errors during slower adjustments.

## Compensate Iterm in Controller Tuning

The I-term Relax in Betaflight is quite complex and not easy to understand just by looking at the settings. To get a simpler but still accurate idea of how it works, we first calculated the analytical PI controller response based on the original PID values. Then we compared this analytical response with the measured frequency response to see how the I-term Relax changes it.
By comparing the two responses, we were able to determine factors at all frequencies that represent the influence of the I-term Relax. This allows us to describe the effect of the I-term Relax without needing to explicitly recreate its internal implementation.

$$
C_{PI,com} = \frac{C_{PI}}{C_{PI,ana}}
$$

This factor is then applied to both the original analytical PI controller and the new PI controller design.
$$
C_{PI,ana,new} = C_{PI,ana} \cdot C_{Pi,com} \qquad
C_{PI,new,new} = C_{PI,new} \cdot C_{Pi,com}
$$
These newly calculated controllers must then be included again in the closed-loop system. This is done as described in the document [Calculate Closed Loop](/CalculateClosedLoop.md).
