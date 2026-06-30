# Linearized Lotka-Volterra Near Equilibrium

This model implements the small-deviation linearization of the predator-prey
Lotka-Volterra equations near their coexistence equilibrium.

For the nonlinear system:

```text
dprey/dt     = alpha*prey - beta*prey*predator
dpredator/dt = delta*prey*predator - gamma*predator
```

the equilibrium is `prey_star = gamma/delta`, `predator_star = alpha/beta`.
For relative deviations `P` and `Q`, the linearized system is:

```text
dP/dt = -alpha*Q
dQ/dt =  gamma*P
```

with `omega = sqrt(alpha*gamma)`. The solution is sinusoidal:

```text
P(t) = P0*cos(omega*t) - alpha*Q0/omega*sin(omega*t)
Q(t) = Q0*cos(omega*t) + gamma*P0/omega*sin(omega*t)
```

## Files

- `linearized_lotka_volterra.bngl` - BNGL ODE model.
- `verify_linearized_lotka_volterra.py` - BioNetGen-vs-analytical check.
