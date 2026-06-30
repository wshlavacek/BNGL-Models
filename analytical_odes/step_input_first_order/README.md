# Step-Input First-Order System

This model implements a first-order linear system driven by a step input:

```text
dX/dt = J(t) - k*X
J(t) = J_base              for t < tau
J(t) = J_base + J_step     for t >= tau
```

The analytical solution is two exponentials joined continuously at `tau`:

```text
X(t) = X_ss0 + (X0 - X_ss0)*exp(-k*t)                         t < tau
X(t) = X_ss1 + (X_tau - X_ss1)*exp(-k*(t - tau))               t >= tau
```

where `X_ss0 = J_base/k`, `X_ss1 = (J_base + J_step)/k`, and
`X_tau = X_ss0 + (X0 - X_ss0)*exp(-k*tau)`.

## Files

- `step_input_first_order.bngl` - BNGL ODE model.
- `verify_step_input_first_order.py` - BioNetGen-vs-analytical check.
