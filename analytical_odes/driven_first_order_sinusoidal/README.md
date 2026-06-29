# Sinusoidally Driven First-Order System

This model represents a first-order decay process driven by a sinusoidal input:

```text
dX/dt = J0 + J1*sin(omega*t) - k*X
```

With:

```text
denom = k^2 + omega^2
X_forced(t) = J0/k + J1*(k*sin(omega*t) - omega*cos(omega*t))/denom
X_forced_0 = J0/k - J1*omega/denom
```

the analytical solution is:

```text
X(t) = X_forced(t) + (X0 - X_forced_0)*exp(-k*t)
```

The transient term decays, leaving a phase-lagged sinusoidal steady response.

The BNGL file reports the closed-form trajectory as `Analytical_X()` in the
`begin functions` block.

## Files

- `driven_first_order_sinusoidal.bngl` - BNGL model for the forced first-order
  ODE system.
- `verify_driven_first_order_sinusoidal.py` - comparison of BioNetGen ODE output
  against the analytical solution.
