# Three-State Irreversible Cycle

This model represents a first-order irreversible cycle:

```text
X1 -> X2 -> X3 -> X1
```

with common rate constant `k`:

```text
dX1/dt = k*X3 - k*X1
dX2/dt = k*X1 - k*X2
dX3/dt = k*X2 - k*X3
```

The total concentration is conserved:

```text
T = X1 + X2 + X3 = X1_0 + X2_0 + X3_0
```

The steady state is:

```text
X_ss = T/3
```

Define initial deviations from steady state:

```text
u1 = X1_0 - X_ss
u2 = X2_0 - X_ss
u3 = X3_0 - X_ss
```

The nonzero eigenvalues are a complex-conjugate pair:

```text
lambda = -3*k/2 +/- i*sqrt(3)*k/2
```

so the real-valued solution is:

```text
gamma = 3*k/2
omega = sqrt(3)*k/2

X1(t) = X_ss + exp(-gamma*t) * [u1*cos(omega*t) + ((u3-u2)/sqrt(3))*sin(omega*t)]
X2(t) = X_ss + exp(-gamma*t) * [u2*cos(omega*t) + ((u1-u3)/sqrt(3))*sin(omega*t)]
X3(t) = X_ss + exp(-gamma*t) * [u3*cos(omega*t) + ((u2-u1)/sqrt(3))*sin(omega*t)]
```

The BNGL file reports these closed-form trajectories as `Analytical_X1()`,
`Analytical_X2()`, and `Analytical_X3()` in the `begin functions` block.

## Files

- `three_state_irreversible_cycle.bngl` - BNGL model for the irreversible cycle.
- `verify_three_state_irreversible_cycle.py` - comparison of BioNetGen ODE output
  against the analytical solution.
