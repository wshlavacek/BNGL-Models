# Piecewise-Linear Gene Switch

This model implements the two-gene piecewise-linear regulatory network from
Casey, de Jong, and Gouze (2006), "Piecewise-linear Models of Genetic Regulatory
Networks: Equilibria and their Stability."

The full model is:

```text
dA/dt = kappa_A*s_plus(B, theta_B_low)*s_minus(A, theta_A_high) - gamma_A*A
dB/dt = kappa_B*s_plus(A, theta_A_low)*s_minus(B, theta_B_high) - gamma_B*B
```

with step functions:

```text
s_plus(X, theta)  = 1 if X > theta, otherwise 0
s_minus(X, theta) = 1 if X < theta, otherwise 0
```

The default initial condition is chosen in the regulatory domain
`A0 < theta_A_low` and `theta_B_low < B0 < theta_B_high`, with `B` crossing
`theta_B_low` before `A` reaches `theta_A_low`. This gives one known switch time:

```text
t_switch = ln(B0/theta_B_low)/gamma_B
```

For `0 <= t < t_switch`:

```text
A(t) = A_on + (A0 - A_on)*exp(-gamma_A*t)
B(t) = B0*exp(-gamma_B*t)
A_on = kappa_A/gamma_A
```

For `t >= t_switch`:

```text
A(t) = A_switch*exp(-gamma_A*(t - t_switch))
B(t) = theta_B_low*exp(-gamma_B*(t - t_switch))
```

where `A_switch = A_on + (A0 - A_on)*exp(-gamma_A*t_switch)`.

The paper's main stability results concern Filippov solutions and singular
equilibrium sets on threshold hyperplanes. Those set-valued differential
inclusion dynamics are not represented here. This catalog entry uses BioNetGen
global functions to encode the PLDE away from ambiguous sliding modes.

## Files

- `piecewise_linear_gene_switch.bngl` - BNGL ODE model for the two-gene
  piecewise-linear regulatory network.
- `verify_piecewise_linear_gene_switch.py` - comparison of BioNetGen ODE output
  against the explicit piecewise solution for the default trajectory.
