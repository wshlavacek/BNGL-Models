# Quadratic Finite-Time Growth

This model implements the scalar nonlinear ODE used as a simple example by
Caravelli and Lin (2023), "On the combinatorics of Lotka-Volterra equations."

The ODE is:

```text
dX/dt = k*X^2
```

For `X(0)=X0`, the analytical solution is:

```text
X(t) = X0 / (1 - k*X0*t)
```

For `k > 0` and `X0 > 0`, the solution blows up at:

```text
t_blowup = 1/(k*X0)
```

The BNGL model uses a nonlinear production flux:

```text
0 -> X      rate k*X^2
```

The simulation interval is kept below `t_blowup` so that the analytical expression
remains finite.

## Files

- `quadratic_finite_time_growth.bngl` - BNGL ODE model for quadratic finite-time
  growth.
- `verify_quadratic_finite_time_growth.py` - comparison of BioNetGen ODE output
  against the analytical solution.
