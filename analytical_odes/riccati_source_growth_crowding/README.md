# Riccati Source-Growth-Crowding

This model represents a scalar Riccati equation with constant source, linear
growth, and quadratic crowding loss:

```text
dX/dt = s + r*X - q*X^2
```

where `s > 0`, `r > 0`, and `q > 0`.

The two roots of the quadratic right-hand side are:

```text
D       = sqrt(r^2 + 4*q*s)
X_plus  = (r + D)/(2*q)
X_minus = (r - D)/(2*q)
```

The ODE can be written as:

```text
dX/dt = -q*(X - X_plus)*(X - X_minus)
```

Since `q*(X_plus - X_minus) = D`, the transformed ratio evolves exponentially:

```text
R(t) = ((X0 - X_plus)/(X0 - X_minus))*exp(-D*t)
```

and the explicit solution is:

```text
X(t) = (X_plus - R(t)*X_minus)/(1 - R(t))
```

The BNGL file reports this closed-form trajectory as `Analytical_X()` in the
`begin functions` block.

## Files

- `riccati_source_growth_crowding.bngl` - BNGL model for the Riccati ODE.
- `verify_riccati_source_growth_crowding.py` - comparison of BioNetGen ODE output
  against the analytical solution.
