# Logistic Growth With Constant Harvesting

This model implements logistic growth with a constant harvest:

```text
dX/dt = r*X*(1 - X/K) - h
```

For `0 < h < r*K/4`, the quadratic right-hand side has two real roots:

```text
D       = sqrt(r^2 - 4*(r/K)*h)
X_plus  = (r + D)/(2*r/K)
X_minus = (r - D)/(2*r/K)
```

The transformed root ratio decays exponentially:

```text
R(t) = ((X0 - X_plus)/(X0 - X_minus))*exp(-D*t)
X(t) = (X_plus - R(t)*X_minus)/(1 - R(t))
```

The BNGL model uses a signed global total-flux rule and is ODE-only.

## Files

- `logistic_harvesting.bngl` - BNGL ODE model.
- `verify_logistic_harvesting.py` - BioNetGen-vs-analytical check.
