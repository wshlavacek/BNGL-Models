# Gompertz Growth

This model represents Gompertz growth:

```text
dX/dt = r*X*ln(K/X)
```

where `K` is the carrying capacity and `r` is the growth-rate parameter. For
`0 < X0 < K`, growth is positive and slows as `X` approaches `K`.

Use:

```text
Y = ln(X/K)
```

Then:

```text
dY/dt = -r*Y
```

so:

```text
Y(t) = ln(X0/K)*exp(-r*t)
X(t) = K*exp(ln(X0/K)*exp(-r*t))
```

The BNGL file reports this closed-form trajectory as `Analytical_X()` in the
`begin functions` block.

## Files

- `gompertz_growth.bngl` - BNGL model for the nonlinear growth ODE.
- `verify_gompertz_growth.py` - comparison of BioNetGen ODE output against the
  analytical solution.
