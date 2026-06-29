# Von Bertalanffy Growth

This model represents surface-area-limited growth with volume-proportional loss:

```text
dW/dt = a*W^(2/3) - b*W
```

where `W` is biomass, mass, or volume. The anabolic term scales with surface area
and the catabolic term scales with mass.

Use the cube-root transform:

```text
Y = W^(1/3)
```

Then:

```text
dY/dt = (a - b*Y)/3
```

so:

```text
Y(t) = a/b - (a/b - W0^(1/3))*exp(-b*t/3)
W(t) = Y(t)^3
```

The BNGL file reports this closed-form trajectory as `Analytical_W()` in the
`begin functions` block.

## Files

- `von_bertalanffy_growth.bngl` - BNGL model for the nonlinear growth ODE.
- `verify_von_bertalanffy_growth.py` - comparison of BioNetGen ODE output against
  the analytical solution.
