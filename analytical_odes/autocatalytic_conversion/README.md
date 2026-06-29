# Autocatalytic Conversion

This model represents irreversible autocatalytic conversion:

```text
A + B -> 2B
```

with mass-action rate constant `k_auto`:

```text
dA/dt = -k_auto A B
dB/dt =  k_auto A B
```

The total concentration is conserved:

```text
T = A + B = A0 + B0
```

so the model reduces to a logistic equation for `B(t)`:

```text
dB/dt = k_auto (T - B) B
```

For `B0 > 0`, the analytical solution is:

```text
B(t) = T / (1 + (A0/B0) exp(-k_auto*T*t))
A(t) = T - B(t)
```

The BNGL file reports these closed-form trajectories as `Analytical_A()` and
`Analytical_B()` in the `begin functions` block.

## Files

- `autocatalytic_conversion.bngl` - BNGL model for the autocatalytic conversion
  ODE system.
- `verify_autocatalytic_conversion.py` - comparison of BioNetGen ODE output
  against the analytical solution.
