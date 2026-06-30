# Reversible First-Order Conversion

This model implements a two-state reversible first-order conversion:

```text
A <-> B
dA/dt = -kf*A + kr*B
dB/dt =  kf*A - kr*B
```

The total `T = A + B` is conserved. With `s = kf + kr`, the equilibrium is:

```text
A_eq = kr*T/s
B_eq = kf*T/s
```

and the analytical solution is:

```text
A(t) = A_eq + (A0 - A_eq)*exp(-s*t)
B(t) = B_eq + (B0 - B_eq)*exp(-s*t)
```

## Files

- `reversible_first_order_conversion.bngl` - BNGL ODE model.
- `verify_reversible_first_order_conversion.py` - BioNetGen-vs-analytical check.
