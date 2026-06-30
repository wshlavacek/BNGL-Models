# Reversible Linear Chain, Three States

This model implements a symmetric reversible first-order chain:

```text
A <-> B <-> C
```

with the same rate `k` on every edge. The total `T=A+B+C` is conserved, and two
orthogonal modes decay with rates `k` and `3*k`:

```text
D(t) = (A0 - C0)*exp(-k*t)
E(t) = (A0 - 2*B0 + C0)*exp(-3*k*t)
B(t) = (T - E(t))/3
A(t) = (2*T + E(t) + 3*D(t))/6
C(t) = (2*T + E(t) - 3*D(t))/6
```

## Files

- `reversible_linear_chain_three_state.bngl` - BNGL ODE model.
- `verify_reversible_linear_chain_three_state.py` - BioNetGen-vs-analytical check.
