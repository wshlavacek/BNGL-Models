# Bateman Chain With Distinct Rates

This model implements a three-state irreversible first-order chain with distinct
rates:

```text
A -> B -> C
```

For `k1 != k2`, `A(0)=A0`, `B(0)=0`, and `C(0)=0`, the Bateman solution is:

```text
A(t) = A0*exp(-k1*t)
B(t) = k1*A0/(k2-k1)*(exp(-k1*t) - exp(-k2*t))
C(t) = A0 - A(t) - B(t)
```

This complements the equal-rate Erlang-chain catalog entry.

## Files

- `bateman_chain_distinct_rates.bngl` - BNGL ODE model.
- `verify_bateman_chain_distinct_rates.py` - BioNetGen-vs-analytical check.
