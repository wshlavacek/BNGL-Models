# One-Compartment Infusion-Washout PK

This model implements one-compartment pharmacokinetics with a constant infusion
that stops at `t_infusion`:

```text
dA/dt = R_in - kel*A      for t < t_infusion
dA/dt = -kel*A            for t >= t_infusion
```

The analytical solution is:

```text
A(t) = A_ss + (A0 - A_ss)*exp(-kel*t)                         t < t_infusion
A(t) = A_stop*exp(-kel*(t - t_infusion))                       t >= t_infusion
```

where `A_ss = R_in/kel` and `A_stop` is the amount at infusion stop.

## Files

- `infusion_washout_pk.bngl` - BNGL ODE model.
- `verify_infusion_washout_pk.py` - BioNetGen-vs-analytical check.
