# Reversible Bimolecular Association

This model represents reversible mass-action association:

```text
A + B <-> C
```

with forward rate constant `kf` and reverse rate constant `kr`:

```text
dA/dt = -kf A B + kr C
dB/dt = -kf A B + kr C
dC/dt =  kf A B - kr C
```

Initial conditions are specified directly:

```text
A(0) = A0
B(0) = B0
C(0) = C0
```

The conserved totals are:

```text
Atot = A0 + C0
Btot = B0 + C0
```

so `A(t) = Atot - C(t)` and `B(t) = Btot - C(t)`. The equation for `C(t)` is:

```text
dC/dt = kf (Atot - C)(Btot - C) - kr C
```

Define:

```text
S      = kf (Atot + Btot) + kr
D      = sqrt(S^2 - 4 kf^2 Atot Btot)
C_low  = (S - D) / (2 kf)
C_high = (S + D) / (2 kf)
R0     = (C0 - C_low) / (C0 - C_high)
R(t)   = R0 exp(-D t)
```

Then:

```text
C(t) = (C_low - R(t) C_high) / (1 - R(t))
A(t) = Atot - C(t)
B(t) = Btot - C(t)
```

The BNGL file reports these closed-form trajectories as `Analytical_A()`,
`Analytical_B()`, and `Analytical_C()` in the `begin functions` block.

## Files

- `reversible_bimolecular_association.bngl` - BNGL model for the reversible
  association ODE system.
- `verify_reversible_bimolecular_association.py` - comparison of BioNetGen ODE
  output against the analytical solution.
