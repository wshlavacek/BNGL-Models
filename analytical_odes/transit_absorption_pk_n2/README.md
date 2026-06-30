# Transit-Compartment Absorption Pharmacokinetics, N=2

This model is a small, elementary specialization of the transit-compartment
pharmacokinetic model in Hof and Bridge (2021), "Exact solutions and
equi-dosing regimen regions for multi-dose pharmacokinetics models with
transit compartments."

The source paper gives exact formulas for a general `n`-transit compartment
cascade and multi-dose input. Those formulas use lower incomplete gamma
functions. For catalog compatibility, this entry fixes `n = 2` and uses a
single bolus dose, giving a closed form with only exponentials and a
`t*exp(-k*t)` repeated-pole term.

The model is:

```text
Transit1 -> Transit2    rate k_transit
Transit2 -> Absorption  rate k_transit
Absorption -> Central   rate k_abs
Central -> 0            rate k_elim
```

with initial condition:

```text
Transit1(0) = F_bio*Dose
Transit2(0) = Absorption(0) = Central(0) = 0
```

The default parameters use the paper's example values:

```text
Dose = 3.5 mg
F_bio = 0.69
k_transit = 12.76 /hour
k_abs = 9.11 /hour
k_elim = 0.96 /hour
```

For `D = F_bio*Dose`, the first two compartments follow the Erlang-chain
solution:

```text
Transit1(t) = D*exp(-k_transit*t)
Transit2(t) = D*k_transit*t*exp(-k_transit*t)
```

The absorption and central solutions are included in the BNGL `functions`
block as `Analytical_Absorption()` and `Analytical_Central()`. The central
solution is the inverse Laplace transform of:

```text
D*k_transit^2*k_abs /
  ((s + k_transit)^2*(s + k_abs)*(s + k_elim))
```

## Files

- `transit_absorption_pk_n2.bngl` - BNGL ODE model for the single-bolus
  two-transit-compartment absorption PK system.
- `verify_transit_absorption_pk_n2.py` - comparison of BioNetGen ODE output
  against the analytical solution.
