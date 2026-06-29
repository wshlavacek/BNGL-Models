# HIV ART Viral Decay

This model represents plasma HIV viral load after antiretroviral therapy begins at
`t = 0`.

Before treatment, viral production is nonzero and constant:

```text
dV/dt = P_pre - c V
```

The initial condition is the pre-treatment steady state:

```text
V(0) = P_pre / c
```

At `t = 0`, therapy is assumed to be 100% efficient, so new viral production is fully
blocked:

```text
P_art = 0
dV/dt = -c V
```

The analytical solution for the simulated post-treatment interval is:

```text
V(t) = (P_pre / c) exp(-c t)
```

## Files

- `hiv_art_viral_decay.bngl` - BNGL model for the post-treatment ODE.
- `verify_hiv_art_viral_decay.py` - comparison of BioNetGen ODE output against the
  analytical solution. BioNetGen writes both `Obs_V` and the BNGL function
  `Analytical_V()` to `.gdat` when `print_functions=>1` is enabled.
