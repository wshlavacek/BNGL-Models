# Logistic Growth

This model implements the Verhulst-Pearl logistic growth equation discussed by
Samoletov and Vasiev (2023), "A statistical interpretation of biologically
inspired growth models."

The ODE is:

```text
dN/dt = r*N*(1 - N/K)
```

where `N` is the population density, `r` is the intrinsic per-capita growth
rate, and `K` is the carrying capacity.

For `N(0)=N0`, the analytical solution is:

```text
N(t) = K / (1 + A*exp(-r*t))
A = (K - N0)/N0
```

The abundance at the inflection point is:

```text
N_inflection = K/2
```

The BNGL model decomposes the ODE into an exponential production term and a
density-dependent crowding loss term:

```text
0 -> N      rate r*N
N -> 0      per-reactant rate r*N/K
```

The paper also derives stochastic population-environment analogues of logistic
and Gompertz growth. Those systems are not included here because BioNetGen's
deterministic ODE simulator does not represent Gaussian white-noise forcing.

## Files

- `logistic_growth.bngl` - BNGL ODE model for Verhulst-Pearl logistic growth.
- `verify_logistic_growth.py` - comparison of BioNetGen ODE output against the
  analytical solution.
