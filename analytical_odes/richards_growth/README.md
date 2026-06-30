# Richards Growth

This model implements the Richards growth equation reviewed by Tsoularis and
Wallace (2002), "Analysis of logistic growth models."

The ODE is:

```text
dN/dt = r*N*(1 - (N/K)^b)
```

where `r` is the intrinsic per-capita growth rate, `K` is the carrying
capacity, and `b` controls the curve shape and inflection point. The model is
the `a = c = 1` special case of the paper's generalized logistic equation.

For `N(0)=N0`, the analytical solution is:

```text
N(t) = K / (1 + A*exp(-b*r*t))^(1/b)
A = (K/N0)^b - 1
```

The abundance at the inflection point is:

```text
N_inflection = K / (1 + b)^(1/b)
```

The BNGL model decomposes the ODE into an exponential production term and a
density-dependent loss term:

```text
0 -> N      rate r*N
N -> 0      per-reactant rate r*(N/K)^b
```

## Files

- `richards_growth.bngl` - BNGL ODE model for Richards growth.
- `verify_richards_growth.py` - comparison of BioNetGen ODE output against
  the analytical solution.
