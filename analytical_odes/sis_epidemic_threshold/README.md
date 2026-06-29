# SIS Epidemic Threshold

This model represents a deterministic SIS epidemic with mass-action infection
and recovery:

```text
S + I -> 2I
I -> S
```

The total population is conserved:

```text
N = S + I = S0 + I0
```

The ODEs are:

```text
dS/dt = -beta*S*I/N + gamma*I
dI/dt =  beta*S*I/N - gamma*I
```

Using `S = N - I`, the infected population follows:

```text
dI/dt = (beta - gamma)*I - (beta/N)*I^2
```

For `beta > gamma` and `I0 > 0`, define:

```text
r = beta - gamma
I_star = N*(beta - gamma)/beta
R0 = (I_star - I0)/I0
```

Then:

```text
I(t) = I_star/(1 + R0*exp(-r*t))
S(t) = N - I(t)
```

The endemic threshold is visible in the steady state: if `beta > gamma`, the
infection approaches `I_star`; if `beta <= gamma`, infection decays to zero
with a different limiting formula.

The BNGL file reports the closed-form trajectories as `Analytical_S()` and
`Analytical_I()` for the supercritical case `beta > gamma`.

## Files

- `sis_epidemic_threshold.bngl` - BNGL model for deterministic SIS dynamics.
- `verify_sis_epidemic_threshold.py` - comparison of BioNetGen ODE output
  against the analytical solution.
