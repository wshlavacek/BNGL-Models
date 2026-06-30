# Two-Species Harmonic Oscillator

This model implements the two species harmonic oscillator (2SHO) from
Hellerstein (2023), "An oscillating reaction network with an exact closed form
solution in the time domain."

The reaction network is:

```text
S1 -> S2       rate k1*S1
S2 -> S1       rate k2*S2
S1 -> 2S1      rate k3*S1
S1 -> 0        total rate k4
S2 -> 0        total rate k5*S1
0  -> S2       rate k6
```

The oscillator constraints are:

```text
k3 = k1 + k2
k5 = k3 + kd
kd > 0
```

With those constraints, the ODE is:

```text
dS1/dt =  k2*S1 + k2*S2 - k4
dS2/dt = -(k2 + kd)*S1 - k2*S2 + k6
```

or `x' = A*x + u`, with:

```text
A = [[k2, k2], [-(k2 + kd), -k2]]
u = [-k4, k6]
theta = sqrt(k2*kd)
```

Since `A^2 = -theta^2 I`, the forced solution is a sinusoid around the offset
`x_offset = -A^-1 u`:

```text
S1(t) = S1_offset + S1_delta0*cos(theta*t) + S1_sin_coeff*sin(theta*t)
S2(t) = S2_offset + S2_delta0*cos(theta*t) + S2_sin_coeff*sin(theta*t)
```

The default parameters use the paper's sensitivity-analysis example:

```text
k1 = 1.00
k2 = 3.91
k3 = 4.91
k4 = 41.77
k5 = 15.00
k6 = 92.20
S1(0) = 5.00
S2(0) = 10.49
```

These values give approximately:

```text
theta = 6.28
Amplitude_S1 = 3.00
Amplitude_S2 = 5.67
Offset_S1 = 5.00
Offset_S2 = 5.68
```

`R4` and `R5` are encoded as signed zero-reactant ODE contribution rules because
their published fluxes are total rates. If they were written as ordinary
degradation rules, BioNetGen would multiply the functional rate law by the
consumed species concentration.

## Files

- `two_species_harmonic_oscillator.bngl` - BNGL ODE model for the 2SHO reaction
  network.
- `verify_two_species_harmonic_oscillator.py` - comparison of BioNetGen ODE
  output against the analytical solution.
