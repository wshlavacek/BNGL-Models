# Quantum Two-Level Rabi Oscillation

This model represents resonant coherent dynamics of a closed quantum two-level
system in Bloch-coordinate form. It is an ODE-only model with signed continuous
state variables, not a stochastic chemical population model.

For resonant drive with Rabi frequency `Omega`, the relevant Bloch equations are:

```text
dU/dt = 0
dV/dt = -Omega*W
dW/dt =  Omega*V
```

where:

```text
W = P_excited - P_ground
P_excited = (1 + W)/2
P_ground  = (1 - W)/2
```

For general initial conditions `V(0)=V0` and `W(0)=W0`:

```text
V(t) = V0*cos(Omega*t) - W0*sin(Omega*t)
W(t) = V0*sin(Omega*t) + W0*cos(Omega*t)
U(t) = U0
```

The default initial condition is the ground state:

```text
U0 = 0
V0 = 0
W0 = -1
```

so:

```text
P_excited(t) = sin^2(Omega*t/2)
P_ground(t)  = cos^2(Omega*t/2)
```

The BNGL file reports the closed-form trajectories as `Analytical_U()`,
`Analytical_V()`, `Analytical_W()`, `Analytical_P_Excited()`, and
`Analytical_P_Ground()`.

## Files

- `quantum_rabi_two_level.bngl` - BNGL ODE model for resonant two-level Rabi
  oscillation.
- `verify_quantum_rabi_two_level.py` - comparison of BioNetGen ODE output against
  the analytical solution.
