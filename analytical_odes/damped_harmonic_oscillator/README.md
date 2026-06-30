# Damped Harmonic Oscillator

This model implements the underdamped second-order oscillator as a two-variable
first-order ODE system:

```text
dx/dt = v
dv/dt = -omega0^2*x - 2*zeta*omega0*v
```

For `0 < zeta < 1`, define `a = zeta*omega0` and
`omega_d = omega0*sqrt(1-zeta^2)`. The analytical solution is:

```text
x(t) = exp(-a*t)*(x0*cos(omega_d*t) + ((v0 + a*x0)/omega_d)*sin(omega_d*t))
v(t) = exp(-a*t)*(v0*cos(omega_d*t) - ((omega0^2*x0 + a*v0)/omega_d)*sin(omega_d*t))
```

The BNGL file uses signed ODE contribution rules and is intended only for ODE
simulation.

## Files

- `damped_harmonic_oscillator.bngl` - BNGL ODE model.
- `verify_damped_harmonic_oscillator.py` - BioNetGen-vs-analytical check.
