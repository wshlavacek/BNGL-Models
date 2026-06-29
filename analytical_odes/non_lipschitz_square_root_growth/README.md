# Non-Lipschitz Square-Root Growth

This model represents the scalar ODE:

```text
dX/dt = a*sqrt(X)
```

For `X0 > 0`, the analytical solution is unique:

```text
X(t) = (sqrt(X0) + a*t/2)^2
```

The interesting case is `X0 = 0`. The right-hand side is continuous but not
Lipschitz at zero, so uniqueness fails. There is a stationary solution:

```text
X(t) = 0
```

but for any delay time `tau >= 0`, there is also a delayed-start solution:

```text
X(t) = 0                         for t <= tau
X(t) = (a*(t - tau)/2)^2          for t > tau
```

This means the same initial condition can have infinitely many valid analytical
solutions. A typical ODE solver initialized exactly at `X0 = 0` will follow the
stationary branch unless perturbed.

The BNGL model uses `X0 > 0` by default so the simulated trajectory can be
verified against the unique nontrivial branch. It also reports two illustrative
zero-initial-condition analytical functions:

```text
ZeroInitial_Stationary_X()
ZeroInitial_Delayed_X()
```

## Files

- `non_lipschitz_square_root_growth.bngl` - BNGL ODE model for square-root growth.
- `verify_non_lipschitz_square_root_growth.py` - comparison of BioNetGen ODE
  output against the analytical solution.
