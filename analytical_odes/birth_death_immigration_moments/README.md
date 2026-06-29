# Birth-Death-Immigration Moments

This model represents the first two moments of a linear birth-death-immigration
branching process:

```text
0 -> X      immigration, rate s
X -> 2X     birth, rate b per individual
X -> 0      death, rate d per individual
```

Let:

```text
m(t) = E[X(t)]
v(t) = Var[X(t)]
```

The moment equations close exactly:

```text
dm/dt = s + (b - d)*m
dv/dt = s + (b + d)*m + 2*(b - d)*v
```

Define:

```text
a        = b - d
c        = b + d
m_ss     = -s/a
m_delta  = M0 - m_ss
v_const  = -(s + c*m_ss)/(2*a)
v_exp1   = -(c*m_delta)/a
v_exp2   = Var0 - v_const - v_exp1
```

For `a != 0`, the analytical solution is:

```text
m(t) = m_ss + m_delta*exp(a*t)
v(t) = v_const + v_exp1*exp(a*t) + v_exp2*exp(2*a*t)
```

The default parameters use `b < d`, so the process is subcritical and the moments
relax to finite steady values. The critical case `b = d` is also analytically
solvable, but uses a separate polynomial limit and is not encoded in this BNGL file.

The BNGL file reports the closed-form trajectories as `Analytical_Mean()` and
`Analytical_Variance()`.

## Files

- `birth_death_immigration_moments.bngl` - BNGL model for the closed moment ODEs.
- `verify_birth_death_immigration_moments.py` - comparison of BioNetGen ODE output
  against the analytical solution.
