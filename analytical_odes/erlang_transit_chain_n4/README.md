# Erlang Transit-Chain Delay, N=4

This model represents a four-stage first-order transit chain:

```text
T0 -> T1 -> T2 -> T3 -> Product
```

with common rate constant `k` and initial dose in `T0`:

```text
T0(0) = Dose
T1(0) = T2(0) = T3(0) = Product(0) = 0
```

The ODE system is:

```text
dT0/dt      = -k*T0
dT1/dt      =  k*T0 - k*T1
dT2/dt      =  k*T1 - k*T2
dT3/dt      =  k*T2 - k*T3
dProduct/dt =  k*T3
```

Let:

```text
tau = k*t
```

Then:

```text
T0(t)      = Dose*exp(-tau)
T1(t)      = Dose*exp(-tau)*tau
T2(t)      = Dose*exp(-tau)*tau^2/2
T3(t)      = Dose*exp(-tau)*tau^3/6
Product(t) = Dose*(1 - exp(-tau)*(1 + tau + tau^2/2 + tau^3/6))
```

The product trajectory is the cumulative distribution function of an Erlang
waiting-time distribution with shape 4 and rate `k`. This makes the model a
compact analytical representation of a distributed maturation or transport delay.

The BNGL file reports the closed-form trajectories as `Analytical_T0()`,
`Analytical_T1()`, `Analytical_T2()`, `Analytical_T3()`, and
`Analytical_Product()` in the `begin functions` block.

## Files

- `erlang_transit_chain_n4.bngl` - BNGL model for the four-stage transit chain.
- `verify_erlang_transit_chain_n4.py` - comparison of BioNetGen ODE output
  against the analytical solution.
