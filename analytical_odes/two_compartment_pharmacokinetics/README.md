# Two-Compartment Pharmacokinetics

This model represents a standard linear two-compartment pharmacokinetic system:

```text
Central <-> Peripheral
Central -> 0
```

with first-order transfer and elimination:

```text
dCentral/dt    = -(k12 + ke)*Central + k21*Peripheral
dPeripheral/dt =  k12*Central - k21*Peripheral
```

The coefficient matrix has two eigenvalues:

```text
s          = k12 + k21 + ke
d          = sqrt(s^2 - 4*ke*k21)
lambda_1   = (-s + d)/2
lambda_2   = (-s - d)/2
```

For initial values `Central_0` and `Peripheral_0`, define:

```text
Central_dot_0 = -(k12 + ke)*Central_0 + k21*Peripheral_0
coef_1 = (Central_dot_0 - lambda_2*Central_0)/(lambda_1 - lambda_2)
coef_2 = (lambda_1*Central_0 - Central_dot_0)/(lambda_1 - lambda_2)
```

Then:

```text
Central(t) = coef_1*exp(lambda_1*t) + coef_2*exp(lambda_2*t)
Peripheral(t) =
  ((lambda_1 + k12 + ke)/k21)*coef_1*exp(lambda_1*t)
  + ((lambda_2 + k12 + ke)/k21)*coef_2*exp(lambda_2*t)
```

The BNGL file reports these closed-form trajectories as `Analytical_Central()`
and `Analytical_Peripheral()`.

## Files

- `two_compartment_pharmacokinetics.bngl` - BNGL model for the biexponential
  two-compartment system.
- `verify_two_compartment_pharmacokinetics.py` - comparison of BioNetGen ODE
  output against the analytical solution.
