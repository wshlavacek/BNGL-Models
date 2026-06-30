# Early Linearized SEIR

This model implements the early-epidemic linearization of an SEIR model with
approximately constant susceptible population:

```text
dE/dt = beta_eff*I - sigma*E
dI/dt = sigma*E - gamma*I
dR/dt = gamma*I
```

The exposed and infectious compartments are a two-mode linear system. If
`lambda_1` and `lambda_2` are the two eigenvalues, then:

```text
E(t) = c1*exp(lambda_1*t) + c2*exp(lambda_2*t)
I(t) = i1*exp(lambda_1*t) + i2*exp(lambda_2*t)
R(t) = R0 + gamma*(i1*(exp(lambda_1*t)-1)/lambda_1
                 + i2*(exp(lambda_2*t)-1)/lambda_2)
```

## Files

- `early_seir_linearized.bngl` - BNGL ODE model.
- `verify_early_seir_linearized.py` - BioNetGen-vs-analytical check.
