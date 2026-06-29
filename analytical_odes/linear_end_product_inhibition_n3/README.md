# Linear End-Product Inhibition Pathway, N=3

This model represents a three-stage biochemical pathway with linear end-product
feedback:

```text
x0 -> x1 -> x2 -> x3 -> 0
```

The upstream substrate `x0` is constant. The input flux into `x1` is inhibited
linearly by the end product `x3`:

```text
v0 = k0*x0 - h*x3
```

The downstream pathway uses a common first-order rate constant `k`:

```text
dx1/dt = k0*x0 - h*x3 - k*x1
dx2/dt = k*x1 - k*x2
dx3/dt = k*x2 - k*x3
```

The steady state is:

```text
x_ss = k0*x0 / (k + h)
x1_ss = x2_ss = x3_ss = x_ss
```

Define the deviations `y_i = x_i - x_ss`, and let:

```text
alpha = (h*k^2)^(1/3)
omega = sqrt(3)*alpha/2
```

For `y3`, write:

```text
y3(t) = exp(-k*t) W(t)
W(t) = A exp(-alpha*t)
     + exp(alpha*t/2) [B cos(omega*t) + C sin(omega*t)]
```

The coefficients are fixed by the initial conditions:

```text
u1 = X1_0 - x_ss
u2 = X2_0 - x_ss
u3 = X3_0 - x_ss

A = (k^2*u1 - alpha*k*u2 + alpha^2*u3) / (3*alpha^2)
B = u3 - A
C = (2*k*u2/alpha - u3 + 3*A) / sqrt(3)
```

Then:

```text
x3(t) = x_ss + exp(-k*t) W(t)
x2(t) = x_ss + exp(-k*t) W'(t)/k
x1(t) = x_ss + exp(-k*t) W''(t)/k^2
```

The BNGL file reports these closed-form trajectories as `Analytical_X1()`,
`Analytical_X2()`, and `Analytical_X3()` in the `begin functions` block.

## Files

- `linear_end_product_inhibition_n3.bngl` - BNGL model for the N=3 feedback
  pathway.
- `verify_linear_end_product_inhibition_n3.py` - comparison of BioNetGen ODE
  output against the analytical solution.
