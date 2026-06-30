# Analytical ODE BNGL Models

This directory is a catalog of small BNGL models whose ODE semantics have closed-form
analytical solutions. The intent is to make each model easy to inspect, simulate, and
verify against an independently implemented formula.

## Catalog

| Model | Governing equation | Analytical solution | BNGL file |
| --- | --- | --- | --- |
| HIV viral decay under fully suppressive antiretroviral therapy | `dV/dt = P_art - c V`, with `P_art = 0` for `t >= 0` and `V(0) = P_pre/c` | `V(t) = (P_pre/c) exp(-c t)` | `hiv_art_viral_decay/hiv_art_viral_decay.bngl` |
| Reversible bimolecular association | `A + B <-> C` with `dC/dt = kf A B - kr C` | Riccati solution for `C(t)` using conserved totals `A+C` and `B+C` | `reversible_bimolecular_association/reversible_bimolecular_association.bngl` |
| Linear end-product inhibition pathway, N=3 | `x0 -> x1 -> x2 -> x3 -> 0`, with input `k0*x0 - h*x3` | Matrix-exponential solution specialized to equal downstream rate `k` | `linear_end_product_inhibition_n3/linear_end_product_inhibition_n3.bngl` |
| Autocatalytic conversion | `A + B -> 2B` with conserved total `A+B` | Logistic solution for `B(t)` | `autocatalytic_conversion/autocatalytic_conversion.bngl` |
| Three-state irreversible cycle | `X1 -> X2 -> X3 -> X1` with equal first-order rate `k` | Conserved total plus damped sine/cosine relaxation modes | `three_state_irreversible_cycle/three_state_irreversible_cycle.bngl` |
| Von Bertalanffy growth | `dW/dt = a*W^(2/3) - b*W` | Cube-root transform solution | `von_bertalanffy_growth/von_bertalanffy_growth.bngl` |
| Logistic growth | `dN/dt = r*N*(1 - N/K)` | Verhulst-Pearl rational-exponential solution | `logistic_growth/logistic_growth.bngl` |
| Piecewise-linear gene switch | Two-gene threshold-regulated synthesis with first-order degradation | Explicit piecewise exponentials for a one-switch trajectory | `piecewise_linear_gene_switch/piecewise_linear_gene_switch.bngl` |
| Two-compartment pharmacokinetics | `Central <-> Peripheral`, `Central -> 0` | Biexponential solution from the two eigenmodes | `two_compartment_pharmacokinetics/two_compartment_pharmacokinetics.bngl` |
| Gompertz growth | `dX/dt = r*X*ln(K/X)` | Log-transform solution | `gompertz_growth/gompertz_growth.bngl` |
| Sinusoidally driven first-order system | `dX/dt = J0 + J1*sin(omega*t) - k*X` | Transient plus phase-lagged sinusoidal steady response | `driven_first_order_sinusoidal/driven_first_order_sinusoidal.bngl` |
| Riccati source-growth-crowding | `dX/dt = s + r*X - q*X^2` | Two-root Riccati solution | `riccati_source_growth_crowding/riccati_source_growth_crowding.bngl` |
| Erlang transit-chain delay, N=4 | `T0 -> T1 -> T2 -> T3 -> Product` with equal rate `k` | Finite polynomial times exponential; product follows an Erlang CDF | `erlang_transit_chain_n4/erlang_transit_chain_n4.bngl` |
| Birth-death-immigration moments | `0 -> X`, `X -> 2X`, `X -> 0`; mean and variance ODEs | Closed moment equations for mean and variance | `birth_death_immigration_moments/birth_death_immigration_moments.bngl` |
| SIS epidemic threshold | `S + I -> 2I`, `I -> S` with conserved `S+I` | Logistic endemic-threshold solution for `I(t)` | `sis_epidemic_threshold/sis_epidemic_threshold.bngl` |
| Quantum two-level Rabi oscillation | Resonant Bloch equations `dV/dt=-Omega*W`, `dW/dt=Omega*V` | Coherent sine/cosine population oscillations | `quantum_rabi_two_level/quantum_rabi_two_level.bngl` |
| Non-Lipschitz square-root growth | `dX/dt = a*sqrt(X)` | Square-law branch; zero initial condition is nonunique | `non_lipschitz_square_root_growth/non_lipschitz_square_root_growth.bngl` |
| Quadratic finite-time growth | `dX/dt = k*X^2` | Rational solution with finite-time pole | `quadratic_finite_time_growth/quadratic_finite_time_growth.bngl` |
| Two-species harmonic oscillator | Linear 2SHO network with constrained `k3=k1+k2`, `k5=k3+kd` | Matrix-exponential sinusoid with `theta=sqrt(k2*kd)` | `two_species_harmonic_oscillator/two_species_harmonic_oscillator.bngl` |
| Transit-compartment absorption pharmacokinetics, N=2 | `Transit1 -> Transit2 -> Absorption -> Central -> 0` with equal transit rates | Elementary specialization of the transit-compartment PK solution | `transit_absorption_pk_n2/transit_absorption_pk_n2.bngl` |
| Richards growth | `dN/dt = r*N*(1 - (N/K)^b)` | Generalized logistic power-transform solution | `richards_growth/richards_growth.bngl` |

## Conventions

- Time is measured in days unless a model README states otherwise.
- Each BNGL file encodes the modeled ODE system over `t >= 0`.
- Parameters that define pre-intervention steady states may appear even when they are
  not active in the post-intervention reaction rules.
- Verification scripts compare BioNetGen ODE output to a closed-form analytical solution.
