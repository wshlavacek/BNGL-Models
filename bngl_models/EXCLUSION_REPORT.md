# Antimony → BNGL Exclusion Report

Per §7.2 of the BNGL house style, models failing the SSA feasibility gate are excluded.
Per §7.9, this report documents each excluded model with its reason category.

## Summary

- **Total Antimony models evaluated:** 117
- **Converted to BNGL:** 55
- **Excluded (SSA-infeasible):** 62

## Excluded Models

### test_models1 (12 excluded out of 29)

| Model | Reason Category |
|-------|----------------|
| m05_three_term_sum | Fractional exponents (Z^(2/3), Z^(-1/2)) on state variables; non-biochemical |
| m10_van_der_Pol | Van der Pol oscillator; continuous amplitudes, not populations |
| m13_composite_func_decomp | Exponential function exp(k*X) on state variable |
| m16_normal_dens | Statistical density function; not a population variable |
| m17_central_t | Student-t density function; statistical distribution |
| m24_scalarExpLog | Composite function exp[(ln Z)^2] on state variable |
| m25_CSTR | Exponential exp(c*Z2/(c+Z2)) on state variable |
| m26_cos_growth | Trigonometric function cos(t) in rate law |
| m27_spiral | Ratio (Z-t)/(Z+t) creates algebraic constraints |
| m28_torus | Fractional power (Z1^2+Z2^2)^(-0.5); geometric trajectory |
| m29_time_varying_beta | Time-dependent sigmoid/tanh branching logic |

### test_models2 (20 excluded out of 28)

| Model | Reason Category |
|-------|----------------|
| S1987_4C_exp_composition | exp[(ln Z)^2]; logarithmic/exponential composition |
| S1987_4D_sum_reduction | Fractional exponents (Z^(2/3), Z^(1/2)) on state variables |
| S1987_5_CSTR | Exponential of rational function of state variables |
| S1987_A2_riccati | Power-law with fractional exponent (Z^(3/2)) |
| S1987_A3_cos_growth | Trigonometric function cos(t) in rate law |
| S1987_A5_spiral | Time-dependent terms (Z-time)/(Z+time) |
| S1987_B_binary | Fractional exponents (Z1^(-1/3), Z2^(1/2)) |
| S1987_B4_torus | Radical term (Z1^2+Z2^2)^(-1/2); geometric |
| S1987_B5_rigid_body | Angular velocities that can go negative |
| S1987_C_implicit_de | Implicit algebraic constraints; DAE system |
| S1987_D1_orbit_e0.1 | Orbital mechanics; positions/velocities, (Z1^2+Z2^2)^(-3/2) |
| S1987_D2_orbit_e0.3 | Orbital mechanics; same as D1 |
| S1987_D3_orbit_e0.5 | Orbital mechanics; same as D1 |
| S1987_D4_orbit_e0.7 | Orbital mechanics; same as D1 |
| S1987_D5_orbit_e0.9 | Orbital mechanics; same as D1 |
| S1987_E1_bessel | Bessel equation with time-dependent coefficients |
| S1987_E2_van_der_pol | Van der Pol oscillator; state can go negative |
| S1987_E3_duffing | Duffing equation with sinusoidal forcing |
| S1987_E4_falling_body | Position/velocity (negative during fall) |
| S1987_E5_implicit | Radical (1+Z2^2)^(1/2) and time-dependent terms |

### test_models3 (24 excluded out of 37)

| Model | Reason Category |
|-------|----------------|
| A2013_power_system_2machine | Trigonometric (sin, cos) on state variables; angles |
| A2013_power_system_3machine | Trigonometric (sin, cos) on state variables; angles |
| A2017_gompertz_recasting | Logarithm ln(y) on state variable |
| DN2015_planetary_motion | Velocities; negative power laws r^(-3) |
| DN2015b_sinx_recasting | Trigonometric sin(x1) on state variable |
| HBF1998_morse_oscillator | Exponential exp(-x) on state; position/momentum |
| HBF1998_three_wave | Wave amplitudes (physical); negative exponents |
| KPW2024_exponential_lifting | Exponential exp(-x) on state variable |
| MS2007_fermentation_yeast | Negative exponents on state vars (X2^(-0.2344)) |
| PP2005_CSTR_arrhenius | Exponential Arrhenius; temperature state variable |
| PP2005_exponential_exact | Exponential exp(-alpha*x) on state |
| RV1990_central_F | Statistical distribution functions |
| RV1990_central_chisquared | Statistical density/cumulative distributions |
| RV1990_central_t_density | Statistical t-distribution density |
| S1993_mixed_terms | Fractional exponent sqrt on state variable |
| S1993_sum_radical | Square root of state variable |
| V1988a_sin_exp_system | sin(Y) and exp(Y^2) on state variables |
| V1988a_weibull_growth | Negative fractional power with time as state |
| V1988b_exponential_ode | Exponential exp(X) on state variable |
| V1990_endemic_exp_infection | Exponential exp(b*A) on state variable |
| V1992_blasius_equation | Fluid dynamics; physical coordinates |
| V1992_log_ode | Logarithm ln(tau) on state |
| V1993_blue_sky | Trigonometric sin(1.5*time); cubic bifurcation |
| V1993_forced_oscillator | Time-dependent sin forcing; mechanical oscillation |
| V1993_rossler_band | Chaotic system; abstract phase-space coordinates |
| Z2022_pll_converter | Trigonometric on phase angle; power electronics |

### test_models4 (5 excluded out of 20)

| Model | Reason Category |
|-------|----------------|
| De_Young1992 | Algebraic constraint (DAE); state conservation |
| Fink2000 | Smooth ReLU approximation of thresholds; time-dependent exp(-k*t) |
| Lev_Bar-Or2000 | Threshold-dependent switching; signal-dependent function |
| McMillen2002 | DAE-like auxiliary variables with manifold drift |
| Weber2018 | Time-dependent tanh sigmoids and cos forcing |

## Reason Category Summary

| Category | Count |
|----------|-------|
| Fractional/negative exponents on state variables | 10 |
| Trigonometric functions on state variables | 10 |
| Exponential/logarithmic functions on state variables | 10 |
| Non-population state variables (positions, velocities, angles) | 10 |
| Orbital/mechanical dynamics | 7 |
| Statistical distributions (not populations) | 4 |
| Algebraic constraints / DAEs | 4 |
| Time-dependent forcing / branching logic | 5 |
| Other non-mass-action kinetics | 2 |

## NFsim Exclusions

The following 12 models were converted to BNGL but NFsim (`nfr`) was excluded because
NFsim failed on them (per §1.1 feasibility rules):

| Model | NFsim Failure Reason |
|-------|---------------------|
| Bergman1989 | Negative propensity in function evaluation |
| Dreisigmeyer2008 | Function evaluation error |
| Goldbeter1996 | Function evaluation error |
| Kholodenko2000 | Composite functions (functions calling other functions) not supported |
| Lander2009 | Function evaluation error |
| MS2007_tryptophan_operon | Memory overflow (>262144 molecules created) |
| Mueller2006_RepLeaky_n3 | Composite functions not supported |
| P2011_branched_SC | Function evaluation error |
| V2005_bistable_gene | Composite functions not supported |
| Xiong2018_ContactKilling_Inhibitor | Function evaluation error |
| m07_SIR | Memory overflow (>262144 molecules created) |
| m21_Haldane_chemostat | Composite functions not supported |
