# Validation — `Cheemalavagu-2024/jak_stat_fig2b`

## Source and scope

- Paper: Cheemalavagu et al., *Cell Systems* 15:37-48.e4 (2024),
  DOI 10.1016/j.cels.2023.12.002.
- Author repository: `ncheemalavagu/STAT_models`, commit
  `28874f0801077f479526157300f70fdccb672013`.
- Fit target: Figure 2B; pSTAT1 and pSTAT3 under IL-6 1/10 ng/mL, IL-10 1/10
  ng/mL, and matched 1+1/10+10 ng/mL combinations, 0–90 min.
- Point of contact: William S. Hlavacek, hlavacek@lanl.gov.

The authors' pooled workbooks were processed exactly as stated in the STAR Methods and
implemented in their `kmeans_clustering_trajs.ipynb`: subtract each experiment's t=0,
clip negative values to zero, divide each pSTAT channel by that experiment's IL-6
10 ng/mL value at 20 min, interpolate missing common times, then average. The fitting SD
is 15% of the mean; this job uses a 0.02 floor at zero/near-zero points to avoid a zero
Gaussian variance and to reflect the plotted precision.

## Gate evidence

1. **Model fidelity — PASS.** The house-style rewrite produces the same 53 species and
   122 unidirectional reactions as the source. Under IL-6 10 ng/mL it matches the authors'
   rounded trajectory with maximum relative error `5.7e-4` (rounding-limited).
2. **Representative published fit — PASS.** Scoring all 2,067 supplied ensemble members
   against the pooled six-condition target selects member 1817 (rows 10902–10907).
   It gives χ² `74.4845` over 84 observations and median relative error `8.778%`.
   Regenerate with `python make_reproduction.py`.
3. **Tier-1 config — PASS.** `check_conf.py` reports edition 2, resolved `job_type=de`,
   bound experiment data, 48 free parameters bound by model ID, and no `__FREE` aliases.
4. **PEtab v2 — PASS.** `petab_roundtrip.py --job-type de` exports six conditions and
   experiments, lints with no errors, and imports successfully. The setup uses only
   arithmetic measurement functions and `chi_sq`; no normalization directive or BPSL
   constraint is present.
5. **Verification notebook — PASS.** The executed notebook independently integrates the
   generated 53-species network with SciPy and agrees with BioNetGen to `4.71e-7` maximum
   scaled error. All notebook assertions pass; the six-condition comparison gives χ²
   `74.4823` and median relative error `8.781%`.
6. **Bounded fit — PASS.** PyBNF v1.6.0 completed the prescribed smoke fit with
   `population_size = 6`, `max_iterations = 2`, and `refine = 0`. It stopped normally
   with objective `959.9372`; this gate checks execution, not convergence in 48 dimensions.

## Adaptations

- Source molecule names were expanded (`L1/L2` → `IL6/IL10`, `S1/S3` →
  `STAT1/STAT3`) without changing the generated network.
- The paper's proposal-dependent channel normalization was replaced by two fitted
  measurement scales. This is the minimal PEtab-v2-compatible representation and adds
  two nuisance parameters to the paper's 46 biological unknowns.
- The job uses differential evolution; the paper used three-chain pTempEst parallel
  tempering. The data, model, free biological parameters, and likelihood uncertainty are
  retained; the optimizer is intentionally a PyBNF-native global method.
