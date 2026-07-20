# Rijal-2025 — exact-SSA fits of two-state promoter noise

Two compact PyBNF edition-2 jobs derived from Rijal and Mehta's two-state promoter model and
the experimental *E. coli* mRNA-noise data of Jones et al. These jobs are the first local step
toward [PyBNF issue #472](https://github.com/lanl/PyBNF/issues/472): a small, real-world,
end-to-end stochastic fit using `experiment: ... method: ssa` through BNGsim.

> Rijal K, Mehta P. "A differentiable Gillespie algorithm for simulating chemical kinetics,
> parameter estimation, and designing synthetic biological circuits." *eLife* 2025;
> 14:RP103877. [doi:10.7554/eLife.103877](https://doi.org/10.7554/eLife.103877)
>
> Jones DL, Brewster RC, Phillips R. "Promoter architecture dictates cell-to-cell variability
> in gene expression." *Science* 2014; 346:1533–1536.
> [doi:10.1126/science.1255301](https://doi.org/10.1126/science.1255301)

## Jobs

- [`lacud5_ssa`](lacud5_ssa/): Jones Fig. 3A `lacUV5`, called `lacUD5` by Rijal;
  fits `r_over_gamma` and `gamma`; exact SSA with 200 trajectories/evaluation;
  reproduction PASS and bounded-fit PARTIAL; confidence 82/100.
- [`five_dl1_ssa`](five_dl1_ssa/): Jones Fig. 3A `5DL1`; fits `r_over_gamma` and
  `gamma`; exact SSA with 200 trajectories/evaluation; reproduction and bounded-fit PARTIAL;
  confidence 70/100.

Both are tiny finite networks: three species and four reactions. The two folders are independent,
self-contained fits so each promoter's published parameter set can be used as the nominal point.

## Model and stochastic measurement

The promoter switches between active and repressor-bound states,

```text
active --kon_R--> inactive       active --r--> active + mRNA
inactive --koff_R--> active      mRNA --gamma--> 0
```

with `koff_R = 1`, as in Rijal and Mehta. Each BNGsim SSA trajectory reports the terminal mRNA
count `m` and `m^2`. PyBNF's `smoothing` averages these over independent trajectories, after which
edition-2 measurement formulas construct

```text
Mean_mRNA = E[m]
mRNA_SD   = sqrt(E[m^2] - E[m]^2).
```

The ordinary `sos` objective is therefore Rijal Eq. 14: squared error in the experimental mean
plus squared error in the experimental mRNA standard deviation. This uses PyBNF's existing
objective and exact BNGsim SSA; no differentiability or custom optimizer is required.

## Data and the necessary reduction

Jones Fig. 3A reports mean mRNA and Fano factor. The exact numeric arrays staged in each job come
from the [authors' DGA repository][dga].
The `.exp` standard deviation is the lossless transformation
`sqrt(mean * Fano)`. The authors' notebook removes one duplicate-mean row per promoter with
`np.unique(..., return_index=True)`, leaving 17 fit rows from each 18-row source array.

Rijal and Mehta fit one latent `kon_R` for every experimental row, but those fitted on-rates were
not tabulated or released. An edition-2 parameter scan also cannot introduce a distinct fitted
parameter for every row. These jobs eliminate that latent input using the exact steady-state mean:

```text
kon_R = koff_R * (r_over_gamma / mean_target - 1).
```

Thus the measured mean indexes the scan and the optimizer estimates `r_over_gamma` and `gamma`;
the transcription rate is derived as `r = r_over_gamma * gamma`. This preserves the four-reaction
stochastic model and makes the real-data fit runnable, but it is a reduced exact-SSA adaptation,
not a byte-for-byte reproduction of the paper's transient differentiable-Gillespie fit.

## What is verified

- Both edition-2 configurations parse and bind the two fitted parameters without legacy
  `__FREE` aliases.
- Both export, lint, and import through PEtab v2 structurally. PEtab does not encode PyBNF's
  SSA-replicate smoothing semantics, so this is a schema round-trip, not a claim that a
  deterministic PEtab runner computes the same ensemble moments.
- A bounded DE smoke fit reaches a finite objective through `BngsimModel`, `method=ssa`, and
  200-trajectory smoothing in about nine seconds per job on the validation machine.
- The committed reproduction figures use 2,000 exact SSA trajectories at the published Fig. 7
  parameter values. See each `VALIDATION.md` for metrics and limitations.

The `lacUD5` spelling follows Rijal's paper and released array. Jones calls the same promoter
series `lacUV5`; the local source table retains the Rijal spelling and documents this discrepancy.

## Run

```bash
export BNGPATH=/path/to/BioNetGen

cd pybnf-jobs/Rijal-2025/lacud5_ssa
pybnf -c lacud5_ssa.conf

cd ../five_dl1_ssa
pybnf -c five_dl1_ssa.conf
```

Each job's `make_data.py` regenerates its `.exp` file from the staged source table, and
`make_reproduction.py` regenerates its exact-SSA comparison figure.

[dga]: https://github.com/Emergent-Behaviors-in-Biology/Differentiable-Gillespie-Algorithm
