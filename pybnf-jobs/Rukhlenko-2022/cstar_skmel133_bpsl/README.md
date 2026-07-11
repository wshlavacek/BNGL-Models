# cstar_skmel133_bpsl — cSTAR DPD / Waddington-landscape BPSL job (PyBNF, native-only)

The constraint-bearing sibling of [`../cstar_skmel133`](../cstar_skmel133), and the one
example in this set that exercises PyBNF's **signature capability — BPSL** (the Biological
Property Specification Language). Instead of fitting numbers, it fits the paper's
**qualitative** claims about the cell-state landscape.

> Rukhlenko OS, Halasz M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

## The idea

The paper's central construct is the **DPD** (Dynamic Phenotype Descriptor) **S** — a signed
cell-state coordinate on a **bistable, Waddington-like landscape** (Sd = death/arrest, S<0;
Sp = proliferation, S>0). The cSTAR headline for the RAF-inhibitor-resistant SKMEL-133 line
is a **threshold switch**: combined **IRS + AKT** inhibition pushes proliferating cells
across the landscape barrier into the death basin — a switch **neither single agent
achieves**. This model (`Null -> S Seq()`, with the piecewise-linear force `Seq()` and the
driving force `betaVal()`) reproduces exactly that, at the published parameters:

| condition | steady-state `Sval` | basin |
|---|---:|---|
| no drug | **+10.0** | proliferation |
| AKT inhibitor alone | +1.5 | proliferation |
| IRS inhibitor alone | +5.7 | proliferation |
| **IRS + AKT combined** | **−3.1** | **death/arrest — switched** |

Those qualitative facts are written as **BPSL `.prop` constraints on `Sval`** (see
`cstar_skmel133_bpsl_switch.png`). The paper itself fit the DPD quantitatively; expressing
the claims as BPSL constraints is this example's contribution.

## The constraints (`.prop` files)

```
baseline.prop:  Sval > 2 at 100000 weight 10                 # no drug -> stays proliferative
akt.prop:       Sval > 0 at 100000 weight 10                 # AKT alone -> no flip
                Sval < baseline.Sval at 100000 weight 5      #   ...but pushed below baseline
irs.prop:       Sval > 0 at 100000 weight 10                 # IRS alone -> no flip
                Sval < baseline.Sval at 100000 weight 5
combo.prop:     Sval < 0 at 100000 weight 20                 # IRS+AKT -> SWITCH to death basin
                Sval < akt_hi.Sval at 100000 weight 10       #   synergy: below AKT alone...
                Sval < irs_hi.Sval at 100000 weight 10       #   ...and below IRS alone
```

`baseline.Sval`, `akt_hi.Sval`, `irs_hi.Sval` are **cross-experiment (dotted-suffix)**
references — the suffix is the experiment name. `at 100000` reads the settled steady state.

## Structure — one model per condition (Miller2025 pattern)

A BPSL constraint-only experiment cannot apply a **non-free-parameter condition**:
`preequilibrate:` needs an `.exp`, and the no-preequilibration condition path only mutates
*free* parameters. So — exactly like the shipped `Miller2025_MEK_Isoforms` example — each
inhibitor **dose is baked into its own model copy**, and each experiment binds that model:

| file | dose | experiment | constraints |
|---|---|---|---|
| `cstar_skmel133_bpsl.bngl` | none (baseline) | `baseline` | `baseline.prop` |
| `cstar_skmel133_bpsl_akt.bngl` | `I_AKT_conc = 3` | `akt_hi` | `akt.prop` |
| `cstar_skmel133_bpsl_irs.bngl` | `I_IRS_conc = 4` | `irs_hi` | `irs.prop` |
| `cstar_skmel133_bpsl_combo.bngl` | `I_IRS_conc = 2, I_AKT_conc = 1.5` | `combo` | `combo.prop` |

The three dose models are generated from the baseline by **`make_models.sh`**. The four free
landscape params bind by name across all models. (As in the quantitative twin, the IRS/AKT
inhibitors act only through `1/(1+I_<x>)` with no binding, so baking the dose into
`I_<x>_conc` is exact.)

## Free parameters (4)

Landscape / driving-force connection coefficients (all feed `pmTOR`, the dominant modulated
term of the driving force), `loguniform` over ≈±½ decade around the published value:
`g_AKTmTOR` (30), `g_IRSAKT` (54.29), `g_SRCmTOR` (35), `g_CDKmTOR` (32).

## ⚠️ Native-only (not PEtab-v2-exportable)

Attaching **any** `.prop`/`.con` makes the job native-only: `export_job` raises
`NotImplementedError` (verified — BPSL has no PEtab v2 representation). Verify with
`job_type = check`, **not** the PEtab round-trip.

## Verification (all pass)

- **Tier-1** (`check_conf.py`): edition 2, `job_type=de` resolves, **4 BPSL constraint sets
  bound**, 4 free params bind by id, no `__FREE`; correctly reported native-only. **PASS.**
- **Model build** (BNG2.pl): each model generates 44 species / 56 reactions.
- **`job_type = check` at published params** (`cstar_skmel133_bpsl_check.conf`):
  **`Satisfied 8 out of 8 constraints`**, objective 0.0 — the model satisfies every
  qualitative claim (including the cross-experiment synergy comparisons).
- **Bounded bngsim `de` fit** (constraints as soft penalties): finite objective (0.0) — the
  simulate→score→propose loop runs across all four models.
- **Native-only assertion**: `pybnf.petab.export_job(...)` raises `NotImplementedError`.
- **Landscape switch** (`cstar_skmel133_bpsl_switch.png`): the DPD reaches +10 / +1.5 / +5.7
  under no-drug / AKT / IRS, and crosses to **−3.1** only under IRS+AKT.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_skmel133_bpsl
./make_models.sh                                     # (re)generate the dose models
pybnf -c cstar_skmel133_bpsl.conf                    # fit the landscape to the qualitative facts
pybnf -c cstar_skmel133_bpsl_check.conf              # or: check satisfaction at published params
```

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_skmel133_bpsl', conf='cstar_skmel133_bpsl.conf', simulator='ode',
    observables=('Sval',),
    system='cSTAR SKMEL-133 DPD / Waddington-landscape bistable switch (Rukhlenko 2022, '
           'PMC9644236); ODE; BPSL .prop constraints on Sval -> NATIVE-ONLY (not '
           'PEtab-exportable); IRS+AKT combination crosses to the death basin'),
    # BPSL: assert export_job raises NotImplementedError, and/or a `check` run reports
    # "Satisfied 8 out of 8" -- do NOT add a PEtab-lint test (it will correctly refuse).
```
