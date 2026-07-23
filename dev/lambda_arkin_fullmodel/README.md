# lambda_arkin_fullmodel — generator for the network-free phage-λ full circuit

Reconstructed generator for the two committed network-free models of the
Arkin, Ross & McAdams (1998) phage-λ lysis/lysogeny decision circuit:

- `models/lambda_switch_arkin1998/lambda_switch_arkin1998_fullcircuit.bngl` (`base`)
- `models/lambda_switch_arkin1998/lambda_switch_arkin1998_fullcircuit_exact.bngl` (`exact`)

The original generator was local-only and never committed (lost on a machine change).
This is a faithful reconstruction: it regenerates **both** committed models
**byte-for-byte** (guarded by `tests/test_lambda_generator.py`).

> The committed `.bngl` headers cite the original files `plus_gen.py`, `or_system.py`,
> `pre_pl.py`, and `plus4d_moi.py` (the assembler). This reconstruction keeps the first
> three names; the assembler is `build_fullcircuit.py` (the equivalent of the original
> `plus4d_moi.py`).

## Why a generator

The models are network-free with a state-increment ("`~PLUS` collapse") RNAP that carries
a `loc~0..2283` position component — a molecule type with ~2300 enumerated states and a
40-configuration Shea-Ackers `O_R` partition function. Editing that by hand does not scale.
The generator makes the structured, hard-to-hand-edit content **spec- and code-driven** so
future refinements are small edits, not surgery on thousands of enumerated states.

## Usage

```bash
# regenerate a model (byte-identical to the committed file)
python dev/lambda_arkin_fullmodel/build_fullcircuit.py base   > models/lambda_switch_arkin1998/lambda_switch_arkin1998_fullcircuit.bngl
python dev/lambda_arkin_fullmodel/build_fullcircuit.py exact  > models/lambda_switch_arkin1998/lambda_switch_arkin1998_fullcircuit_exact.bngl

# verify both variants still match the committed files
python dev/lambda_arkin_fullmodel/build_fullcircuit.py --check
pytest tests/test_lambda_generator.py
```

## Layout

| file | role |
|------|------|
| `build_fullcircuit.py` | assembler — splices generated content into the static prose; `base`/`exact`/`--check` |
| `plus_gen.py`   | `~PLUS` RNAP molecule type (loc `0..2283`) + per-operon initiation/elongation/landmark rules, from the `OPERONS` spec |
| `or_system.py`  | 40-config Shea-Ackers `O_R` partition function + `A_PR()`/`A_PRM()` (Table 1/2) |
| `pre_pl.py`     | `P_RE` (CII-activated) and `P_L` (Cro2/CI2-repressed) functions (Table 1) |
| `_sections.py`  | verbatim static prose (header, Table-3 parameters, R1–R5 turnover, growth/dilution) — model **data**, carried not generated |

## Extending toward the full published model (issue follow-up)

The deferred fidelity items live in the generated layer — add them here rather than by hand:

- **Minor terminators** `tR0` (in *cro*), `tR2` (after *P*), `tL2a` (before *cIII*): add a
  `term`-like landmark to the relevant `OPERONS` entry in `plus_gen.py` and a fall-off
  `k_fo_*` parameter in `_sections.py`. `tL2a` (trims *cIII*) is the highest-impact for the
  Fig-6a tails.
- **Explicit ribosome / RNase-E translation** and **per-nt elongation**: extend
  `plus_gen.py` (new landmark kinds / a per-nt stepping mode) — these restore the
  establishment-timing variance that sharpens the Fig-6a onset and tails.
- **Promoter occlusion**: a post-initiation refractory landmark in `plus_gen.py`.

After any change, re-run `--check` (it will report a diff against the committed files) and
update the committed `.bngl` + verification with `run_fullcircuit.py`.
