# `Borghans_BiophysChem1997/campaign/`

The search campaign that produced the numbers in this slug's `README.md` — the paired
`ms`-vs-`gntr` tables, the multistart tallies, the box sweeps. It is kept here, in a
subdirectory, because the slug's top level follows this collection's convention: one
`Borghans_BiophysChem1997.conf`, the model, the data, and the scoring artifacts. A search
campaign is not one of those, but it is also not disposable — without it the tallies in
the README are unreproducible assertions.

## What is here, and what is not

Committed: the **drivers** (hand-written) and the **config templates** they expand.

Not committed: everything a run emits — one config per run, the status tables, the
progress logs, the PyBNF console captures. Those are ignored by the
`Driven-campaign fan-out` block in the repo's root `.gitignore`, and they are
recoverable: a generated config comes back by re-running the driver that wrote it, and
the start-point JSONs regenerate exactly, since `make_radius_starts.py` and
`find_oscillating_starts.py` both default to `--seed 20260814` through
`np.random.default_rng`.

## Running

Every script cd's to the **job directory** (one level up) before doing anything, because
the confs name their model and `.exp` data relative to it. Run them from anywhere:

```bash
campaign/<script>.sh
```

PyBNF lives in a separate checkout with its own venv, so two variables matter. Both are
exported by `.envrc.local` (untracked; direnv sources it from `.envrc`), and both fall
back to PATH if unset:

| variable | what it is | who needs it |
|---|---|---|
| `PYBNF_BIN` | the `pybnf` entry point | anything that runs a fit |
| `PYBNF_PY` | PyBNF's **interpreter** | the scripts that `import pybnf` |

`PYBNF_PY` is not interchangeable with this repo's `python3`: `pybnf.pset` and
`pybnf.parse` pull in `roadrunner` and `distributed`, which are installed in PyBNF's
venv and not in BNGL-Models'. Under plain `python3` the top-level import succeeds and the
failure surfaces later, inside the run.

Each script's own header says which of the two it needs.
