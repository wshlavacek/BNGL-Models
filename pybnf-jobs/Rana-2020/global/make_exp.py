"""Build the .exp files for the Rana-2020 PyBNF jobs from the digitized Fig. 1 CSVs.

One script, four jobs: run it from inside any slug folder and it writes that slug's tables.

Source of the data
------------------
Rana et al. (2020) tabulate nothing. The ThT traces of Fig. 1 are recovered by
`models/amyloid_beta_competing_aggregation_pathways_rana2020/digitize_rana2020.py`
(600 dpi rasterization of PDF page 4, ink-vector colour separation, 0.5 h bins); the
committed CSVs live in that model's `reference/` directory and are read here.

Panel assignment
----------------
The published caption labels Fig. 1b "micelle addition event" and Fig. 1c "micelle removal
event", but the body text cites Fig. 1c for the micelle-addition fit ("the three
experiments", SSE 4.12, Sec. IV-B) and Fig. 1b for the micelle-removal fit (SSE 1.22,
Sec. IV-C), and the kinetics agree with the body text: every trace of panel (b) rises
without a lag, as expected when fatty acid is present from t=0, while panel (c) holds two
traces that stay flat until 3 h and 24 h. The body text is followed here.

Time base
---------
A switching experiment is expressed in the conf as a fixed-length `preequilibrate` phase
under one condition followed by the measured phase under the other, and the measured phase
restarts the clock at zero. Points before the event are therefore dropped and the remaining
times are shifted by the event time. Nothing is lost: for the addition experiments the
pre-event trace is the on-pathway experiment, and for the removal experiments it is the
fatty-acid control -- both are fit in their own right.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SLUG = HERE.name
CSV = (HERE.parents[2] / "models" / "amyloid_beta_competing_aggregation_pathways_rana2020"
       / "reference")

# slug -> {exp file stem: (panel, curve, event time in h, observable column)}
JOBS = {
    "on_pathway": {
        "on_pathway": ("fig1a", "on_pathway", 0.0, "ThT"),
    },
    "micelle_addition": {
        "fa_control": ("fig1c", "fa_control", 0.0, "ThT"),
        "addition_3h": ("fig1c", "addition_3h", 3.0, "ThT"),
        "addition_24h": ("fig1c", "addition_24h", 24.0, "ThT"),
    },
    "micelle_removal": {
        "fa_control": ("fig1b", "fa_control", 0.0, "ThT"),
        "removal_5h": ("fig1b", "removal_5h", 5.0, "ThT"),
        "removal_24h": ("fig1b", "removal_24h", 24.0, "ThT"),
    },
    # The paper's global fit is over "all the five curves: on pathway data, micelle
    # addition at 3h and 24 hour, micelle removal at 5h and 24 hour" (Sec. IV-D). The two
    # fatty-acid controls are not among them and are left to the individual jobs.
    # Table I gives the global fit two off-pathway mapping constants, k_off1 and k_off2;
    # ThT_add and ThT_rem are those two measurement models.
    "global": {
        "on_pathway": ("fig1a", "on_pathway", 0.0, "ThT_add"),
        "addition_3h": ("fig1c", "addition_3h", 3.0, "ThT_add"),
        "addition_24h": ("fig1c", "addition_24h", 24.0, "ThT_add"),
        "removal_5h": ("fig1b", "removal_5h", 5.0, "ThT_rem"),
        "removal_24h": ("fig1b", "removal_24h", 24.0, "ThT_rem"),
    },
}


def write(slug=None, outdir=None):
    slug = slug or SLUG
    outdir = Path(outdir or HERE)
    total = 0
    for stem, (panel, curve, t_event, obs) in JOBS[slug].items():
        df = pd.read_csv(CSV / f"rana2020_{panel}_digitized.csv", comment="#")
        df = df[df.curve == curve].copy()
        n_all = len(df)
        df = df[df.time_h >= t_event]
        df["time"] = (df.time_h - t_event).round(4)
        # PyBNF reads the column header from line 1 of the .exp and parses every later
        # line as data, so provenance cannot be carried in the file; it is recorded here
        # and in the job README instead.
        path = outdir / f"{stem}.exp"
        with open(path, "w") as fh:
            fh.write(f"# time {obs}\n")
            for t, v in zip(df.time, df.ThT_au):
                fh.write(f"{t:>10.4f} {v:>9.5f}\n")
        dropped = f", {n_all - len(df)} pre-event points dropped" if t_event else ""
        total += len(df)
        print(f"  {path.name:20s} {len(df):4d} points  "
              f"t = {df.time.min():.2f} .. {df.time.max():.2f} h{dropped}")
    print(f"  {slug}: {total} points total")


if __name__ == "__main__":
    write(*sys.argv[1:])
