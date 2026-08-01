#!/usr/bin/env python
"""Reproduction figure for the Rohrs-2018 cd3zeta_competitive_inhibition job.

Simulates cd3zeta_competitive_inhibition.bngl over the nine fitted conditions -- wild-type
CD3zeta on 10%, 0% and 45% POPS liposomes and the six tyrosine-to-phenylalanine ITAM point
mutants -- and overlays the six site-specific percent-phosphorylation curves on the fit
data (the *.exp files, digitized from Fig. S5 of Rohrs et al. 2018).

By default it simulates at the .bngl nominals, which are the authors' reported best fit
(Data S4, mmc5.pdf). Pass --params-file <PyBNF Results/sorted_params.txt> to overlay a
fitted parameter set instead; the script then prints both parameter sets side by side and
the sum of squared error of each, which is the quantity the paper scores with (Methods,
"Comparison of model structures"; 3.47e4 reported over all ten data sets, of which this job
fits nine -- Fig. 1C is raster and was not digitized).

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy.
Usage: BNGPATH=... python make_reproduction.py [--params-file output/Results/sorted_params.txt]
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
MODEL = HERE / "cd3zeta_competitive_inhibition.bngl"
SITES = ["A1", "A2", "B1", "B2", "C1", "C2"]
COLORS = {"A1": "#0072BD", "A2": "#D95319", "B1": "#EDB120",
          "B2": "#7E2F8E", "C1": "#77AC30", "C2": "#4DBEEE"}
FITTED = ["LCK_T", "KmA1", "KmA2", "KmB2", "KmC1", "KmC2", "Xi"]

# (exp file, title, the live_<site> gate this condition switches off)
CONDITIONS = [
    ("wt_10pops_rep2", "wild type, 10% POPS", None),
    ("wt_0pops", "wild type, 0% POPS", None),
    ("wt_45pops", "wild type, 45% POPS", None),
    ("mut_A1", "A1F mutant", "A1"),
    ("mut_A2", "A2F mutant", "A2"),
    ("mut_B1", "B1F mutant", "B1"),
    ("mut_B2", "B2F mutant", "B2"),
    ("mut_C1", "C1F mutant", "C1"),
    ("mut_C2", "C2F mutant", "C2"),
]
SAMPLE_TIMES = np.logspace(-1, 3, 121)


def read_exp(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path)
    return header, np.atleast_2d(data)


def read_gdat(path):
    with open(path) as fh:
        cols = fh.readline().lstrip("#").split()
    return cols, np.atleast_2d(np.loadtxt(path))


def simulate(params, dead_site, workdir):
    """Run the model once at `params`, with `dead_site`'s gate switched off."""
    text = MODEL.read_text()
    for name, value in params.items():
        text = re.sub(rf"^(\s*{name}\s+)\S+(\s*#.*)$", rf"\g<1>{value!r}\g<2>",
                      text, flags=re.MULTILINE)
    if dead_site:
        text = re.sub(rf"^(\s*live_{dead_site}\s+)\S+(\s*#.*)$", r"\g<1>0\g<2>",
                      text, flags=re.MULTILINE)
    times = ",".join(f"{t:.6g}" for t in SAMPLE_TIMES)
    text += ("\nbegin actions\n  generate_network({overwrite=>1})\n"
             '  simulate({method=>"ode",suffix=>"ode",t_start=>0,print_functions=>1,'
             f"sample_times=>[{times}]}})\n"
             "end actions\n")
    run = Path(workdir) / "run.bngl"
    run.write_text(text)
    r = subprocess.run([str(Path(os.environ["BNGPATH"]) / "BNG2.pl"), run.name],
                       cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0 or "ABORT" in r.stdout:
        sys.exit(r.stdout[-3000:])
    cols, data = read_gdat(Path(workdir) / "run_ode.gdat")
    keep = data[:, cols.index("time")] > 0        # drop t = 0; the axis is logarithmic
    out = {"time": data[keep, cols.index("time")]}
    for s in SITES:
        out[s] = data[keep, cols.index(f"pct_{s}")]
    return out


def nominal_params():
    out = {}
    for line in MODEL.read_text().splitlines():
        m = re.match(r"\s*(\w+)\s+([0-9.eE+-]+)\s*#", line)
        if m and m.group(1) in FITTED:
            out[m.group(1)] = float(m.group(2))
    return out


def read_best_fit(path):
    """Read a fitted parameter set.

    Accepts PyBNF's ``Results/sorted_params.txt`` (a tab-separated table whose
    header names the parameters and whose first data row is the best fit) or a
    plain ``<name> <value>`` file.
    """
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    header = lines[0].lstrip("#").split()
    if set(FITTED) <= set(header):
        row = lines[1].split()
        row = row[-len(header):]                  # the leading rank column may be blank
        return {name: float(row[header.index(name)]) for name in FITTED}
    out = {}
    for line in lines:
        parts = line.split()
        if len(parts) == 2 and parts[0] in FITTED:
            out[parts[0]] = float(parts[1])
    missing = set(FITTED) - set(out)
    if missing:
        sys.exit(f"{path}: could not find values for {sorted(missing)}")
    return out


def sse(params, workdir):
    total, n = 0.0, 0
    for stem, _, dead in CONDITIONS:
        header, data = read_exp(HERE / f"{stem}.exp")
        sim = simulate(params, dead, workdir)
        for j, col in enumerate(header[1:], start=1):
            site = col.replace("pct_", "").replace("()", "")
            pred = np.interp(np.log10(data[:, 0]), np.log10(sim["time"]), sim[site])
            resid = pred - data[:, j]
            ok = np.isfinite(resid)
            total += float((resid[ok] ** 2).sum())
            n += int(ok.sum())
    return total, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params-file",
                    help="PyBNF Results/sorted_params.txt to overlay instead of nominals")
    ap.add_argument("--out", default=str(HERE / "cd3zeta_competitive_inhibition_reproduction.png"))
    args = ap.parse_args()

    if "BNGPATH" not in os.environ:
        sys.exit("set BNGPATH to the folder containing BNG2.pl")

    published = nominal_params()
    fitted = read_best_fit(args.params_file) if args.params_file else None
    workdir = tempfile.mkdtemp(prefix="rohrs_repro_")
    try:
        print(f"{'parameter':10s} {'published (Data S4)':>20s}"
              + (f" {'PyBNF fit':>14s} {'ratio':>8s}" if fitted else ""))
        for name in FITTED:
            line = f"{name:10s} {published[name]:20.4g}"
            if fitted:
                line += f" {fitted[name]:14.4g} {fitted[name] / published[name]:8.2f}"
            print(line)

        sse_pub, n = sse(published, workdir)
        print(f"\nSSE at the published parameters: {sse_pub:.4g} over {n} points "
              f"(rms {np.sqrt(sse_pub / n):.2f} percentage points)")
        if fitted:
            sse_fit, _ = sse(fitted, workdir)
            print(f"SSE at the PyBNF best fit:       {sse_fit:.4g} over {n} points "
                  f"(rms {np.sqrt(sse_fit / n):.2f} percentage points)")
        print("Rohrs et al. report SSE 3.47e4 for this mechanism over all ten data sets; "
              "this job fits nine of them.")

        fig, axes = plt.subplots(3, 3, figsize=(13.5, 10.5), sharex=True, sharey=True)
        for ax, (stem, title, dead) in zip(axes.ravel(), CONDITIONS, strict=True):
            header, data = read_exp(HERE / f"{stem}.exp")
            sim_pub = simulate(published, dead, workdir)
            sim_fit = simulate(fitted, dead, workdir) if fitted else None
            for j, col in enumerate(header[1:], start=1):
                site = col.replace("pct_", "").replace("()", "")
                ax.plot(sim_pub["time"], sim_pub[site], "-", lw=1.4, color=COLORS[site],
                        label=site)
                if sim_fit is not None:
                    ax.plot(sim_fit["time"], sim_fit[site], "--", lw=1.1,
                            color=COLORS[site], alpha=0.9)
                ok = np.isfinite(data[:, j])
                ax.plot(data[ok, 0], data[ok, j], "o", ms=4.5, mfc="none",
                        mec=COLORS[site])
            ax.set_xscale("log")
            ax.set_xlim(0.08, 1200)
            ax.set_ylim(-4, 104)
            ax.set_title(title, fontsize=10)
        for ax in axes[-1]:
            ax.set_xlabel("time (min)")
        for ax in axes[:, 0]:
            ax.set_ylabel("% phosphorylation")
        axes[0, 0].legend(ncol=3, fontsize=8, frameon=False, loc="upper left")
        subtitle = ("solid = published parameters (Data S4), open circles = Fig. S5 data"
                    + (", dashed = PyBNF best fit" if fitted else ""))
        fig.suptitle("Rohrs et al. (2018) CD3ζ competitive inhibition fit\n" + subtitle,
                     fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(args.out, dpi=140)
        print(f"\nwrote {args.out}")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
