#!/usr/bin/env python
"""Compute the six model-predicted dose responses of Fig. 3 in Imoto et al. (2026).

Each panel of Fig. 3 plots steady-state ppERK against inhibitor concentration in
units of that inhibitor's Kd, at five KSR1 abundances, normalized to the no-drug
steady state at the same KSR1 abundance. Reproducing all six panels needs ~725
steady states of a 2709-species network, which is why this is a driver script
rather than a notebook cell; `verify_imoto2026.ipynb` reads the CSV it writes.

Conditions come from Tables S2 and S3 of Imoto et al. (2026):

    cell line     RAS-GTP 300 nM (PSN1, KRAS-G12R/WT) or 35 nM (MCF7, KRAS WT)
    Type II RAFi  f 0.005, g1 0.044, g2 4.346, g3 1, fK 0.025
    Type I-1/2    f 0.01,  g1 0.429, g2 102.95, g3 1, fK 50
    Cobimetinib   fM 1 (see below), fKR 1e-4

Table S3 lists fM = 0.1 for Cobimetinib, but two published artifacts say fM = 1:
the deposited SBML stores 1.0, and both MEKi panels of Fig. 3 reproduce at 1.0
(median |delta| 0.003 and 0.007 in normalized ppERK) and not at 0.1 (0.114 and
0.010, with maxima of 0.185 and 0.119). The figure is what this campaign exists
to reproduce, so fM = 1 is used here and `--fm 0.1` regenerates the Table S3
counterfactual that the verification notebook reports alongside it.

All other parameters keep the values in the `.bngl`, including fKR = 1e-4, and
only one inhibitor is present at a time. Doses enter as the clamped concentration
of the drug species, so the dose ladder needs no network change; each dose starts
from the previous dose's steady state, which is what makes the campaign feasible.

`--fm` re-runs the two MEKi panels at another value of fM; `--fm 0.1` is the
Table S3 counterfactual the notebook reports.

Usage
-----
    python run_imoto2026.py [--jobs N] [--net reference/<model>.net] [--fm X]

Writes `reference/imoto2026_fig3_model_doseresponse.csv`, or
`..._fM<X>.csv` covering only the two MEKi panels when `--fm` is given.
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from independent_imoto2026 import Network

HERE = Path(__file__).resolve().parent
DEFAULT_NET = HERE / "reference/ksr1_raf_mek_inhibitor_response_imoto2026.net"
OUT = HERE / "reference/imoto2026_fig3_model_doseresponse.csv"

KSR_LEVELS = (0.0, 15.0, 30.0, 50.0, 100.0)
ATOL = RTOL = 1e-7
T_END = 5.0e3  # the network is at steady state well before this (see the notebook)

RAFI_DOSES = np.unique(np.concatenate([
    np.linspace(0.0, 2.0, 11),
    np.array([2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0,
              45.0, 50.0]),
]))
COBI_DOSES_10 = np.unique(np.concatenate([
    np.linspace(0.0, 1.0, 6),
    np.array([1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]),
]))
COBI_DOSES_20 = np.unique(np.concatenate([
    np.linspace(0.0, 1.0, 6),
    np.array([1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]),
]))

TYPE_II = dict(fa=0.005, g1a=0.044, g2a=4.346, g3a=1.0, fKa=0.025)
TYPE_I_HALF = dict(fa=0.01, g1a=0.429, g2a=102.95, g3a=1.0, fKa=50.0)
# fM = 1 rather than Table S3's 0.1; see the module docstring.
COBIMETINIB = dict(fM=1.0)

PANELS = {
    "A": dict(cell="MCF7", ras=35.0, drug="type_II_RAFi", factors=TYPE_II,
              species="RAFi1", doses=RAFI_DOSES),
    "B": dict(cell="PSN1", ras=300.0, drug="type_II_RAFi", factors=TYPE_II,
              species="RAFi1", doses=RAFI_DOSES),
    "C": dict(cell="MCF7", ras=35.0, drug="type_I_half_RAFi", factors=TYPE_I_HALF,
              species="RAFi1", doses=RAFI_DOSES),
    "D": dict(cell="PSN1", ras=300.0, drug="type_I_half_RAFi", factors=TYPE_I_HALF,
              species="RAFi1", doses=RAFI_DOSES),
    "E": dict(cell="MCF7", ras=35.0, drug="cobimetinib", factors=COBIMETINIB,
              species="Cobimetinib", doses=COBI_DOSES_10),
    "F": dict(cell="PSN1", ras=300.0, drug="cobimetinib", factors=COBIMETINIB,
              species="Cobimetinib", doses=COBI_DOSES_20),
}

_NET: Network | None = None


def _load(path):
    global _NET
    if _NET is None:
        _NET = Network(path)
    return _NET


def drug_index(net, molecule):
    """Index of the clamped free-drug species.

    The species must be the lone molecule of its complex: `Cobimetinib(mek)`, not
    the many complexes whose alphabetically first molecule is Cobimetinib.
    """
    hits = [i for i, n in enumerate(net.names)
            if n.split("::")[-1].lstrip("$").startswith(molecule + "(") and "." not in n]
    if len(hits) != 1:
        raise RuntimeError(f"{molecule}: {len(hits)} candidate species")
    return hits[0]


def dose_response(job):
    """One (panel, KSR1) curve: steady-state ppERK across the dose ladder."""
    # The thermodynamic factors travel in the job rather than being read from
    # module state: workers are spawned, not forked, so they re-import this
    # module and would not see a parent-process override.
    net_path, key, ksr, factors = job
    panel = PANELS[key]
    net = _load(net_path)
    net.set_param(RAS_0=panel["ras"], KSR_0=ksr, **factors)
    j = drug_index(net, panel["species"])
    # Kd of the drug being titrated: doses are reported in units of it.
    kd = net.env["Kd_MEK_Cobimetinib" if panel["species"] == "Cobimetinib"
                 else "Kd_BRAF_RAFi1"]

    rows, x, base = [], None, None
    t0 = time.time()
    for dose in panel["doses"]:
        start = net.x0.copy() if x is None else x.copy()
        start[j] = dose * kd
        x = net.steady_state(x0=start, t_end=T_END, atol=ATOL, rtol=RTOL)
        pperk = float(net.observe("Obs_Tot_ppERK", x)[0])
        if base is None:
            base = pperk
        rows.append((key, panel["cell"], panel["drug"], ksr, float(dose),
                     dose * kd, pperk, pperk / base))
    print(f"  panel {key} KSR={ksr:5.0f} nM  {len(rows):3d} doses  "
          f"{time.time() - t0:6.0f} s  ppERK(0) = {base:8.2f} nM", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--net", type=Path, default=DEFAULT_NET)
    ap.add_argument("--fm", type=float, default=None,
                    help="override the Cobimetinib fM factor and run only panels E and F")
    args = ap.parse_args()

    out_path = OUT
    keys = list(PANELS)
    factors = {k: dict(p["factors"]) for k, p in PANELS.items()}
    if args.fm is not None:
        keys = ["E", "F"]
        for k in keys:
            factors[k]["fM"] = args.fm
        out_path = OUT.with_name(OUT.stem + f"_fM{args.fm:g}.csv")
    jobs = [(str(args.net), key, ksr, factors[key]) for key in keys for ksr in KSR_LEVELS]
    print(f"{len(jobs)} dose-response curves on {args.jobs} workers", flush=True)
    t0 = time.time()
    with Pool(args.jobs) as pool:
        out = pool.map(dose_response, jobs)
    print(f"campaign finished in {(time.time() - t0) / 60:.1f} min")

    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# Steady-state ppERK dose responses of the curated model, computed to"])
        w.writerow(["# reproduce the six model panels of Fig. 3 in Imoto et al. (2026)."])
        w.writerow([f"# Written by run_imoto2026.py; integrator BDF, atol=rtol={ATOL:g},"])
        w.writerow([f"# t_end={T_END:g} s, continuation along the dose ladder."])
        fm = args.fm if args.fm is not None else COBIMETINIB["fM"]
        w.writerow([f"# Cobimetinib fM = {fm:g}."])
        w.writerow(["panel", "cell", "drug", "KSR_nM", "dose_over_Kd", "dose_nM",
                    "ppERK_nM", "ppERK_norm"])
        for rows in out:
            for r in rows:
                w.writerow([r[0], r[1], r[2], f"{r[3]:g}", f"{r[4]:g}", f"{r[5]:g}",
                            f"{r[6]:.6f}", f"{r[7]:.6f}"])
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
