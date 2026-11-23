#!/usr/bin/env python
"""Check reduced_onoff.bngl against the authors' own on-off model file, with BNG2.pl only.

The authors ship their reduced model at the Table 1 fitted values as S1 Code of the article
(fitted_reduced_model_on_off_protocol.bngl), vendored byte-identical here under reference/ and
pinned by URL and sha256 in upstream.json. Its actions block runs the on-off protocol of Table 2
as three CVODE phases with TR set by setParameter:

    TR = 0 on [0, 3600) s   ->   TR = 1 on [3600, 10800) s   ->   TR = 0 on [10800, 46800] s

for WT cells, then again for A20KO cells (c_deg = 0). reduced_onoff.bngl is the same reduced model
written for an edition-2 conf: it has no actions block, it starts at TNF onset, and TR is a
function of an inert clock, TRfunc() = if(simtime < t_tnf_off, 1, 0) with t_tnf_off = 7200 s, so
the 2 h TNF step is a discontinuity inside one integration rather than a phase boundary.

This script runs both files through BNG2.pl and asserts, against tolerances fixed before the run:

  (a) the authors' TR = 0 phase leaves the seed state (IKKn = 1, everything else 0) unchanged:
      every species stays within 1e-12 of its seed over [0, 3600] s, in both genotypes. That is the
      ground for reduced_onoff.bngl starting at TNF onset with no pre-equilibration;
  (b) the generated networks are the same: the authors' 6 species and 14 reactions, plus only
      the Clock() species and its 0 -> Clock() reaction in reduced_onoff.bngl;
  (c) reduced_onoff.bngl at t equals the authors' file at t + 3600 s, for all six species in WT
      and A20KO, on the 145 points of the authors' 300 s output grid over the 12 h after onset:
      max |difference| <= 1e-5 x the observable's peak over the protocol (A20 in A20KO is
      identically zero in both, so it must agree to 1e-12 absolute).

Why 1e-5: the two files encode the same six ODEs with the same 13 parameter values, so any
difference comes from the integrator: the authors' file asks CVODE for rtol = atol = 1e-16, the
exact file runs at BNG2.pl's defaults (1e-8), and the TNF-off step is crossed by CVODE's error
control rather than by a restart. 1e-5 is a thousand times the looser of those tolerances and two
orders below the smallest change a parameter edit could make, since the Table 1 values carry
three significant digits (a one-unit change in the third digit is at least 1e-3 relative).

Nothing here uses bngsim. Run from the repository root or this folder, with BNGPATH set to the
folder holding BNG2.pl:

    BNGPATH=... uv run python pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_onoff/check_author_file.py

It writes reference/author_file_check.json and reduced_onoff_author_file_check.png, and exits 0
only if every check passes.
"""

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
UPSTREAM = os.path.join(HERE, "upstream.json")
EXACT = os.path.join(HERE, "reduced_onoff.bngl")
AUTHOR_REL = "reference/upstream_fitted_reduced_model_on_off_protocol.bngl"
OUT_JSON = os.path.join(HERE, "reference", "author_file_check.json")
OUT_PNG = os.path.join(HERE, "reduced_onoff_author_file_check.png")

SPECIES = ["IKK_n", "IKK_a", "NFkB_n", "A20", "IkBa", "tIkBa"]
T_ONSET = 3600.0  # TNF onset in the authors' clock
T_SPAN = 43200.0  # 12 h after onset: 2 h on, 10 h off
DT = 300.0  # the authors' output spacing after onset
TOL_PEAK = 1e-5  # (c), fraction of the observable's peak
TOL_ZERO = 1e-12  # (a), and (c) for an observable that is zero throughout


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def bng2pl():
    bngpath = os.environ.get("BNGPATH")
    if not bngpath or not os.path.isfile(os.path.join(bngpath, "BNG2.pl")):
        sys.exit("set BNGPATH to the folder that holds BNG2.pl")
    return os.path.join(bngpath, "BNG2.pl")


def run_bng(bngl_text, name, workdir):
    path = os.path.join(workdir, name + ".bngl")
    with open(path, "w") as fh:
        fh.write(bngl_text)
    r = subprocess.run(
        ["perl", bng2pl(), path], cwd=workdir, capture_output=True, text=True, check=False
    )
    if r.returncode != 0:
        sys.exit(f"BNG2.pl failed on {name}:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    return r.stdout


def read_gdat(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, comments="#", ndmin=2)
    return {h: data[:, i] for i, h in enumerate(header)}


def read_net(path):
    """Species strings and reactions (reactant names, product names, rate law) of a .net."""
    text = open(path).read()

    def block(name):
        m = re.search(rf"^begin {name}\n(.*?)^end {name}", text, re.S | re.M)
        return [ln.split("#")[0].split() for ln in m.group(1).splitlines() if ln.split("#")[0]]

    species = {int(c[0]): c[1] for c in block("species")}
    reactions = []
    for c in block("reactions"):
        reac = tuple(sorted(species[int(i)] for i in c[1].split(",") if i != "0"))
        prod = tuple(sorted(species[int(i)] for i in c[2].split(",") if i != "0"))
        reactions.append((reac, prod, c[3]))
    return set(species.values()), reactions


def exact_with_protocol(a20ko):
    """reduced_onoff.bngl, unedited, followed by one ODE simulation over the 12 h on-off span."""
    src = open(EXACT).read()
    knockout = 'setParameter("c_deg",0)\n' if a20ko else ""
    n = int(round(T_SPAN / DT))
    return (
        src
        + "\nbegin actions\ngenerate_network({overwrite=>1})\n"
        + knockout
        + f'simulate({{method=>"ode",t_start=>0,t_end=>{T_SPAN:g},n_steps=>{n}}})\n'
        + "end actions\n"
    )


def main():
    pin = json.load(open(UPSTREAM))["files"][AUTHOR_REL]
    author = os.path.join(HERE, AUTHOR_REL)
    got = sha256(author)
    if got != pin["sha256"]:
        sys.exit(f"{AUTHOR_REL}: sha256 {got} does not match upstream.json {pin['sha256']}")

    with tempfile.TemporaryDirectory() as wd:
        # The authors' file, verbatim, with its own actions block.
        shutil.copy(author, os.path.join(wd, "authors.bngl"))
        r = subprocess.run(
            ["perl", bng2pl(), "authors.bngl"], cwd=wd, capture_output=True, text=True, check=False
        )
        if r.returncode != 0:
            sys.exit(f"BNG2.pl failed on the authors' file:\n{r.stdout[-3000:]}")
        auth = {
            "WT": read_gdat(os.path.join(wd, "authors_on_off.gdat")),
            "A20KO": read_gdat(os.path.join(wd, "authors_on_off_A20KO.gdat")),
        }
        net_auth = read_net(os.path.join(wd, "authors.net"))
        exact = {}
        for geno in ("WT", "A20KO"):
            name = f"exact_{geno}"
            run_bng(exact_with_protocol(geno == "A20KO"), name, wd)
            exact[geno] = read_gdat(os.path.join(wd, name + ".gdat"))
        net_exact = read_net(os.path.join(wd, "exact_WT.net"))

    report = {
        "exact_file": "reduced_onoff.bngl",
        "exact_file_sha256": sha256(EXACT),
        "author_file": AUTHOR_REL,
        "author_file_sha256": got,
        "engine": "BioNetGen BNG2.pl (run_network, CVODE); no bngsim",
        "tolerances": {
            "pre_stimulus_constancy_abs": TOL_ZERO,
            "trajectory_max_abs_over_peak": TOL_PEAK,
            "identically_zero_observable_abs": TOL_ZERO,
        },
        "checks": {},
    }
    ok = True

    # (a) the authors' TR = 0 phase does not move the seed state.
    pre = {}
    for geno, g in auth.items():
        m = g["time"] <= T_ONSET + 1e-9
        dev = max(float(np.max(np.abs(g[s][m] - g[s][0]))) for s in SPECIES)
        seed_ok = g["IKK_n"][0] == 1.0 and all(g[s][0] == 0.0 for s in SPECIES if s != "IKK_n")
        pre[geno] = {
            "n_points": int(m.sum()),
            "max_abs_change": dev,
            "seed_is_resting_state": seed_ok,
        }
        ok &= dev <= TOL_ZERO and seed_ok
    report["checks"]["a_pre_stimulus_phase_is_steady"] = pre

    # (b) networks: authors' species and reactions, plus only the clock.
    sp_a, rx_a = net_auth
    sp_e, rx_e = net_exact
    extra_sp = sorted(sp_e - sp_a)
    missing_sp = sorted(sp_a - sp_e)
    topo = lambda rx: sorted((r, p) for r, p, _ in rx)  # noqa: E731
    clock_rx = [(r, p) for r, p in topo(rx_e) if "Clock()" in r + p]
    same_rx = topo(rx_a) == [x for x in topo(rx_e) if x not in clock_rx]
    net_ok = (
        extra_sp == ["Clock()"] and not missing_sp and clock_rx == [((), ("Clock()",))] and same_rx
    )
    report["checks"]["b_network"] = {
        "authors_species": len(sp_a),
        "authors_reactions": len(rx_a),
        "exact_species": len(sp_e),
        "exact_reactions": len(rx_e),
        "species_only_in_exact": extra_sp,
        "species_only_in_authors": missing_sp,
        "same_reactions_apart_from_clock": same_rx,
        "pass": net_ok,
    }
    ok &= net_ok

    # (c) trajectories after onset: exact(t) against authors(t + 3600).
    traj = {}
    grid = np.arange(0.0, T_SPAN + DT / 2, DT)
    for geno in ("WT", "A20KO"):
        a, e = auth[geno], exact[geno]
        ia = [int(np.argmin(np.abs(a["time"] - (T_ONSET + t)))) for t in grid]
        ie = [int(np.argmin(np.abs(e["time"] - t))) for t in grid]
        if not (
            np.allclose(a["time"][ia], T_ONSET + grid, atol=1e-6)
            and np.allclose(e["time"][ie], grid, atol=1e-6)
        ):
            sys.exit(f"{geno}: output grids do not line up with the 300 s comparison grid")
        rows = {}
        for s in SPECIES:
            ya, ye = a[s][ia], e[s][ie]
            peak = float(np.max(np.abs(ya)))
            diff = float(np.max(np.abs(ye - ya)))
            if peak == 0.0:
                metric, tol = diff, TOL_ZERO
            else:
                metric, tol = diff / peak, TOL_PEAK
            big = np.abs(ya) > 1e-3 * peak if peak else np.zeros_like(ya, bool)
            rel = float(np.max(np.abs(ye[big] - ya[big]) / np.abs(ya[big]))) if big.any() else 0.0
            rows[s] = {
                "peak_authors": peak,
                "max_abs_diff": diff,
                "metric": "max_abs_diff/peak" if peak else "max_abs_diff (zero observable)",
                "value": metric,
                "tol": tol,
                "pass": metric <= tol,
                "max_rel_diff_where_above_1e-3_of_peak": rel,
            }
            ok &= metric <= tol
        traj[geno] = {"n_points": len(grid), "observables": rows}
    report["checks"]["c_trajectories_after_onset"] = traj
    report["verdict"] = "PASS" if ok else "FAIL"

    with open(OUT_JSON, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")

    print("reduced_onoff.bngl vs the authors' fitted_reduced_model_on_off_protocol.bngl (BNG2.pl)")
    for geno, v in pre.items():
        n, dev = v["n_points"], v["max_abs_change"]
        print(f"  (a) {geno:5s} TR=0 phase, {n} points: max change {dev:.1e}")
    nb = report["checks"]["b_network"]
    print(
        f"  (b) network: authors {nb['authors_species']} species / {nb['authors_reactions']} rxns,"
        f" exact {nb['exact_species']} / {nb['exact_reactions']} (+ {extra_sp}):"
        f" {'PASS' if net_ok else 'FAIL'}"
    )
    for geno, v in traj.items():
        for s, row in v["observables"].items():
            print(
                f"  (c) {geno:5s} {s:7s} {row['metric']:31s} {row['value']:.2e}"
                f" (tol {row['tol']:.0e}) {'PASS' if row['pass'] else 'FAIL'}"
            )
    print(f"verdict: {report['verdict']}  -> {os.path.relpath(OUT_JSON, HERE)}")

    plot(auth, exact, grid, traj, report["verdict"])
    return 0 if ok else 1


def plot(auth, exact, grid, traj, verdict):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(SPECIES, plt.rcParams["axes.prop_cycle"].by_key()["color"], strict=False))
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    for col, geno in enumerate(("WT", "A20KO")):
        a, e = auth[geno], exact[geno]
        ax, rx = axes[0, col], axes[1, col]
        m = a["time"] >= T_ONSET - 1e-9
        for s in SPECIES:
            ax.plot(e["time"] / 60, e[s], "-", color=colors[s], lw=1.5, label=f"{s} (exact file)")
            ax.plot(
                (a["time"][m][::4] - T_ONSET) / 60,
                a[s][m][::4],
                "o",
                mfc="none",
                mec=colors[s],
                ms=5,
            )
            ia = [int(np.argmin(np.abs(a["time"] - (T_ONSET + t)))) for t in grid]
            ie = [int(np.argmin(np.abs(e["time"] - t))) for t in grid]
            peak = np.max(np.abs(a[s][ia]))
            if peak > 0:
                r = np.abs(e[s][ie] - a[s][ia]) / peak
                rx.semilogy(grid / 60, np.maximum(r, 1e-16), "-", color=colors[s], lw=1.2, label=s)
        ax.axvline(120, color="0.6", ls=":", lw=1)
        ax.set_title(f"{geno}: lines = reduced_onoff.bngl, open markers = authors' file")
        ax.set_ylabel("non-dimensional amount")
        rx.axhline(TOL_PEAK, color="k", ls="--", lw=1, label=f"tolerance {TOL_PEAK:.0e}")
        rx.axvline(120, color="0.6", ls=":", lw=1)
        worst = max(v["value"] for v in traj[geno]["observables"].values())
        rx.set_title(f"{geno}: |difference| / peak (worst {worst:.1e})")
        rx.set_xlabel("time after TNF onset (min); dotted: TNF off at 120 min")
        rx.set_ylim(1e-16, 1e-3)
    axes[0, 0].legend(fontsize=7, frameon=False)
    axes[1, 0].legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle(
        "Jaruszewicz-Blonska et al. 2023, on-off protocol: reduced_onoff.bngl vs the authors'"
        " S1 Code file, both through BNG2.pl\n"
        f"pre-stimulus phase steady, networks equal bar the clock, 6 species x 2 genotypes x 145"
        f" points within {TOL_PEAK:.0e} of peak: {verdict}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=110)
    print("wrote", os.path.relpath(OUT_PNG, HERE))


if __name__ == "__main__":
    sys.exit(main())
