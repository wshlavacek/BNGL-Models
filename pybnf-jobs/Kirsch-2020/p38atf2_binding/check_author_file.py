#!/usr/bin/env python
"""Check p38atf2_binding.bngl against the authors' own model file, with BNG2.pl only.

The authors' model is Supplementary Software 1 of the article (Bionetgen_JNK-p38-ATF2_model.bngl),
vendored byte-identical under reference/ and pinned by URL and sha256 in upstream.json. Its actions
block is the published stimulation protocol: 20000 s of equilibration with the MAPK activation
constants set by setParameter to their basal values (kstim7 = keq7, kstim6 = keq6), then 2400 s of
anisomycin stimulation at kstim7 and kstim6. p38atf2_binding.bngl is the same model written for an
edition-2 conf: no actions block, the basal and stimulated rates as parameters, and a gate Stim
(0 basal, 1 stimulated) selecting between them through two functions, k7_act() and k6_act().

Asserted, against tolerances fixed before the run:

  (a) the generated networks are the same: 48 species and 152 reactions with the same reactants
      and products, and the same rate law on every reaction except the two activation steps, where
      the authors' kstim7 / kstim6 is k7_act() / k6_act() here;
  (b) the exact file driven through the same protocol with Stim (0 for -20000..0 s, 1 from 0 s)
      reproduces all four observables of the authors' run on its 2401-point grid (0-2400 s, 1 s
      apart), the equilibrated state at t = 0 included: max |difference| <= 1e-6 x the observable's
      peak. Both runs use the same network and BNG2.pl's default CVODE settings, so the only source
      of difference is rounding in keq + Stim*(kstim - keq); 1e-6 is a hundred times CVODE's rtol
      and three orders below what any edit of a published value (three significant digits) moves.

Reported, NOT asserted: the exact file against the Fig. 4a NanoBit data (Source Data, sheet
Figure_4) at the published point, under the same stimulated protocol, as treated/untreated =
p38ATF2all(t) / p38ATF2all(0):

  - WT control, the curve this job fits (p38atf2_binding.exp);
  - WT + JNK-IN-8 (reference/fig4a_wt_jnkin8.exp), simulated with k1 = k2 = 0 as the Methods say.

The Methods state that kstim6, dp2 and dp4 were fit with PyBioNetFit to BOTH curves, so the
single-curve comparison is not the authors' objective; the joint 1/2 chi^2 is printed beside it. The
data's error bars are propagated standard errors of n = 3 (Source Data cell C27), and at that
precision the published point does not fit either curve: see the README.

Nothing here uses bngsim. Run from the repository root or this folder, with BNGPATH set:

    BNGPATH=... uv run python pybnf-jobs/Kirsch-2020/p38atf2_binding/check_author_file.py

It writes reference/author_file_check.json and p38atf2_binding_author_file_check.png, and exits 0
only if (a) and (b) pass.
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
EXACT = os.path.join(HERE, "p38atf2_binding.bngl")
AUTHOR_REL = "reference/upstream_Bionetgen_JNK-p38-ATF2_model.bngl"
DATA = {
    "WT control": os.path.join(HERE, "p38atf2_binding.exp"),
    "WT + JNK-IN-8": os.path.join(HERE, "reference", "fig4a_wt_jnkin8.exp"),
}
OUT_JSON = os.path.join(HERE, "reference", "author_file_check.json")
OUT_PNG = os.path.join(HERE, "p38atf2_binding_author_file_check.png")

OBS = ["pT69pT71", "p38ATF2all", "ppp38", "ppJNK"]
GATED = {"k7_act": "kstim7", "k6_act": "kstim6"}  # exact-file rate law -> authors' rate law
TOL_PEAK = 1e-6


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def bng2pl():
    bngpath = os.environ.get("BNGPATH")
    if not bngpath or not os.path.isfile(os.path.join(bngpath, "BNG2.pl")):
        sys.exit("set BNGPATH to the folder that holds BNG2.pl")
    return os.path.join(bngpath, "BNG2.pl")


def run_bng(stem, workdir):
    r = subprocess.run(
        ["perl", bng2pl(), stem + ".bngl"], cwd=workdir, capture_output=True, text=True
    )
    if r.returncode != 0:
        sys.exit(f"BNG2.pl failed on {stem}:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")


def read_gdat(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, comments="#", ndmin=2)
    return {h: data[:, i] for i, h in enumerate(header)}


def read_net(path):
    text = open(path).read()

    def block(name):
        m = re.search(rf"^begin {name}\n(.*?)^end {name}", text, re.S | re.M)
        return [ln.split("#")[0].split() for ln in m.group(1).splitlines() if ln.split("#")[0]]

    species = {int(c[0]): c[1] for c in block("species")}

    def side(field):
        return tuple(sorted(species[int(i)] for i in field.split(",") if i != "0"))

    reactions = {(side(c[1]), side(c[2])): c[3] for c in block("reactions")}
    return set(species.values()), reactions


def exact_with_protocol(stem_suffix, knockout_jnk=False):
    """p38atf2_binding.bngl, unedited, then the published equilibrate -> stimulate protocol."""
    inhibitor = 'setParameter("k1",0)\nsetParameter("k2",0)\n' if knockout_jnk else ""
    return (
        open(EXACT).read()
        + "\nbegin actions\ngenerate_network({overwrite=>1})\n"
        + inhibitor
        + 'setParameter("Stim",0)\n'
        + f'simulate({{prefix=>"eq{stem_suffix}",method=>"ode",t_start=>-20000,t_end=>0,'
        + "n_steps=>10})\n"
        + 'setParameter("Stim",1)\n'
        + f'simulate({{prefix=>"st{stem_suffix}",method=>"ode",t_start=>0,t_end=>2466,'
        + "n_steps=>2466})\nend actions\n"
    )


def read_exp(path):
    a = np.loadtxt(path, comments="#", ndmin=2)
    return a[:, 0], a[:, 1], a[:, 2]


def data_comparison(sim, exp_path):
    t, y, se = read_exp(exp_path)
    fold = sim["p38ATF2all"] / sim["p38ATF2all"][0]
    idx = [int(np.flatnonzero(sim["time"] == ti)[0]) for ti in t]
    ym = fold[idx]
    z = (ym - y) / se
    rel = (ym - y) / y
    k_data, k_model = int(np.argmax(y)), int(np.argmax(fold))
    return {
        "n_points": len(t),
        "half_chi2": float(0.5 * np.sum(z**2)),
        "median_abs_rel_err": float(np.median(np.abs(rel))),
        "max_abs_rel_err": float(np.max(np.abs(rel))),
        "median_abs_z": float(np.median(np.abs(z))),
        "max_abs_z": float(np.max(np.abs(z))),
        "z_min": float(z.min()),
        "z_max": float(z.max()),
        "n_beyond_3_se": int(np.sum(np.abs(z) > 3)),
        "data_peak": [float(y[k_data]), float(t[k_data])],
        "model_peak": [float(fold[k_model]), float(sim["time"][k_model])],
        "t0_data_vs_1_in_se": float((y[0] - 1.0) / se[0]),
        "points": [
            {"t_s": float(a), "data": float(b), "se": float(c), "model": float(d), "z": float(e)}
            for a, b, c, d, e in zip(t, y, se, ym, z, strict=True)
        ],
    }


def main():
    pin = json.load(open(os.path.join(HERE, "upstream.json")))["files"][AUTHOR_REL]
    author = os.path.join(HERE, AUTHOR_REL)
    got = sha256(author)
    if got != pin["sha256"]:
        sys.exit(f"{AUTHOR_REL}: sha256 {got} does not match upstream.json {pin['sha256']}")

    with tempfile.TemporaryDirectory() as wd:
        shutil.copy(author, os.path.join(wd, "authors.bngl"))
        run_bng("authors", wd)
        auth = read_gdat(os.path.join(wd, "01.gdat"))  # the stimulation phase, 0-2400 s
        net_a = read_net(os.path.join(wd, "authors.net"))
        for stem, ko in (("exact", False), ("exact_jnkin8", True)):
            open(os.path.join(wd, stem + ".bngl"), "w").write(exact_with_protocol(stem, ko))
            run_bng(stem, wd)
        exact = read_gdat(os.path.join(wd, "stexact.gdat"))
        exact_ko = read_gdat(os.path.join(wd, "stexact_jnkin8.gdat"))
        net_e = read_net(os.path.join(wd, "exact.net"))

    ok = True
    # (a) networks
    sp_a, rx_a = net_a
    sp_e, rx_e = net_e
    law_diff = sorted({(rx_e[k], rx_a[k]) for k in rx_a if k in rx_e and rx_e[k] != rx_a[k]})
    net_ok = (
        sp_a == sp_e
        and len(sp_a) == 48
        and set(rx_a) == set(rx_e)
        and len(rx_a) == 152
        and law_diff == sorted(GATED.items())
    )
    ok &= net_ok

    # (b) trajectories on the authors' grid
    n = len(auth["time"])
    grid_ok = np.array_equal(auth["time"], exact["time"][:n])
    traj = {}
    for s in OBS:
        peak = float(np.max(np.abs(auth[s])))
        diff = float(np.max(np.abs(exact[s][:n] - auth[s])))
        traj[s] = {
            "peak_authors": peak,
            "max_abs_diff": diff,
            "value": diff / peak,
            "tol": TOL_PEAK,
        }
        traj[s]["pass"] = bool(grid_ok and diff / peak <= TOL_PEAK)
        ok &= traj[s]["pass"]

    # (c) reported only: the published point against the two NanoBit curves
    data = {
        "WT control": data_comparison(exact, DATA["WT control"]),
        "WT + JNK-IN-8": data_comparison(exact_ko, DATA["WT + JNK-IN-8"]),
    }
    joint = sum(v["half_chi2"] for v in data.values())

    report = {
        "exact_file": "p38atf2_binding.bngl",
        "exact_file_sha256": sha256(EXACT),
        "author_file": AUTHOR_REL,
        "author_file_sha256": got,
        "engine": "BioNetGen BNG2.pl (run_network, CVODE); no bngsim",
        "asserted": {
            "a_network": {
                "species": [len(sp_a), len(sp_e)],
                "reactions": [len(rx_a), len(rx_e)],
                "rate_laws_that_differ_exact_vs_authors": law_diff,
                "pass": bool(net_ok),
            },
            "b_trajectories_on_authors_grid": {
                "n_points": n,
                "t_range_s": [float(auth["time"][0]), float(auth["time"][-1])],
                "observables": traj,
            },
        },
        "reported_not_asserted": {
            "what": "p38ATF2all(t)/p38ATF2all(0) at the published point vs Source Data Figure_4 "
            "(Fig. 4a); error bars are propagated standard errors, n = 3",
            "curves": data,
            "joint_half_chi2_both_curves": joint,
            "note": "kstim6, dp2 and dp4 were fit by the authors to both curves (Methods); this "
            "job fits the WT control alone, so its single-curve objective is not the authors'.",
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")

    print("p38atf2_binding.bngl vs the authors' Bionetgen_JNK-p38-ATF2_model.bngl (BNG2.pl)")
    print(
        f"  (a) network: {len(sp_a)}/{len(sp_e)} species, {len(rx_a)}/{len(rx_e)} reactions,"
        f" rate laws differ only at {law_diff}: {'PASS' if net_ok else 'FAIL'}"
    )
    for s, row in traj.items():
        print(
            f"  (b) {s:10s} {n} points, max |diff|/peak {row['value']:.1e}"
            f" (tol {TOL_PEAK:.0e}) {'PASS' if row['pass'] else 'FAIL'}"
        )
    print("  reported, not asserted (published point vs Fig. 4a NanoBit, SE n = 3):")
    for name, v in data.items():
        print(
            f"    {name:14s} 1/2 chi^2 {v['half_chi2']:6.2f};"
            f" median |rel| {v['median_abs_rel_err']:.3f}, max {v['max_abs_rel_err']:.3f};"
            f" z {v['z_min']:+.2f}..{v['z_max']:+.2f};"
            f" {v['n_beyond_3_se']}/{v['n_points']} beyond 3 SE; peak data {v['data_peak'][0]:.3f}"
            f" @ {v['data_peak'][1]:.0f} s, model {v['model_peak'][0]:.3f}"
            f" @ {v['model_peak'][1]:.0f} s"
        )
    print(f"    joint 1/2 chi^2 over both curves: {joint:.2f}")
    print(f"verdict (a, b): {report['verdict']}  -> {os.path.relpath(OUT_JSON, HERE)}")
    plot(auth, exact, exact_ko, traj, data, joint, report["verdict"])
    return 0 if ok else 1


def plot(auth, exact, exact_ko, traj, data, joint, verdict):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(OBS, plt.rcParams["axes.prop_cycle"].by_key()["color"], strict=False))
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    n = len(auth["time"])
    for s in OBS:
        peak = traj[s]["peak_authors"]
        ax[0].plot(exact["time"] / 60, exact[s] / peak, "-", color=colors[s], lw=1.5, label=s)
        ax[0].plot(
            auth["time"][::120] / 60, auth[s][::120] / peak, "o", mfc="none", mec=colors[s], ms=5
        )
        r = np.abs(exact[s][:n] - auth[s]) / peak
        ax[1].semilogy(auth["time"] / 60, np.maximum(r, 1e-18), "-", color=colors[s], lw=1, label=s)
    ax[0].set(
        xlabel="time after anisomycin (min)",
        ylabel="observable / its peak",
        title="lines = p38atf2_binding.bngl, open markers = authors' file",
    )
    ax[0].legend(frameon=False, fontsize=8)
    worst = max(v["value"] for v in traj.values())
    ax[1].axhline(TOL_PEAK, color="k", ls="--", lw=1, label=f"tolerance {TOL_PEAK:.0e}")
    ax[1].set_ylim(1e-18, 1e-3)
    if worst == 0.0:
        ax[1].text(
            0.5,
            0.3,
            "every difference is exactly 0.0\n(4 observables x 2401 points)",
            transform=ax[1].transAxes,
            ha="center",
            fontsize=11,
        )
    ax[1].set(
        xlabel="time after anisomycin (min)",
        ylabel="|difference| / peak (zeros drawn at 1e-18)",
        title=f"author-file check, asserted: worst {worst:.1e}",
    )
    ax[1].legend(frameon=False, fontsize=8)
    for (name, v), sim, col in zip(
        data.items(), (exact, exact_ko), ("#222222", "#1f77b4"), strict=True
    ):
        pts = v["points"]
        t = np.array([p["t_s"] for p in pts]) / 60
        ax[2].errorbar(
            t,
            [p["data"] for p in pts],
            yerr=[p["se"] for p in pts],
            fmt="o",
            ms=4,
            color=col,
            capsize=2,
            label=f"{name}, Fig. 4a (mean +/- SE, n = 3)",
        )
        fold = sim["p38ATF2all"] / sim["p38ATF2all"][0]
        ax[2].plot(
            sim["time"] / 60,
            fold,
            "-",
            color=col,
            lw=1.6,
            label=f"model, published point: 1/2 chi^2 {v['half_chi2']:.1f},"
            f" {v['n_beyond_3_se']}/{v['n_points']} beyond 3 SE",
        )
    ax[2].set(
        xlabel="time after anisomycin (min)",
        ylabel="p38-ATF2 NanoBit, treated/untreated",
        title=f"data, REPORTED NOT ASSERTED (joint 1/2 chi^2 {joint:.1f})",
    )
    ax[2].legend(frameon=False, fontsize=7.5, loc="upper left")
    ax[2].grid(alpha=0.25)
    fig.suptitle(
        "Kirsch et al. 2020: p38atf2_binding.bngl vs the authors' Supplementary Software 1 model"
        f" (BNG2.pl). Network identical and trajectories within {TOL_PEAK:.0e} of peak: {verdict}."
        "\nThe published point was fit jointly to both NanoBit curves; at their stated SE it fits"
        " neither (residuals of opposite sign).",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT_PNG, dpi=110)
    print("wrote", os.path.relpath(OUT_PNG, HERE))


if __name__ == "__main__":
    sys.exit(main())
