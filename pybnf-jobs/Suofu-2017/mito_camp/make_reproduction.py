#!/usr/bin/env python
"""Reproduction figure for the Suofu-2017 mito_camp job.

Runs the two stimulations of Suofu et al. (2017) Fig. 4J -- melatonin (MT1 on the plasma
membrane AND the outer mitochondrial membrane, the paper's "model 2") and DAMGO (mu-opioid
receptor on the plasma membrane only, "model 1") -- from a common receptor-silent steady
state, and compares the result with

  * the fit targets: the averaged cAMP time courses of Fig. 4H and Fig. 4I, digitized into
    `melatonin.exp` and `damgo.exp`;
  * the authors' own fit of the same model to the same data, the four curves plotted in
    Fig. 4J, digitized into the library model's
    `reference/suofu2017_fig4J_model_curves_digitized.csv`.

Metrics per curve: RMSE in percentage points of the pre-agonist maximum (the units of both
the data and the observables) and the median absolute error. Percentage points, not
relative error, because three of the four curves decay into the 0-30% band where a relative
error divides by a small number.

Simulation goes through BNG2.pl, so the figure reproduces without the bngsim toolchain. The
basal state is reached by a long unstimulated integration rather than by the steady-state
solve PyBNF uses for `preequilibrate:`; the script asserts that it converged.

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy.
Usage: BNGPATH=... python make_reproduction.py
"""
import os
import subprocess
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "mito_camp.bngl")
FIG4J = os.path.join(HERE, "..", "..", "..", "models",
                     "mitochondrial_mt1_camp_signaling_suofu2017", "reference",
                     "suofu2017_fig4J_model_curves_digitized.csv")

T_EQUIL = 1.0e7     # s, unstimulated; convergence is asserted below
T_END = 600.0       # s, the span of Fig. 4H/4I/4J
# BioNetGen writes function observables into the .gdat without parentheses;
# the .exp headers carry them, as PyBNF requires.
OBS = ("camp_CY_pct", "camp_ML_pct")

# The two stimulations, and which digitized Fig. 4J columns they correspond to
ARMS = {
    "melatonin": dict(switch="MT1_isActive", label="Melatonin (MT1: PM + OMM)",
                      paper=("MT1_PM", "MT1_OMM")),
    "damgo": dict(switch="OR_isActive", label="DAMGO (mu-OR: PM only)",
                  paper=("uOR_PM", "uOR_OMM")),
}
CURVES = ((OBS[0], "PM sensor / cytosol", "tab:blue"),
          (OBS[1], "OMM sensor / lumen", "k"))


def simulate(overrides=None):
    """Equilibrate mito_camp.bngl unstimulated, then run both stimulations.

    Returns {arm: (t, {obs: y})} plus the basal free-cAMP amounts.
    """
    text = open(MODEL).read()
    text = text[:text.index("end model") + len("end model")]
    setp = "".join(f'  setParameter("{k}",{v!r})\n' for k, v in (overrides or {}).items())
    acts = ["begin actions", "  generate_network({overwrite=>1})", setp.rstrip("\n"),
            f'  simulate({{method=>"ode",suffix=>"eq",t_start=>0,t_end=>{T_EQUIL},'
            "n_steps=>20,atol=>1e-10,rtol=>1e-10,print_functions=>1})",
            '  saveConcentrations("ceq")']
    for arm, spec in ARMS.items():
        acts += [f'  resetConcentrations("ceq")',
                 *[f'  setParameter("{s["switch"]}",{1 if a == arm else 0})'
                   for a, s in ARMS.items()],
                 f'  simulate({{method=>"ode",suffix=>"{arm}",t_start=>0,t_end=>{T_END},'
                 "n_steps=>600,atol=>1e-10,rtol=>1e-10,print_functions=>1})"]
    acts.append("end actions")
    text += "\n\n" + "\n".join(a for a in acts if a.strip()) + "\n"

    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        open(os.path.join(tmp, "m.bngl"), "w").write(text)
        res = subprocess.run(["perl", BNG, "m.bngl"], cwd=tmp,
                             capture_output=True, text=True)
        assert res.returncode == 0, res.stderr or res.stdout[-3000:]
        for tag in ("eq", *ARMS):
            path = os.path.join(tmp, f"m_{tag}.gdat")
            with open(path) as fh:
                header = fh.readline().lstrip("#").split()
            d = np.loadtxt(path)
            out[tag] = (d[:, 0], {o: d[:, header.index(o)] for o in header[1:]})

    t, y = out["eq"]
    drift = max(abs(y[o][-1] / y[o][-2] - 1.0) for o in OBS)
    assert drift < 1e-4, f"basal state not converged at t = {T_EQUIL:g} s (drift {drift:.2e})"
    return out


def load_exp(arm):
    d = np.atleast_2d(np.loadtxt(os.path.join(HERE, f"{arm}.exp")))
    return {OBS[0]: d[:, [0, 1]], OBS[1]: d[:, [0, 2]]}


def main():
    sim = simulate()
    paper = np.genfromtxt(FIG4J, delimiter=",", names=True)

    print(f"basal free cAMP: cytosol {sim['eq'][1]['cAMP_CY'][-1]:.4g} molecules "
          f"(Table 2: 1.0e7), lumen {sim['eq'][1]['cAMP_ML'][-1]:.4g} "
          f"(Table 2: 1.131e5)\n")

    print(f"{'arm':10s} {'curve':22s} {'n':>3s} {'RMSE (pp)':>10s} "
          f"{'RMSE Fig. 4J':>13s} {'model 600 s':>12s} {'Fig. 4J 600 s':>14s}")
    sse_tot = sse_paper = n_tot = 0
    for arm, spec in ARMS.items():
        exp = load_exp(arm)
        t, y = sim[arm]
        for (obs, label, _), pcol in zip(CURVES, spec["paper"]):
            d = exp[obs]
            d = d[~np.isnan(d[:, 1])]
            resid = np.interp(d[:, 0], t, y[obs]) - d[:, 1]
            # the authors' own fit, scored on exactly the same points
            rp = 100 * np.interp(d[:, 0], paper["time_s"], paper[pcol]) - d[:, 1]
            sse_tot += float((resid ** 2).sum())
            sse_paper += float((rp ** 2).sum())
            n_tot += len(d)
            print(f"{arm:10s} {label:22s} {len(d):3d} "
                  f"{np.sqrt((resid**2).mean()):10.2f} {np.sqrt((rp**2).mean()):13.2f} "
                  f"{y[obs][-1]:11.1f}% {100*paper[pcol][-1]:13.1f}%")
    print(f"\ntotal over {n_tot} points -- this fit: SSE {sse_tot:.1f}, "
          f"RMSE {np.sqrt(sse_tot/n_tot):.2f} pp, pybnf sos {sse_tot/2:.1f}")
    print(f"{'':21s}Suofu Fig. 4J: SSE {sse_paper:.1f}, "
          f"RMSE {np.sqrt(sse_paper/n_tot):.2f} pp, pybnf sos {sse_paper/2:.1f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, (arm, spec) in zip(axes, ARMS.items()):
        exp = load_exp(arm)
        t, y = sim[arm]
        for (obs, label, colour), pcol in zip(CURVES, spec["paper"]):
            ax.plot(t, y[obs], color=colour, lw=1.8, label=f"this fit -- {label}")
            ax.plot(paper["time_s"], 100 * paper[pcol], color=colour, lw=1.2, ls="--",
                    alpha=0.8, label=f"Suofu Fig. 4J -- {label}")
            d = exp[obs]
            d = d[~np.isnan(d[:, 1])]
            ax.plot(d[:, 0], d[:, 1], "o", mfc="none", mec=colour, ms=5,
                    label=f"Fig. 4{'H' if arm == 'melatonin' else 'I'} data -- {label}")
        ax.set_title(spec["label"], fontsize=11)
        ax.set_xlabel("time after agonist (s)")
        ax.set_xlim(0, T_END)
        ax.set_ylim(0, 115)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("cAMP (% of pre-agonist basal)")
    axes[0].legend(frameon=False, fontsize=7, loc="upper right", ncol=1)
    fig.suptitle("Mitochondrial MT1 cAMP signaling (Suofu et al. 2017) -- PyBNF refit vs. "
                 "the Fig. 4H/4I data and the authors' Fig. 4J curves", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(HERE, "mito_camp_reproduction.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {os.path.basename(out)}")


if __name__ == "__main__":
    main()
