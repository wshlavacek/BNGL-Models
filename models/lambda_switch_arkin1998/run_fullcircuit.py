"""Driver + verification for lambda_switch_arkin1998_fullcircuit.bngl -- the full,
closed-loop, network-free bacteriophage-lambda decision circuit of Arkin, Ross &
McAdams (1998), Genetics 149:1633-1648.

Self-contained: reads the committed .bngl next to it, overrides the MOI parameter by
text substitution, generates BNG-XML with the `bionetgen` CLI, and simulates the XML
network-free with NFsim or RuleMonkey via `bngsim`. Reproduces Arkin Figs 3 and 6.

Requirements: `bionetgen` on PATH (BNG >= 2.9.3) and the `bngsim` package (which wraps
NFsim + RuleMonkey). Run with a Python that has bngsim + numpy, e.g.
    python run_fullcircuit.py fig3   60          # Fig 3a: MOI=6 avg Cro2/CI2 +/-1sigma
    python run_fullcircuit.py fig6   40           # Fig 6: f_lysogeny(MOI) + Poisson->API
    python run_fullcircuit.py traj   40  6  rm    # one MOI mean trajectory (engine parity)

Units: V_ref=1.66054e-15 L makes 1 molecule ~ 1 nM, so dimer COUNTS are the nM
concentrations on Arkin's axes.

Lysogeny classifier (Arkin footnote 12, p.1639-1640): a cell commits to lysogeny when
free-CII production activated P_RE and [CI2] > [Cro2] at the end of the 35-min cell
cycle. Operationally per run: LYSOGENIC <=> CI2_dimer > Cro2_dimer at t=35 min (the
footnote states this end-state test is itself the indication that P_RE activation
locked in the CI feedback loop). A stricter "regime" cross-check also requires CI2 to
be in the lysogenic band (CI2 >= 50 nM).
"""
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# Which .bngl to drive. Default is the fast decision-level full circuit. Set env var
# FULLCIRCUIT_MODEL=exact to drive the paper-exact companion (issue #6: cell-growth
# dilution + Table-2 antitermination saturation); or pass an explicit .bngl path.
_MODEL_ALIASES = {
    "base": "lambda_switch_arkin1998_fullcircuit.bngl",
    "exact": "lambda_switch_arkin1998_fullcircuit_exact.bngl",
}
_sel = os.environ.get("FULLCIRCUIT_MODEL", "base")
MODEL = os.path.join(HERE, _MODEL_ALIASES.get(_sel, _sel))
T_END = 2100.0                      # 35-min cell cycle

# observables carried by the model (see the .bngl observables block)
TRAJ = ["CI2_dimer", "Cro2_dimer", "Obs_CII", "CII_total", "CIII_tot", "Obs_N",
        "CImon", "Crmon"]


def build_xml(moi, tag=None):
    """Write a MOI-substituted copy of the model to a temp dir, run bionetgen
    writeXML, and return the XML path. Uses a system temp dir so running the driver
    never litters the model directory."""
    tag = tag or f"moi{moi}"
    workdir = tempfile.mkdtemp(prefix=f"fullcircuit_{tag}_")
    text = open(MODEL).read()
    text, n = re.subn(r"(?m)^(\s*MOI\s+)\d+", rf"\g<1>{moi}", text)
    assert n == 1, f"expected exactly one MOI parameter line, found {n}"
    bngl = os.path.join(workdir, f"fullcircuit_{tag}.bngl")
    open(bngl, "w").write(text)
    out = os.path.join(workdir, "out")
    r = subprocess.run(["bionetgen", "run", "-i", bngl, "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    return os.path.join(out, f"fullcircuit_{tag}.xml")


def ensemble(xml, method, n_seeds, n_points=8, seed0=30000):
    """Return (tvec, {obs: array[n_seeds, n_points]}) of RAW per-seed trajectories.
    method: 'rm' (RuleMonkey, the fast workhorse) or 'nf' (NFsim)."""
    from bngsim import NfsimSession, RuleMonkeySession
    Sess = NfsimSession if method == "nf" else RuleMonkeySession
    stacks = {o: [] for o in TRAJ}
    tvec = None
    for s in range(n_seeds):
        with Sess(xml) as sess:
            sess.initialize(seed=seed0 + s)
            res = sess.simulate(0, T_END, n_points=n_points)
            names = list(sess.get_observable_names())
            Y = np.asarray(res.observables)
            tvec = np.asarray(res.time)
            for o in TRAJ:
                stacks[o].append(Y[:, names.index(o)])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


def classify(stacks):
    """Footnote-12 per-seed lysogeny classification from final-time (t=35min) counts."""
    ci2 = stacks["CI2_dimer"][:, -1]
    cro2 = stacks["Cro2_dimer"][:, -1]
    lyso = ci2 > cro2
    return {"lyso": lyso, "regime": lyso & (ci2 >= 50.0), "ci2": ci2, "cro2": cro2}


# --------------------------------------------------------------------- Fig 3 (MOI 6)
def fig3(n_seeds=60, method="rm", moi=6):
    print(f"=== Fig 3a: MOI={moi} average Cro2/CI2 with +/-1sigma (16/84 pctile) band, "
          f"{n_seeds} {method} seeds ===")
    xml = build_xml(moi)
    t0 = time.time()
    tvec, st = ensemble(xml, method, n_seeds, n_points=8, seed0=30000)
    print(f"  [{time.time()-t0:.0f}s]  t(min):" + "".join(f"{t/60:7.0f}" for t in tvec))
    for o, lab in (("Cro2_dimer", "Cro2"), ("CI2_dimer", "CI2")):
        print(f"  {lab} avg :" + "".join(f"{v:7.1f}" for v in st[o].mean(0)))
        print(f"  {lab} 16pc:" + "".join(f"{v:7.1f}" for v in np.percentile(st[o], 16, 0)))
        print(f"  {lab} 84pc:" + "".join(f"{v:7.1f}" for v in np.percentile(st[o], 84, 0)))
    c = classify(st)
    print(f"  => lysogenic fraction (CI2>Cro2 @35min): {c['lyso'].mean():.2f}  "
          f"(regime CI2>=50: {c['regime'].mean():.2f}) of {n_seeds}")
    print("  paper Fig 3a: Cro2 avg ~55-60 nM plateau, CI2 lower with a broad band.")
    return tvec, st


# ------------------------------------------------------------------- Fig 6 (MOI, API)
def lyso_fraction(moi, n_seeds, method="rm"):
    xml = build_xml(moi)
    t0 = time.time()
    _, st = ensemble(xml, method, n_seeds, n_points=2, seed0=20000)
    c = classify(st)
    return c["lyso"].mean(), c["regime"].mean(), c["ci2"].mean(), c["cro2"].mean(), time.time() - t0


def poisson_weight(f_of_moi, api_values, m_max=40):
    """F(API) = sum_{M>=1} Poisson(M|API) * f(M) (Arkin Eq 1-2). Below the smallest
    simulated MOI use f=0 (paper: ~0 at low MOI); above the largest, hold the last
    value (saturating); between simulated MOIs, linear-interpolate."""
    from math import exp, lgamma, log
    mois = sorted(f_of_moi)

    def fval(m):
        if m < mois[0]:
            return 0.0
        if m >= mois[-1]:
            return f_of_moi[mois[-1]]
        if m in f_of_moi:
            return f_of_moi[m]
        lo = max(x for x in mois if x <= m)
        hi = min(x for x in mois if x >= m)
        w = (m - lo) / (hi - lo)
        return (1 - w) * f_of_moi[lo] + w * f_of_moi[hi]

    out = {}
    for api in api_values:
        out[api] = sum(exp(-api + m * log(api) - lgamma(m + 1)) * fval(m)
                       for m in range(1, m_max + 1))
    return out


def fig6(n_seeds=40, method="rm", mois=(1, 2, 3, 4, 5, 6, 8, 10)):
    print(f"=== Fig 6a: f_lysogeny(MOI), {n_seeds} {method} seeds/MOI ===")
    print("  paper 'Full' curve (Table-3 proteolysis): ~0 for MOI<3, rapid rise MOI>3.")
    f = {}
    for moi in mois:
        fl, fr, mci2, mcro2, secs = lyso_fraction(moi, n_seeds, method)
        f[moi] = fl
        print(f"  MOI={moi:2d}: f_lyso={fl:.2f}  (regime={fr:.2f})  "
              f"<CI2>={mci2:5.1f} <Cro2>={mcro2:5.1f}  [{secs:.0f}s]")
    print("\n=== Fig 6b: Poisson-weighted fraction lysogens vs API (Eq 1-2) ===")
    apis = [0.5, 1, 2, 3, 5, 7, 10, 15]
    F = poisson_weight(f, apis)
    for api in apis:
        print(f"  API={api:5.1f}: F(lyso)={F[api]:.4f}")
    return f, F


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "fig3"
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else (60 if mode == "fig3" else 40)

    if mode == "fig3":
        fig3(n_seeds, sys.argv[3] if len(sys.argv) > 3 else "rm")
    elif mode == "fig6":
        fig6(n_seeds, sys.argv[3] if len(sys.argv) > 3 else "rm")
    elif mode == "traj":
        moi = int(sys.argv[3]) if len(sys.argv) > 3 else 6
        method = sys.argv[4] if len(sys.argv) > 4 else "rm"
        xml = build_xml(moi)
        t0 = time.time()
        tvec, st = ensemble(xml, method, n_seeds)
        print(f"MOI={moi} {method} {n_seeds} seeds [{time.time()-t0:.0f}s]  t(min):" +
              "".join(f"{t/60:6.0f}" for t in tvec))
        for o in ("Cro2_dimer", "CI2_dimer", "Obs_N", "Obs_CII", "CIII_tot"):
            print(f"  {o:>10s}:" + "".join(f"{v:6.1f}" for v in st[o].mean(0)))
    else:
        print(f"unknown mode {mode!r}; use fig3 | fig6 | traj")
