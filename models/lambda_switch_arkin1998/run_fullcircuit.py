"""Driver + verification for lambda_switch_arkin1998_fullcircuit.bngl -- the full,
closed-loop, network-free bacteriophage-lambda decision circuit of Arkin, Ross &
McAdams (1998), Genetics 149:1633-1648.

Self-contained: reads the committed .bngl next to it, overrides the MOI parameter by
text substitution, generates BNG-XML with the `bionetgen` CLI, and simulates the XML
network-free with NFsim or RuleMonkey via `bngsim`. Reproduces Arkin Figs 3 and 6.

Requirements: `bionetgen` on PATH (BNG >= 2.9.3) and the `bngsim` package (which wraps
NFsim + RuleMonkey). Run with a Python that has bngsim + numpy, e.g.
    python run_fullcircuit.py fig3       60       # Fig 3a: MOI=6 avg Cro2/CI2 +/-1sigma
    python run_fullcircuit.py fig6       48       # Fig 6: f_lysogeny(MOI) + Poisson->API
    python run_fullcircuit.py agreement 100       # NFsim vs RuleMonkey, 100 seeds/engine
    python run_fullcircuit.py traj       40  6 rm # one MOI mean trajectory
Every mode caches to reference/<prefix>_*.npz (prefix = fullcircuit | fullcircuit_exact),
records the engine that produced it, and runs with the seed count the committed cache was
produced with when none is given -- so a bare `fig3` / `fig6` / `agreement` regenerates the
committed artifact rather than overwriting it with a differently-sized one.

VERIFICATION LEVEL 1 -- `agreement`. There is no reaction network to integrate here, so the
independent implementation is a different network-free ALGORITHM on the same rules, not a
different build of the same one: bngsim's RuleMonkey (exact -- computes event times over
every particle, rejects nothing) standing witness for bngsim's NFsim (rejection-based). Both
sample the same CTMC, so their ensemble means must agree; a disagreement can only be in the
sampler. Both engines read the SAME XML with MOI already baked in by build_xml() -- never
set_param() after initialize(), which NFsim does not propagate into rule rates and RuleMonkey
does not propagate into seed-species counts (lanl/bngsim#44).

Two arms, not the usual three: BNG2.pl's bundled NFsim binary is not an arm because the
committed .bngl is writeXML-only by design (its promoter rates read Species-type dimer
observables, which NFsim freezes -- RuleWorld/nfsim#86 -- and the Molecules-type workaround
plus complex bookkeeping is what the driver supplies), so there is no legacy-binary run a
reader reproduces from the model file. -bscb is inert for this model and is left at the
bngsim default. That follows from the rules -- every bond-forming rule (R_dim_CI, R_dim_Cro,
and the four protease-binding rules) joins two FREE monomers, the largest species is a dimer,
and no rule can bind two sites already in a common complex, so there is nothing for
same-complex blocking to act on and no traversal depth beyond 2 -- and it was measured rather
than left as an argument: 50 seeds at MOI 6 with block_same_complex_binding True vs False give
BIT-IDENTICAL ensembles on every observable, and the offset against RuleMonkey below is
unchanged by the flag. So the flag is not the explanation for that offset.

The headline statistic is the LYSOGENY FRACTION, which is binomial: its noise is
sqrt(p(1-p)/n) and does not shrink by averaging over observables, so `agreement` runs 100
seeds per engine (cf. lambda_switch_cortes2017) rather than the 12-24 that suit a continuous
observable. Ensemble means are compared two ways: pairwise z over 8 observables x 36 time
points x each MOI, judged on max|z| against the expected maximum of |N| standard normals
together with the fraction below 3 (never on max|z| < 3 alone -- see
skills/curate-model/references/network-free-verification.md Sec. 5); and a per-observable
ENDPOINT z combined across MOIs. The second test exists because the first one missed a real
defect: hundreds of early-time comparisons where both engines sit near zero dilute a
systematic late-time offset, and a difference that shows up as |z| ~ 2.2 at each of two MOIs
IN THE SAME DIRECTION never trips a max|z| threshold, though combined it is |z| = 3.1.

RESULT (100 seeds/engine, MOI 4 and 6, bngsim 0.12.2+dee86da42547):

  * lambda_switch_arkin1998_fullcircuit.bngl -- the engines AGREE. max|z| = 2.81 over 558
    trajectory comparisons against an expected maximum of 3.12, 100% below 3; largest
    combined endpoint |z| = 1.5; lysogeny fractions differ by <= 0.07 against a 0.05
    binomial noise scale. Cached in reference/fullcircuit_agreement.npz.

  * lambda_switch_arkin1998_fullcircuit_exact.bngl -- the engines DO NOT AGREE, and the
    cause is an NFsim defect, not the model. NFsim runs the free N / CIII / CII pools ~30%
    high (combined endpoint z = +3.4 for Obs_N) and the CI2 / Cro2 dimers ~15% low
    (z = -2.9), the same direction at both MOIs. A same-engine null (RuleMonkey against
    itself on a disjoint seed block, identical statistic) gives |z| <= 1.9 with no sign
    structure, so the error bars are sized correctly and this is not heavy-tailed noise.

    Localized to a minimal model and filed as lanl/bngsim#195: for a rule whose reactant
    pattern is SYMMETRIC and whose
    rate is a FUNCTION, NFsim does not divide out the pattern's automorphism and fires at
    twice the intended rate. A(d!1).A(d!1) -> 0 at a functional dilution rate decays to
    497.5 where the analytic answer and RuleMonkey both give ~1000. The same rule with a
    CONSTANT rate is correct in both engines, and an asymmetric dimer C(d!1).D(d!1) at the
    same functional rate is correct in both -- so it is the combination that breaks.

    In this model that is Rdil_CI2 and Rdil_Cro2, the growth-dilution rules the exact
    variant adds (`CI(d!1).CI(d!1) -> 0 mu_dil()`). Under NFsim the dimers are diluted at
    2x, which weakens Shea-Ackers repression at P_R/P_L and raises transcription -- every
    observed sign follows. The base model has no functional-rate dilution rules at all,
    which is exactly why it agrees.

    CONSEQUENCE: run the exact variant with RuleMonkey. Every committed artifact in
    reference/ is RuleMonkey and is unaffected; RuleMonkey reproduces the analytic dilution
    solution to 0.5%. Re-run `agreement` after any bngsim upgrade -- the cache records the
    build it was produced with, and `rescore` recomputes every statistic from the cached
    per-seed ensembles without re-simulating.

Units: V_ref=1.66054e-15 L makes 1 molecule ~ 1 nM, so dimer COUNTS are the nM
concentrations on Arkin's axes.

Lysogeny classifier (Arkin footnote 12, p.1639-1640): a cell commits to lysogeny iff
BOTH (i) free-CII production activated P_RE -- operationally, the total P_RE open-complex
rate A_PRE()*MOI reached >= 1 open complex / 2 min over a contiguous 4-min window -- AND
(ii) [CI2] > [Cro2] at the end of the 35-min cell cycle. Both conditions are evaluated
per seed (see classify()); condition (i) needs the free-CII trajectory at ~1-min
resolution, so the Fig-6 sweep samples 36 points. Requiring (i) matters at low MOI,
where it removes cells that reach CI2>Cro2 without a productive P_RE burst -- it lowers
the Fig-6a onset (e.g. MOI 3 ~0.2 -> ~0.08, matching the paper) while leaving MOI>=4
essentially unchanged. A stricter "regime" cross-check also requires CI2 >= 50 nM.
NOTE: at high MOI this model saturates to ~100% lysogeny -- higher than Arkin's Fig 6a
"Full" plateau (~82%). The cause of Arkin's lower plateau is unidentified, and his code
could not be found (as of 2026-07-23), so the difference is not diagnosable against the
source. However, ~100% high-MOI saturation is NOT unreasonable: the modern CII-threshold
model of Cortes, Trinh, Zeng & Balazsi (2017, Biophys J 113:2110-2120) has a lysogeny
probability that likewise "saturates at 1," and reports agreement with experiment. So
this saturation is a defensible behavior, and Arkin's plateau may be the outlier (its
criteria were "ad hoc near the decision boundary"). Documented on
fullcircuit_exact_verification.png (issue #15).
"""
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REFDIR = os.path.join(HERE, "reference")
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

# reference/ cache prefix for whichever model is selected: fullcircuit | fullcircuit_exact
PREFIX = os.path.basename(MODEL)[len("lambda_switch_arkin1998_"):-len(".bngl")]
IS_EXACT = PREFIX.endswith("_exact")

# Protocol constants of the COMMITTED caches. A bare `fig3`/`fig6`/`agreement` reruns
# exactly these, so the cache is regenerable; pass a seed count to run something else.
FIG3_SEEDS, FIG3_MOI, FIG3_POINTS, FIG3_SEED0 = 60, 6, 8, 30000
FIG6_SEEDS, FIG6_POINTS, FIG6_SEED0 = (48 if IS_EXACT else 50), 36, 20000
FIG6_MOIS = (1, 2, 3, 4, 5, 6, 7, 8, 10, 12) if IS_EXACT else (1, 2, 3, 4, 5, 6, 8, 10)
FIG6_APIS = ((0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20) if IS_EXACT else
             (0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 7, 10, 15))
# `agreement`: MOI 4 and 6 bracket the Fig-6a onset the exact variant is claimed to restore,
# and are where f_lyso is far enough from 0 and 1 for a binomial comparison to have power.
AGREE_SEEDS, AGREE_MOIS, AGREE_POINTS, AGREE_SEED0 = 100, (4, 6), 36, 10000

# observables carried by the model (see the .bngl observables block)
TRAJ = ["CI2_dimer", "Cro2_dimer", "Obs_CII", "CII_total", "CIII_tot", "Obs_N",
        "CImon", "Crmon"]

# --- Arkin footnote-12 lysogeny classifier ------------------------------------------
# A cell is committed to lysogeny iff (i) free-CII production ACTIVATED P_RE -- defined
# as the total P_RE open-complex rate A_PRE()*MOI reaching >= 1 open complex / 2 min over
# a contiguous 4-min window -- AND (ii) [CI2] > [Cro2] at the end of the 35-min cycle.
# A_PRE() is the Table-1 P_RE partition function (same expression as in the .bngl); the
# total rate over the cell's MOI genome copies is A_PRE()*MOI.
_RT = 1.9872e-3 * 310.15
_M2M = 1e-9                              # 1 molecule ~ 1 nM at V_ref -> M/molecule
_RNAP_M = 30 * _M2M                      # buffered free RNAP = 30 nM
_PRE_RATE = 1.0 / 120.0                  # 1 open complex / 2 min, in 1/s
_PRE_WIN_MIN = 4.0                       # contiguous activation window (min)


def _A_PRE(cii):
    """P_RE open-complex rate (1/s) vs free-CII (nM), Table-1 P_RE partition function."""
    ci = np.asarray(cii) * _M2M
    w2 = np.exp(9.9 / _RT) * _RNAP_M     # exp(-dG_RE_2/RT), dG_RE_2 = -9.9
    w3 = np.exp(9.7 / _RT) * ci          # dG_RE_3 = -9.7
    w4 = np.exp(21.5 / _RT) * ci * _RNAP_M  # dG_RE_4 = -21.5
    return (w2 * 4e-5 + w4 * 0.015) / (1 + w2 + w3 + w4)


def _pre_activated(cii_traj, tvec, moi):
    """footnote-12 (i): does total P_RE activation A_PRE()*MOI stay >= 1 OC / 2 min over
    some contiguous 4-min window of the trajectory?"""
    apre = _A_PRE(cii_traj) * moi
    dt = (tvec[1] - tvec[0]) / 60.0 if len(tvec) > 1 else _PRE_WIN_MIN
    w = max(1, int(round(_PRE_WIN_MIN / dt)))
    if len(apre) < w:
        return apre.mean() >= _PRE_RATE
    cs = np.concatenate([[0.0], np.cumsum(apre)])
    return bool(((cs[w:] - cs[:-w]) / w >= _PRE_RATE).any())


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


def classify(stacks, tvec, moi):
    """Arkin footnote-12 per-seed lysogeny classification: LYSOGENIC iff P_RE activated
    (i) AND [CI2]>[Cro2] at t=35 min (ii). Needs the CII trajectory (for (i)) at fine
    time resolution and the MOI (total P_RE rate = A_PRE()*MOI). Also returns the
    end-state-only fraction (lyso_end) for reference."""
    ci2 = stacks["CI2_dimer"][:, -1]
    cro2 = stacks["Cro2_dimer"][:, -1]
    end_switch = ci2 > cro2
    cii = stacks["Obs_CII"]
    pre = np.array([_pre_activated(cii[i], tvec, moi) for i in range(cii.shape[0])])
    lyso = end_switch & pre
    return {"lyso": lyso, "lyso_end": end_switch, "regime": lyso & (ci2 >= 50.0),
            "ci2": ci2, "cro2": cro2}


def _bngsim_build():
    """Version + build id of the bngsim core in use. An engine-agreement number is a
    statement about one build, so the build belongs in the artifact (lanl/bngsim#125)."""
    try:
        import bngsim
        v = getattr(bngsim, "__version__", "?")
        from bngsim import _bngsim_core as core
        return f"{v}+{getattr(core, '__build_commit__', '?')}"
    except Exception:
        return "?"


def _cache(name, **arrays):
    """Write reference/<PREFIX>_<name>.npz. Every cache records `engine` and `model`, so
    which sampler produced a committed number is a recorded fact and not an inference from
    a default argument (issue #39)."""
    os.makedirs(REFDIR, exist_ok=True)
    path = os.path.join(REFDIR, f"{PREFIX}_{name}.npz")
    np.savez_compressed(path, model=os.path.basename(MODEL), **arrays)
    print(f"  cached -> reference/{os.path.basename(path)}")
    return path


# --------------------------------------------------------------------- Fig 3 (MOI 6)
def fig3(n_seeds=FIG3_SEEDS, method="rm", moi=FIG3_MOI, save=True):
    print(f"=== Fig 3a: MOI={moi} average Cro2/CI2 with +/-1sigma (16/84 pctile) band, "
          f"{n_seeds} {method} seeds ===")
    xml = build_xml(moi)
    t0 = time.time()
    tvec, st = ensemble(xml, method, n_seeds, n_points=FIG3_POINTS, seed0=FIG3_SEED0)
    print(f"  [{time.time()-t0:.0f}s]  t(min):" + "".join(f"{t/60:7.0f}" for t in tvec))
    for o, lab in (("Cro2_dimer", "Cro2"), ("CI2_dimer", "CI2")):
        print(f"  {lab} avg :" + "".join(f"{v:7.1f}" for v in st[o].mean(0)))
        print(f"  {lab} 16pc:" + "".join(f"{v:7.1f}" for v in np.percentile(st[o], 16, 0)))
        print(f"  {lab} 84pc:" + "".join(f"{v:7.1f}" for v in np.percentile(st[o], 84, 0)))
    c = classify(st, tvec, moi)
    print(f"  => lysogenic fraction (footnote-12): {c['lyso'].mean():.2f}  "
          f"(CI2>Cro2 only: {c['lyso_end'].mean():.2f}; "
          f"regime CI2>=50: {c['regime'].mean():.2f}) of {n_seeds}")
    print("  paper Fig 3a: Cro2 avg ~55-60 nM plateau, CI2 lower with a broad band.")
    if save:
        _cache("fig3", tvec=tvec, engine=method, moi=moi, n_seeds=n_seeds,
               seed0=FIG3_SEED0,
               cro2_mean=st["Cro2_dimer"].mean(0),
               cro2_p16=np.percentile(st["Cro2_dimer"], 16, 0),
               cro2_p84=np.percentile(st["Cro2_dimer"], 84, 0),
               ci2_mean=st["CI2_dimer"].mean(0),
               ci2_p16=np.percentile(st["CI2_dimer"], 16, 0),
               ci2_p84=np.percentile(st["CI2_dimer"], 84, 0))
    return tvec, st


# ------------------------------------------------------------------- Fig 6 (MOI, API)
def lyso_fraction(moi, n_seeds, method="rm"):
    xml = build_xml(moi)
    t0 = time.time()
    # ~1-min sampling so the footnote-12 4-min activation window is resolved
    tvec, st = ensemble(xml, method, n_seeds, n_points=FIG6_POINTS, seed0=FIG6_SEED0)
    c = classify(st, tvec, moi)
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


def fig6(n_seeds=FIG6_SEEDS, method="rm", mois=FIG6_MOIS, apis=FIG6_APIS, save=True):
    print(f"=== Fig 6a: f_lysogeny(MOI), {n_seeds} {method} seeds/MOI ===")
    print("  paper 'Full' curve (Table-3 proteolysis): ~0 for MOI<3, rapid rise MOI>3.")
    f = {}
    for moi in mois:
        fl, fr, mci2, mcro2, secs = lyso_fraction(moi, n_seeds, method)
        f[moi] = fl
        print(f"  MOI={moi:2d}: f_lyso={fl:.2f}  (regime={fr:.2f})  "
              f"<CI2>={mci2:5.1f} <Cro2>={mcro2:5.1f}  [{secs:.0f}s]")
    print("\n=== Fig 6b: Poisson-weighted fraction lysogens vs API (Eq 1-2) ===")
    F = poisson_weight(f, apis)
    for api in apis:
        print(f"  API={api:5.1f}: F(lyso)={F[api]:.4f}")
    if save:
        _cache("fig6", mois=np.array(mois, float), f=np.array([f[m] for m in mois]),
               apis=np.array(apis, float), F=np.array([F[a] for a in apis]),
               n_seeds=n_seeds, engine=method, seed0=FIG6_SEED0)
    return f, F


# ------------------------------------------------- NFsim / RuleMonkey agreement (level 1)
def _mean_se(A):
    """Ensemble mean and standard error of the mean over seeds (axis 0)."""
    return A.mean(0), A.std(0, ddof=1) / np.sqrt(A.shape[0])


def _endpoint_z(res, mois, thresh=1.0):
    """Per-observable z of the ENDPOINT ensemble means, and the same z combined across
    MOIs. The trajectory summary in _traj_z() averages a systematic offset away among
    hundreds of early-time comparisons where both engines sit near zero; a per-observable
    endpoint test does not, and it is what exposed the symmetric-pattern defect described
    in the module docstring. Combining independent MOIs as sum(z)/sqrt(k) turns a
    consistent sign into significance instead of letting it hide below |z|=3 twice over."""
    out = {}
    for o in TRAJ:
        zs = []
        for moi in mois:
            a = res[("nf", moi)][o][:, -1]
            b = res[("rm", moi)][o][:, -1]
            if max(a.mean(), b.mean()) <= thresh:
                continue
            sa = a.std(ddof=1) / np.sqrt(a.size)
            sb = b.std(ddof=1) / np.sqrt(b.size)
            den = np.sqrt(sa**2 + sb**2)
            zs.append((a.mean() - b.mean()) / den if den > 0 else 0.0)
        if zs:
            out[o] = (zs, sum(zs) / np.sqrt(len(zs)))
    return out


def _traj_z(res, mois, thresh=1.0):
    """Pairwise z of NFsim vs RuleMonkey ensemble means over active observables x time
    points x MOI. An observable is active at a given MOI if either engine's mean ever
    exceeds `thresh` molecules; points where both s.e.m. vanish carry no information and
    are dropped rather than scored as z=0."""
    z_all, active = [], []
    for moi in mois:
        for o in TRAJ:
            a, b = res[("nf", moi)][o], res[("rm", moi)][o]
            if max(a.mean(0).max(), b.mean(0).max()) <= thresh:
                continue
            active.append((moi, o))
            ma, sa = _mean_se(a)
            mb, sb = _mean_se(b)
            den = np.sqrt(sa**2 + sb**2)
            z_all.append(np.abs(ma - mb)[den > 0] / den[den > 0])
    z = np.concatenate(z_all) if z_all else np.array([0.0])
    return z, active


def agreement(n_seeds=AGREE_SEEDS, mois=AGREE_MOIS, n_points=AGREE_POINTS, save=True):
    """Level-1 cross-check: bngsim NFsim (rejection-based) against bngsim RuleMonkey
    (exact) on the SAME XML, the SAME seeds and the SAME protocol -- see the module
    docstring for why that pair is the check and why there is no third arm."""
    print(f"=== agreement: bngsim-NFsim vs bngsim-RuleMonkey, {n_seeds} seeds/engine ===")
    print(f"    {os.path.basename(MODEL)}, MOI {list(mois)}, t_end={T_END/60:.0f} min, "
          f"{n_points} points, seeds {AGREE_SEED0}..{AGREE_SEED0+n_seeds-1} (shared)")
    res, cls, tvec = {}, {}, None
    for moi in mois:
        xml = build_xml(moi)                       # one XML, both engines
        for eng in ("nf", "rm"):
            t0 = time.time()
            tvec, st = ensemble(xml, eng, n_seeds, n_points=n_points, seed0=AGREE_SEED0)
            res[(eng, moi)] = st
            cls[(eng, moi)] = classify(st, tvec, moi)
            print(f"  [{time.time()-t0:5.0f}s] MOI={moi} {eng} done")

    print("\n  lysogeny fraction (footnote-12; binomial, s.e. = sqrt(p(1-p)/n)):")
    print("    MOI    f_NF +/- se     f_RM +/- se     |diff|     z")
    dmax = 0.0
    for moi in mois:
        p = {e: cls[(e, moi)]["lyso"].mean() for e in ("nf", "rm")}
        se = {e: np.sqrt(p[e] * (1 - p[e]) / n_seeds) for e in ("nf", "rm")}
        d = abs(p["nf"] - p["rm"])
        dmax = max(dmax, d)
        z = d / np.sqrt(se["nf"] ** 2 + se["rm"] ** 2)
        print(f"    {moi:3d}  {p['nf']:6.2f} +/- {se['nf']:.3f}  {p['rm']:6.2f} +/- "
              f"{se['rm']:.3f}   {d:6.3f}  {z:5.2f}")

    z, active = _traj_z(res, mois)
    n_cmp = z.size
    exp_max = statistics.NormalDist().inv_cdf(1 - 1 / (2 * n_cmp))
    zmax, frac3 = float(z.max()), float((z < 3).mean())
    print(f"\n  trajectory means: {len(active)} active (MOI, observable) pairs, "
          f"{n_cmp} comparisons")
    print(f"    max|z| = {zmax:.2f}   (expected max of {n_cmp} standard normals = "
          f"{exp_max:.2f})")
    print(f"    fraction |z| < 3 = {frac3:.3f}   (expect ~0.997 if both engines sample "
          f"the same process)")

    ez = _endpoint_z(res, mois)
    print("\n  endpoint means per observable (+ = NFsim high); a consistent sign across "
          "MOIs\n  is what a systematic engine defect looks like:")
    print("    observable   " + "".join(f"  MOI{m:>2d}" for m in mois) + "   combined")
    for o, (zs, comb) in ez.items():
        print(f"    {o:>12s} " + "".join(f"{v:7.2f}" for v in zs) + f"   {comb:8.2f}")
    ez_max = max((abs(c) for _, c in ez.values()), default=0.0)
    worst = max(ez, key=lambda o: abs(ez[o][1])) if ez else "-"

    ok = zmax < 1.5 * exp_max and frac3 > 0.99 and ez_max < 3.0
    verdict = "AGREE" if ok else "DISAGREE"
    print(f"\n  => {verdict}: max |f_lyso(NF) - f_lyso(RM)| = {dmax:.3f} over MOI "
          f"{list(mois)}, against a per-engine binomial noise scale of "
          f"{np.sqrt(0.25/n_seeds):.3f} at n={n_seeds};")
    print(f"     largest combined endpoint |z| = {ez_max:.2f} ({worst}).")
    if not ok:
        print("     A combined endpoint |z| >= 3 is a SYSTEMATIC difference, not sampling "
              "noise.\n     See the module docstring: NFsim doubles the rate of a rule "
              "whose reactant\n     pattern is symmetric AND whose rate is a function "
              "(Rdil_CI2 / Rdil_Cro2).")

    if save:
        arrays = {f"{e}__moi{m}__{o}": res[(e, m)][o]
                  for (e, m) in res for o in TRAJ}
        arrays.update({f"lyso__{e}__moi{m}": cls[(e, m)]["lyso"] for (e, m) in cls})
        _cache("agreement", tvec=tvec, engines=np.array(["nf", "rm"]),
               mois=np.array(mois, float), observables=np.array(TRAJ),
               n_seeds=n_seeds, seed0=AGREE_SEED0, t_end=T_END,
               zmax=zmax, frac_z_lt3=frac3, expected_zmax=exp_max,
               n_comparisons=n_cmp, max_df_lyso=dmax, endpoint_zmax=ez_max,
               verdict=verdict, bngsim=_bngsim_build(), **arrays)
    return res, dict(zmax=zmax, frac_z_lt3=frac3, expected_zmax=exp_max,
                     max_df_lyso=dmax, endpoint_zmax=ez_max, verdict=verdict)


# ------------------------------------------------- the engine defect, as its own artifact
_SYM_BNGL = """begin model
begin parameters
  k0 {k0}
  Vol0 {vol0}
  k_grow k0*Vol0
  X0 {x0}
end parameters
begin molecule types
  Src()
  Vol()
  X()
  A(d)
  C(d)
  D(d)
end molecule types
begin seed species
  Src() 1
  Vol() Vol0
  X() X0
  A(d!1).A(d!1) X0
  C(d!1).D(d!1) X0
end seed species
begin observables
  Molecules Obs_Vol Vol()
  Molecules Obs_X   X()
  Species   Sym     A(d!1).A(d!1)
  Species   Asym    C(d!1).D(d!1)
end observables
begin functions
  mu_dil() = k_grow/Obs_Vol
end functions
begin reaction rules
  Rgrow:     Src()           -> Src() + Vol() k_grow
  Rdil_X:    X()             -> 0 mu_dil() DeleteMolecules
  Rdil_Sym:  A(d!1).A(d!1)   -> 0 mu_dil() DeleteMolecules
  Rdil_Asym: C(d!1).D(d!1)   -> 0 mu_dil() DeleteMolecules
end reaction rules
end model
begin actions
  writeXML()
end actions
"""


def symcheck(n_seeds=60):
    """Minimal reproducer for the NFsim defect that makes the exact variant's engine
    agreement fail (issue #39). Reproduces this model's growth-dilution mechanism in
    isolation: a clock species drives Obs_Vol = Vol0*(1+k0*t), and three species are
    diluted at the FUNCTIONAL rate mu_dil() = k_grow/Obs_Vol, whose exact solution is
    X(t) = X0/(1+k0*t) -- an exact halving over the 35-min doubling time.

    A plain molecule and an ASYMMETRIC dimer are correct in both engines. The SYMMETRIC
    dimer A(d!1).A(d!1), whose pattern matches each dimer twice, is halved AGAIN by NFsim:
    the automorphism is not divided out when the rate is a function, so the rule fires at
    2x. `symcheck` is the re-check to run after any bngsim upgrade -- if the Sym row comes
    back at ~1.00x for NFsim, the defect is fixed and `agreement` should be re-run."""
    k0, vol0, x0 = 4.76e-4, 1000, 2000
    print(f"=== symcheck: symmetric reactant pattern at a functional rate, "
          f"{n_seeds} seeds/engine ===")
    print(f"    exact: every diluted species halves over {T_END/60:.0f} min "
          f"(X0={x0} -> {x0/2:.0f}); Vol doubles")
    workdir = tempfile.mkdtemp(prefix="symcheck_")
    bngl = os.path.join(workdir, "symcheck.bngl")
    open(bngl, "w").write(_SYM_BNGL.format(k0=k0, vol0=vol0, x0=x0))
    out = os.path.join(workdir, "out")
    r = subprocess.run(["bionetgen", "run", "-i", bngl, "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    xml = os.path.join(out, "symcheck.xml")

    from bngsim import NfsimSession, RuleMonkeySession
    obs = ["Obs_Vol", "Obs_X", "Sym", "Asym"]
    got = {}
    for eng, Sess in (("rm", RuleMonkeySession), ("nf", NfsimSession)):
        rows = []
        for s in range(n_seeds):
            with Sess(xml) as sess:
                sess.initialize(seed=8000 + s)
                res = sess.simulate(0, T_END, n_points=2)
                names = list(sess.get_observable_names())
                rows.append(np.asarray(res.observables)[-1])
        A = np.vstack(rows)
        got[eng] = {o: A[:, names.index(o)] for o in obs}

    exact = {"Obs_Vol": vol0 * (1 + k0 * T_END), "Obs_X": x0 / (1 + k0 * T_END)}
    exact["Sym"] = exact["Asym"] = exact["Obs_X"]
    print(f"\n  {'observable':>10s} {'exact':>9s} {'RuleMonkey':>20s} {'NFsim':>20s}")
    bad = []
    for o in obs:
        e = exact[o]
        cells = []
        for eng in ("rm", "nf"):
            m = got[eng][o].mean()
            cells.append(f"{m:11.1f} ({m/e:.2f}x)")
            if abs(m / e - 1) > 0.05:
                bad.append((o, eng, m / e))
        print(f"  {o:>10s} {e:9.1f} {cells[0]:>20s} {cells[1]:>20s}")
    if bad:
        print("\n  => DEFECT PRESENT: " + "; ".join(
            f"{eng.upper()} {o} at {r:.2f}x the exact value" for o, eng, r in bad))
        print("     A rule whose reactant pattern is symmetric AND whose rate is a function"
              "\n     fires at twice the intended rate. The exact variant's Rdil_CI2 /"
              "\n     Rdil_Cro2 are exactly this; run that model with RuleMonkey.")
    else:
        print("\n  => no defect: both engines reproduce the analytic dilution to 5%."
              "\n     If this previously failed, re-run `agreement` for both models.")
    return got


def rescore():
    """Recompute every agreement statistic from the cached PER-SEED ensembles, without
    re-running a single simulation. This is what caching per seed rather than per summary
    buys: the endpoint test below was added after the campaign had already been run, and
    it is the test that found the defect."""
    path = os.path.join(REFDIR, f"{PREFIX}_agreement.npz")
    d = np.load(path, allow_pickle=True)
    mois = [int(m) for m in d["mois"]]
    res = {(e, m): {o: d[f"{e}__moi{m}__{o}"] for o in TRAJ}
           for e in ("nf", "rm") for m in mois}
    n_seeds = int(d["n_seeds"])
    print(f"=== rescore {os.path.basename(path)}: {n_seeds} seeds/engine, MOI {mois} ===")
    print(f"    produced by bngsim {str(d['bngsim']) if 'bngsim' in d.files else '(unrecorded)'}")
    for moi in mois:
        p = {e: d[f"lyso__{e}__moi{moi}"].mean() for e in ("nf", "rm")}
        print(f"  MOI {moi}: f_lyso NF={p['nf']:.2f} RM={p['rm']:.2f} "
              f"|diff|={abs(p['nf']-p['rm']):.3f}")
    z, active = _traj_z(res, mois)
    exp_max = statistics.NormalDist().inv_cdf(1 - 1 / (2 * z.size))
    print(f"  trajectory: {z.size} comparisons, max|z| = {z.max():.2f} "
          f"(expected {exp_max:.2f}), fraction |z|<3 = {(z < 3).mean():.3f}")
    print("  endpoint per observable (+ = NFsim high):")
    for o, (zs, comb) in _endpoint_z(res, mois).items():
        print(f"    {o:>12s} " + "".join(f"{v:7.2f}" for v in zs) + f"   combined {comb:7.2f}")
    return res


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "fig3"
    _dflt = {"fig3": FIG3_SEEDS, "fig6": FIG6_SEEDS, "agreement": AGREE_SEEDS}.get(mode, 40)
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else _dflt

    if mode == "fig3":
        fig3(n_seeds, sys.argv[3] if len(sys.argv) > 3 else "rm")
    elif mode == "fig6":
        fig6(n_seeds, sys.argv[3] if len(sys.argv) > 3 else "rm")
    elif mode == "agreement":
        mois = tuple(int(m) for m in sys.argv[3].split(",")) if len(sys.argv) > 3 \
            else AGREE_MOIS
        agreement(n_seeds, mois)
    elif mode == "rescore":
        rescore()
    elif mode == "symcheck":
        symcheck(n_seeds if len(sys.argv) > 2 else 60)
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
        print(f"unknown mode {mode!r}; "
              f"use fig3 | fig6 | agreement | rescore | symcheck | traj")
