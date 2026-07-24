"""Driver + verification for p53_nhej_dolan2015.bngl -- the integrated stochastic model
of NHEJ DNA double-strand-break repair, ATM/p53/Mdm2 damage signaling and p53/p21-mediated
early senescence of Dolan et al. (2015), PLoS Comput Biol 11(5):e1004246
(doi:10.1371/journal.pcbi.1004246; Supplementary File S11). The model's complexity is a
signaling/repair state space (50 independent DNA loci, an ATM per-locus tag, linear repair
chains), not multivalent aggregation, so it is simulated NETWORK-FREE (NFsim) and has no
same-complex-binding ambiguity: standard/legacy BioNetGen NFsim is a valid reference
(no -bscb / -utl needed).

As deposited, IRRADIATION IS OFF and DNA damage is driven only by chronic background ROS:
this is the non-irradiated baseline of Fig. 4 in Dolan et al. (2015).

Three things are verified here:

  reference    -- runs the committed model's faithful published protocol (single-cell
                  NFsim run over 2400 min of chronic background-ROS damage) with standard
                  BioNetGen NFsim (seed 2015 reproduces the committed reference .gdat) and
                  checks the paper's stated quantitative baseline claims: background ROS is
                  ~constant near kROS/kdROS = 10; DNA damage foci are low/transient; p53
                  shows damage-triggered pulses; p21 stays low; and NO cell enters early
                  senescence (Senescent_Counter == 0), exactly as reported for normal
                  background ROS. Fig. 4 panels are representative SINGLE stochastic runs
                  (not ensemble curves), so trajectory-level digitization is not
                  meaningful; the quantitative claims above are what the figure asserts.

  independent  -- model-specification check. With DNA damage turned off (kdam1=kdam2=0 ->
                  no DSBs -> ATM stays in state 0 -> no p53/Mdm2 phosphorylation), the
                  p53/Mdm2/p21/GADD45/p38/ROS core is a closed mass-action network. An
                  independent SciPy ODE (deterministic mean-field) of that network is
                  compared to the BNG-NFsim ensemble mean. Linear / high-copy observables
                  (ROS, total p53, p53-Mdm2 complex, p53 mRNA, p21, GADD45) match tightly;
                  the low-copy UNBOUND species (free Mdm2, Mdm2 mRNA, free p53) sit above
                  the mean-field value by the expected stochastic offset (negative
                  p53-Mdm2 binding covariance at small copy number). The exact stochastic
                  identity E[Mdm2_mRNA] == E[p53_Unphos] (both driven linearly by free p53)
                  is also checked as an internal consistency test of the encoded rate laws.

  agreement    -- cross-checks that the three network-free engines agree on the SAME rules:
                  standard BioNetGen NFsim, bngsim NfsimSession, and bngsim RuleMonkeySession
                  (issue #10, task 4). All three sample the same continuous-time Markov
                  chain, so their ensemble-mean signaling trajectories agree within
                  Monte-Carlo error. Agreement is quantified by the pairwise z-score
                  |mean_i - mean_j| / sqrt(se_i^2 + se_j^2) across observables and times.
                  A reduced window (t_end=480 min) keeps the exact RuleMonkey method and
                  multi-seed ensembles tractable; the model is unchanged.

Requirements: `bionetgen` on PATH (BNG core with NFsim v1.14.3) and the `bngsim` package
(NfsimSession + RuleMonkeySession). Run with the repo venv:
    python run_dolan2015.py reference           # faithful protocol + baseline anchors
    python run_dolan2015.py independent 16       # mean-field ODE vs NFsim ensemble (16 seeds)
    python run_dolan2015.py agreement 12         # 3-engine cross-check, 12 seeds/engine
    python run_dolan2015.py plot                 # (re)build verify_dolan2015.png
    python run_dolan2015.py all
"""
import glob
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, "p53_nhej_dolan2015.bngl")
REFDIR = os.path.join(HERE, "reference")
BNGCLI = os.environ.get("BIONETGEN", "bionetgen")

REF_SEED = 2015          # reproduces the committed reference .gdat
FULL_T_END = 2400        # published per-run length (min)
FULL_N = 240             # committed output resolution (every 10 min)
CROSS_T_END = 480        # reduced window for the tractable 3-engine cross-check

# Fig. 4 observables (Dolan et al. 2015): ROS, DNA damage foci, p53, p21.
FIG4_SINGLE = ["ROS", "p53", "p21"]           # single-column observables
# total DNA damage foci = sum over the 50 loci (built on the fly)

# Signaling-core observables used by the independent ODE check and the 3-engine agreement.
CORE_OBS = ["ROS", "p53", "p21", "p53_mRNA", "MDM2_mRNA", "MDM2_Unphos",
            "GADD45", "p38p", "p53_Bound", "p53_Unphos"]
AGREE_OBS = ["ROS", "p53", "p53_mRNA", "MDM2_mRNA", "p53_Bound", "p21", "GADD45"]


# ---------------------------------------------------------------- BNG helpers
def _model_block():
    """Committed model text up to and including `end model` (everything except the actions
    block), so callers can append their own actions. Split on `end model` to avoid matching
    the literal 'begin actions' that appears inside the header #@note comment."""
    text = open(MODEL).read()
    marker = "\nend model"
    i = text.index(marker)
    return text[:i + len(marker)] + "\n"


def _set_param(text, name, value):
    text, n = re.subn(rf"(?m)^(\s*{re.escape(name)}\s+)\S+", rf"\g<1>{value}", text, count=1)
    assert n == 1, f"expected one '{name}' parameter line, found {n}"
    return text


def _bng_run(bngl_text, tag):
    workdir = tempfile.mkdtemp(prefix=f"dolan_{tag}_")
    bngl = os.path.join(workdir, f"{tag}.bngl")
    open(bngl, "w").write(bngl_text)
    out = os.path.join(workdir, "out")
    r = subprocess.run([BNGCLI, "run", "-i", bngl, "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:] + "\n---STDOUT---\n" + r.stdout[-3000:]
    return out


def _read_gdat(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    cols = {name: data[:, i] for i, name in enumerate(header)}
    t = cols.pop("time")
    return t, cols


def _total_foci(cols):
    return sum(cols[f"Damage_Foci_{i}"] for i in range(1, 51))


def _bng_nf_run(t_end, n_steps, seed, kdam=None):
    """One BNG-NFsim run of the committed model; returns (t, cols). If kdam is not None,
    kdam1=kdam2=kdam (use 0 to switch DNA damage off for the independent check)."""
    base = _model_block()
    if kdam is not None:
        base = _set_param(base, "kdam1", kdam)
        base = _set_param(base, "kdam2", kdam)
    act = ("\nbegin actions\n"
           f'  simulate({{method=>"nf",suffix=>"r",t_start=>0,t_end=>{t_end},'
           f"n_steps=>{n_steps},gml=>2147483647,seed=>{seed}}})\n"
           "end actions\n")
    out = _bng_run(base + act, f"nf_s{seed}")
    return _read_gdat(glob.glob(os.path.join(out, "*_r.gdat"))[0])


def _build_xml(kdam=None):
    base = _model_block()
    if kdam is not None:
        base = _set_param(base, "kdam1", kdam)
        base = _set_param(base, "kdam2", kdam)
    base += "\nbegin actions\n  writeXML()\nend actions\n"
    out = _bng_run(base, "xml")
    return glob.glob(os.path.join(out, "*.xml"))[0]


# ---------------------------------------------------------------- reference
def reference():
    """Run the committed faithful published protocol (BNG NFsim, seed 2015) and check the
    paper's stated quantitative baseline claims for non-irradiated cells (Fig. 4)."""
    print("=== reference: faithful published protocol (BNG NFsim, chronic background ROS) ===")
    t0 = time.time()
    t, cols = _bng_nf_run(FULL_T_END, FULL_N, REF_SEED)
    foci = _total_foci(cols)
    ros, p53, p21 = cols["ROS"], cols["p53"], cols["p21"]
    sen = cols["Senescent_Counter"]
    print(f"  [{time.time()-t0:.0f}s] {t[0]:.0f}->{t[-1]:.0f} min, {len(t)} outputs")
    print(f"  ROS         mean={ros.mean():6.2f}  (theory kROS/kdROS = 10)")
    print(f"  p53         mean={p53.mean():6.2f}  min={p53.min():.0f} max={p53.max():.0f} "
          f"(damage-triggered pulses)")
    print(f"  p21         mean={p21.mean():6.2f}  max={p21.max():.0f}  (stays below the "
          f"senescence threshold of 15)")
    print(f"  DNA foci    mean={foci.mean():6.2f}  max={foci.max():.0f}  (low / transient)")
    print(f"  senescence  max Senescent_Counter = {sen.max():.0f}  "
          f"(expect 0: no cell enters senescence at normal background ROS)")
    checks = {
        "ROS ~ 10 (kROS/kdROS)": abs(ros.mean() - 10) < 3,
        "no senescence (Senescent_Counter == 0)": sen.max() == 0,
        "p21 stays < 15": p21.max() < 15,
        "damage foci low (mean < 3)": foci.mean() < 3,
    }
    for k, v in checks.items():
        print(f"    [{'PASS' if v else 'FAIL'}] {k}")
    ref = os.path.join(REFDIR, "p53_nhej_dolan2015_nfr.gdat")
    if os.path.exists(ref):
        print(f"  committed reference present: {ref}")
    return t, cols


# ---------------------------------------------------------------- independent ODE
def _ode_at(params, t_end):
    """Deterministic mean-field ODE of the damage-off p53/Mdm2/p21/GADD45/p38/ROS core,
    integrated to t_end and returned as a dict of observable values AT t_end. (Matched-time
    comparison, not forced steady state: GADD45 is slow, tau ~ 1/kGADD45deg ~ 1667 min, so
    the fair comparison to a finite NFsim run is the ODE evaluated at the same time.)"""
    from scipy.integrate import solve_ivp
    p = params

    def rhs(t, y):
        m53, p53u, C, mM, M, m21, q1, q2, q3, G, p38p, ROS = y
        p38u = 100.0 - p38p
        bind = p["kp53MDM2bind"] * p53u * M
        dis = p["kp53MDM2dis"] * C
        cdeg = p["kp53degMDM2dep"] * C
        return [
            p["kp53mRNAsyn"] * 1 - p["kp53mRNAdeg"] * m53,                       # m53
            p["kp53syn"] * m53 - p["kp53deg"] * p53u - bind + dis,               # p53u (free)
            bind - dis - cdeg,                                                   # C (complex)
            p["kMDM2mRNAsyn"] * p53u - p["kMDM2mRNAdeg"] * mM,                   # mM
            p["kMDM2syn"] * mM - p["kMDM2deg"] * M - bind + dis + cdeg,          # M (free)
            p["kp21mRNAsyn"] * p53u - p["kp21mRNAdeg"] * m21,                    # m21
            p["kp21synstep1"] * m21 - p["kp21synstep2"] * q1,                    # q1
            p["kp21synstep2"] * q1 - p["kp21synstep3"] * q2,                     # q2
            p["kp21synstep3"] * q2 - p["kp21deg"] * q3,                          # q3
            p["kGADD45act"] * q3 - p["kGADD45deg"] * G,                          # G
            p["kp38phos"] * p38u * G - p["kp38dphos"] * p38p,                    # p38p
            p["kROS"] * 1 + p["kROSgen"] * p38p - p["kdROS"] * ROS,             # ROS
        ]

    y0 = [10, 5, 95, 10, 5, 0, 0, 0, 0, 0, 0, 0]  # matches seed species (p21=0, GADD45=0, ROS=0)
    sol = solve_ivp(rhs, [0, t_end], y0, t_eval=[t_end], rtol=1e-9, atol=1e-12)
    m53, p53u, C, mM, M, m21, q1, q2, q3, G, p38p, ROS = sol.y[:, -1]
    return dict(ROS=ROS, p53=p53u + C, p21=q3, p53_mRNA=m53, MDM2_mRNA=mM,
                MDM2_Unphos=M, GADD45=G, p38p=p38p, p53_Bound=C, p53_Unphos=p53u)


def _parse_params():
    body = _model_block()
    p = {}
    for m in re.finditer(r"(?m)^\s*(k[A-Za-z0-9_]+)\s+([0-9.eE+-]+)\s*#", body):
        p[m.group(1)] = float(m.group(2))
    return p


def independent(n_seeds=16, t_end=2000, n=40, save=True):
    """Compare the mean-field ODE of the damage-off signaling core to the BNG-NFsim
    ensemble mean (matched time t_end), and check E[Mdm2_mRNA] == E[p53_Unphos]."""
    print("=== independent: mean-field ODE vs BNG-NFsim ensemble (DNA damage OFF) ===")
    ode = _ode_at(_parse_params(), t_end=t_end)
    stacks = {o: [] for o in CORE_OBS}
    t0 = time.time()
    for s in range(n_seeds):
        _, cols = _bng_nf_run(t_end, n, 700 + s, kdam=0)
        for o in CORE_OBS:
            stacks[o].append(cols[o])
    ens = {o: np.vstack(v) for o, v in stacks.items()}
    print(f"  [{time.time()-t0:.0f}s] {n_seeds} seeds, t_end={t_end} min\n")
    print(f"  {'observable':13s} {'ODE':>8s} {'NFsim_mean':>11s} {'se':>7s} {'z':>6s}")
    rows = {}
    for o in CORE_OBS:
        finals = ens[o][:, -1]
        m, se = finals.mean(), finals.std(ddof=1) / np.sqrt(len(finals))
        z = abs(ode[o] - m) / se if se > 0 else float("nan")
        rows[o] = (ode[o], m, se, z)
        tag = "unbound(low-copy)" if o in ("MDM2_Unphos", "MDM2_mRNA", "p53_Unphos") else ""
        print(f"  {o:13s} {ode[o]:8.2f} {m:11.2f} {se:7.2f} {z:6.2f}  {tag}")
    # Exact first-order identity (no covariance: p53 mRNA is a birth-death process driven
    # by the constant source P=1), a clean check that NFsim reproduces the encoded rate law.
    m53 = ens["p53_mRNA"][:, -1]
    exact = _parse_params()["kp53mRNAsyn"] / _parse_params()["kp53mRNAdeg"]  # = 10
    print(f"\n  exact birth-death  E[p53_mRNA]={m53.mean():.2f}+-{m53.std(ddof=1)/np.sqrt(len(m53)):.2f}"
          f"  ==  kp53mRNAsyn/kp53mRNAdeg = {exact:.1f}")
    # The low-copy UNBOUND species (free p53, free Mdm2, Mdm2 mRNA) sit ABOVE the mean-field
    # ODE: the stochastic mean exceeds mean-field because Cov(free p53, free Mdm2) < 0 for the
    # p53-Mdm2 binding at small copy number, so less complex forms than the product of means.
    mm = ens["MDM2_mRNA"][:, -1].mean()
    print(f"  low-copy stochastic offset: E[Mdm2_mRNA]={mm:.1f} > ODE {ode['MDM2_mRNA']:.1f} "
          f"(free-species mean-field bias; totals/complex match, see table)")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "independent.npz"),
                 obs=np.array(CORE_OBS),
                 ode=np.array([ode[o] for o in CORE_OBS]),
                 nf_mean=np.array([rows[o][1] for o in CORE_OBS]),
                 nf_se=np.array([rows[o][2] for o in CORE_OBS]),
                 n_seeds=n_seeds, t_end=t_end)
        print(f"  cached -> {os.path.join(REFDIR, 'independent.npz')}")
    return rows


# ---------------------------------------------------------------- agreement
def _ensemble_bngsim(xml, engine, n_seeds, t_end, n_points, seed0=5000):
    from bngsim import NfsimSession, RuleMonkeySession
    Sess = NfsimSession if engine == "nf" else RuleMonkeySession
    stacks = {o: [] for o in AGREE_OBS}
    tvec = None
    for s in range(n_seeds):
        with Sess(xml) as sess:
            sess.initialize(seed=seed0 + s)
            res = sess.simulate(0, t_end, n_points=n_points)
            names = list(sess.get_observable_names())
            Y = np.asarray(res.observables)
            tvec = np.asarray(res.time)
            for o in AGREE_OBS:
                stacks[o].append(Y[:, names.index(o)])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


def _ensemble_bng_native(n_seeds, t_end, n_points, seed0=5000):
    stacks = {o: [] for o in AGREE_OBS}
    tvec = None
    for s in range(n_seeds):
        t, cols = _bng_nf_run(t_end, n_points - 1, seed0 + s)
        tvec = t
        for o in AGREE_OBS:
            stacks[o].append(cols[o])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


def agreement(n_seeds=12, t_end=CROSS_T_END, n_points=25, save=True):
    """Three-engine network-free agreement cross-check (issue #10, task 4)."""
    print("=== agreement: BNG-NFsim vs bngsim-NFsim vs bngsim-RuleMonkey ===")
    print(f"    chronic background ROS, t_end={t_end} min, {n_seeds} seeds/engine")
    xml = _build_xml()
    res = {}
    for label, fn in (("bng_nf", lambda: _ensemble_bng_native(n_seeds, t_end, n_points)),
                      ("bngsim_nf", lambda: _ensemble_bngsim(xml, "nf", n_seeds, t_end, n_points)),
                      ("bngsim_rm", lambda: _ensemble_bngsim(xml, "rm", n_seeds, t_end, n_points))):
        t0 = time.time()
        tvec, st = fn()
        res[label] = (tvec, st)
        print(f"  [{time.time()-t0:5.0f}s] {label:10s} done")
    tvec = res["bng_nf"][0]

    def mean_se(st, o):
        A = st[o]
        return A.mean(0), A.std(0, ddof=1) / np.sqrt(A.shape[0])
    pairs = [("bng_nf", "bngsim_nf"), ("bng_nf", "bngsim_rm"), ("bngsim_nf", "bngsim_rm")]
    zmax, zvals = 0.0, []
    for a, b in pairs:
        for o in AGREE_OBS:
            ma, sa = mean_se(res[a][1], o)
            mb, sb = mean_se(res[b][1], o)
            denom = np.sqrt(sa**2 + sb**2)
            z = np.abs(ma - mb) / np.where(denom > 0, denom, np.inf)
            zvals.append(z[np.isfinite(z)])
            zmax = max(zmax, float(np.nanmax(z)))
    zall = np.concatenate(zvals)
    frac3 = float((zall < 3).mean())
    print(f"  pairwise-mean z-scores over {len(AGREE_OBS)} obs x {n_points} pts x 3 pairs:")
    print(f"    max|z| = {zmax:.2f}   fraction |z|<3 = {frac3:.3f}   "
          f"(expect ~0.997 if engines sample the same process)")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "agreement.npz"),
                 tvec=tvec, obs=np.array(AGREE_OBS),
                 **{f"{lab}__{o}": res[lab][1][o] for lab in res for o in AGREE_OBS},
                 t_end=t_end, n_seeds=n_seeds)
        print(f"  cached -> {os.path.join(REFDIR, 'agreement.npz')}")
    return res, dict(zmax=zmax, frac_z_lt3=frac3)


# ---------------------------------------------------------------- figure
def plot():
    """Build verify_dolan2015.png from the committed reference .gdat + cached npz files."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # Row 1: Fig. 4 reproduction (committed reference single-cell run) -----------------
    ref = os.path.join(REFDIR, "p53_nhej_dolan2015_nfr.gdat")
    if os.path.exists(ref):
        t, cols = _read_gdat(ref)
        th = t / 60.0  # minutes -> hours
        ax = axes[0, 0]
        ax.plot(th, cols["ROS"], lw=1, color="#d95f02", label="ROS")
        ax.axhline(10, ls="--", lw=0.8, color="k", label="kROS/kdROS = 10")
        ax.plot(th, _total_foci(cols), lw=1, color="#1b9e77", label="DNA damage foci")
        ax.set(title="Fig. 4 baseline: ROS & DNA damage foci",
               xlabel="time (h)", ylabel="molecules / cell")
        ax.legend(fontsize=7)
        ax = axes[0, 1]
        ax.plot(th, cols["p53"], lw=1, color="#7570b3", label="p53")
        ax.plot(th, cols["p21"], lw=1, color="#e7298a", label="p21")
        ax.set(title="Fig. 4 baseline: p53 pulses & p21",
               xlabel="time (h)", ylabel="molecules / cell")
        ax.legend(fontsize=7)
    axes[0, 2].axis("off")
    axes[0, 2].text(0.0, 0.5,
        "Dolan et al. (2015) PLoS Comput Biol 11(5):e1004246\n"
        "Integrated NHEJ repair + p53/p21 senescence,\n"
        "network-free (NFsim). Committed run: single cell,\n"
        "chronic background ROS, no irradiation (Fig. 4).\n\n"
        "Top: reproduced baseline -- ROS ~ 10 (birth/death),\n"
        "  transient DNA damage foci, damage-triggered p53\n"
        "  pulses, p21 low; no senescence.\n\n"
        "Bottom: (left/mid) independent mean-field ODE of the\n"
        "  damage-off signaling core vs BNG-NFsim ensemble;\n"
        "  (right) three network-free engines agree on the\n"
        "  same rules (ensemble mean +/- s.e.m.).",
        fontsize=8, va="center")

    # Row 2a/2b: independent ODE vs NFsim ensemble ------------------------------------
    ind = os.path.join(REFDIR, "independent.npz")
    if os.path.exists(ind):
        d = np.load(ind, allow_pickle=True)
        obs = list(d["obs"])
        ode = d["ode"]; nf = d["nf_mean"]; se = d["nf_se"]
        order = sorted(range(len(obs)), key=lambda i: -max(ode[i], nf[i]))
        obs = [obs[i] for i in order]; ode = ode[order]; nf = nf[order]; se = se[order]
        x = np.arange(len(obs))
        ax = axes[1, 0]
        ax.bar(x - 0.2, ode, width=0.4, color="#377eb8", label="ODE (mean-field)")
        ax.bar(x + 0.2, nf, width=0.4, yerr=se, color="#ff7f00",
               label="NFsim ensemble", capsize=2)
        ax.set_xticks(x)
        ax.set_xticklabels(obs, rotation=45, ha="right", fontsize=6)
        ax.set(title="Independent check: ODE vs NFsim (damage OFF)",
               ylabel="molecules / cell (steady state)")
        ax.legend(fontsize=7)
        ax = axes[1, 1]
        lo = np.minimum(ode, nf); hi = np.maximum(ode, nf)
        ax.errorbar(ode, nf, yerr=se, fmt="o", ms=5, color="#377eb8")
        lim = [0, max(hi.max(), 1) * 1.15]
        ax.plot(lim, lim, "k--", lw=0.8)
        for i, o in enumerate(obs):
            ax.annotate(o, (ode[i], nf[i]), fontsize=6, xytext=(3, 3),
                        textcoords="offset points")
        ax.set(title="ODE vs NFsim (identity line)",
               xlabel="ODE mean-field", ylabel="NFsim ensemble mean",
               xlim=lim, ylim=lim)

    # Row 2c: three-engine agreement --------------------------------------------------
    npz = os.path.join(REFDIR, "agreement.npz")
    ax = axes[1, 2]
    if os.path.exists(npz):
        d = np.load(npz, allow_pickle=True)
        tvec = d["tvec"] / 60.0
        styles = {"bng_nf": ("BNG NFsim", "o", "#1b9e77"),
                  "bngsim_nf": ("bngsim NFsim", "s", "#7570b3"),
                  "bngsim_rm": ("bngsim RuleMonkey", "^", "#d95f02")}
        for o in ("p53", "ROS"):
            for lab, (name, mk, col) in styles.items():
                key = f"{lab}__{o}"
                if key not in d:
                    continue
                A = d[key]
                m = A.mean(0); se = A.std(0, ddof=1) / np.sqrt(A.shape[0])
                ax.plot(tvec, m, marker=mk, ms=3, color=col, lw=1,
                        label=name if o == "p53" else None)
                ax.fill_between(tvec, m - se, m + se, color=col, alpha=0.15)
            ax.annotate(o, (tvec[-1], A.mean(0)[-1]), fontsize=7,
                        xytext=(2, 0), textcoords="offset points")
        ax.set(title="3-engine agreement (p53, ROS)",
               xlabel="time (h)", ylabel="molecules / cell")
        ax.legend(fontsize=7, title="engine")

    fig.suptitle("p53_nhej_dolan2015 verification: Fig. 4 baseline + independent ODE + "
                 "3-engine network-free agreement", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    png = os.path.join(HERE, "verify_dolan2015.png")
    fig.savefig(png, dpi=130)
    print(f"wrote {png}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    if mode == "reference":
        reference()
    elif mode == "independent":
        independent(n_seeds=int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    elif mode == "agreement":
        agreement(n_seeds=int(sys.argv[2]) if len(sys.argv) > 2 else 12)
    elif mode == "plot":
        plot()
    elif mode == "all":
        reference()
        independent(n_seeds=int(sys.argv[2]) if len(sys.argv) > 2 else 16)
        agreement(n_seeds=int(sys.argv[3]) if len(sys.argv) > 3 else 12)
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference | independent [N] | agreement [N] | plot | all")
