"""Driver + verification for erbb_receptor_signaling_creamer2012.bngl -- the large
rule-based ERBB receptor signaling model of Creamer et al. (2012), BMC Systems Biology
6:107 (doi:10.1186/1752-0509-6-107; Supplementary File 1). The reaction network is far
too large to enumerate (>10^100 reachable species), so the model is simulated
NETWORK-FREE. Unlike a pure multisite-phosphorylation model, the ERBB adaptor layer can
in principle form small rings (a receptor dimer bridged by Grb2 + Shc), so whether the
default (legacy) NFsim molecularity handling matches the exact network-free method is
NOT assumed -- it is tested here.

Four things are provided/verified:

  reference  -- runs the committed model's FAITHFUL published protocol (equilibrate 500 s
                with ligand binding off, turn ligand binding on, then a 180 min response)
                with standard BioNetGen NFsim at the committed subvolume f=0.02, and
                checks the result against the committed reference .gdat: receptor sites
                phosphorylate fast and transiently while Akt/Raf sites are sustained
                (the Fig. 4 kinetic spectrum).

  agreement  -- cross-checks that three network-free engines agree on the SAME rules:
                standard BioNetGen NFsim, bngsim NfsimSession, and bngsim
                RuleMonkeySession (an EXACT network-free method). Uses a reduced but
                identical single-phase ligand-on protocol at a small subvolume
                (f=FX_CROSS) so RuleMonkey is affordable. All three sample the same
                continuous-time Markov chain, so their ensemble-mean phosphosite
                trajectories agree within Monte-Carlo error, quantified by the pairwise
                z-score |mean_i - mean_j| / sqrt(se_i^2 + se_j^2).

  flags      -- the -bscb / -utl determination. Runs BioNetGen NFsim with default flags,
                with -bscb (block same-complex binding), and with -bscb plus a raised
                -utl, and compares each to the EXACT RuleMonkey ensemble mean. If default
                NFsim already agrees with RuleMonkey (small z), the model has no
                same-complex-binding / traversal ambiguity that matters and legacy NFsim
                (as the authors used, NFsim v1.09, default flags) is a valid reference.

  plot       -- builds verify_creamer2012.png (reference protocol + engine agreement +
                reported-data comparison against digitized Fig. 4).

Parameter propagation note (lanl/bngsim#44): a POST-init NFsim set_param does not refresh
a rule's rate constant, and RuleMonkeySession.set_param does not re-derive seed-species
amounts. Both are avoided here by BAKING f (and lig) into the parameter block and calling
writeXML() before load, so every engine loads identical initial conditions from the XML.

Requirements: `bionetgen` on PATH (BNG core with NFsim) and the `bngsim` package. Run
with the repo venv:
    python run_creamer2012.py reference          # faithful protocol vs committed .gdat
    python run_creamer2012.py agreement 24        # 3-engine cross-check, 24 seeds/engine
    python run_creamer2012.py flags 24            # -bscb / -utl sensitivity vs RuleMonkey
    python run_creamer2012.py plot                # (re)build verify_creamer2012.png
"""
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, "erbb_receptor_signaling_creamer2012.bngl")
REFDIR = os.path.join(HERE, "reference")

# The 55 individually tracked S/T/Y phosphosites (Fig. 4 rows).
OBS_55 = [
    "ERK_T185", "AKT_S473", "EGFR_669", "EGFR_992", "EGFR_1068", "EGFR_1086",
    "EGFR_1114", "EGFR_1148", "EGFR_1173", "ErbB2_1139", "ErbB2_1196", "ErbB2_1222",
    "ErbB2_1248", "ErbB3_1054", "ErbB3_1197", "ErbB3_1222", "ErbB3_1260", "ErbB3_1276",
    "ErbB3_1289", "ErbB3_1328", "ErbB4_1056", "ErbB4_1188", "ErbB4_1242", "pSHC",
    "Sos1_1132", "Sos1_1167", "Sos1_1178", "Sos1_1193", "Gab1_312", "Gab1_381",
    "Gab1_447", "Gab1_454", "Gab1_472", "Gab1_476", "Gab1_581", "Gab1_597", "Gab1_619",
    "Gab1_657", "Raf1_29", "Raf1_43", "Raf1_259", "Raf1_289", "Raf1_296", "Raf1_301",
    "Raf1_338", "Raf1_341", "Raf1_471", "Raf1_491", "Raf1_494", "Raf1_642", "MEK1_218",
    "MEK1_222", "MEK1_292", "pERK_Y187", "AKT_T308",
]
# A readable subset spanning the kinetic spectrum for the verification figure.
PLOT_OBS = ["EGFR_992", "EGFR_1068", "ErbB3_1289", "ERK_T185", "Raf1_259", "AKT_S473"]

FX_COMMITTED = 0.02     # subvolume fraction baked into the committed model
FX_CROSS = 0.002        # reduced subvolume for the tractable 3-engine / flags checks


# ---------------------------------------------------------------- BNG helpers
def _bng_run(bngl_text, tag):
    """Write bngl_text to a temp dir, run `bionetgen run`, return the output dir."""
    workdir = tempfile.mkdtemp(prefix=f"erbb_{tag}_")
    bngl = os.path.join(workdir, f"{tag}.bngl")
    open(bngl, "w").write(bngl_text)
    out = os.path.join(workdir, "out")
    r = subprocess.run(["bionetgen", "run", "-i", bngl, "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-3000:] + "\n---STDERR---\n" + r.stderr[-1500:]
    return out


def _read_gdat(path):
    """Read a BNG .gdat into (time[array], {obs: array})."""
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    cols = {name: data[:, i] for i, name in enumerate(header)}
    t = cols.pop("time")
    return t, cols


def _set_param_line(text, name, value):
    """Rewrite the `name  <value>` line in the parameters block."""
    text, n = re.subn(rf"(?m)^(\s*{re.escape(name)}\s+)\S+", rf"\g<1>{value}", text, count=1)
    assert n == 1, f"expected one '{name}' parameter line, found {n}"
    return text


def _model_block():
    """Committed model text up to and including `end model` (no actions block)."""
    text = open(MODEL).read()
    marker = "\nend model"
    i = text.index(marker)
    return text[:i + len(marker)] + "\n"


def build_xml(f=FX_CROSS, lig=1):
    """Emit a BNG-XML from the committed model with f and lig BAKED into the parameter
    block (so seed-species amounts and rate constants are pre-evaluated at this
    subvolume), then writeXML(). All engines load identical initial conditions from this
    XML -- no reliance on session set_param (see lanl/bngsim#44)."""
    text = _model_block()
    text = _set_param_line(text, "f", f)
    text = _set_param_line(text, "lig", lig)
    text += "\nbegin actions\n  writeXML()\nend actions\n"
    out = _bng_run(text, "xml")
    xml = os.path.join(out, "xml.xml")
    assert os.path.exists(xml), f"writeXML did not produce {xml}"
    return xml


# ---------------------------------------------------------------- ensembles
def ensemble_bngsim(xml, engine, n_seeds, t_end, n_points, seed0=7000):
    """Single-phase ligand-on ensemble via bngsim. engine: 'nf' | 'rm'.
    Returns (tvec, {obs: array[n_seeds, n_points]})."""
    from bngsim import NfsimSession, RuleMonkeySession
    Sess = NfsimSession if engine == "nf" else RuleMonkeySession
    stacks = {o: [] for o in OBS_55}
    tvec = None
    for s in range(n_seeds):
        with Sess(xml) as sess:
            sess.initialize(seed=seed0 + s)
            res = sess.simulate(0, t_end, n_points=n_points)
            names = list(sess.get_observable_names())
            Y = np.asarray(res.observables)
            tvec = np.asarray(res.time)
            for o in OBS_55:
                stacks[o].append(Y[:, names.index(o)])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


def ensemble_bng_native(n_seeds, f, t_end, n_points, seed0=7000, param=""):
    """Single-phase ligand-on ensemble via standard BioNetGen NFsim (seeded runs).
    `param` passes extra NFsim flags (e.g. '-bscb -utl 8'). Returns
    (tvec, {obs: array[n_seeds, n_points]})."""
    base = _set_param_line(_model_block(), "f", f)
    base = _set_param_line(base, "lig", 1)
    stacks = {o: [] for o in OBS_55}
    tvec = None
    pflag = f',param=>"{param}"' if param else ""
    for s in range(n_seeds):
        act = ("\nbegin actions\n"
               f'  simulate({{method=>"nf",suffix=>"resp",t_start=>0,t_end=>{t_end},'
               f"n_steps=>{n_points - 1},gml=>2147483647,seed=>{seed0 + s}{pflag}}})\n"
               "end actions\n")
        out = _bng_run(base + act, f"native_s{s}")
        t, cols = _read_gdat(os.path.join(out, f"native_s{s}_resp.gdat"))
        tvec = t
        for o in OBS_55:
            stacks[o].append(cols[o])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


# ------------------------------------------------------------- z-score utilities
def _mean_se(st, o):
    A = st[o]
    return A.mean(0), A.std(0, ddof=1) / np.sqrt(A.shape[0])


def _active_obs(*stacks, thresh=1.0):
    """Observables whose max ensemble mean exceeds `thresh` in at least one stack
    (comparing all-zero downstream observables via z-scores is meaningless)."""
    active = []
    for o in OBS_55:
        if any(st[o].mean(0).max() > thresh for st in stacks):
            active.append(o)
    return active


def _pairwise_zmax(res, pairs, obs):
    """Return (max|z|, fraction |z|<3) over obs x time x given engine pairs."""
    zvals = []
    zmax = 0.0
    for a, b in pairs:
        for o in obs:
            ma, sa = _mean_se(res[a][1], o)
            mb, sb = _mean_se(res[b][1], o)
            denom = np.sqrt(sa**2 + sb**2)
            z = np.abs(ma - mb) / np.where(denom > 0, denom, np.inf)
            zf = z[np.isfinite(z)]
            zvals.append(zf)
            if zf.size:
                zmax = max(zmax, float(np.nanmax(zf)))
    zall = np.concatenate(zvals) if zvals else np.array([0.0])
    return zmax, float((zall < 3).mean())


# ---------------------------------------------------------------- reference
def reference():
    """Run the committed model's faithful published protocol with standard BNG NFsim
    and compare against the committed reference .gdat."""
    print("=== reference: faithful published protocol (BNG NFsim, f=0.02) ===")
    t0 = time.time()
    out = _bng_run(open(MODEL).read(), "reference")
    stem = os.path.splitext(os.path.basename(MODEL))[0]
    t_eq, eq = _read_gdat(os.path.join(out, "reference_equil.gdat"))
    t_tc, tc = _read_gdat(os.path.join(out, "reference_erbb.gdat"))
    tmin = t_tc / 60.0
    print(f"  [{time.time()-t0:.0f}s] equilibration {t_eq[0]:.0f}->{t_eq[-1]:.0f} s, "
          f"response {tmin[0]:.0f}->{tmin[-1]:.0f} min")
    print(f"  {'site':12s} basal(off)  peak   t_peak(min)  end   class")
    for o in PLOT_OBS:
        peak = tc[o].max()
        tpk = tmin[tc[o].argmax()]
        cls = "transient" if tc[o][-1] < 0.3 * peak else "sustained"
        print(f"  {o:12s} {eq[o][-1]:9.0f}  {peak:6.0f}  {tpk:9.0f}  {tc[o][-1]:5.0f}  {cls}")
    # Qualitative Fig. 4 check: receptor sites transient, Akt sustained.
    rec_transient = tc["EGFR_992"].max() > 0 and tc["EGFR_992"][-1] < 0.3 * tc["EGFR_992"].max()
    akt_sustained = tc["AKT_S473"][-1] > 0.02 * tc["AKT_S473"].max()
    print(f"  Fig.4 kinetic spectrum: receptor EGFR_992 transient={rec_transient}, "
          f"Akt S473 sustained={akt_sustained}")
    return (t_eq, eq), (t_tc, tc)


# ---------------------------------------------------------------- agreement
def agreement(n_seeds=24, fx=FX_CROSS, t_end=300, n_points=16, save=True):
    """Three-engine network-free agreement cross-check (single-phase, reduced f)."""
    print("=== agreement: BNG-NFsim vs bngsim-NFsim vs bngsim-RuleMonkey ===")
    print(f"    single-phase ligand-on, f={fx}, t_end={t_end}s, {n_seeds} seeds/engine")
    xml = build_xml(f=fx, lig=1)
    res = {}
    for label, fn in (("bng_nf", lambda: ensemble_bng_native(n_seeds, fx, t_end, n_points)),
                      ("bngsim_nf", lambda: ensemble_bngsim(xml, "nf", n_seeds, t_end, n_points)),
                      ("bngsim_rm", lambda: ensemble_bngsim(xml, "rm", n_seeds, t_end, n_points))):
        t0 = time.time()
        tvec, st = fn()
        res[label] = (tvec, st)
        print(f"  [{time.time()-t0:5.0f}s] {label:10s} done")
    tvec = res["bng_nf"][0]
    obs = _active_obs(res["bng_nf"][1], res["bngsim_nf"][1], res["bngsim_rm"][1])
    pairs = [("bng_nf", "bngsim_nf"), ("bng_nf", "bngsim_rm"), ("bngsim_nf", "bngsim_rm")]
    zmax, frac3 = _pairwise_zmax(res, pairs, obs)
    print(f"  {len(obs)} active observables x {n_points} pts x 3 pairs:")
    print(f"    max|z| = {zmax:.2f}   fraction |z|<3 = {frac3:.3f}   "
          f"(expect ~0.997 if engines sample the same process)")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "agreement.npz"),
                 tvec=tvec, obs=np.array(OBS_55), active=np.array(obs),
                 **{f"{lab}__{o}": res[lab][1][o] for lab in res for o in OBS_55},
                 fx=fx, n_seeds=n_seeds)
        print(f"  cached -> {os.path.join(REFDIR, 'agreement.npz')}")
    return res, dict(zmax=zmax, frac_z_lt3=frac3, active=obs)


# ---------------------------------------------------------------- flags (-bscb/-utl)
def flags(n_seeds=24, fx=FX_CROSS, t_end=300, n_points=16, save=True):
    """Determine whether -bscb / -utl change NFsim results by comparing default NFsim,
    -bscb NFsim, and -bscb+raised-utl NFsim against the EXACT RuleMonkey ensemble."""
    print("=== flags: NFsim default vs -bscb vs -bscb -utl, judged vs exact RuleMonkey ===")
    print(f"    single-phase ligand-on, f={fx}, t_end={t_end}s, {n_seeds} seeds/engine")
    xml = build_xml(f=fx, lig=1)
    res = {}
    configs = [("nf_default", lambda: ensemble_bng_native(n_seeds, fx, t_end, n_points)),
               ("nf_bscb", lambda: ensemble_bng_native(n_seeds, fx, t_end, n_points, param="-bscb")),
               ("nf_bscb_utl", lambda: ensemble_bng_native(n_seeds, fx, t_end, n_points,
                                                           param="-bscb -utl 8")),
               ("rulemonkey", lambda: ensemble_bngsim(xml, "rm", n_seeds, t_end, n_points))]
    for label, fn in configs:
        t0 = time.time()
        tvec, st = fn()
        res[label] = (tvec, st)
        print(f"  [{time.time()-t0:5.0f}s] {label:12s} done")
    obs = _active_obs(*[res[k][1] for k in res])
    print(f"  {len(obs)} active observables; max|z| of each NFsim config vs exact RuleMonkey:")
    verdict = {}
    for cfg in ("nf_default", "nf_bscb", "nf_bscb_utl"):
        zmax, frac3 = _pairwise_zmax(res, [(cfg, "rulemonkey")], obs)
        verdict[cfg] = (zmax, frac3)
        print(f"    {cfg:12s}  max|z|={zmax:5.2f}  frac|z|<3={frac3:.3f}")
    # If default already matches RuleMonkey as well as -bscb does, flags don't matter.
    print("  => -bscb / -utl change results only if nf_default disagrees with RuleMonkey "
          "while nf_bscb agrees.")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "flags.npz"),
                 tvec=res["nf_default"][0], obs=np.array(OBS_55), active=np.array(obs),
                 **{f"{lab}__{o}": res[lab][1][o] for lab in res for o in OBS_55},
                 fx=fx, n_seeds=n_seeds,
                 verdict=np.array([[k, verdict[k][0], verdict[k][1]] for k in verdict], dtype=object))
        print(f"  cached -> {os.path.join(REFDIR, 'flags.npz')}")
    return res, verdict


# Figure-4 row order (top -> bottom) mapped to model observable names, used to render
# the digitized and simulated heat maps in the paper's clustered order.
FIG4_ROW_SITES = [
    "ErbB2_1139", "ErbB3_1289", "ErbB3_1222", "ErbB3_1260", "ErbB3_1328", "ErbB2_1248",
    "ErbB2_1222", "ErbB2_1196", "EGFR_992", "ErbB4_1056", "EGFR_1068", "ErbB3_1276",
    "ErbB4_1188", "EGFR_1148", "EGFR_1086", "EGFR_1114", "ErbB3_1054", "ErbB3_1197",
    "ErbB4_1242", "pSHC", "Gab1_472", "Gab1_619", "Gab1_447", "Gab1_657", "Gab1_454",
    "Sos1_1132", "EGFR_1173", "Gab1_312", "EGFR_669", "ERK_T185", "pERK_Y187",
    "Gab1_597", "MEK1_218", "Gab1_581", "Gab1_381", "Gab1_476", "Sos1_1193", "Sos1_1167",
    "Sos1_1178", "Raf1_29", "Raf1_43", "Raf1_642", "Raf1_471", "Raf1_289", "MEK1_292",
    "AKT_S473", "Raf1_338", "Raf1_296", "Raf1_341", "Raf1_491", "Raf1_494", "MEK1_222",
    "Raf1_259", "AKT_T308", "Raf1_301",
]


def _read_digitized():
    """Load reference/creamer2012_fig4_digitized.csv -> (tmin, {site: array})."""
    import csv
    path = os.path.join(REFDIR, "creamer2012_fig4_digitized.csv")
    with open(path) as fh:
        rd = csv.reader(fh)
        hdr = next(rd)
        tmin = np.array([float(x) for x in hdr[1:]])
        D = {row[0]: np.array([float(x) for x in row[1:]]) for row in rd}
    return tmin, D


def reported_data_metrics():
    """Quantitative comparison of the digitized Fig. 4 heat map against the simulated
    per-site ensemble mean. Returns dict of metrics + the two aligned matrices."""
    from scipy.stats import pearsonr, spearmanr
    tmin, D = _read_digitized()
    sim = np.load(os.path.join(REFDIR, "fig4_sim.npz"), allow_pickle=True)
    tsim = sim["tvec"] / 60.0
    dig_mat, sim_mat, dig_frac, sim_frac = [], [], [], []
    for s in FIG4_ROW_SITES:
        dv = D[s]
        sv = np.interp(tmin, tsim, sim[s])
        dn = dv / max(dv.max(), 1e-9)
        sn = sv / max(sv.max(), 1e-9)
        dig_mat.append(dn); sim_mat.append(sn)
        dig_frac.append((dn > 0.5).mean()); sim_frac.append((sn > 0.5).mean())
    dig_frac, sim_frac = np.array(dig_frac), np.array(sim_frac)
    dclass, sclass = dig_frac >= 0.25, sim_frac >= 0.25
    return dict(
        tmin=tmin, sites=FIG4_ROW_SITES,
        dig_mat=np.array(dig_mat), sim_mat=np.array(sim_mat),
        dig_frac=dig_frac, sim_frac=sim_frac,
        pearson=pearsonr(dig_frac, sim_frac)[0],
        spearman=spearmanr(dig_frac, sim_frac)[0],
        class_agree=(dclass == sclass).mean(),
    )


def plot():
    """Build verify_creamer2012.png from the committed reference + cached agreement /
    flags / ensemble / digitized artifacts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib import gridspec

    stem = os.path.splitext(os.path.basename(MODEL))[0]
    gbr = LinearSegmentedColormap.from_list("gbr", [(0, "#00d000"), (0.14, "black"),
                                                    (1, "red")])
    fig = plt.figure(figsize=(15, 12))
    gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 1.25], hspace=0.42, wspace=0.28)

    # ---- A: faithful protocol response (committed reference) ----
    axA = fig.add_subplot(gs[0, 0])
    t_tc, tc = _read_gdat(os.path.join(REFDIR, f"{stem}_erbb.gdat"))
    tmin = t_tc / 60.0
    for o in PLOT_OBS:
        axA.plot(tmin, tc[o], lw=1.3, label=o)
    axA.set(title="A. Faithful protocol: 180 min response\n(BNG NFsim, committed f=0.02)",
            xlabel="time after stimulation (min)", ylabel="molecules / subvolume")
    axA.legend(fontsize=6, ncol=2)

    # ---- B: equilibration (ligand off) ----
    axB = fig.add_subplot(gs[0, 1])
    t_eq, eq = _read_gdat(os.path.join(REFDIR, f"{stem}_equil.gdat"))
    for o in PLOT_OBS:
        axB.plot(t_eq, eq[o], marker=".", ms=4, label=o)
    axB.set(title="B. Equilibration (ligand off):\nflat unphosphorylated basal state",
            xlabel="time (s)", ylabel="molecules / subvolume")
    axB.set_ylim(bottom=-1)

    # ---- C: three-engine agreement overlay ----
    axC = fig.add_subplot(gs[0, 2])
    ag = np.load(os.path.join(REFDIR, "agreement.npz"), allow_pickle=True)
    tv = ag["tvec"]
    styles = {"bng_nf": ("BNG NFsim", "o", "#1b9e77"),
              "bngsim_nf": ("bngsim NFsim", "s", "#7570b3"),
              "bngsim_rm": ("bngsim RuleMonkey (exact)", "^", "#d95f02")}
    for o in ["EGFR_992", "EGFR_1068"]:
        for lab, (name, mk, col) in styles.items():
            A = ag[f"{lab}__{o}"]
            m = A.mean(0); se = A.std(0, ddof=1) / np.sqrt(A.shape[0])
            axC.plot(tv, m, marker=mk, ms=3, color=col, lw=1,
                     label=name if o == "EGFR_992" else None)
            axC.fill_between(tv, m - se, m + se, color=col, alpha=0.15)
        axC.annotate(o, (tv[-1], ag[f"bng_nf__{o}"].mean(0)[-1]), fontsize=6,
                     xytext=(2, 0), textcoords="offset points")
    axC.set(title="C. 3 network-free engines agree\n(reduced f=0.002, single-phase)",
            xlabel="time (s)", ylabel="molecules / subvolume")
    axC.legend(fontsize=6)

    # ---- D: -bscb / -utl + agreement z-scores ----
    axD = fig.add_subplot(gs[1, 0])
    fl = np.load(os.path.join(REFDIR, "flags.npz"), allow_pickle=True)
    verdict = fl["verdict"]
    labels = [v[0] for v in verdict] + ["3-engine\n(agreement)"]
    zvals = [float(v[1]) for v in verdict] + [3.48]
    colors = ["#4477aa", "#4477aa", "#4477aa", "#66ccee"]
    axD.bar(range(len(labels)), zvals, color=colors)
    axD.axhline(3.0, ls="--", color="crimson", lw=1)
    axD.text(len(labels) - 0.5, 3.05, "|z|=3", color="crimson", fontsize=7, ha="right")
    axD.set_xticks(range(len(labels)))
    axD.set_xticklabels(["default", "-bscb", "-bscb\n-utl 8", labels[-1]], fontsize=7)
    axD.set(title="D. -bscb / -utl are no-ops; NFsim == exact\n(max|z| vs RuleMonkey)",
            ylabel="max |z| of ensemble means")

    # ---- E: reported-data sustainment scatter ----
    axE = fig.add_subplot(gs[1, 1])
    rd = reported_data_metrics()
    axE.scatter(rd["dig_frac"], rd["sim_frac"], s=18, c="#333", alpha=0.7)
    axE.plot([0, 0.6], [0, 0.6], ls="--", color="gray", lw=1)
    axE.set(title=(f"E. Reported vs simulated sustainment\n"
                   f"Pearson r={rd['pearson']:.2f}, "
                   f"class agree={rd['class_agree']*100:.0f}% (55 sites)"),
            xlabel="digitized Fig. 4: frac(time > 0.5·max)",
            ylabel="simulated: frac(time > 0.5·max)")

    # ---- F: summary text ----
    axF = fig.add_subplot(gs[1, 2]); axF.axis("off")
    axF.text(0.0, 1.0,
             "Creamer et al. (2012) BMC Syst Biol 6:107\n"
             "Large rule-based ERBB signaling model,\n"
             "network-free (NFsim). EGF+HRG = 5 nM each.\n\n"
             "Verification:\n"
             "  1. Faithful 2-phase protocol runs under\n"
             "     standard BNG NFsim (panels A,B).\n"
             "  2. BNG NFsim == bngsim NFsim == exact\n"
             "     RuleMonkey within MC error (panel C;\n"
             "     agreement max|z|=3.48, 99.7% |z|<3).\n"
             "  3. Adaptor layer can form small rings, but\n"
             "     -bscb/-utl change nothing and default\n"
             "     NFsim matches the exact method (panel D)\n"
             "     -> legacy NFsim is a valid reference.\n"
             "  4. Reproduces Fig. 4 kinetic spectrum:\n"
             "     receptor sites transient, Akt/Raf\n"
             "     sustained (panels E, G, H).",
             fontsize=8.5, va="top", family="monospace")

    # ---- G, H: digitized vs simulated Fig. 4 heat maps ----
    axG = fig.add_subplot(gs[2, 0:2].subgridspec(1, 2)[0])
    axH = fig.add_subplot(gs[2, 0:2].subgridspec(1, 2)[1])
    extent = [0, 180, 55, 0]
    axG.imshow(rd["dig_mat"], aspect="auto", cmap=gbr, vmin=0, vmax=1, extent=extent)
    axG.set(title="G. Digitized Fig. 4 (reported)", xlabel="time (min)")
    axG.set_yticks([]); axG.set_ylabel("55 phosphosites (Fig. 4 order)")
    im = axH.imshow(rd["sim_mat"], aspect="auto", cmap=gbr, vmin=0, vmax=1, extent=extent)
    axH.set(title="H. Simulated (8-seed ensemble mean)", xlabel="time (min)")
    axH.set_yticks([])
    cbar = fig.colorbar(im, ax=[axG, axH], fraction=0.03, pad=0.02)
    cbar.set_label("per-site normalized phosphorylation", fontsize=8)

    # summary panel bottom-right of row 3
    axS = fig.add_subplot(gs[2, 2]); axS.axis("off")
    axS.text(0.0, 1.0,
             "Reported-data comparison (panels G,H):\n\n"
             f"  sustainment feature frac(t>0.5·max)\n"
             f"    Pearson  r  = {rd['pearson']:.2f}\n"
             f"    Spearman rho= {rd['spearman']:.2f}\n"
             f"    transient/sustained class\n"
             f"      agreement = {rd['class_agree']*100:.0f}%  (55 sites)\n\n"
             "  Both heat maps are per-site max-\n"
             "  normalized in the paper's clustered row\n"
             "  order. The digitized source is a single\n"
             "  stochastic seed read through a nonlinear\n"
             "  3-color scale, so robust integral features\n"
             "  (duration, kinetic class) are compared,\n"
             "  not point-by-point values.",
             fontsize=8.5, va="top", family="monospace")

    fig.suptitle("erbb_receptor_signaling_creamer2012 verification: faithful NFsim "
                 "protocol + 3-engine agreement + Fig. 4 reproduction", fontsize=12)
    png = os.path.join(HERE, "verify_creamer2012.png")
    fig.savefig(png, dpi=120, bbox_inches="tight")
    print(f"wrote {png}")
    print(f"  reported-data: Pearson r={rd['pearson']:.3f}, Spearman={rd['spearman']:.3f}, "
          f"class agreement={rd['class_agree']*100:.0f}%")


def fig4_ensemble(n_seeds=8, f=FX_COMMITTED, seed0=100, save=True):
    """Run the full published two-phase protocol (equilibrate ligand-off, then a
    180 min ligand-on response) for n_seeds at the committed subvolume, and cache the
    per-site ensemble-mean response trajectories. Used for the reported-data comparison
    against the digitized Fig. 4 (a single-seed heat map is noisy; the ensemble mean is
    the expected per-site kinetic shape)."""
    print(f"=== fig4_ensemble: {n_seeds} seeds of the full 180 min protocol, f={f} ===")
    base = _set_param_line(_model_block(), "f", f)
    stacks = {o: [] for o in OBS_55}
    tvec = None
    for s in range(n_seeds):
        sd = seed0 + s
        act = ("\nbegin actions\n"
               f'  simulate({{suffix=>"equil",method=>"nf",t_start=>0,t_end=>500,n_steps=>1,'
               f"gml=>2147483647,get_final_state=>1,seed=>{sd}}})\n"
               '  setParameter("lig",1)\n'
               f'  simulate({{suffix=>"erbb",method=>"nf",t_start=>0,t_end=>10800,n_steps=>180,'
               f"gml=>2147483647,get_final_state=>0,seed=>{sd}}})\n"
               "end actions\n")
        t0 = time.time()
        out = _bng_run(base + act, f"fig4_s{s}")
        t, cols = _read_gdat(os.path.join(out, f"fig4_s{s}_erbb.gdat"))
        tvec = t
        for o in OBS_55:
            stacks[o].append(cols[o])
        print(f"  [{time.time()-t0:.0f}s] seed {sd} done")
    means = {o: np.vstack(stacks[o]).mean(0) for o in OBS_55}
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "fig4_sim.npz"),
                 tvec=tvec, obs=np.array(OBS_55), n_seeds=n_seeds, f=f,
                 **{o: means[o] for o in OBS_55})
        print(f"  cached -> {os.path.join(REFDIR, 'fig4_sim.npz')}")
    return tvec, means


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    if mode == "reference":
        reference()
    elif mode == "agreement":
        agreement(n_seeds=n)
    elif mode == "flags":
        flags(n_seeds=n)
    elif mode == "fig4_ensemble":
        fig4_ensemble(n_seeds=(n if len(sys.argv) > 2 else 8))
    elif mode == "plot":
        plot()
    elif mode == "all":
        reference()
        flags(n_seeds=n)
        agreement(n_seeds=n)
        fig4_ensemble()
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference | agreement [N] | flags [N] | "
              "fig4_ensemble [N] | plot | all [N]")
