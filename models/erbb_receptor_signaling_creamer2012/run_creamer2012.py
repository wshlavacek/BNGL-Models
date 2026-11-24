"""Driver + verification for erbb_receptor_signaling_creamer2012.bngl -- the large
rule-based ERBB receptor signaling model of Creamer et al. (2012), BMC Systems Biology
6:107 (doi:10.1186/1752-0509-6-107; Supplementary File 1). The reaction network is far
too large to enumerate (>10^100 reachable species), so the model is simulated
NETWORK-FREE. Unlike a pure multisite-phosphorylation model, the ERBB adaptor layer can
in principle form small rings (a receptor dimer bridged by Grb2 + Shc), so whether the
default (legacy) NFsim molecularity handling matches the exact network-free method is
NOT assumed -- it is tested here.

Verified here, through BioNetGen's own NFsim (BNG2.pl) unless stated:

  fig4_ensemble -- runs the committed model's FAITHFUL published protocol -- the file's own
                actions block (equilibrate 500 s with ligand binding off, setParameter
                lig=1, 180 min response), unchanged except for pinned seeds (1000+r
                equilibration, 5000+r response) and a 10-s output grid (the Fig. 4 bitmap's
                own column spacing; the trajectory is unchanged) -- N=144 times in parallel
                at the committed f=0.02, caching per-run features (F_hi, F_half, auc,
                t_peak), per-site sums for the ensemble mean and its s.e.m., and replicate
                1's trajectories to reference/fig4_ensemble_f0.02.npz. `fig4_negctl` does
                the same with the pre-fix rules (ten association constants -- Grb2kp1-4,
                Gab1kp156-161 -- not divided by f), the negative control.

  data       -- level 2 against Fig. 4 (one NFsim run at the paper's compartment scale 0.2,
                per-site max-normalized; digitize_creamer2012.py). Per site: F_hi, the
                fraction of the 180 min above 0.14 of the site's maximum (the colorbar's
                black stop: red vs green in the figure, no colour-blending assumption). z =
                (published - ensemble mean) / sd, sd^2 = digitization^2 + (f=0.02 single-run
                SD)^2 x 0.02/0.2 (the published run's own scatter: a count's relative
                fluctuation goes as 1/sqrt(system size)) + s.e.m.^2, over the sites whose
                ensemble-mean peak is >= 20 molecules. Asserted: max |z| <= Phi^-1(1 -
                0.01/2k) (Bonferroni, 1% family-wise, valid under the strong dependence
                between sites of one cascade) and rms z <= sqrt(chi2_0.999(k)/k); the
                negative control must fail both.

  reference  -- runs the committed file verbatim (`reference save` rewrites reference/*.gdat);
                `xml` rewrites reference/<stem>.xml at the shipped nominals.

  agreement / flags -- the older three-engine cross-check (standard BioNetGen NFsim, bngsim
                NfsimSession, bngsim RuleMonkeySession) and -bscb/-utl determination, on a
                reduced single-phase ligand-on protocol at f=FX_CROSS. Two arms run through
                bngsim; the cached results predate the scaling fix (they compare engines on
                the pre-fix rules).

  plot       -- builds verify_creamer2012.png from committed artifacts only.

Parameter propagation note (lanl/bngsim#44): a POST-init NFsim set_param does not refresh
a rule's rate constant, and RuleMonkeySession.set_param does not re-derive seed-species
amounts. The bngsim arms avoid both by BAKING f (and lig) into the parameter block and
calling writeXML() before load.

Requirements: BioNetGen 2.9.3 (BNG2.pl via $BNG2_PL / $BNGPATH, or the bionetgen CLI);
`agreement` and `flags` alone need the `bngsim` package. Run with the repo venv:
    python run_creamer2012.py fig4_ensemble 144 0.02 [workers] [cache]   # level 2 (hours)
    python run_creamer2012.py fig4_negctl 72 0.02 [workers] [cache]      # negative control
    python run_creamer2012.py data                # level 2, from the caches
    python run_creamer2012.py reference [save]    # faithful protocol, verbatim
    python run_creamer2012.py agreement 24        # 3-engine cross-check (bngsim)
    python run_creamer2012.py flags 24            # -bscb / -utl vs RuleMonkey (bngsim)
    python run_creamer2012.py plot                # (re)build verify_creamer2012.png
"""
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from functools import partial

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
def find_bionetgen():
    """BNG2.pl from $BNG2_PL or $BNGPATH, then BNG2.pl on PATH, then the `bionetgen` CLI.
    Returns (path, is_cli)."""
    for cand in (os.environ.get("BNG2_PL"),
                 os.path.join(os.environ["BNGPATH"], "BNG2.pl") if os.environ.get("BNGPATH")
                 else None,
                 shutil.which("BNG2.pl")):
        if cand and os.path.isfile(cand):
            return os.path.abspath(cand), False
    cli = shutil.which("bionetgen")
    if cli:
        return cli, True
    raise FileNotFoundError("BioNetGen not found: set BNG2_PL (or BNGPATH) to a BioNetGen "
                            "2.9.3 BNG2.pl, or put the bionetgen CLI on PATH.")


@contextmanager
def _scratch(tag):
    """A per-run temporary directory inside this folder, removed afterwards."""
    d = tempfile.mkdtemp(prefix=f".tmp_{tag}_", dir=HERE)
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _bng_run(bngl_text, tag, workdir):
    """Write bngl_text to workdir/<tag>.bngl, run BioNetGen there, return the output dir."""
    bng, cli = find_bionetgen()
    bngl = os.path.join(workdir, f"{tag}.bngl")
    open(bngl, "w").write(bngl_text)
    cmd = ([bng, "run", "-i", bngl, "-o", workdir] if cli
           else ["perl", bng, "--outdir", workdir, bngl])
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    assert r.returncode == 0, r.stdout[-3000:] + "\n---STDERR---\n" + r.stderr[-1500:]
    return workdir


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
    out = _bng_run(text, "xml", tempfile.mkdtemp(prefix=".tmp_xml_", dir=HERE))
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
        with _scratch(f"native_s{s}") as wd:
            out = _bng_run(base + act, f"native_s{s}", wd)
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
def write_reference_xml():
    """reference/<stem>.xml: writeXML() of the committed model at its shipped nominals
    (f=0.02, lig=0). The temporary file carries the model's stem, which BioNetGen writes into
    the XML as the model id."""
    stem = os.path.splitext(os.path.basename(MODEL))[0]
    with _scratch("xml") as wd:
        _bng_run(_model_block() + "\nbegin actions\n  writeXML({overwrite=>1})\nend actions\n",
                 stem, wd)
        shutil.copy(os.path.join(wd, f"{stem}.xml"), os.path.join(REFDIR, f"{stem}.xml"))
    print(f"  wrote reference/{stem}.xml")


def reference(save=False):
    """Run the committed model's faithful published protocol (the file verbatim) with
    standard BNG NFsim; with save=True the two .gdat become the committed reference/."""
    print("=== reference: faithful published protocol (BNG NFsim, f=0.02) ===")
    t0 = time.time()
    stem = os.path.splitext(os.path.basename(MODEL))[0]
    with _scratch("reference") as wd:
        out = _bng_run(open(MODEL).read(), "reference", wd)
        t_eq, eq = _read_gdat(os.path.join(out, "reference_equil.gdat"))
        t_tc, tc = _read_gdat(os.path.join(out, "reference_erbb.gdat"))
        if save:
            for suf in ("equil", "erbb"):
                shutil.copy(os.path.join(out, f"reference_{suf}.gdat"),
                            os.path.join(REFDIR, f"{stem}_{suf}.gdat"))
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
    shutil.rmtree(os.path.dirname(xml), ignore_errors=True)
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
    shutil.rmtree(os.path.dirname(xml), ignore_errors=True)
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


# ---------------------------------------------------------------- Fig. 4 (level 2)
# The heat map's own grid (digitize_creamer2012.py): 1081 columns = 0..10800 s every 10 s.
FIG4_DT = 10.0
FIG4_NSTEPS = 1080
BLACK = 0.14                 # the colorbar's black stop: 0.14 of each site's maximum
DIG_NOISE = (0.005, 0.03)    # JPEG noise in value units, green / red side (digitizer header)
N_FIG4 = 144                 # committed replicate count at f=0.02
F_PAPER = 0.2                # Creamer et al. Methods: "compartment sizes were scaled by 0.2"
DATA_FWER = 0.01             # family-wise false-alarm rate of the per-site bound
HIGH_COUNT = 20.0            # a site enters the data test if its ensemble-mean peak at
                             # f=0.02 is >= 20 molecules (>= 200 in the paper's f=0.2 run)
# The ten association constants the pre-fix file (sha256 6672c0a8...) left unscaled by f.
PREFIX_UNSCALED = (["Grb2kp1", "Grb2kp2", "Grb2kp3", "Grb2kp4"]
                   + [f"Gab1kp{k}" for k in range(156, 162)])


def fig4_npz(f=FX_COMMITTED, negctl=False):
    return os.path.join(REFDIR, f"fig4_ensemble_f{f:g}{'_negctl' if negctl else ''}.npz")


def _fig4_actions(rep):
    """The committed actions block, verbatim, except: the response is written every 10 s
    (n_steps 180 -> 1080, the Fig. 4 column grid; the trajectory is unchanged) and seeds are
    pinned (1000+rep for the equilibration, 5000+rep for the response)."""
    text = open(MODEL).read()
    acts = text[text.index("\nend model") + len("\nend model"):]
    acts, n1 = re.subn(r'(suffix=>"erbb".*?n_steps=>)180\b', rf"\g<1>{FIG4_NSTEPS}", acts,
                       flags=re.S)
    assert n1 == 1, "response simulate() not found"
    seeds = iter((1000 + rep, 5000 + rep))
    acts, n = re.subn(r"\}\)", lambda _: f",seed=>{next(seeds)}}})", acts)
    assert n == 2, f"expected two simulate() calls, found {n}"
    return acts


def _prefix_block():
    """The model block with the ten constants of PREFIX_UNSCALED un-divided by f: the rules
    of the pre-fix file, for the negative control."""
    blk = _model_block()
    for name in PREFIX_UNSCALED:
        blk, n = re.subn(rf"(?m)^({name}\s+)\(([^)]*)\)/f", r"\g<1>\2", blk)
        assert n == 1, name
    return blk


def _fig4_rep(rep, f=FX_COMMITTED, negctl=False, cache=None):
    """One replicate of the published protocol; returns (t, X[55 sites, 1081])."""
    if cache:
        fn = os.path.join(cache, f"rep_{rep:04d}.npz")
        if os.path.exists(fn):
            d = np.load(fn)
            return d["t"], d["X"]
    tag = f"fig4_{rep:04d}"
    blk = _prefix_block() if negctl else _model_block()
    if f != FX_COMMITTED:
        blk = _set_param_line(blk, "f", f)
    with _scratch(tag) as wd:
        out = _bng_run(blk + _fig4_actions(rep), tag, wd)
        t, cols = _read_gdat(os.path.join(out, f"{tag}_erbb.gdat"))
    X = np.array([cols[s_] for s_ in FIG4_ROW_SITES])
    assert X.shape[1] == FIG4_NSTEPS + 1
    if cache:
        np.savez(fn + ".part.npz", t=t, X=X)
        os.replace(fn + ".part.npz", fn)
    return t, X


def normalize(X):
    """Per-site max normalization along the last (time) axis, as in Fig. 4."""
    X = np.asarray(X, float)
    mx = X.max(-1, keepdims=True)
    return np.divide(X, mx, out=np.zeros_like(X), where=mx > 0)


def fig4_features(Y):
    """Per-site features of max-normalized time courses Y[..., 1081]. F_hi: fraction of the
    180 min spent above the colorbar's black stop (0.14 of the maximum) -- red vs green in
    the figure, so it needs no assumption about how colours are blended between stops."""
    Y = np.asarray(Y, float)
    return dict(F_hi=(Y > BLACK).mean(-1), F_half=(Y > 0.5).mean(-1), auc=Y.mean(-1),
                t_peak=Y.argmax(-1) * FIG4_DT / 60.0)


def _bng_version():
    bng, cli = find_bionetgen()
    r = subprocess.run([bng, "-v"] if cli else ["perl", bng, "-v"], capture_output=True,
                       text=True)
    line = [ln for ln in r.stdout.splitlines() if "BioNetGen" in ln]
    return line[0].strip() if line else "unknown"


def fig4_ensemble(n=N_FIG4, f=FX_COMMITTED, workers=None, negctl=False, cache=None,
                  reps=None):
    """n replicates of the published two-phase protocol (the committed file's own actions;
    10-s output; pinned seeds) at subvolume f, in parallel; caches per-replicate features,
    per-site sums for the ensemble mean and its s.e.m., and replicate 1's trajectories to
    reference/fig4_ensemble_f<f>[_negctl].npz. `cache` (a directory) keeps finished
    replicates so an interrupted campaign resumes; `reps` (a list) assembles exactly those
    replicates, e.g. the ones an interrupted campaign finished."""
    from multiprocessing import Pool
    workers = workers or min(n, os.cpu_count() or 1)
    if cache:
        os.makedirs(cache, exist_ok=True)
    what = "pre-fix rules (negative control)" if negctl else "committed rules"
    print(f"=== fig4_ensemble: {n} replicates, f={f}, {what}, {workers} workers ===")
    t0 = time.time()
    reps = np.array(sorted(reps) if reps is not None else range(1, n + 1))
    n = len(reps)
    with Pool(workers) as pool:
        runs = pool.map(partial(_fig4_rep, f=f, negctl=negctl, cache=cache), [int(r) for r in reps])
    t = runs[0][0]
    X = np.array([r[1] for r in runs], dtype=float)
    F = fig4_features(normalize(X))
    path = fig4_npz(f, negctl)
    np.savez_compressed(path, t=t, sites=np.array(FIG4_ROW_SITES), f=f, n=n,
                        reps=reps, seed_equil=1000 + reps, seed_resp=5000 + reps,
                        **{k: v.astype(np.float32) for k, v in F.items()},
                        runmax=X.max(-1).astype(np.float32), sum=X.sum(0), sumsq=(X**2).sum(0),
                        run1=X[0].astype(np.float32), negctl=negctl,
                        bionetgen=np.array(_bng_version()),
                        nfsim_flags=np.array("-cb (BNG2.pl default); no -bscb; default -utl; "
                                             "gml 2147483647"))
    print(f"  [{time.time()-t0:.0f}s] cached -> {path}")


def read_digitized():
    """reference/creamer2012_fig4_digitized.csv -> (t in min, V[55, 1081]), Fig. 4 row order."""
    path = os.path.join(REFDIR, "creamer2012_fig4_digitized.csv")
    with open(path) as fh:
        rows = list(csv.reader(ln for ln in fh if not ln.startswith("#")))
    tmin = np.array([float(x) for x in rows[0][1:]])
    sites = [r[0] for r in rows[1:]]
    assert sites == FIG4_ROW_SITES, "digitized rows are not in Fig. 4 order"
    return tmin, np.array([[float(x) for x in r[1:]] for r in rows[1:]])


def digitized_features():
    """Features of the published run, and the digitization error of F_hi: half the change
    when the black threshold moves by the JPEG noise on either side (floor: one column)."""
    _, V = read_digitized()
    F = fig4_features(V)
    lo = (V > BLACK - DIG_NOISE[0]).mean(-1)
    hi = (V > BLACK + DIG_NOISE[1]).mean(-1)
    return F, np.maximum(0.5 * np.abs(lo - hi), 1.0 / V.shape[1])


def load_fig4(path):
    d = np.load(path)
    n = int(d["n"])
    mean = d["sum"] / n
    var = np.maximum(d["sumsq"] / n - mean**2, 0.0) * n / (n - 1)
    return d, mean, np.sqrt(var / n)


def data_test(path=None):
    """Level 2 at the committed subvolume. The published Fig. 4 is ONE run at the paper's
    f=0.2, the ensemble here is n runs of the committed file at f=0.02, so the test compares
    the published run's F_hi, site by site, with the ensemble's mean single-run F_hi, over the
    sites whose single-run F_hi is a faithful estimate of the mean (ensemble-mean peak >=
    HIGH_COUNT molecules). Error budget per site: digitization (threshold sensitivity), the
    published run's own run-to-run scatter -- the f=0.02 single-run SD scaled by
    sqrt(0.02/0.2) (fluctuations of a count go as 1/sqrt(system size)) -- and the s.e.m. of
    the ensemble mean. Returns per-site z and the pooled statistics."""
    path = path or fig4_npz()
    d, mean, _ = load_fig4(path)
    Fp, sd_dig = digitized_features()
    Fs = d["F_hi"].astype(float)
    n = Fs.shape[0]
    mu, sd1 = Fs.mean(0), Fs.std(0, ddof=1)
    f = float(d["f"])
    sd_pub = sd1 * np.sqrt(f / F_PAPER)
    den = np.sqrt(sd_dig**2 + sd_pub**2 + sd1**2 / n + (1.0 / 1081) ** 2)
    z = (Fp["F_hi"] - mu) / den
    sel = mean.max(-1) >= HIGH_COUNT
    k = int(sel.sum())
    from scipy.stats import chi2, norm
    D = float(np.sqrt(np.mean(z[sel] ** 2)))
    D_tol = float(np.sqrt(chi2.ppf(0.999, k) / k))
    # Sites move together (the five Raf1 activation-loop sites and MEK1 S222 are one
    # cascade), so the per-site bound is Bonferroni, valid under any dependence: 1%
    # family-wise over the k sites.
    zmax_tol = float(norm.ppf(1.0 - DATA_FWER / (2 * k)))
    return dict(z=z, sel=sel, k=k, D=D, D_tol=D_tol, zmax=float(np.abs(z[sel]).max()),
                zmax_tol=zmax_tol,
                F_pub=Fp["F_hi"], F_sim=mu, sd_den=den, sd_dig=sd_dig, n=n, f=f,
                peak=mean.max(-1))


def f_comparison(fa=FX_COMMITTED, fb=F_PAPER):
    """Committed f=0.02 against the paper's f=0.2: ensemble-mean concentrations (count / f)
    per site over the 180 min, as z = difference / combined s.e.m., over the sites with a
    mean peak >= HIGH_COUNT at fa."""
    A, ma, sa = load_fig4(fig4_npz(fa))
    B, mb, sb = load_fig4(fig4_npz(fb))
    ca, cb = ma / fa, mb / fb
    se = np.sqrt((sa / fa) ** 2 + (sb / fb) ** 2)
    z = np.divide(ca - cb, se, out=np.zeros_like(ca), where=se > 0)
    sel = ma.max(-1) >= HIGH_COUNT
    rel = np.abs(ca - cb).max(-1) / np.maximum(cb.max(-1), 1e-300)
    return dict(z=z, sel=sel, zmax=float(np.abs(z[sel]).max()),
                frac_lt3=float((np.abs(z[sel]) < 3).mean()), rel=rel, na=int(A["n"]),
                nb=int(B["n"]), ca=ca, cb=cb)


def plot():
    """Build verify_creamer2012.png from committed artifacts only."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import LinearSegmentedColormap

    stem = os.path.splitext(os.path.basename(MODEL))[0]
    gbr = LinearSegmentedColormap.from_list(
        "gbr", [(0, (51 / 255, 254 / 255, 50 / 255)), (BLACK, "black"), (1, "red")])
    fig = plt.figure(figsize=(16, 13))
    gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 1.3], hspace=0.45, wspace=0.28)

    # A, B: the faithful protocol (committed reference run).
    axA = fig.add_subplot(gs[0, 0])
    t_tc, tc = _read_gdat(os.path.join(REFDIR, f"{stem}_erbb.gdat"))
    for o in PLOT_OBS:
        axA.plot(t_tc / 60.0, tc[o], lw=1.3, label=o)
    axA.set(title="A. Faithful protocol: 180 min response\n(BNG NFsim, committed f=0.02)",
            xlabel="time after stimulation (min)", ylabel="molecules / subvolume")
    axA.legend(fontsize=6, ncol=2)
    axB = fig.add_subplot(gs[0, 1])
    t_eq, eq = _read_gdat(os.path.join(REFDIR, f"{stem}_equil.gdat"))
    for o in PLOT_OBS:
        axB.plot(t_eq, eq[o], marker=".", ms=4, label=o)
    axB.set(title="B. Equilibration (ligand off):\nflat unphosphorylated basal state",
            xlabel="time (s)", ylabel="molecules / subvolume")
    axB.set_ylim(bottom=-1)

    # C: level 2 -- per-site F_hi, published run vs the f=0.02 ensemble.
    r = data_test()
    neg = data_test(fig4_npz(negctl=True))
    axC = fig.add_subplot(gs[0, 2])
    sel = r["sel"]
    axC.plot([0, 1], [0, 1], color="gray", lw=0.8)
    axC.errorbar(r["F_pub"][sel], r["F_sim"][sel], yerr=r["sd_den"][sel], fmt="o", ms=3.5,
                 color="#1b9e77", label=f"committed file ({r['k']} high-count sites)")
    axC.plot(r["F_pub"][~sel], r["F_sim"][~sel], "o", ms=3.5, mfc="none", color="#1b9e77",
             label="low-count sites (not tested)")
    axC.plot(neg["F_pub"][sel], neg["F_sim"][sel], "x", ms=4, color="crimson",
             label="pre-fix rules (negative control)")
    axC.set(xlim=(0, 1), ylim=(0, 1), xlabel="published Fig. 4: time above 0.14 max (fraction)",
            ylabel=f"BNG NFsim, {r['n']} runs: mean single-run value")
    axC.set_title(f"C. Level 2 vs Fig. 4: rms z {r['D']:.2f} (tol {r['D_tol']:.2f}), "
                  f"max|z| {r['zmax']:.2f} (tol {r['zmax_tol']:.2f})\npre-fix rules: rms z "
                  f"{neg['D']:.1f}, max|z| {neg['zmax']:.0f} -> rejected", fontsize=9)
    axC.legend(fontsize=6, loc="lower right")

    # D: per-site z.
    axD = fig.add_subplot(gs[1, 0:2])
    x = np.arange(len(FIG4_ROW_SITES))
    axD.bar(x - 0.2, np.where(sel, r["z"], 0), width=0.4, color="#1b9e77", label="committed")
    axD.bar(x + 0.2, np.where(sel, neg["z"], 0), width=0.4, color="crimson", alpha=0.6,
            label="pre-fix rules")
    axD.axhline(0, color="k", lw=0.6)
    axD.set_xticks(x)
    axD.set_xticklabels(FIG4_ROW_SITES, rotation=90, fontsize=6)
    axD.set_ylim(-15, 15)
    axD.set(ylabel="z = (sim - published) / sd", title="D. Per-site z of the time above "
            "0.14 max (high-count sites; low-count sites blank; clipped at +/-15)")
    axD.legend(fontsize=7)

    # E: summary text.
    axE = fig.add_subplot(gs[1, 2])
    axE.axis("off")
    lines = ["Creamer et al. (2012) BMC Syst Biol 6:107",
             "BNGsim-independent checks (BioNetGen NFsim):",
             f"  level 2 vs digitized Fig. 4, {r['n']} runs at f={r['f']:g}:",
             f"    rms z {r['D']:.2f} <= {r['D_tol']:.2f} over {r['k']} sites,",
             f"    max |z| {r['zmax']:.2f} <= {r['zmax_tol']:.2f}  "
             + ("PASS" if r["D"] <= r["D_tol"] and r["zmax"] <= r["zmax_tol"] else "FAIL"),
             f"  negative control (pre-fix rules): rms z {neg['D']:.1f} -> rejected"]
    if os.path.exists(fig4_npz(F_PAPER)):
        fc = f_comparison()
        lines += [f"  f=0.02 vs the paper's f=0.2 ({fc['na']} vs {fc['nb']} runs):",
                  f"    mean concentration max|z| {fc['zmax']:.1f}, "
                  f"frac|z|<3 {fc['frac_lt3']:.3f}"]
    lines += ["", "Panels F, G: digitized Fig. 4 and one",
              "  simulated f=0.02 run, same colour map.",
              "The older three-engine and -bscb/-utl",
              "  caches ran the pre-fix rules (bngsim arms)",
              "  and are not shown."]
    axE.text(0, 1, "\n".join(lines), fontsize=8.5, va="top", family="monospace")

    # F, G: heat maps.
    tmin, V = read_digitized()
    d = np.load(fig4_npz())
    extent = [0, 180, 55, 0]
    for k_, (M, ttl) in enumerate(((V, "F. Digitized Fig. 4 (one run, f=0.2)"),
                                   (normalize(d["run1"]), "G. BNG NFsim, committed file, "
                                    "replicate 1 (f=0.02)"))):
        ax = fig.add_subplot(gs[2, k_])
        im = ax.imshow(M, aspect="auto", cmap=gbr, vmin=0, vmax=1, extent=extent)
        ax.set(title=ttl, xlabel="time (min)")
        ax.set_yticks(np.arange(55) + 0.5)
        ax.set_yticklabels(FIG4_ROW_SITES if k_ == 0 else [], fontsize=4.5)
    ax = fig.add_subplot(gs[2, 2])
    mean = d["sum"] / int(d["n"])
    ax.imshow(normalize(mean), aspect="auto", cmap=gbr, vmin=0, vmax=1, extent=extent)
    ax.set(title=f"H. Ensemble mean of {int(d['n'])} runs (f=0.02)", xlabel="time (min)")
    ax.set_yticks([])
    fig.colorbar(im, ax=fig.axes[-3:], fraction=0.02, pad=0.01).set_label(
        "per-site normalized level", fontsize=8)

    fig.suptitle("erbb_receptor_signaling_creamer2012: published two-phase protocol under "
                 "BioNetGen NFsim vs Creamer et al. (2012) Fig. 4", fontsize=12)
    png = os.path.join(HERE, "verify_creamer2012.png")
    fig.savefig(png, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {png}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    n = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 24
    if mode == "reference":
        reference(save=len(sys.argv) > 2 and sys.argv[2] == "save")
    elif mode == "xml":
        write_reference_xml()
    elif mode == "agreement":
        agreement(n_seeds=n)
    elif mode == "flags":
        flags(n_seeds=n)
    elif mode in ("fig4_ensemble", "fig4_negctl"):
        # fig4_ensemble [N] [f] [workers] [cache_dir] [first_rep]
        n = int(sys.argv[2]) if len(sys.argv) > 2 else N_FIG4
        f = float(sys.argv[3]) if len(sys.argv) > 3 else FX_COMMITTED
        w = int(sys.argv[4]) if len(sys.argv) > 4 else None
        c = sys.argv[5] if len(sys.argv) > 5 else None
        r0 = int(sys.argv[6]) if len(sys.argv) > 6 else 1
        fig4_ensemble(n, f, w, negctl=(mode == "fig4_negctl"), cache=c,
                      reps=range(r0, r0 + n))
    elif mode == "fig4_assemble":
        # fig4_assemble f cache_dir [negctl]: the npz from whatever replicates are cached
        f, c = float(sys.argv[2]), sys.argv[3]
        reps = sorted(int(x[4:8]) for x in os.listdir(c) if re.fullmatch(r"rep_\d{4}\.npz", x))
        fig4_ensemble(f=f, workers=2, negctl=len(sys.argv) > 4, cache=c, reps=reps)
    elif mode == "data":
        r = data_test()
        print(f"  {r['k']} high-count sites; rms z = {r['D']:.2f} (tol {r['D_tol']:.2f}), "
              f"max |z| = {r['zmax']:.2f} (tol {r['zmax_tol']:.2f})")
        assert r["D"] <= r["D_tol"] and r["zmax"] <= r["zmax_tol"]
    elif mode == "plot":
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference | agreement [N] | flags [N] | "
              "fig4_ensemble [N] [f] [workers] [cache] | fig4_negctl [...] | data | plot")
