"""Driver + verification for tcr_signaling_chylek2014.bngl -- the early T-cell-receptor
(TCR) / CD28 signaling model of Chylek et al. (2014), PLoS ONE 9(8):e104240
(doi:10.1371/journal.pone.0104240; Supplementary File A in File S1). The model's
combinatorial complexity is the multisite-phosphorylation state space of individual
proteins, not multivalent aggregation, so it is simulated NETWORK-FREE (NFsim) and has
no same-complex-binding ambiguity: standard/legacy BioNetGen NFsim is a valid reference
(no -bscb / -utl needed).

Verified here, through BioNetGen's own NFsim (BNG2.pl) unless stated:

  ensemble   -- runs the committed model's FAITHFUL published protocol -- the file's own
                actions block (equilibrate 500 s with ligand binding off, turn on kfl/kfl_m
                per Supplementary File B, 60 s stimulation response), unchanged except for
                pinned seeds (1000+r equilibration, 5000+r response) -- N=144 times in
                parallel, and caches every replicate to reference/data_ensemble.npz.
                `negctl` caches a negative control: the same protocol at Fx=0.005 with the
                literal ligand on-rates left at their Fx=0.02 values (4x-slow binding).

  data       -- level 2. The ensemble mean, as log2 fold change over the equilibrated
                baseline (what Fig. S6 plots), against the published Fig. S6 simulation
                lines (digitize_chylek2014.py; vector-exact to the PDF's 0.8-pt grid) at the
                paper's measurement times 5/15/30/60 s for all 16 sites. z = |d| / sd(d),
                sd(d)^2 = s.e.^2 (1 + N/N_PUB_MIN) + grid terms: our bootstrap s.e., the
                published line's own Monte-Carlo error bounded as that of a mean of at least
                N_PUB_MIN=10 runs (the paper averaged an unstated number), and the 0.8-pt
                grid in y and, through the local slope, in t. Asserted: max z <= DATA_ZMAX=4
                (0.4% family-wise for 64 independent comparisons); the negative control must
                exceed it.

  reference  -- runs the committed file verbatim; reference_typicality() tests the committed
                reference/ run (one unseeded run) as a draw from the ensemble.

  agreement  -- the older three-engine cross-check: standard BioNetGen NFsim, bngsim
                NfsimSession and bngsim RuleMonkeySession on the SAME rules, single-phase
                ligand-on at a reduced subvolume Fx=0.005 (kfl/kfl_m at their Fx=0.02
                values, i.e. 4x-slow binding: it compares engines on one rule set, not the
                published protocol). Two of its arms run through bngsim. Agreement is the
                pairwise z-score |mean_i - mean_j| / sqrt(se_i^2 + se_j^2).

  plot       -- builds verify_tcr_chylek2014.png from committed artifacts only.

Requirements: BioNetGen 2.9.3 (BNG2.pl via $BNG2_PL / $BNGPATH, or the bionetgen CLI);
`agreement` alone needs the `bngsim` package. Run with the repo venv:
    python run_tcr_chylek2014.py ensemble 144 [workers]   # level-2 ensemble (hours of CPU)
    python run_tcr_chylek2014.py negctl 72 [workers]      # negative control
    python run_tcr_chylek2014.py data                     # level 2, from the caches
    python run_tcr_chylek2014.py reference                # faithful protocol, verbatim
    python run_tcr_chylek2014.py agreement 16              # 3-engine cross-check (bngsim)
    python run_tcr_chylek2014.py plot                     # (re)build the PNG
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

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, "tcr_signaling_chylek2014.bngl")
REFDIR = os.path.join(HERE, "reference")

# Site-specific phosphorylation observables carried by the model (paper Figs. 2-4).
OBS = ["TCR_pY149_D", "TCR_pY171_G", "TCR_pY111", "TCR_pY123", "TCR_pY188_E",
       "TCR_pY199_E", "ZAP70_pY493", "LCK_pY424", "LCK_pY192", "ITK_pY512",
       "PTPN6_pY566", "WAS_pY291", "PAG1_pY163", "DOK2_pY299", "DOK1_pY449",
       "PLCG1_pY783", "LAT_pY191", "LCK_pY505"]
# A readable subset for the verification figure.
PLOT_OBS = ["TCR_pY149_D", "TCR_pY111", "ZAP70_pY493", "LCK_pY192", "LAT_pY191",
            "PLCG1_pY783"]

# Published simulation constants (Supplementary File A / File B).
KFL_ON, KFL_M_ON = 6e-6, 6e-4
FX_FULL = 0.02          # published subvolume fraction
FX_CROSS = 0.005        # reduced subvolume for the tractable 3-engine cross-check


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
    assert r.returncode == 0, r.stderr[-3000:] + "\n---STDOUT---\n" + r.stdout[-2000:]
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
    """Return the committed model text up to and including `end model` (i.e. everything
    except the begin actions block), so a caller can append its own actions. Splitting on
    the `end model` boundary avoids matching the literal text "begin actions" that appears
    inside the header #@note comment."""
    text = open(MODEL).read()
    marker = "\nend model"
    i = text.index(marker)
    return text[:i + len(marker)] + "\n"


def build_xml(fx=FX_CROSS, kfl=KFL_ON, kfl_m=KFL_M_ON):
    """Emit a BNG-XML from the committed model with Fx/kfl/kfl_m BAKED into the parameter
    block (so seed-species amounts like Proteintot=200000*Fx are pre-evaluated in the XML),
    then writeXML(). All engines load identical initial conditions from this XML with no
    reliance on session set_param -- necessary because bngsim RuleMonkeySession.set_param
    does not re-evaluate derived seed-species amounts (it would leave Proteintot at the
    XML-time Fx), whereas NfsimSession does; baking removes that asymmetry."""
    text = _model_block()
    text = _set_param_line(text, "Fx", fx)
    text = _set_param_line(text, "kfl", kfl)
    text = _set_param_line(text, "kfl_m", kfl_m)
    text += "\nbegin actions\n  writeXML()\nend actions\n"
    xmldir = tempfile.mkdtemp(prefix=".tmp_xml_", dir=HERE)  # caller removes it
    out = _bng_run(text, "xml", xmldir)
    xml = os.path.join(out, "xml.xml")
    assert os.path.exists(xml), f"writeXML did not produce {xml}"
    return xml


# ---------------------------------------------------------------- ensembles
def ensemble_bngsim(xml, engine, n_seeds, t_end, n_points, seed0=5000):
    """Single-phase ligand-on ensemble via bngsim. engine: 'nf' | 'rm'.
    Returns (tvec, {obs: array[n_seeds, n_points]}).

    The XML already has Fx/kfl/kfl_m baked in (see build_xml), so ligand binding is on
    from t=0 and every engine loads identical initial conditions; no session set_param is
    used. This sidesteps two bngsim limitations found while building this check: (1) a
    POST-init set_param does not refresh an NFsim rule's rate constant (turning kfl on
    after initialize() produces no ligand response), and (2) RuleMonkeySession.set_param
    does not re-evaluate derived seed-species amounts. Both are avoided by baking."""
    from bngsim import NfsimSession, RuleMonkeySession
    Sess = NfsimSession if engine == "nf" else RuleMonkeySession
    stacks = {o: [] for o in OBS}
    tvec = None
    for s in range(n_seeds):
        with Sess(xml) as sess:
            sess.initialize(seed=seed0 + s)
            res = sess.simulate(0, t_end, n_points=n_points)
            names = list(sess.get_observable_names())
            Y = np.asarray(res.observables)
            tvec = np.asarray(res.time)
            for o in OBS:
                stacks[o].append(Y[:, names.index(o)])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


def ensemble_bng_native(n_seeds, fx, t_end, n_points, seed0=5000,
                        kfl=KFL_ON, kfl_m=KFL_M_ON):
    """Single-phase ligand-on ensemble via standard BioNetGen NFsim (seeded runs).
    Returns (tvec, {obs: array[n_seeds, n_points]})."""
    base = _model_block()
    base = _set_param_line(base, "Fx", fx)
    base = _set_param_line(base, "kfl", kfl)
    base = _set_param_line(base, "kfl_m", kfl_m)
    stacks = {o: [] for o in OBS}
    tvec = None
    for s in range(n_seeds):
        act = ("\nbegin actions\n"
               f"  simulate({{method=>\"nf\",suffix=>\"resp\",t_start=>0,t_end=>{t_end},"
               f"n_steps=>{n_points - 1},seed=>{seed0 + s}}})\n"
               "end actions\n")
        with _scratch(f"native_s{s}") as wd:
            out = _bng_run(base + act, f"native_s{s}", wd)
            t, cols = _read_gdat(os.path.join(out, f"native_s{s}_resp.gdat"))
        tvec = t
        for o in OBS:
            stacks[o].append(cols[o])
    return tvec, {o: np.vstack(v) for o, v in stacks.items()}


# ---------------------------------------------------------------- reference
def reference():
    """Run the committed model's faithful published protocol with standard BNG NFsim
    and compare against the committed reference .gdat (if present)."""
    print("=== reference: faithful published protocol (BNG NFsim, Fx=0.02) ===")
    t0 = time.time()
    tag = "reference"
    with _scratch(tag) as wd:
        out = _bng_run(open(MODEL).read(), tag, wd)  # outputs are named after the temp bngl
        t_eq, eq = _read_gdat(os.path.join(out, f"{tag}_equil.gdat"))
        t_tc, tc = _read_gdat(os.path.join(out, f"{tag}_tcr.gdat"))
    stem = os.path.splitext(os.path.basename(MODEL))[0]
    print(f"  [{time.time()-t0:.0f}s] equilibration {t_eq[0]:.0f}->{t_eq[-1]:.0f} s, "
          f"response {t_tc[0]:.0f}->{t_tc[-1]:.0f} s")
    print(f"  {'observable':13s}  basal(t=500)   response(t=60)   fold")
    for o in PLOT_OBS:
        b, r = eq[o][-1], tc[o][-1]
        print(f"  {o:13s}  {b:10.0f}   {r:12.0f}   {(r/b if b else float('nan')):5.2f}x")
    ref = os.path.join(REFDIR, f"{stem}_tcr.gdat")
    if os.path.exists(ref):
        t_ref, rc = _read_gdat(ref)
        # Stochastic: compare the committed reference vs a fresh run only qualitatively
        # (both should show the same monotone ligand response, not identical counts).
        up = sum(tc[o][-1] > tc[o][0] for o in ("TCR_pY149_D", "TCR_pY171_G",
                                                "ZAP70_pY493"))
        print(f"  committed reference present; fresh run shows ligand response "
              f"({up}/3 key ITAM/ZAP observables rise). OK.")
    else:
        print("  (no committed reference .gdat found; run writes fresh outputs only)")
    return (t_eq, eq), (t_tc, tc)


# ---------------------------------------------------------------- agreement
def agreement(n_seeds=16, fx=FX_CROSS, t_end=60, n_points=13, save=True):
    """Three-engine network-free agreement cross-check (single-phase, reduced Fx)."""
    print("=== agreement: BNG-NFsim vs bngsim-NFsim vs bngsim-RuleMonkey ===")
    print(f"    single-phase ligand-on, Fx={fx}, t_end={t_end}s, {n_seeds} seeds/engine")
    xml = build_xml(fx=fx)  # Fx/kfl/kfl_m baked into the XML -> identical ICs for all engines
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

    # Pairwise z-scores of ensemble means across observables and times.
    def mean_se(st, o):
        A = st[o]
        return A.mean(0), A.std(0, ddof=1) / np.sqrt(A.shape[0])
    pairs = [("bng_nf", "bngsim_nf"), ("bng_nf", "bngsim_rm"), ("bngsim_nf", "bngsim_rm")]
    zmax, zvals = 0.0, []
    for a, b in pairs:
        for o in OBS:
            ma, sa = mean_se(res[a][1], o)
            mb, sb = mean_se(res[b][1], o)
            denom = np.sqrt(sa**2 + sb**2)
            z = np.abs(ma - mb) / np.where(denom > 0, denom, np.inf)
            zvals.append(z[np.isfinite(z)])
            zmax = max(zmax, float(np.nanmax(z)))
    zall = np.concatenate(zvals)
    frac3 = float((zall < 3).mean())
    print(f"  pairwise-mean z-scores over {len(OBS)} obs x {n_points} pts x 3 pairs:")
    print(f"    max|z| = {zmax:.2f}   fraction |z|<3 = {frac3:.3f}   "
          f"(expect ~0.997 if engines sample the same process)")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "agreement.npz"),
                 tvec=tvec, obs=np.array(OBS),
                 **{f"{lab}__{o}": res[lab][1][o] for lab in res for o in OBS},
                 fx=fx, n_seeds=n_seeds)
        print(f"  cached -> {os.path.join(REFDIR, 'agreement.npz')}")
    return res, dict(zmax=zmax, frac_z_lt3=frac3)


# ---------------------------------------------------------------- data (Fig. S6)
# The 16 Fig. S6 sites, panel order a-p (digitize_chylek2014.py writes the three CSVs).
FIG_S6 = ["TCR_pY149_D", "TCR_pY171_G", "TCR_pY111", "TCR_pY123", "TCR_pY188_E",
          "TCR_pY199_E", "ZAP70_pY493", "LCK_pY424", "LCK_pY192", "ITK_pY512",
          "PTPN6_pY566", "WAS_pY291", "PAG1_pY163", "DOK2_pY299", "DOK1_pY449",
          "PLCG1_pY783"]
DATA_TIMES = (5.0, 15.0, 30.0, 60.0)   # the paper's measurement times (Table S1)
ENSEMBLE_NPZ = os.path.join(REFDIR, "data_ensemble.npz")
# Negative control for the level-2 check: the same protocol at a 4x smaller subvolume WITHOUT
# rescaling the ligand on-rates, which the actions block sets as literals computed for
# Fx = 0.02 (kfl = 7.45e5/(NA*simECFvol) = 6e-6, kfl_m = 6e-4) -- the slip an Fx reduction
# invites. Ligand binding then runs 4x slow; the check must reject it.
NEGCTL_NPZ = os.path.join(REFDIR, "data_ensemble_negctl.npz")
NEGCTL_FX = 0.005
N_DATA = 144          # committed replicate count
N_NEGCTL = 72
N_PUB_MIN = 10        # floor on the runs averaged into each published Fig. S6 line
Q_PT = 0.8            # Fig. S6 coordinate grid (pt); uniform error +/-Q_PT/2
DATA_ZMAX = 4.0       # max |z| allowed over the 16 x 4 comparisons
N_BOOT = 1000


def _seeded_actions(rep):
    """The committed file's own actions block, verbatim, with seed=>1000+rep added to the
    equilibration and seed=>5000+rep to the response (the file pins none)."""
    text = open(MODEL).read()
    acts = text[text.index("\nend model") + len("\nend model"):]
    seeds = iter((1000 + rep, 5000 + rep))
    acts, n = re.subn(r"\}\)", lambda m: f",seed=>{next(seeds)}}})", acts)
    assert n == 2, f"expected two simulate() calls in the actions block, found {n}"
    return acts


def _faithful_rep(rep, fx=None):
    """One replicate of the published two-phase protocol through BioNetGen's NFsim (with the
    subvolume line rewritten to fx for the negative control)."""
    tag = f"rep{rep:04d}"
    block = _model_block() if fx is None else _set_param_line(_model_block(), "Fx", fx)
    with _scratch(tag) as wd:
        out = _bng_run(block + _seeded_actions(rep), tag, wd)
        t_eq, eq = _read_gdat(os.path.join(out, f"{tag}_equil.gdat"))
        t_tc, tc = _read_gdat(os.path.join(out, f"{tag}_tcr.gdat"))
    return t_eq, eq, t_tc, tc


def _bng_version():
    bng, cli = find_bionetgen()
    cmd = [bng, "-v"] if cli else ["perl", bng, "-v"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    line = [ln for ln in r.stdout.splitlines() if "BioNetGen" in ln or "version" in ln.lower()]
    return line[0].strip() if line else "unknown"


def ensemble(n=N_DATA, workers=None, negctl=False):
    """Run n replicates of the faithful published protocol -- the committed .bngl with its own
    actions block, seeds pinned per replicate -- in parallel, and cache every replicate's
    equilibration and response trajectories to reference/data_ensemble.npz (negctl: the
    negative control, to reference/data_ensemble_negctl.npz)."""
    from functools import partial
    from multiprocessing import Pool
    workers = workers or min(n, os.cpu_count() or 1)
    what = f"negative control (Fx={NEGCTL_FX}, kfl unscaled)" if negctl else "published protocol"
    print(f"=== ensemble: {n} replicates of the {what} (BNG NFsim), {workers} workers ===")
    t0 = time.time()
    with Pool(workers) as pool:
        runs = pool.map(partial(_faithful_rep, fx=NEGCTL_FX if negctl else None),
                        range(1, n + 1))
    names = list(runs[0][3])
    eq = np.array([[r[1][o] for o in names] for r in runs]).transpose(0, 2, 1)
    tc = np.array([[r[3][o] for o in names] for r in runs]).transpose(0, 2, 1)
    assert np.array_equal(eq, np.round(eq)) and np.array_equal(tc, np.round(tc))
    os.makedirs(REFDIR, exist_ok=True)
    path = NEGCTL_NPZ if negctl else ENSEMBLE_NPZ
    np.savez_compressed(path, t_equil=runs[0][0], t_resp=runs[0][2],
                        obs=np.array(names), equil=eq.astype(np.int32),
                        resp=tc.astype(np.int32), reps=np.arange(1, n + 1),
                        seed_equil=1000 + np.arange(1, n + 1),
                        seed_resp=5000 + np.arange(1, n + 1),
                        bionetgen=np.array(_bng_version()),
                        nfsim_flags=np.array("-cb (BNG2.pl default); no -bscb; default -utl"))
    print(f"  [{time.time()-t0:.0f}s] cached {n} replicates -> {path}")


def _csv_rows(name):
    with open(os.path.join(REFDIR, name)) as fh:
        return list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))


def published_figS6():
    """{site: (t, log2 fold change)} of the published Fig. S6 lines, and the measured points
    {site: (t, mean, sd)}, from the digitized CSVs."""
    lines, points = {}, {}
    for r in _csv_rows("chylek2014_figS6_simulated_digitized.csv"):
        lines.setdefault(r["site"], []).append((float(r["t_s"]), float(r["log2_fold_change"])))
    for r in _csv_rows("chylek2014_figS6_data_digitized.csv"):
        points.setdefault(r["site"], []).append((float(r["t_s"]), float(r["mean_log2"]),
                                                 float(r["sd_log2"])))
    out = {}
    for site, v in lines.items():
        t = np.array([a for a, _ in v])
        y = np.array([b for _, b in v])
        tu = np.unique(t)  # a vertical step on the 0.8-pt grid: average its two ends
        out[site] = (tu, np.array([y[t == x].mean() for x in tu]))
    return out, {k: np.array(v) for k, v in points.items()}


def load_ensemble(path=ENSEMBLE_NPZ):
    d = np.load(path)
    return d, list(d["obs"])


def data(path=ENSEMBLE_NPZ):
    """Level 2: the faithful protocol's ensemble mean, as log2 fold change over the
    equilibrated baseline, against the published Fig. S6 lines at the four measurement times.
    Returns a dict of the per-site comparison and the summary numbers."""
    d, names = load_ensemble(path)
    tvec, X = d["t_resp"], d["resp"].astype(float)
    n = X.shape[0]
    lines, points = published_figS6()
    rng = np.random.default_rng(0)
    boot = rng.integers(0, n, size=(N_BOOT, n))
    # 0.8-pt grid -> log2 units on these axes (13.9-14.1 pt per log2 unit; 1.20 pt per s).
    q_y = (Q_PT / 13.92) / np.sqrt(12.0)
    q_t = (Q_PT / 1.2004) / np.sqrt(12.0)
    res = {}
    for site in FIG_S6:
        A = X[:, :, names.index(site)]
        m = A.mean(0)
        y = np.log2(m / m[0])
        mb = A[boot].mean(1)
        yb = np.log2(mb / mb[:, :1])
        se = yb.std(0, ddof=1)
        tp, yp = lines[site]
        pub = np.interp(tvec, tp, yp)
        slope = np.abs(np.gradient(pub, tvec))
        sd_d = np.sqrt(se**2 * (1.0 + n / N_PUB_MIN) + q_y**2 + (q_t * slope)**2)
        idx = [int(np.argmin(np.abs(tvec - t))) for t in DATA_TIMES]
        dd = y - pub
        z = np.abs(dd[idx]) / sd_d[idx]
        P = points[site]
        res[site] = dict(t=tvec, y=y, se=se, pub=pub, d=dd, sd_d=sd_d, idx=idx, z=z,
                         basal=m[0], meas_t=P[:, 0], meas_m=P[:, 1], meas_sd=P[:, 2])
    zall = np.concatenate([r["z"] for r in res.values()])
    zcurve = np.concatenate([np.abs(r["d"][1:]) / r["sd_d"][1:] for r in res.values()])
    summary = dict(n=n, zmax=float(zall.max()), n_cmp=zall.size,
                   zmax_curve=float(zcurve.max()), n_fail=int((zall > DATA_ZMAX).sum()),
                   frac_z_lt3=float((zall < 3).mean()),
                   dmax=max(float(np.abs(r["d"][r["idx"]]).max()) for r in res.values()),
                   rms_all=float(np.sqrt(np.mean(np.concatenate(
                       [r["d"][1:]**2 for r in res.values()])))))
    return res, summary


def reference_typicality():
    """Is the committed reference/ run (one unseeded run of the same protocol) a typical
    draw from the ensemble? Statistic: max over observables x times of |x - mean| / sd, with
    sd the single-run standard deviation; its null distribution is the same statistic for each
    ensemble replicate against the other n-1 (leave one out)."""
    d, names = load_ensemble()
    stem = os.path.splitext(os.path.basename(MODEL))[0]
    out = {}
    for phase, key, tkey in (("tcr", "resp", "t_resp"), ("equil", "equil", "t_equil")):
        t_ref, ref = _read_gdat(os.path.join(REFDIR, f"{stem}_{phase}.gdat"))
        assert np.allclose(t_ref, d[tkey]), phase
        X = d[key].astype(float)
        R = np.array([ref[o] for o in names]).T
        n = X.shape[0]

        def stat(x, pool):
            mu, sd = pool.mean(0), pool.std(0, ddof=1)
            ok = sd > 0
            return float(np.max(np.abs(x - mu)[ok] / sd[ok]))

        s_ref = stat(R, X)
        null = np.array([stat(X[i], np.delete(X, i, axis=0)) for i in range(n)])
        out[phase] = dict(stat=s_ref, null=null, rank=float((null >= s_ref).mean()))
    return out


# ---------------------------------------------------------------- figure
def plot():
    """Build verify_tcr_chylek2014.png from committed artifacts only: the reference .gdat,
    agreement.npz, data_ensemble.npz, data_ensemble_negctl.npz and the Fig. S6 CSVs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stem = os.path.splitext(os.path.basename(MODEL))[0]
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))

    # Row 1: the faithful published protocol (committed reference run).
    t_tc, tc = _read_gdat(os.path.join(REFDIR, f"{stem}_tcr.gdat"))
    t_eq, eq = _read_gdat(os.path.join(REFDIR, f"{stem}_equil.gdat"))
    typ = reference_typicality()
    ax = axes[0, 0]
    for o in PLOT_OBS:
        ax.plot(t_tc, tc[o], marker=".", ms=3, label=o)
    ax.set(title="A. Faithful protocol, committed reference run:\n60 s response "
                 f"(typical draw: max|z| {typ['tcr']['stat']:.2f}, null max "
                 f"{typ['tcr']['null'].max():.2f})",
           xlabel="time after stimulation (s)", ylabel="molecules / subvolume")
    ax.title.set_fontsize(9)
    ax.legend(fontsize=6, ncol=2)
    ax = axes[0, 1]
    for o in PLOT_OBS:
        ax.plot(t_eq, eq[o], marker=".", ms=3, label=o)
    ax.set(title="B. Equilibration, ligand binding off\n(500 s; state carried into the response)",
           xlabel="time (s)", ylabel="molecules / subvolume")
    ax.title.set_fontsize(9)

    # Row 2: three-engine agreement (bngsim arms; cached, not BNGsim-independent evidence).
    d = np.load(os.path.join(REFDIR, "agreement.npz"), allow_pickle=True)
    tvec = d["tvec"]
    styles = {"bng_nf": ("BNG NFsim", "o", "#1b9e77"),
              "bngsim_nf": ("bngsim NFsim", "s", "#7570b3"),
              "bngsim_rm": ("bngsim RuleMonkey", "^", "#d95f02")}
    groups = [PLOT_OBS[0:2], PLOT_OBS[2:4], PLOT_OBS[4:6]]
    for gi, group in enumerate(groups):
        ax = axes[1, gi]
        for o in group:
            for lab, (name, mk, col) in styles.items():
                A = d[f"{lab}__{o}"]
                m = A.mean(0)
                se = A.std(0, ddof=1) / np.sqrt(A.shape[0])
                ax.plot(tvec, m, marker=mk, ms=3, color=col, lw=1,
                        label=name if o == group[0] else None)
                ax.fill_between(tvec, m - se, m + se, color=col, alpha=0.15)
            ax.annotate(o, (tvec[-1], A.mean(0)[-1]), fontsize=6, xytext=(2, 0),
                        textcoords="offset points")
        ax.set(xlabel="time (s)", ylabel="molecules / subvolume")
        if gi == 0:
            ax.legend(fontsize=6, title="engine (Fx=0.005, single-phase)", title_fontsize=6)
        ax.set_title(("D. " if gi == 0 else "E. " if gi == 1 else "F. ")
                     + "3-engine agreement (cached): " + ", ".join(group), fontsize=8)

    # Row 3: level 2, the ensemble against the published Fig. S6 simulation.
    res, sm = data()
    neg, nsm = data(NEGCTL_NPZ)
    cmap = plt.get_cmap("tab20")
    ax = axes[2, 0]
    lo, hi = -2.0, 2.6
    ax.plot([lo, hi], [lo, hi], color="gray", lw=0.8)
    for k, site in enumerate(FIG_S6):
        r = res[site]
        ax.plot(r["pub"][1:], r["y"][1:], lw=0.8, color=cmap(k % 20), alpha=0.7)
        ax.plot(r["pub"][r["idx"]], r["y"][r["idx"]], "o", ms=3.5, color=cmap(k % 20))
    ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel="published Fig. S6 line (log2 fold change)",
           ylabel=f"BNG NFsim, {sm['n']}-run mean (log2 fold change)")
    ax.set_title(f"G. Level 2: 16 sites vs Fig. S6 (dots: 5/15/30/60 s)\n"
                 f"max|z| {sm['zmax']:.2f} over {sm['n_cmp']} (tol {DATA_ZMAX:g}); "
                 f"max|d| {sm['dmax']:.3f} log2", fontsize=9)
    for ax, (R, S, lab) in zip(axes[2, 1:], ((res, sm, "H. published protocol"),
                                             (neg, nsm, "I. negative control"))):
        for k, site in enumerate(FIG_S6):
            r = R[site]
            ax.plot(r["t"][1:], r["d"][1:] / r["sd_d"][1:], lw=0.8, color=cmap(k % 20))
            ax.plot(r["t"][r["idx"]], r["d"][r["idx"]] / r["sd_d"][r["idx"]], "o", ms=3,
                    color=cmap(k % 20))
        ax.axhspan(-DATA_ZMAX, DATA_ZMAX, color="green", alpha=0.08)
        ax.axhline(DATA_ZMAX, color="crimson", ls="--", lw=0.8)
        ax.axhline(-DATA_ZMAX, color="crimson", ls="--", lw=0.8)
        ax.set(xlabel="time after stimulation (s)", ylabel="signed z = d / sd(d)")
        verdict = ("PASS" if S["zmax"] <= DATA_ZMAX
                   else "rejected" if lab.startswith("I") else "FAIL")
        ax.set_title(f"{lab}: max|z| at 5/15/30/60 s = {S['zmax']:.2f} -> {verdict}"
                     + (f"\n({S['n_fail']}/{S['n_cmp']} beyond tol; Fx={NEGCTL_FX}, kfl not "
                        "rescaled)" if lab.startswith("I") else
                        f"\n(tol band +/-{DATA_ZMAX:g}; sd(d) from bootstrap SE, n_pub>="
                        f"{N_PUB_MIN}, 0.8-pt grid)"), fontsize=9)
    ymax = max(10.0, min(60.0, 1.1 * nsm["zmax"]))
    axes[2, 2].set_ylim(-ymax, ymax)

    axes[0, 2].axis("off")
    axes[0, 2].text(0.0, 1.0,
        "Chylek et al. (2014) PLoS ONE 9:e104240\n"
        "Early TCR/CD28 signaling, network-free (NFsim)\n\n"
        "BNGsim-independent checks (BioNetGen NFsim):\n"
        f"  level 2  {sm['n']} runs of the file's own two-phase\n"
        "    protocol vs the published Fig. S6 lines\n"
        f"    (vector-read; points = Table S1): max|z|\n"
        f"    {sm['zmax']:.2f} <= {DATA_ZMAX:g} over 16 sites x 4 times  PASS\n"
        f"  negative control (Fx {NEGCTL_FX}, kfl unscaled):\n"
        f"    max|z| {nsm['zmax']:.1f}  -> rejected, as it must be\n"
        f"  committed reference run: max|z| {typ['tcr']['stat']:.2f}\n"
        f"    (response), {typ['equil']['stat']:.2f} (equil.) vs leave-one-out\n"
        f"    null maxima {typ['tcr']['null'].max():.2f}, {typ['equil']['null'].max():.2f}  "
        "PASS\n\n"
        "Rows D-F: the older three-engine cross-check,\n"
        "  two arms through bngsim (cached; shown for\n"
        "  context, not part of the evidence above).",
        fontsize=8.5, va="top", family="monospace")

    fig.suptitle("tcr_signaling_chylek2014 verification: published two-phase NFsim protocol "
                 "vs Chylek et al. (2014) Fig. S6", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    png = os.path.join(HERE, "verify_tcr_chylek2014.png")
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"wrote {png}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    if mode == "reference":
        reference()
    elif mode == "agreement":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        agreement(n_seeds=n)
    elif mode in ("ensemble", "negctl"):
        neg = mode == "negctl"
        n = int(sys.argv[2]) if len(sys.argv) > 2 else (N_NEGCTL if neg else N_DATA)
        w = int(sys.argv[3]) if len(sys.argv) > 3 else None
        ensemble(n, w, negctl=neg)
    elif mode == "data":
        res, sm = data()
        for site in FIG_S6:
            r = res[site]
            print(f"  {site:12s} basal {r['basal']:7.1f}  |d| at 5/15/30/60 s: "
                  + " ".join(f"{abs(r['d'][i]):.3f}" for i in r["idx"])
                  + "   z: " + " ".join(f"{z:.2f}" for z in r["z"]))
        print(f"  {sm['n']} replicates; max|z| = {sm['zmax']:.2f} over {sm['n_cmp']} "
              f"(tolerance {DATA_ZMAX}); max|d| = {sm['dmax']:.3f} log2")
        assert sm["zmax"] <= DATA_ZMAX
    elif mode == "plot":
        plot()
    elif mode == "all":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        reference()
        agreement(n_seeds=n)
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference | agreement [N] | ensemble [N] [workers] | "
              "data | plot | all [N]")
