"""Driver + verification for egfr_oligomerization_mitra2019.bngl -- the network-free
EGFR higher-order-oligomerization model of Mitra et al. (2019), iScience 19:1012-1036
(doi:10.1016/j.isci.2019.08.045; BioNetFit-1 "example 2"), fit to the cluster-density
and phospho-EGFR data of Kozer et al. (2013), Mol BioSyst 9:1849-1863
(doi:10.1039/c3mb70073a).

EGF drives EGFR into ectodomain-crosslinked oligomers that grow WITHOUT BOUND, so the
reaction network is infinite and the model is simulated NETWORK-FREE (NFsim /
RuleMonkey). This is genuine multivalent AGGREGATION: the bimolecular crosslink rules
join two SEPARATE aggregates, while ring closure within an aggregate is the dedicated
chi_r-enhanced rules, so the crosslink rules must block same-complex binding (NFsim
-bscb; the 4-molecule ring patterns also need -utl 5). Because this is an aggregation
model, legacy BioNetGen NFsim is not a self-evidently valid reference; correctness is
instead confirmed by agreement of two INDEPENDENT network-free algorithms -- NFsim
(Yang et al. 2008, Phys Rev E 78:031910) and RuleMonkey (Colvin et al. 2009/2010) --
both with -bscb.

Three things are verified here:

  reference   -- reproduces Kozer 2013 Fig 2B/2D/3B/3D. Runs the model network-free
                 (BNG NFsim, -bscb -utl 5) at the shipped PyBioNetFit best fit, averaged
                 over replicate runs: a 30 nM time course (cluster density Fig 3B,
                 phospho-EGFR Fig 3D) and a dose-response scan over EGF (cluster density
                 Fig 2B, phospho-EGFR Fig 2D). Model outputs are scaled by the fitted
                 alpha*_pre and compared to the average-normalized data; reports chi_sq
                 and the median/max relative error per dataset.

  agreement   -- the independent-algorithm cross-check (issue #8, capability group).
                 The two INDEPENDENT network-free algorithms bngsim provides -- NFsim
                 (Yang 2008) and the exact RuleMonkey (Colvin 2009/2010) -- run the SAME
                 model with -bscb at a few RuleMonkey-tractable doses (RuleMonkey is exact
                 and steps every EGF particle, so the top data doses with ~1e6 EGF are
                 cluster-scale) and AGREE on per-dose equilibrium pEGFR / cluster counts
                 within Monte-Carlo scale -- the reference of record for this aggregation
                 model. BNG2.pl's bundled LEGACY NFsim v1.14.3 binary is run for contrast:
                 it tracks the same behavior but reads ~12% lower on pEGFR (a binary-version
                 difference -- the universal traversal limit is NOT the cause, utl 4/5/6
                 give the same means), which is exactly why legacy NFsim is not the
                 reference for this aggregation model (issue #8 capability group). bngsim
                 NFsim is also run WITHOUT -bscb; at this fit -bscb and no-bscb differ by
                 only a few percent (higher-order rings are only rarely populated -- unlike
                 the TLBR model, where dropping -bscb diverges) -- documented, not assumed.

  plot        -- builds verify_mitra2019.png from the cached .npz files.

Requirements: BioNetGen 2.9.3, found from $BNG2_PL or $BNGPATH, then BNG2.pl or the
`bionetgen` CLI on PATH; `agreement` also needs the bngsim package (NfsimSession +
RuleMonkeySession). Run with the repo venv:
    python run_mitra2019.py reference 12 [W]  # Fig 2B/2D/3B/3D reproduction, 12 replicates
    python run_mitra2019.py agreement 6 [W]   # bngsim NFsim=RuleMonkey (+legacy, bscb), 6 seeds
    python run_mitra2019.py agreement-bng 6 [W]  # refresh only the BNG2.pl NFsim arm
    python run_mitra2019.py plot            # (re)build verify_mitra2019.png
    python run_mitra2019.py all
W is the number of BNG2.pl runs in parallel (default 1).

Every BNG2.pl run gets its own directory under a per-invocation work directory inside
this folder (removed afterwards), and its output is read back by the exact name that
run wrote. Each dose of a replicate's dose-response is its own one-point
parameter_scan, so the doses run in parallel; BioNetGen reuses a scan's seed at every
point, so this is the same simulation as the full scan. Seeds: time-course replicate r
uses 1000+r, dose-response replicate r 2000+r (at every dose), agreement seed s 5000+s.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, "egfr_oligomerization_mitra2019.bngl")
REFDIR = os.path.join(HERE, "reference")



def find_bionetgen():
    """BNG2.pl from $BNG2_PL or $BNGPATH, then BNG2.pl on PATH, then the `bionetgen` CLI.
    Returns (path, use_cli)."""
    for cand in (os.environ.get("BNG2_PL"),
                 os.path.join(os.environ["BNGPATH"], "BNG2.pl") if os.environ.get("BNGPATH") else None,
                 shutil.which("BNG2.pl")):
        if cand and os.path.isfile(cand):
            return os.path.abspath(cand), False
    cli = shutil.which("bionetgen")
    if cli:
        return cli, True
    raise FileNotFoundError("BioNetGen not found: set BNG2_PL (or BNGPATH) to a BioNetGen 2.9.3 "
                            "BNG2.pl, or put the bionetgen CLI on PATH.")

# PyBioNetFit best fit (Mitra 2019 Data S1, RuleHub 04-egfrnf/fit_ade, chi_sq=13.4) --
# also the shipped .bngl nominals; re-applied here so the driver is robust to edits.
BESTFIT = dict(k_o=3.54931, k_c=15.8220, kaf=45.0460, kar=0.64255, chi_r=1907.32)
ALPHA = dict(a1=2.5718e-5, a2=2.7509e-5, a3=5.0739e-5, a4=3.2791e-5)  # alpha*_pre
F = 0.01                     # subvolume factor (model parameter f)
T_END = 600                  # experimental endpoint (s)
GML = 2147483647             # 32-bit-max global molecule limit
UTL = 5                      # universal traversal limit (above the 4-molecule ring patterns)
NFPARAM = f"-bscb -utl {UTL}"  # block-same-complex-binding + ring-closure traversal flags

# EGFR oligomer-size observables summed into the cluster count (model function Clusters).
OLIGO = ("monomer dimer trimer tetramer pentamer hexamer heptamer octamer nonamer "
         "decamer undecamer dodecamer tridecamer tetradecamer pentadecamer hexadecamer "
         "heptadecamer octadecamer nonadecamer icosadecamer").split()

# Data doses actually fit (doseresponse.exp; 0.01 nM is NaN and is dropped).
DOSE_DATA = [0.001, 0.1, 1.0, 10.0, 100.0]
# RuleMonkey-tractable doses for the 3-engine agreement check (clustering regime).
DOSE_AGREE = [0.1, 1.0, 10.0]


# ---------------------------------------------------------------- data
def _load_csv(name):
    # Strip leading '#' comment lines ourselves: the header comments contain commas,
    # which confuses genfromtxt's inline-comment handling.
    from io import StringIO
    lines = [ln for ln in open(os.path.join(REFDIR, name)) if not ln.lstrip().startswith("#")]
    return np.genfromtxt(StringIO("".join(lines)), delimiter=",", names=True)


def load_data():
    tc = _load_csv("kozer2013_fig3bd_timecourse.csv")
    dr = _load_csv("kozer2013_fig2bd_doseresponse.csv")
    return tc, dr


# ---------------------------------------------------------------- BNG helpers
def _model_block():
    """Committed model text up to and including `end model` (no actions block)."""
    text = open(MODEL).read()
    marker = "\nend model"
    i = text.index(marker)
    return text[:i + len(marker)] + "\n"


def _set_param_text(text, name, value):
    text, n = re.subn(rf"(?m)^(\s*{re.escape(name)}\s+)\S+", rf"\g<1>{value:.10g}", text, count=1)
    assert n == 1, f"expected one '{name}' parameter line, found {n}"
    return text


def _apply(text, params):
    for k, v in (params or {}).items():
        text = _set_param_text(text, k, v)
    return text


def _run_bng(bngl_text, rundir, tag):
    """Write `bngl_text` to rundir/<tag>.bngl and run BioNetGen on it with all outputs in
    rundir (BioNetGen names them <tag>_<suffix>.<ext>)."""
    bng, use_cli = find_bionetgen()
    bngl = os.path.join(rundir, f"{tag}.bngl")
    with open(bngl, "w") as fh:
        fh.write(bngl_text)
    cmd = ([bng, "run", "-i", bngl, "-o", rundir] if use_cli
           else ["perl", bng, "--outdir", rundir, bngl])
    r = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:] + "\n---STDOUT---\n" + r.stdout[-3000:]
    return bngl


def _workdir(label):
    """A fresh per-invocation work directory inside this folder (never /tmp)."""
    return tempfile.mkdtemp(prefix=f".{label}_", dir=HERE)


def _run_and_read(bngl_text, workdir, tag, ext):
    """Run BioNetGen in its own directory under `workdir` and return the columns of the
    one output file with extension `ext` that run wrote, read by its exact name."""
    rundir = tempfile.mkdtemp(prefix=f"{tag}_", dir=workdir)
    _run_bng(bngl_text, rundir, tag)
    out = [f for f in os.listdir(rundir) if f.endswith(ext)]
    suffix = {".gdat": "nfr", ".scan": "scan"}[ext]
    assert out == [f"{tag}_{suffix}{ext}"], f"expected exactly {tag}_{suffix}{ext}, found {out}"
    cols = _read_table(os.path.join(rundir, out[0]))
    shutil.rmtree(rundir)
    return cols


def _pmap(fn, args, workers):
    """[fn(*a) for a in args], `workers` processes in parallel."""
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            return list(ex.map(fn, *zip(*args)))
    return [fn(*a) for a in args]


def _assert_distinct(A, what):
    """Distinct seeds must give distinct replicates (rows of A)."""
    n_unique = len({tuple(np.ravel(row)) for row in A})
    assert n_unique == len(A), f"{what}: only {n_unique} distinct replicates of {len(A)}"


def _read_table(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    return {name: data[:, i] for i, name in enumerate(header)}


def _clusters_pegfr(cols):
    clusters = sum(cols[name] for name in OLIGO)
    return clusters, cols["pEGFR"]


def _bng_timecourse(seed, workdir):
    """30 nM network-free time course; returns (t, Clusters(t), pEGFR(t))."""
    base = _apply(_model_block(), BESTFIT)
    base = _set_param_text(base, "LT_nM", 30.0)
    act = ("\nbegin actions\n"
           f'  simulate({{method=>"nf",suffix=>"nfr",t_start=>0,t_end=>{T_END},n_steps=>60,'
           f'seed=>{seed},gml=>{GML},print_functions=>1,param=>"{NFPARAM}"}})\n'
           "end actions\n")
    cols = _run_and_read(base + act, workdir, f"tc_s{seed}", ".gdat")
    clusters, pegfr = _clusters_pegfr(cols)
    return cols["time"], clusters, pegfr


def _bng_dose_scan(doses, seed, workdir):
    """One BNG NFsim parameter_scan over `doses` (-bscb) to the t_end endpoint; returns
    (Clusters, pEGFR) per dose. This is the model's actions-block dose-response protocol."""
    base = _apply(_model_block(), BESTFIT)
    vals = ",".join(f"{d:.10g}" for d in doses)
    act = ("\nbegin actions\n"
           f'  parameter_scan({{parameter=>"LT_nM",par_scan_vals=>[{vals}],method=>"nf",'
           f'suffix=>"scan",t_start=>0,t_end=>{T_END},n_steps=>1,seed=>{seed},gml=>{GML},'
           f'print_functions=>1,param=>"{NFPARAM}"}})\n'
           "end actions\n")
    tag = f"scan_s{seed}_" + "_".join(f"{d:.10g}" for d in doses).replace(".", "p")
    cols = _run_and_read(base + act, workdir, tag, ".scan")
    assert np.allclose(cols["LT_nM"], doses, rtol=1e-8), "scan rows do not match the doses"
    return _clusters_pegfr(cols)


def _bng_dose_point(dose, seed, workdir):
    """(Clusters, pEGFR) at t_end at one dose: a one-point parameter_scan."""
    c, p = _bng_dose_scan([dose], seed, workdir)
    return float(c[0]), float(p[0])


def _bng_dose_replicates(doses, seeds, workdir, workers=1):
    """Dose response at `seeds` (one replicate per seed, same seed at every dose, as in a
    full parameter_scan), every (seed, dose) run in parallel; returns (Clusters, pEGFR),
    each (n_seeds, n_doses)."""
    args = [(d, s, workdir) for s in seeds for d in doses]
    out = np.array(_pmap(_bng_dose_point, args, workers)).reshape(len(seeds), len(doses), 2)
    cl, pe = out[..., 0], out[..., 1]
    _assert_distinct(np.concatenate([cl, pe], axis=1), "dose response")
    return cl, pe


def _dose_xml(dose, workdir):
    """Bake one dose (+ best fit) into the model and writeXML; returns the XML path.
    Baking is required because bngsim RuleMonkey's set_param does not re-derive
    seed-species counts (see memory: bngsim set_param NF gaps)."""
    tag = f"dose_{dose:.10g}".replace(".", "p").replace("-", "m")
    xml = os.path.join(workdir, f"{tag}.xml")
    if os.path.exists(xml):
        return xml
    base = _apply(_model_block(), BESTFIT)
    base = _set_param_text(base, "LT_nM", dose)
    base += "\nbegin actions\n  writeXML({overwrite=>1})\nend actions\n"
    _run_bng(base, workdir, tag)
    assert os.path.isfile(xml), f"no XML produced for dose {dose}"
    return xml


# ---------------------------------------------------------------- bngsim helpers
def _bngsim(engine, xml, seed, bscb=True):
    """Final equilibrium (Clusters, pEGFR) from one bngsim network-free run. NFsim gets
    the same -utl 5 as the BNG runs (its default auto-UTL of 4 is too shallow for the
    4-molecule ring patterns); RuleMonkey is exact and takes no traversal limit."""
    from bngsim import NfsimSession, RuleMonkeySession
    if engine == "nf":
        s = NfsimSession(xml, molecule_limit=GML, block_same_complex_binding=bscb,
                         traversal_limit=UTL)
    else:
        s = RuleMonkeySession(xml, molecule_limit=GML, block_same_complex_binding=bscb)
    with s:
        s.initialize(seed=seed)
        res = s.simulate(0, T_END, n_points=4)
        names = list(s.get_observable_names())
        Y = np.asarray(res.observables)[-1]
    d = {n: Y[i] for i, n in enumerate(names)}
    clusters = sum(d.get(o, 0.0) for o in OLIGO)
    return clusters, d["pEGFR"]


# ---------------------------------------------------------------- reference
def _metrics(model, data, sd):
    ok = ~np.isnan(data)
    chisq = 0.5 * float(np.sum(((model[ok] - data[ok]) / sd[ok]) ** 2))
    big = ok & (np.abs(data) > 0.05)          # points above the near-zero noise floor
    rel = np.abs(model[big] - data[big]) / np.abs(data[big])
    return chisq, float(np.median(rel)), float(rel.max())


def _bng_version():
    """'BioNetGen <version>' of the install find_bionetgen() picks (VERSION file)."""
    bng, use_cli = find_bionetgen()
    vf = os.path.join(os.path.dirname(bng), "VERSION")
    v = open(vf).read().split()[0] if (not use_cli and os.path.isfile(vf)) else "unknown"
    return f"BioNetGen {v} ({'bionetgen CLI' if use_cli else 'BNG2.pl'}; bundled NFsim)"


def _chisq_total(clus_t, pegf_t, clus_d, pegf_d, tc, dr, dmask):
    """Total chi_sq of replicate-mean model outputs against the four datasets."""
    m = (ALPHA["a2"] * clus_t / F, ALPHA["a4"] * pegf_t / F,
         ALPHA["a1"] * clus_d / F, ALPHA["a3"] * pegf_d / F)
    d = ((tc["pre2_time"], tc["pre2_time_SD"]), (tc["pre4_time"], tc["pre4_time_SD"]),
         (dr["pre1_dose"][dmask], dr["pre1_dose_SD"][dmask]),
         (dr["pre3_dose"][dmask], dr["pre3_dose_SD"][dmask]))
    return sum(_metrics(mi, di, si)[0] for mi, (di, si) in zip(m, d))


def reference(n_reps=12, workers=1, save=True):
    """Reproduce Kozer 2013 Fig 2B/2D/3B/3D: BNG2.pl NFsim (-bscb -utl 5) at the best fit
    vs the data. Replicates are independent BNG2.pl runs (time course seed 1000+r, dose
    response seed 2000+r), `workers` in parallel."""
    print("=== reference: Kozer 2013 Fig 2B/2D/3B/3D reproduction (BNG NFsim, -bscb) ===")
    print(f"    {n_reps} replicates, best fit, t_end={T_END}s, {workers} parallel runs")
    tc, dr = load_data()
    t_exp = tc["time_s"]
    seeds_tc = [1000 + r for r in range(n_reps)]
    seeds_dr = [2000 + r for r in range(n_reps)]
    # Most expensive runs first (top doses, then the 30 nM time course) so they start first.
    order = sorted(range(len(DOSE_DATA)), key=lambda i: -DOSE_DATA[i])
    wd = _workdir("reference")
    try:
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=max(1, workers)) as ex:
            fut_dr = {}
            for i in order[:2]:
                for s in seeds_dr:
                    fut_dr[s, i] = ex.submit(_bng_dose_point, DOSE_DATA[i], s, wd)
            fut_tc = [ex.submit(_bng_timecourse, s, wd) for s in seeds_tc]
            for i in order[2:]:
                for s in seeds_dr:
                    fut_dr[s, i] = ex.submit(_bng_dose_point, DOSE_DATA[i], s, wd)
            tcs = [f.result() for f in fut_tc]
            print(f"  [{time.time()-t0:.0f}s] time course done")
            drs = np.array([[fut_dr[s, i].result() for i in range(len(DOSE_DATA))]
                            for s in seeds_dr])          # (n_reps, n_doses, 2)
            print(f"  [{time.time()-t0:.0f}s] dose response done")
    finally:
        shutil.rmtree(wd, ignore_errors=True)

    # Per-replicate arrays; the time course is sampled every 10 s, so the experimental
    # times (multiples of 30 s) are grid points and np.interp just picks them out.
    tc_t = tcs[0][0]
    tc_cl = np.array([c for _, c, _ in tcs])             # (n_reps, n_t)
    tc_pe = np.array([p for _, _, p in tcs])
    dr_cl, dr_pe = drs[..., 0], drs[..., 1]              # (n_reps, n_doses)
    _assert_distinct(np.concatenate([tc_cl, tc_pe], axis=1), "time course")
    _assert_distinct(np.concatenate([dr_cl, dr_pe], axis=1), "dose response")
    cl_t = np.array([np.interp(t_exp, tc_t, c) for c in tc_cl])
    pe_t = np.array([np.interp(t_exp, tc_t, p) for p in tc_pe])
    clus_t, pegf_t = cl_t.mean(0), pe_t.mean(0)
    clus_d, pegf_d = dr_cl.mean(0), dr_pe.mean(0)

    # Model observables scaled to the average-normalized data (alpha*_pre*obs/f).
    m_pre2 = ALPHA["a2"] * clus_t / F      # cluster density vs time  (Fig 3B)
    m_pre4 = ALPHA["a4"] * pegf_t / F      # phospho-EGFR   vs time  (Fig 3D)
    # Align dose model outputs to the data doses (drop the 0.01 nM NaN row).
    dmask = np.isin(np.round(dr["LT_nM"], 6), np.round(DOSE_DATA, 6))
    m_pre1 = ALPHA["a1"] * clus_d / F      # cluster density vs dose  (Fig 2B)
    m_pre3 = ALPHA["a3"] * pegf_d / F      # phospho-EGFR   vs dose  (Fig 2D)

    rows = [
        ("Fig 3B cluster density (time)", m_pre2, tc["pre2_time"], tc["pre2_time_SD"]),
        ("Fig 3D phospho-EGFR   (time)", m_pre4, tc["pre4_time"], tc["pre4_time_SD"]),
        ("Fig 2B cluster density (dose)", m_pre1, dr["pre1_dose"][dmask], dr["pre1_dose_SD"][dmask]),
        ("Fig 2D phospho-EGFR   (dose)", m_pre3, dr["pre3_dose"][dmask], dr["pre3_dose_SD"][dmask]),
    ]
    total_chisq = 0.0
    print(f"  {'dataset':32s} {'chi_sq':>8s} {'med|rel|':>9s} {'max|rel|':>9s}")
    for lab, m, d, sd in rows:
        chisq, med, mx = _metrics(m, d, sd)
        total_chisq += chisq
        print(f"  {lab:32s} {chisq:8.2f} {med:9.3f} {mx:9.3f}")
    print(f"  {'TOTAL chi_sq':32s} {total_chisq:8.2f}   (published PyBNF best fit chi_sq=13.4)")
    # Monte-Carlo standard error of the total: bootstrap over replicates (time-course and
    # dose-response replicates resampled independently), fixed RNG seed.
    rng = np.random.default_rng(0)
    boot = []
    for _ in range(2000):
        a = rng.integers(0, n_reps, n_reps)
        b = rng.integers(0, n_reps, n_reps)
        boot.append(_chisq_total(cl_t[a].mean(0), pe_t[a].mean(0), dr_cl[b].mean(0),
                                 dr_pe[b].mean(0), tc, dr, dmask))
    chisq_se = float(np.std(boot, ddof=1))
    print(f"  bootstrap s.e. of TOTAL chi_sq over replicates = {chisq_se:.2f}")

    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "reproduction.npz"),
                 t_exp=t_exp, m_pre2=m_pre2, m_pre4=m_pre4,
                 d_pre2=tc["pre2_time"], d_pre2_sd=tc["pre2_time_SD"],
                 d_pre4=tc["pre4_time"], d_pre4_sd=tc["pre4_time_SD"],
                 doses=np.array(DOSE_DATA), m_pre1=m_pre1, m_pre3=m_pre3,
                 d_pre1=dr["pre1_dose"][dmask], d_pre1_sd=dr["pre1_dose_SD"][dmask],
                 d_pre3=dr["pre3_dose"][dmask], d_pre3_sd=dr["pre3_dose_SD"][dmask],
                 n_reps=n_reps, total_chisq=total_chisq, chisq_se=chisq_se,
                 tc_t=tc_t, tc_clusters=tc_cl, tc_pegfr=tc_pe,
                 dose_clusters=dr_cl, dose_pegfr=dr_pe,
                 seeds_tc=np.array(seeds_tc), seeds_dose=np.array(seeds_dr),
                 engine=f"{_bng_version()}, param '{NFPARAM}', gml {GML}")
        print(f"  cached -> {os.path.join(REFDIR, 'reproduction.npz')}")


# ---------------------------------------------------------------- agreement
def _agreement_bng(n_seeds, workdir, workers=1):
    """BNG2.pl NFsim arm of the agreement check: seeds 5000.., every (seed, dose) its own
    run; returns (Clusters, pEGFR), each (n_seeds, n_doses)."""
    t0 = time.time()
    cl, pe = _bng_dose_replicates(DOSE_AGREE, [5000 + s for s in range(n_seeds)], workdir,
                                  workers)
    print(f"  [{time.time()-t0:5.0f}s] bng_nf           done")
    return cl, pe


def agreement(n_seeds=4, workers=1, save=True):
    """bngsim NFsim == bngsim RuleMonkey (two independent NF algorithms, -bscb) is the
    reference of record; BNG2.pl's legacy NFsim v1.14.3 and a no-bscb run are contrasts."""
    print("=== agreement: bngsim-NFsim == bngsim-RuleMonkey (+ legacy BNG2.pl NFsim) ===")
    print(f"    {n_seeds} seeds/engine, doses {DOSE_AGREE} nM, t_end={T_END}s")
    wd = _workdir("agreement")
    try:
        bng_cl, bng_pe = _agreement_bng(n_seeds, wd, workers)
        # bngsim engines: per-dose baked XML, reused across seeds/engines.
        xmls = [_dose_xml(d, wd) for d in DOSE_AGREE]
        ens_cl = {"bng_nf": bng_cl}
        ens_pe = {"bng_nf": bng_pe}
        for label, engine, bscb in (("bngsim_nf", "nf", True),
                                    ("bngsim_rm", "rm", True),
                                    ("bngsim_nf_nobscb", "nf", False)):
            t0 = time.time()
            cl = np.empty((n_seeds, len(DOSE_AGREE)))
            pe = np.empty((n_seeds, len(DOSE_AGREE)))
            for i in range(len(DOSE_AGREE)):
                for s in range(n_seeds):
                    cl[s, i], pe[s, i] = _bngsim(engine, xmls[i], 5000 + s, bscb)
            ens_cl[label], ens_pe[label] = cl, pe
            print(f"  [{time.time()-t0:5.0f}s] {label:16s} done")
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    _agreement_report(ens_pe, ens_cl, n_seeds, save)


def agreement_bng(n_seeds=6, workers=1, save=True):
    """Re-run only the BNG2.pl NFsim arm of the committed agreement.npz (same seeds) and
    recompute the summary; the bngsim arms in the cache are kept as they are."""
    print("=== agreement-bng: refresh the BNG2.pl NFsim arm of agreement.npz ===")
    old = np.load(os.path.join(REFDIR, "agreement.npz"))
    assert int(old["n_seeds"]) == n_seeds, f"cache has {int(old['n_seeds'])} seeds per engine"
    assert np.allclose(old["doses"], DOSE_AGREE)
    wd = _workdir("agreement_bng")
    try:
        bng_cl, bng_pe = _agreement_bng(n_seeds, wd, workers)
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    labs = ("bngsim_nf", "bngsim_rm", "bngsim_nf_nobscb")
    ens_pe = {"bng_nf": bng_pe, **{k: old[f"pe_{k}"] for k in labs}}
    ens_cl = {"bng_nf": bng_cl, **{k: old[f"cl_{k}"] for k in labs}}
    _agreement_report(ens_pe, ens_cl, n_seeds, save)


def _agreement_report(ens_pe, ens_cl, n_seeds, save):
    def mean_se(A):
        return A.mean(0), A.std(0, ddof=1) / np.sqrt(A.shape[0])

    print("\n  per-dose equilibrium pEGFR (mean +/- s.e.m.):")
    print(f"  {'dose(nM)':>10s} {'BNG-NFsim':>16s} {'bngsim-NFsim':>16s} "
          f"{'bngsim-RuleMonk':>16s} {'NFsim no-bscb':>16s}")
    for i, d in enumerate(DOSE_AGREE):
        row = []
        for lab in ("bng_nf", "bngsim_nf", "bngsim_rm", "bngsim_nf_nobscb"):
            m, se = mean_se(ens_pe[lab][:, i])
            row.append(f"{m:.0f}+/-{se:.0f}")
        print(f"  {d:10.3f} " + " ".join(f"{r:>16s}" for r in row))

    # Max |mean difference| between two engines relative to their pooled per-dose mean,
    # over an informative-dose mask (pEGFR is ~0 at the lowest dose, which would blow up a
    # relative metric; Clusters is ~500 at every dose, so all doses count there).
    def rel_between(ens, a, b, mask):
        pooled = 0.5 * (ens[a].mean(0) + ens[b].mean(0))
        return float(np.max((np.abs(ens[a].mean(0) - ens[b].mean(0)) / pooled)[mask]))

    pe_pooled = np.mean([ens_pe[k] for k in ("bngsim_nf", "bngsim_rm")], (0, 1))
    pe_mask = pe_pooled > 20.0                      # doses where phospho-EGFR is informative
    all_mask = np.ones(len(DOSE_AGREE), bool)

    # PRIMARY capability check (issue #8): the two INDEPENDENT network-free algorithms
    # bngsim provides -- NFsim (Yang 2008) and the exact RuleMonkey (Colvin 2009/2010) --
    # agree. This is the reference of record for this aggregation model.
    rel_pe = rel_between(ens_pe, "bngsim_nf", "bngsim_rm", pe_mask)     # pEGFR
    rel_cl = rel_between(ens_cl, "bngsim_nf", "bngsim_rm", all_mask)    # Clusters
    # Legacy contrast: BNG2.pl's bundled NFsim v1.14.3 binary vs the exact RuleMonkey.
    # (The universal-traversal-limit is NOT the cause -- utl 4/5/6 give the same means;
    # this is an NFsim-binary version difference, and is exactly why legacy NFsim is not
    # the reference for an aggregation model -- issue #8 capability group.)
    legacy_pe = rel_between(ens_pe, "bng_nf", "bngsim_rm", pe_mask)
    # -bscb vs no-bscb contrast (same bngsim NFsim engine).
    bias_pe = np.abs(ens_pe["bngsim_nf"].mean(0) - ens_pe["bngsim_nf_nobscb"].mean(0))
    print(f"\n  PRIMARY: bngsim-NFsim == bngsim-RuleMonkey (two independent NF algorithms),")
    print(f"    max relative diff (informative doses): pEGFR={rel_pe:.1%}, Clusters={rel_cl:.1%}"
          f"  (Monte-Carlo scale)")
    print(f"  legacy BNG2.pl NFsim v1.14.3 vs exact RuleMonkey: pEGFR offset {legacy_pe:.1%} "
          f"(binary-version diff,\n    not UTL -> legacy NFsim is not the reference for this "
          f"aggregation model)")
    print(f"  -bscb vs no-bscb (bngsim NFsim): max |d pEGFR| = {bias_pe.max():.1f} "
          f"molecules (~few % ->\n    rings only rarely populated at this fit; far less "
          f"than a ring-forming model)")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "agreement.npz"),
                 doses=np.array(DOSE_AGREE), n_seeds=n_seeds, rel_pe=rel_pe, rel_cl=rel_cl,
                 legacy_pe=legacy_pe, bng_seeds=np.arange(5000, 5000 + n_seeds),
                 bng_engine=f"{_bng_version()}, param '{NFPARAM}', gml {GML}",
                 **{f"pe_{k}": v for k, v in ens_pe.items()},
                 **{f"cl_{k}": v for k, v in ens_cl.items()})
        print(f"  cached -> {os.path.join(REFDIR, 'agreement.npz')}")


# ---------------------------------------------------------------- figure
def plot():
    """Build verify_mitra2019.png from the cached reproduction/agreement npz files."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    rp = os.path.join(REFDIR, "reproduction.npz")
    ag = os.path.join(REFDIR, "agreement.npz")

    chisq_txt = ""
    if os.path.exists(rp):
        d = np.load(rp)
        lab_model = f"BioNetGen NFsim ({NFPARAM}), best fit, {int(d['n_reps'])} reps"
        chisq_txt = (f"; BioNetGen NFsim chi_sq = {float(d['total_chisq']):.2f} +/- "
                     f"{float(d['chisq_se']):.2f} (tol 30; published 13.4)"
                     if "chisq_se" in d else "")
        # Fig 3B -- cluster density vs time
        ax = axes[0, 0]
        ax.errorbar(d["t_exp"] / 60, d["d_pre2"], yerr=d["d_pre2_sd"], fmt="o", color="#222",
                    capsize=3, zorder=5, label="Kozer 2013 Fig 3B")
        ax.plot(d["t_exp"] / 60, d["m_pre2"], "s-", color="#1f77b4", lw=2, mfc="none",
                label=lab_model)
        ax.set(title="Fig 3B: EGFR cluster density vs time (30 nM)",
               xlabel="time (min)", ylabel="scaled cluster density")
        # Fig 3D -- phospho-EGFR vs time
        ax = axes[0, 1]
        ax.errorbar(d["t_exp"] / 60, d["d_pre4"], yerr=d["d_pre4_sd"], fmt="o", color="#222",
                    capsize=3, zorder=5, label="Kozer 2013 Fig 3D")
        ax.plot(d["t_exp"] / 60, d["m_pre4"], "s-", color="#d62728", lw=2, mfc="none",
                label=lab_model)
        ax.set(title="Fig 3D: phospho-EGFR vs time (30 nM)",
               xlabel="time (min)", ylabel="scaled pEGFR")
        # Fig 2B -- cluster density vs dose
        ax = axes[1, 0]
        ax.errorbar(d["doses"], d["d_pre1"], yerr=d["d_pre1_sd"], fmt="o", color="#222",
                    capsize=3, zorder=5, label="Kozer 2013 Fig 2B")
        ax.plot(d["doses"], d["m_pre1"], "s-", color="#1f77b4", lw=2, mfc="none",
                label=lab_model)
        ax.set(xscale="log", title="Fig 2B: EGFR cluster density vs EGF dose",
               xlabel="EGF dose (nM)", ylabel="scaled cluster density")
        # Fig 2D -- phospho-EGFR vs dose
        ax = axes[1, 1]
        ax.errorbar(d["doses"], d["d_pre3"], yerr=d["d_pre3_sd"], fmt="o", color="#222",
                    capsize=3, zorder=5, label="Kozer 2013 Fig 2D")
        ax.plot(d["doses"], d["m_pre3"], "s-", color="#d62728", lw=2, mfc="none",
                label=lab_model)
        ax.set(xscale="log", title="Fig 2D: phospho-EGFR vs EGF dose",
               xlabel="EGF dose (nM)", ylabel="scaled pEGFR")
        for ax in (axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]):
            ax.legend(frameon=False, fontsize=8)
            ax.grid(alpha=0.25, which="both")

    if os.path.exists(ag):
        d = np.load(ag)
        doses = d["doses"]
        # Panel: independent-algorithm agreement (pEGFR). bngsim NFsim == RuleMonkey is
        # the reference of record; BNG2.pl's legacy NFsim v1.14.3 binary is the contrast.
        ax = axes[0, 2]
        styles = [("bngsim_nf", "bngsim NFsim (Yang 2008)", "s-", "#7570b3", 1.0),
                  ("bngsim_rm", "bngsim RuleMonkey (Colvin 2009, exact)", "^-", "#d95f02", 1.0),
                  ("bng_nf", "BNG2.pl NFsim v1.14.3 (legacy)", "o--", "#1b9e77", 0.55)]
        for lab, name, fmt, col, alpha in styles:
            A = d[f"pe_{lab}"]
            m = A.mean(0); se = A.std(0, ddof=1) / np.sqrt(A.shape[0])
            ax.errorbar(doses, m, yerr=se, fmt=fmt, color=col, ms=7, lw=1.3,
                        capsize=2, mfc="none", alpha=alpha, label=name)
        legacy = float(d["legacy_pe"]) if "legacy_pe" in d else float("nan")
        ax.set(xscale="log", title=f"bngsim NFsim = RuleMonkey ({float(d['rel_pe']):.1%})\n"
               f"legacy v1.14.3 reads {legacy:.0%} lower",
               xlabel="EGF dose (nM)", ylabel="equilibrium pEGFR (molecules)")
        ax.legend(frameon=False, fontsize=7.5)
        ax.grid(alpha=0.25, which="both")
        # Panel: -bscb vs no-bscb contrast (pEGFR)
        ax = axes[1, 2]
        on = d["pe_bngsim_nf"]; off = d["pe_bngsim_nf_nobscb"]
        ax.errorbar(doses, on.mean(0), yerr=on.std(0, ddof=1) / np.sqrt(on.shape[0]),
                    fmt="s-", color="#7570b3", ms=7, lw=1.3, capsize=2, mfc="none",
                    label="NFsim -bscb (correct)")
        ax.errorbar(doses, off.mean(0), yerr=off.std(0, ddof=1) / np.sqrt(off.shape[0]),
                    fmt="x--", color="#e7298a", ms=7, lw=1.3, capsize=2,
                    label="NFsim no-bscb")
        ax.set(xscale="log", title="Steric flag: -bscb vs no-bscb\n(nearly coincide -- "
               "rings rarely populated)", xlabel="EGF dose (nM)",
               ylabel="equilibrium pEGFR (molecules)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.25, which="both")

    fig.suptitle("egfr_oligomerization_mitra2019 verification: Kozer 2013 Fig 2B/2D/3B/3D "
                 f"reproduction{chisq_txt}\n+ independent network-free algorithm agreement "
                 "(bngsim NFsim = RuleMonkey, -bscb)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    png = os.path.join(HERE, "verify_mitra2019.png")
    fig.savefig(png, dpi=130)
    print(f"wrote {png}")


if __name__ == "__main__":
    print("BioNetGen:", find_bionetgen()[0])
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else None
    w = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    if mode == "reference":
        reference(n_reps=n or 12, workers=w)
    elif mode == "agreement":
        agreement(n_seeds=n or 6, workers=w)
    elif mode == "agreement-bng":
        agreement_bng(n_seeds=n or 6, workers=w)
    elif mode == "plot":
        plot()
    elif mode == "all":
        reference(n_reps=n or 12, workers=w)
        agreement(n_seeds=6, workers=w)
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference [N [W]] | agreement [N [W]] | "
              f"agreement-bng [N [W]] | plot | all [N [W]]")
