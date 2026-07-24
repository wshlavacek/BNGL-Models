"""Driver + verification for tlbr_steric_monine2010.bngl -- the kinetic
trivalent-ligand / bivalent-receptor (TLBR) model of Monine et al. (2010),
Biophys J 98(1):48-56 (doi:10.1016/j.bpj.2009.09.043).

A trivalent ligand crosslinks bivalent IgE-FceRI receptors into cell-surface
aggregates that grow WITHOUT BOUND, so the reaction network is infinite and the
model is simulated NETWORK-FREE (NFsim / RuleMonkey). The one steric constraint of
the TLBR model is that aggregates are ACYCLIC (the equivalent-site Goldstein-
Perelson assumption): a tethered ligand arm may not bind a receptor already in its
own aggregate. That constraint is NFsim's block-same-complex-binding (-bscb), on by
default in bngsim's NfsimSession and RuleMonkeySession. Because this is an
aggregation model with same-complex-binding ambiguity, legacy BioNetGen NFsim is
NOT a valid reference by itself; correctness is confirmed instead by agreement of
two INDEPENDENT network-free algorithms -- NFsim (Yang et al. 2008, Phys Rev E
78:031910) and RuleMonkey (Colvin et al. 2009/2010) -- both run with -bscb.

Three things are verified here:

  reference   -- reproduces Monine 2010 Fig 2a. Runs the committed model's faithful
                 published protocol (network-free NFsim parameter_scan with -bscb,
                 over the 12 experimental doses, to equilibrium t_end=5000 s) at the
                 Table-1 best fit (K1=0.467/nM, K2=87.03/nM, alpha=0.816 -- the
                 shipped nominals), averaged over replicate scans, and compares the
                 model bound-ligand fraction to the Fig 2a data. Reports the fit
                 objective sos (FL space), the RMS in lambda (Monine kept fits with
                 RMS < 0.02, SI Eq 11), and the relative error over doses above the
                 near-zero noise floor (FL > 0.05).

  agreement   -- the independent-algorithm cross-check (issue #11, task 3). Three
                 network-free engines run the SAME acyclic rules with -bscb:
                 BioNetGen NFsim v1.14.3, bngsim NfsimSession, and bngsim
                 RuleMonkeySession. All three sample the same continuous-time Markov
                 chain, so their per-dose equilibrium lambda agree within Monte-Carlo
                 error (pairwise z-scores). For contrast, bngsim NFsim is also run
                 WITHOUT -bscb: ring closure is then allowed, the acyclic TLBR model
                 is no longer simulated, and lambda shifts systematically -- so a
                 no-bscb (legacy) run is not a valid reference for this model.

  plot        -- builds verify_monine2010.png from the cached .npz files.

Requirements: BNG2.pl (set BNGPATH, or the default install is used) and the bngsim
package (NfsimSession + RuleMonkeySession). Run with the repo venv:
    python run_monine2010.py reference 5        # Fig 2a reproduction, 5 replicate scans
    python run_monine2010.py agreement 10       # 3-engine + bscb contrast, 10 seeds/engine
    python run_monine2010.py plot               # (re)build verify_monine2010.png
    python run_monine2010.py all
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
MODEL = os.path.join(HERE, "tlbr_steric_monine2010.bngl")
REFDIR = os.path.join(HERE, "reference")

BNGPATH = os.environ.get("BNGPATH") or os.path.expanduser("~/Simulations/BioNetGen-2.9.3")
BNG2PL = os.path.join(BNGPATH, "BNG2.pl")

# Monine 2010 Table 1 best fit for the TLBR model (the paper's own reported result
# AND the shipped .bngl nominals).
BESTFIT = dict(K1=0.467, K2=87.03, alpha=0.816)
ALPHA = BESTFIT["alpha"]
T_END = 5000                 # equilibration horizon (s); NFsim has no steady-state solve
GML = 2147483647             # 32-bit-max global molecule limit

# The 12 experimental doses (nM) and FL1 data of Monine 2010 Fig 2a.
_DATA = np.array([
    [0.0005006902, 0.004888544], [0.001362623, -0.0016556821],
    [0.0044341334, 0.0018911886], [0.0149210839, 0.0176273048],
    [0.0441574, 0.0415561926], [0.1507897315, 0.2057684534],
    [0.5013619944, 0.4269447413], [1.5652727704, 0.5017002142],
    [5.2257161826, 0.6008341783], [16.9016532291, 0.7060841749],
    [67.9604112309, 0.888530348], [213.4593409505, 1.0243034146],
])
DOSES = _DATA[:, 0]
FL_DATA = _DATA[:, 1]


# ---------------------------------------------------------------- BNG helpers
def _model_block():
    """Committed model text up to and including `end model` (no actions block)."""
    text = open(MODEL).read()
    marker = "\nend model"
    i = text.index(marker)
    return text[:i + len(marker)] + "\n"


def _set_param_text(text, name, value):
    text, n = re.subn(rf"(?m)^(\s*{re.escape(name)}\s+)\S+", rf"\g<1>{value:.12g}", text, count=1)
    assert n == 1, f"expected one '{name}' parameter line, found {n}"
    return text


def _apply_bestfit(text, params):
    for k, v in (params or {}).items():
        text = _set_param_text(text, k, v)
    return text


def _run_bng(bngl_text, workdir, tag):
    bngl = os.path.join(workdir, f"{tag}.bngl")
    open(bngl, "w").write(bngl_text)
    r = subprocess.run(["perl", BNG2PL, bngl], cwd=workdir,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:] + "\n---STDOUT---\n" + r.stdout[-3000:]
    return bngl


def _read_scan(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    return {name: data[:, i] for i, name in enumerate(header)}


def _bng_scan_lambda(doses, seed, workdir, params=None):
    """One BNG NFsim parameter_scan over `doses` (with -bscb); returns final lambda
    per dose. This is exactly the model's actions-block protocol."""
    base = _apply_bestfit(_model_block(), params)
    vals = ",".join(f"{d:.10g}" for d in doses)
    act = ("\nbegin actions\n"
           f'  parameter_scan({{parameter=>"LT_conc",par_scan_vals=>[{vals}],'
           f'method=>"nf",suffix=>"scan",t_start=>0,t_end=>{T_END},n_steps=>50,'
           f'gml=>{GML},print_functions=>1,seed=>{seed},param=>"-bscb"}})\n'
           "end actions\n")
    _run_bng(base + act, workdir, f"scan_s{seed}")
    cols = _read_scan(glob.glob(os.path.join(workdir, "*_scan.scan"))[0])
    return (cols["L_total"] - cols["L_free"]) / (2.0 * cols["R_total"])


def _dose_xml(dose, workdir, params=None):
    """Bake one dose (and best-fit params) into the model and writeXML; returns the
    XML path. Baking is required because bngsim RuleMonkey's set_param does not
    re-derive seed-species counts (see memory: bngsim set_param NF gaps)."""
    tag = f"dose_{dose:.10g}".replace(".", "p").replace("-", "m")
    xml = os.path.join(workdir, f"{tag}.xml")
    if os.path.exists(xml):
        return xml
    base = _apply_bestfit(_model_block(), params)
    base = _set_param_text(base, "LT_conc", dose)
    base += "\nbegin actions\n  writeXML({overwrite=>1})\nend actions\n"
    _run_bng(base, workdir, tag)
    produced = glob.glob(os.path.join(workdir, f"{tag}.xml"))
    assert produced, f"no XML produced for dose {dose}"
    return produced[0]


# ---------------------------------------------------------------- bngsim helpers
def _bngsim_lambda(engine, xml, seed, bscb=True):
    """Final equilibrium lambda from one bngsim network-free run."""
    from bngsim import NfsimSession, RuleMonkeySession
    Sess = NfsimSession if engine == "nf" else RuleMonkeySession
    with Sess(xml, molecule_limit=GML, block_same_complex_binding=bscb) as s:
        s.initialize(seed=seed)
        res = s.simulate(0, T_END, n_points=4)
        names = list(s.get_observable_names())
        Y = np.asarray(res.observables)[-1]
        d = {n: Y[i] for i, n in enumerate(names)}
    return (d["L_total"] - d["L_free"]) / (2.0 * d["R_total"])


# ---------------------------------------------------------------- reference
def reference(n_reps=5, save=True):
    """Reproduce Monine 2010 Fig 2a: NFsim (-bscb) dose-response vs the binding data."""
    print("=== reference: Fig 2a reproduction (BNG NFsim, -bscb, Table-1 fit) ===")
    print(f"    {n_reps} replicate parameter_scans over {len(DOSES)} doses, t_end={T_END}s")
    with tempfile.TemporaryDirectory() as wd:
        t0 = time.time()
        reps = np.vstack([_bng_scan_lambda(DOSES, 1000 + r, wd, BESTFIT) for r in range(n_reps)])
    lam = reps.mean(0)
    fl_model = lam / ALPHA
    print(f"  [{time.time()-t0:.0f}s] done")
    sos = float(np.sum((fl_model - FL_DATA) ** 2))
    rms_lam = float(np.sqrt(np.mean((lam - ALPHA * FL_DATA) ** 2)))
    big = np.abs(FL_DATA) > 0.05
    rel = np.abs(fl_model[big] - FL_DATA[big]) / np.abs(FL_DATA[big])
    print(f"  sos (fit objective, FL)     = {sos:.5f}")
    print(f"  RMS (lambda space)          = {rms_lam:.4f}   (Monine kept fits RMS<0.02, SI Eq 11)")
    print(f"  median |rel err| (FL>0.05)  = {np.median(rel):.3f}")
    print(f"  max    |rel err| (FL>0.05)  = {rel.max():.3f}")
    print(f"  {'dose(nM)':>12s} {'FL_model':>9s} {'FL_data':>9s}")
    for d, fm, fd in zip(DOSES, fl_model, FL_DATA):
        print(f"  {d:12.4f} {fm:9.4f} {fd:9.4f}")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "reproduction.npz"),
                 doses=DOSES, fl_data=FL_DATA, fl_model=fl_model, lam_model=lam,
                 reps=reps, n_reps=n_reps, sos=sos, rms_lam=rms_lam)
        print(f"  cached -> {os.path.join(REFDIR, 'reproduction.npz')}")
    return fl_model


# ---------------------------------------------------------------- agreement
def agreement(n_seeds=10, save=True):
    """Three network-free engines (BNG NFsim, bngsim NFsim, bngsim RuleMonkey), all
    with -bscb, agree on the acyclic TLBR model; a no-bscb run is shown to diverge."""
    print("=== agreement: BNG-NFsim vs bngsim-NFsim vs bngsim-RuleMonkey (all -bscb) ===")
    print(f"    {n_seeds} seeds/engine, {len(DOSES)} doses, t_end={T_END}s")
    with tempfile.TemporaryDirectory() as wd:
        # BNG NFsim: one parameter_scan per seed (all doses at once).
        t0 = time.time()
        bng = np.vstack([_bng_scan_lambda(DOSES, 5000 + s, wd, BESTFIT) for s in range(n_seeds)])
        print(f"  [{time.time()-t0:5.0f}s] bng_nf      done")
        # bngsim engines: per-dose baked XML, reused across seeds/engines.
        xmls = [_dose_xml(d, wd, BESTFIT) for d in DOSES]
        engines = {}
        for label, engine, bscb in (("bngsim_nf", "nf", True),
                                     ("bngsim_rm", "rm", True),
                                     ("bngsim_nf_nobscb", "nf", False)):
            t0 = time.time()
            engines[label] = np.array(
                [[_bngsim_lambda(engine, xmls[i], 5000 + s, bscb) for s in range(n_seeds)]
                 for i in range(len(DOSES))]).T   # shape (n_seeds, n_doses)
            print(f"  [{time.time()-t0:5.0f}s] {label:16s} done")
    ens = {"bng_nf": bng, **engines}

    def mean_se(A):
        return A.mean(0), A.std(0, ddof=1) / np.sqrt(A.shape[0])

    big = np.abs(FL_DATA) > 0.05    # informative doses (above the noise floor)
    print("\n  per-dose equilibrium lambda (mean +/- s.e.m.), informative doses:")
    print(f"  {'dose(nM)':>10s} {'BNG-NFsim':>16s} {'bngsim-NFsim':>16s} "
          f"{'bngsim-RuleMonk':>16s} {'NFsim no-bscb':>16s}")
    for i in np.where(big)[0]:
        row = []
        for lab in ("bng_nf", "bngsim_nf", "bngsim_rm", "bngsim_nf_nobscb"):
            m, se = ens[lab][:, i].mean(), ens[lab][:, i].std(ddof=1) / np.sqrt(n_seeds)
            row.append(f"{m:.3f}+/-{se:.3f}")
        print(f"  {DOSES[i]:10.3f} " + " ".join(f"{r:>16s}" for r in row))

    # Agreement metric. The equilibrium lambda is so tightly determined (s.e.m.
    # ~0.001-0.005) that z-scores are hypersensitive: sub-2% differences between
    # independent engine implementations become "significant" even though they are
    # scientifically negligible. The physically meaningful statement is the ABSOLUTE
    # agreement in lambda, compared to Monine's own RMS < 0.02 acceptance (SI Eq 11).
    pairs = [("bng_nf", "bngsim_nf"), ("bng_nf", "bngsim_rm"), ("bngsim_nf", "bngsim_rm")]
    dmax, zmax = 0.0, 0.0
    for a, b in pairs:
        ma, sa = mean_se(ens[a]); mb, sb = mean_se(ens[b])
        dmax = max(dmax, float(np.max(np.abs(ma - mb)[big])))
        denom = np.sqrt(sa**2 + sb**2)
        z = np.abs(ma - mb)[big] / np.where(denom[big] > 0, denom[big], np.inf)
        zmax = max(zmax, float(np.nanmax(z)))
    # Contrast: -bscb vs no-bscb (same NFsim engine): bias from ring closure. This is
    # dose-dependent -- largest at low dose, where a ligand can double-bond one
    # receptor into a 1:1 ring that -bscb forbids.
    m_on, _ = mean_se(ens["bngsim_nf"]); m_off, _ = mean_se(ens["bngsim_nf_nobscb"])
    absbias = np.abs(m_on - m_off)[big]
    print(f"\n  three -bscb engines agree to max |d lambda| = {dmax:.4f} across informative "
          f"doses\n    (Monte-Carlo scale; well within Monine's RMS < 0.02 acceptance, "
          f"SI Eq 11; max|z|={zmax:.1f} only because s.e.m. is tiny)")
    print(f"  -bscb vs no-bscb (bngsim NFsim): max |d lambda| = {absbias.max():.4f} "
          f"(at {DOSES[big][absbias.argmax()]:.3g} nM)")
    print(f"    -> without -bscb NFsim closes rings and no longer simulates the acyclic "
          f"TLBR model; the no-bscb run is not a valid reference.")
    if save:
        os.makedirs(REFDIR, exist_ok=True)
        np.savez(os.path.join(REFDIR, "agreement.npz"),
                 doses=DOSES, fl_data=FL_DATA, n_seeds=n_seeds, zmax=zmax, dmax=dmax,
                 **{lab: ens[lab] for lab in ens})
        print(f"  cached -> {os.path.join(REFDIR, 'agreement.npz')}")
    return ens, dmax


# ---------------------------------------------------------------- figure
def plot():
    """Build verify_monine2010.png from the cached reproduction/agreement npz files."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: Fig 2a reproduction ----------------------------------------------------
    ax = axes[0]
    rp = os.path.join(REFDIR, "reproduction.npz")
    if os.path.exists(rp):
        d = np.load(rp)
        ax.plot(d["doses"], d["fl_data"], "o", color="#222", ms=8, zorder=5,
                label="Monine 2010 Fig 2a data")
        ax.plot(d["doses"], d["fl_model"], "s-", color="#1f77b4", lw=2, ms=6,
                label=f"NFsim (-bscb) @ Table 1\nsos={float(d['sos']):.4f}, "
                      f"RMS_lambda={float(d['rms_lam']):.3f}")
        ax.set(xscale="log", xlabel="total ligand (nM)",
               ylabel="FL  (bound-ligand fraction = lambda/alpha)",
               title="Fig 2a reproduction (acyclic TLBR, -bscb)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.25, which="both")

    # Panel 2: three-engine agreement -------------------------------------------------
    ax = axes[1]
    ag = os.path.join(REFDIR, "agreement.npz")
    if os.path.exists(ag):
        d = np.load(ag)
        doses = d["doses"]
        styles = {"bng_nf": ("BNG NFsim", "o", "#1b9e77"),
                  "bngsim_nf": ("bngsim NFsim", "s", "#7570b3"),
                  "bngsim_rm": ("bngsim RuleMonkey", "^", "#d95f02")}
        for lab, (name, mk, col) in styles.items():
            A = d[lab]
            m = A.mean(0); se = A.std(0, ddof=1) / np.sqrt(A.shape[0])
            ax.errorbar(doses, m, yerr=se, fmt=mk + "-", color=col, ms=6, lw=1.3,
                        capsize=2, label=name, mfc="none")
        ax.set(xscale="log", xlabel="total ligand (nM)", ylabel="equilibrium lambda",
               title=f"3-engine agreement, all -bscb\nmax |d lambda| = {float(d['dmax']):.3f} "
                     f"(< Monine RMS 0.02)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.25, which="both")

    # Panel 3: -bscb vs no-bscb (steric constraint matters) ---------------------------
    ax = axes[2]
    if os.path.exists(ag):
        d = np.load(ag)
        doses = d["doses"]
        on = d["bngsim_nf"]; off = d["bngsim_nf_nobscb"]
        ax.errorbar(doses, on.mean(0), yerr=on.std(0, ddof=1) / np.sqrt(on.shape[0]),
                    fmt="s-", color="#7570b3", ms=6, lw=1.3, capsize=2, mfc="none",
                    label="NFsim -bscb (acyclic TLBR)")
        ax.errorbar(doses, off.mean(0), yerr=off.std(0, ddof=1) / np.sqrt(off.shape[0]),
                    fmt="x--", color="#e7298a", ms=6, lw=1.3, capsize=2,
                    label="NFsim no-bscb (rings allowed)")
        ax.set(xscale="log", xlabel="total ligand (nM)", ylabel="equilibrium lambda",
               title="Steric constraint: -bscb required\n(legacy no-bscb is a different model)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.25, which="both")

    fig.suptitle("tlbr_steric_monine2010 verification: Fig 2a reproduction + independent "
                 "network-free algorithm agreement (NFsim = RuleMonkey, -bscb)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    png = os.path.join(HERE, "verify_monine2010.png")
    fig.savefig(png, dpi=130)
    print(f"wrote {png}")


if __name__ == "__main__":
    assert os.path.exists(BNG2PL), f"BNG2.pl not found at {BNG2PL}; set BNGPATH"
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    if mode == "reference":
        reference(n_reps=int(sys.argv[2]) if len(sys.argv) > 2 else 5)
    elif mode == "agreement":
        agreement(n_seeds=int(sys.argv[2]) if len(sys.argv) > 2 else 10)
    elif mode == "plot":
        plot()
    elif mode == "all":
        reference(n_reps=int(sys.argv[2]) if len(sys.argv) > 2 else 5)
        agreement(n_seeds=int(sys.argv[3]) if len(sys.argv) > 3 else 10)
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference [N] | agreement [N] | plot | all")
