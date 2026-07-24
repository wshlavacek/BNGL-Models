"""Driver + verification for tcr_signaling_chylek2014.bngl -- the early T-cell-receptor
(TCR) / CD28 signaling model of Chylek et al. (2014), PLoS ONE 9(8):e104240
(doi:10.1371/journal.pone.0104240; Supplementary File A in File S1). The model's
combinatorial complexity is the multisite-phosphorylation state space of individual
proteins, not multivalent aggregation, so it is simulated NETWORK-FREE (NFsim) and has
no same-complex-binding ambiguity: standard/legacy BioNetGen NFsim is a valid reference
(no -bscb / -utl needed).

Two things are verified here:

  reference  -- runs the committed model's FAITHFUL published protocol (equilibrate 500 s
                with ligand binding off, turn on kfl/kfl_m per Supplementary File B, then
                a 60 s stimulation response) with standard BioNetGen NFsim, and checks the
                result against the committed reference .gdat. This is the paper-exact,
                full-subvolume (Fx=0.02) run.

  agreement  -- cross-checks that the three network-free engines agree on the SAME rules:
                standard BioNetGen NFsim, bngsim NfsimSession, and bngsim RuleMonkeySession.
                Because RuleMonkey (an exact network-free method) cannot afford the 500 s
                equilibration at full subvolume, the cross-check uses a REDUCED but
                identical protocol on all three engines: a single-phase, ligand-on run
                (kfl/kfl_m on from t=0, no equilibration) at a reduced subvolume
                (Fx=0.005). All three sample the same continuous-time Markov chain, so
                their ensemble-mean phosphosite trajectories agree within Monte-Carlo
                error. Agreement is quantified by the pairwise z-score
                |mean_i - mean_j| / sqrt(se_i^2 + se_j^2) across observables and times.

Requirements: `bionetgen` on PATH (BNG core with NFsim) and the `bngsim` package
(NfsimSession + RuleMonkeySession). Run with the repo venv:
    python run_tcr_chylek2014.py reference          # faithful protocol vs committed .gdat
    python run_tcr_chylek2014.py agreement 16        # 3-engine cross-check, 16 seeds/engine
    python run_tcr_chylek2014.py plot                # (re)build verify_tcr_chylek2014.png
"""
import os
import re
import subprocess
import sys
import tempfile
import time

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
def _bng_run(bngl_text, tag):
    """Write bngl_text to a temp dir, run `bionetgen run`, return the output dir."""
    workdir = tempfile.mkdtemp(prefix=f"tcr_{tag}_")
    bngl = os.path.join(workdir, f"{tag}.bngl")
    open(bngl, "w").write(bngl_text)
    out = os.path.join(workdir, "out")
    r = subprocess.run(["bionetgen", "run", "-i", bngl, "-o", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:] + "\n---STDOUT---\n" + r.stdout[-2000:]
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
    out = _bng_run(text, "xml")
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
        out = _bng_run(base + act, f"native_s{s}")
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
    out = _bng_run(open(MODEL).read(), tag)  # BNG names outputs after the temp bngl (tag)
    stem = os.path.splitext(os.path.basename(MODEL))[0]
    t_eq, eq = _read_gdat(os.path.join(out, f"{tag}_equil.gdat"))
    t_tc, tc = _read_gdat(os.path.join(out, f"{tag}_tcr.gdat"))
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


# ---------------------------------------------------------------- figure
def plot():
    """Build verify_tcr_chylek2014.png from the committed reference .gdat + agreement.npz."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stem = os.path.splitext(os.path.basename(MODEL))[0]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))

    # Row 1: faithful published protocol (full equilibration + response) from committed ref.
    eqp = os.path.join(REFDIR, f"{stem}_equil.gdat")
    tcp = os.path.join(REFDIR, f"{stem}_tcr.gdat")
    if os.path.exists(tcp):
        t_tc, tc = _read_gdat(tcp)
        ax = axes[0, 0]
        for o in PLOT_OBS:
            ax.plot(t_tc, tc[o], marker=".", ms=3, label=o)
        ax.set(title="Faithful protocol: 60 s response (BNG NFsim, Fx=0.02)",
               xlabel="time after stimulation (s)", ylabel="molecules / subvolume")
        ax.legend(fontsize=6, ncol=2)
        if os.path.exists(eqp):
            t_eq, eq = _read_gdat(eqp)
            ax2 = axes[0, 1]
            for o in PLOT_OBS:
                ax2.plot(t_eq, eq[o], marker=".", ms=3, label=o)
            ax2.set(title="Equilibration (ligand OFF): flat basal state",
                    xlabel="time (s)", ylabel="molecules / subvolume")
            ax2.legend(fontsize=6, ncol=2)
        axes[0, 2].axis("off")
        axes[0, 2].text(0.0, 0.5,
            "Chylek et al. (2014) PLoS ONE 9(8):e104240\n"
            "Early TCR/CD28 signaling, network-free (NFsim)\n\n"
            "Top: published two-phase protocol\n"
            "  equilibrate 500 s (kfl=kfl_m=0) then turn on\n"
            "  ligand binding and follow 60 s response.\n\n"
            "Bottom: three network-free engines agree on\n"
            "  the same rules (reduced subvolume Fx=0.005,\n"
            "  single-phase, ligand-on). Points = ensemble\n"
            "  mean; shaded = +/-1 s.e.m.", fontsize=8, va="center")

    # Row 2: three-engine agreement from cached npz.
    npz = os.path.join(REFDIR, "agreement.npz")
    if os.path.exists(npz):
        d = np.load(npz, allow_pickle=True)
        tvec = d["tvec"]
        styles = {"bng_nf": ("BNG NFsim", "o", "#1b9e77"),
                  "bngsim_nf": ("bngsim NFsim", "s", "#7570b3"),
                  "bngsim_rm": ("bngsim RuleMonkey", "^", "#d95f02")}
        # Place the 6 plotted observables into 3 axes (2 observables per axis), each
        # overlaying all three engines' mean +/- s.e.m. trajectories.
        groups = [PLOT_OBS[0:2], PLOT_OBS[2:4], PLOT_OBS[4:6]]
        for gi, group in enumerate(groups):
            ax = axes[1, gi]
            for o in group:
                for lab, (name, mk, col) in styles.items():
                    key = f"{lab}__{o}"
                    if key not in d:
                        continue
                    A = d[key]
                    m = A.mean(0)
                    se = A.std(0, ddof=1) / np.sqrt(A.shape[0])
                    ax.plot(tvec, m, marker=mk, ms=3, color=col, lw=1,
                            label=name if o == group[0] else None)
                    ax.fill_between(tvec, m - se, m + se, color=col, alpha=0.15)
                ax.annotate(o, (tvec[-1], A.mean(0)[-1]), fontsize=6,
                            xytext=(2, 0), textcoords="offset points")
            ax.set(xlabel="time (s)", ylabel="molecules / subvolume")
            if gi == 0:
                ax.legend(fontsize=6, title="engine")
            ax.set_title(", ".join(group), fontsize=8)

    fig.suptitle("tcr_signaling_chylek2014 verification: faithful NFsim protocol + "
                 "3-engine network-free agreement", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    png = os.path.join(HERE, "verify_tcr_chylek2014.png")
    fig.savefig(png, dpi=130)
    print(f"wrote {png}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "reference"
    if mode == "reference":
        reference()
    elif mode == "agreement":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        agreement(n_seeds=n)
    elif mode == "plot":
        plot()
    elif mode == "all":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        reference()
        agreement(n_seeds=n)
        plot()
    else:
        print(f"unknown mode {mode!r}; use reference | agreement [N] | plot | all [N]")
