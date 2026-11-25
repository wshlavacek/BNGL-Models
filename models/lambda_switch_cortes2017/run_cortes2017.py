"""Driver + verification for the network-free Cortes et al. (2017) phage-lambda
lysis/lysogeny decision models (issue #20).

Cortes MG, Trinh JT, Zeng L, Balazsi G (2017), Biophys J 113:2110-2120. This drives
the QSSA-reduced "detailed" model (lambda_switch_cortes2017.bngl, SI Eqs. 38-55) --
and, optionally, the per-genome operator-state variant -- network-free with NFsim or
RuleMonkey via bngsim. It reproduces:

  * Fig. S3 / Fig. 2A: probability of lysogeny vs MOI for initially small / average /
    large cells (rescaled initial volumes v0 = 0.5 / 1.0 / 1.5), and
  * Fig. S2: example lysogenic and lytic single-cell trajectories.

Self-contained: reads the committed .bngl next to it, overrides the MOI and v0
parameters by text substitution (post-init NFsim rate constants are NOT re-derived
from changed parameters, so parameters must be baked into the XML), generates BNG-XML
with the `bionetgen` CLI, and simulates network-free with `bngsim`.

Decision classifier (SI Sec. 1.10), evaluated per stochastic run from the CI/CII/Q
concentration trajectories ([X] = Obs_X * Vscale / Obs_Vol, in nM):

  Criteria #1 (matches experiment; Fig. S3 solid curves):
    LYSOGENY  if [Q] never crosses Q_T  AND  [CI] crosses CI_T
    LYSIS     if [Q] crosses Q_T at any time
  Criteria #2 (used to derive the simple CII-Q model; Fig. S3 dashed curves):
    LYSOGENY  if [CII] crosses CII_T before [Q] crosses Q_T
    LYSIS     if [Q] crosses Q_T before [CII] crosses CII_T
  Runs where neither regulator crosses its threshold are indeterminate and are
  excluded from the statistics (denominator = lysogenic + lytic).

Requirements: `bionetgen` on PATH (BNG >= 2.9.3) and the `bngsim` package. Examples:
    python run_cortes2017.py sweep 400 rm          # %lysogeny(MOI) for 3 cell sizes
    python run_cortes2017.py traj  1               # example lyso/lytic trajectories
    python run_cortes2017.py parity 100            # NFsim vs RuleMonkey agreement
The `bng` mode uses no bngsim: BNG2.pl only (found from $BNG2_PL or $BNGPATH), with numpy --
    python run_cortes2017.py bng 1000 36           # per-genome SSA + NFsim, detailed SSA
(see the section above bng() for its three arms and what each run is checked for).
"""

import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_MODEL_ALIASES = {
    "detailed": "lambda_switch_cortes2017.bngl",
    "pergenome": "lambda_switch_cortes2017_pergenome.bngl",
}
_sel = os.environ.get("CORTES_MODEL", "detailed")
MODEL = os.path.join(HERE, _MODEL_ALIASES.get(_sel, _sel))

# Decision parameters (Table S2 / Fig. S4); Vscale must match the .bngl.
VSCALE = 1000.0
Q_T = 134.897  # nM   -- Q lysis threshold
CI_T = 61.0073  # nM   -- CI lysogeny threshold (criteria #1)
CII_T = 22.0  # nM   -- CII lysogeny threshold (criteria #2, Fig. S4)
T_END = 60.0  # min  -- decision window (cf. Fig. S2)
N_POINTS = 121  # 0.5-min sampling, fine enough to catch threshold crossings

# initial rescaled volumes for the three experimental cell-size classes (Fig. S3)
CELL_SIZES = {"small": 0.5, "average": 1.0, "large": 1.5}

TRAJ = ["Obs_CI", "Obs_Cro", "Obs_CII", "Obs_Q", "Obs_DNA", "Obs_Vol"]


def build_xml(moi, v0, tag=None):
    """Write an MOI/v0-substituted, writeXML-only copy of MODEL to a temp dir, run
    `bionetgen`, and return the XML path. Uses a system temp dir so the driver never
    litters the model directory."""
    tag = tag or f"moi{moi}_v{v0}"
    workdir = tempfile.mkdtemp(prefix=f"cortes_{tag}_")
    text = open(MODEL).read()
    text, n = re.subn(r"(?m)^(\s*MOI\s+)[0-9.]+", rf"\g<1>{moi}", text)
    assert n == 1, f"expected exactly one MOI parameter line, found {n}"
    text, n = re.subn(r"(?m)^(\s*v0\s+)[0-9.]+", rf"\g<1>{v0}", text)
    assert n == 1, f"expected exactly one v0 parameter line, found {n}"
    # force writeXML-only actions (network-free)
    text = text[: text.index("begin actions")] + "begin actions\n  writeXML()\nend actions\n"
    bngl = os.path.join(workdir, f"cortes_{tag}.bngl")
    open(bngl, "w").write(text)
    out = os.path.join(workdir, "out")
    r = subprocess.run(["bionetgen", "run", "-i", bngl, "-o", out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-3000:]
    return os.path.join(out, f"cortes_{tag}.xml")


def ensemble(xml, method, n_seeds, n_points=N_POINTS, seed0=40000):
    """Return (tvec, {obs: array[n_seeds, n_points]}) of raw per-seed trajectories.
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


def _conc(stacks, obs):
    """Concentration trajectory [X](t) = Obs_X * Vscale / Obs_Vol, in nM."""
    return stacks[obs] * VSCALE / stacks["Obs_Vol"]


def classify(stacks, tvec):
    """Per-seed lysis/lysogeny classification by criteria #1 and #2 (SI Sec. 1.10).
    Returns dict of boolean/def arrays over seeds."""
    ci = _conc(stacks, "Obs_CI")
    cii = _conc(stacks, "Obs_CII")
    q = _conc(stacks, "Obs_Q")
    q_cross = q >= Q_T
    ci_cross = ci >= CI_T
    cii_cross = cii >= CII_T
    q_ever = q_cross.any(1)
    ci_ever = ci_cross.any(1)

    # criteria #1
    lyso1 = (~q_ever) & ci_ever
    lysis1 = q_ever

    # criteria #2: first-passage time comparison (index of first crossing; inf if none)
    def first(idx_bool):
        out = np.full(idx_bool.shape[0], np.inf)
        for i in range(idx_bool.shape[0]):
            w = np.where(idx_bool[i])[0]
            if len(w):
                out[i] = tvec[w[0]]
        return out

    t_cii = first(cii_cross)
    t_q = first(q_cross)
    lyso2 = t_cii < t_q
    lysis2 = t_q < t_cii
    return {
        "lyso1": lyso1,
        "lysis1": lysis1,
        "lyso2": lyso2,
        "lysis2": lysis2,
        "q_ever": q_ever,
        "ci_ever": ci_ever,
    }


def lyso_prob(c, key="1"):
    """Probability of lysogeny excluding indeterminate runs."""
    lyso, lysis = c[f"lyso{key}"], c[f"lysis{key}"]
    denom = lyso.sum() + lysis.sum()
    return 100.0 * lyso.sum() / denom if denom else float("nan")


# --------------------------------------------------------------- %lysogeny vs MOI
def lyso_fraction(moi, v0, n_seeds, method="rm", seed0=40000):
    xml = build_xml(moi, v0)
    tvec, st = ensemble(xml, method, n_seeds, seed0=seed0)
    c = classify(st, tvec)
    return lyso_prob(c, "1"), lyso_prob(c, "2"), st


def sweep(n_seeds=400, method="rm", mois=(1, 2, 3, 4, 5)):
    print(f"=== Fig. S3: %lysogeny vs MOI, {n_seeds} {method} seeds/point ===")
    print("  (criteria #1 = experiment-matching solid curves; #2 = dashed curves)")
    hdr = "  cell(v0)   " + "".join(f"   MOI{m}(1/2) " for m in mois)
    print(hdr)
    result = {}
    for size, v0 in CELL_SIZES.items():
        row = []
        for moi in mois:
            t0 = time.time()
            p1, p2, _ = lyso_fraction(moi, v0, n_seeds, method)
            row.append((p1, p2))
            print(
                f"    {size:7s}({v0}) MOI={moi}: crit#1={p1:5.1f}%  crit#2={p2:5.1f}%"
                f"  [{time.time() - t0:.0f}s]"
            )
        result[size] = row
    return result


# ------------------------------------------------------------- Fig. S2 trajectories
def traj(moi=1, v0=1.0, n_seeds=60, method="rm"):
    """Print mean lysogenic and mean lytic trajectories (criteria #1 split)."""
    print(f"=== Fig. S2: example trajectories, MOI={moi}, v0={v0}, {n_seeds} {method} seeds ===")
    xml = build_xml(moi, v0)
    tvec, st = ensemble(xml, method, n_seeds)
    c = classify(st, tvec)
    for lab, mask in (("LYSOGENY", c["lyso1"]), ("LYSIS", c["lysis1"])):
        if mask.sum() == 0:
            print(f"  {lab}: none in this sample")
            continue
        print(f"  {lab} (n={mask.sum()}):  t(min) " + "".join(f"{t:5.0f}" for t in tvec[::20]))
        for o in ("Obs_CI", "Obs_Cro", "Obs_CII", "Obs_Q", "Obs_DNA"):
            m = _conc(st, o)[mask].mean(0) if o != "Obs_DNA" else st[o][mask].mean(0)
            unit = "cnt" if o == "Obs_DNA" else "nM"
            print(f"    {o[4:]:>4s}({unit}):" + "".join(f"{v:5.0f}" for v in m[::20]))
    return tvec, st, c


# ------------------------------------------------------------- NFsim/RuleMonkey parity
def parity(n_seeds=100, moi=3, v0=1.0):
    print(f"=== NFsim vs RuleMonkey parity: MOI={moi}, v0={v0}, {n_seeds} seeds ===")
    out = {}
    for method in ("nf", "rm"):
        t0 = time.time()
        p1, p2, _ = lyso_fraction(moi, v0, n_seeds, method)
        out[method] = (p1, p2)
        print(f"  {method}: crit#1={p1:.1f}%  crit#2={p2:.1f}%  [{time.time() - t0:.0f}s]")
    d1 = abs(out["nf"][0] - out["rm"][0])
    d2 = abs(out["nf"][1] - out["rm"][1])
    print(
        f"  |NF-RM| crit#1={d1:.1f}%  crit#2={d2:.1f}%  "
        f"(expect within stochastic noise ~ sqrt(p(1-p)/n))"
    )
    return out


# ------------------------------------------------ the same sweep with no bngsim (BNG2.pl only)
# Every mode above simulates through bngsim. `bng` reruns the Fig. S3 sweep with BioNetGen alone,
# three arms on the SAME conditions, seeds and classifier:
#
#   pg_ssa   the per-genome file's generated network (generate_network: 22 species, 82
#            reactions -- the file is finite although its notes call it network-free only),
#            simulated with BNG2.pl's exact SSA (run_network);
#   pg_nf    the per-genome file network-free, BNG2.pl simulate_nf, i.e. the stock NFsim 1.14.3
#            binary at BNG2.pl's default flags (-cb; no -bscb, -utl automatic);
#   det_ssa  the primary (detailed, mean-field promoter) file's generated network, BNG2.pl SSA.
#
# Each block is ONE BNG2.pl run: setParameter("MOI"), setParameter("v0"), then for each seed a
# simulate with resetConcentrations() in between, exactly the action list a reader would write.
# Every trajectory is checked to start at Obs_DNA = MOI and Obs_Vol = v0*Vscale, so a parameter
# that did not take effect cannot pass silently. Per seed, the cache keeps the criteria #1/#2
# calls, the CII and Q first-passage times, and Obs_* at t = 20/40/60 min.
BNG_ARMS = {"pg_ssa": ("pergenome", "ssa"), "pg_nf": ("pergenome", "nf"),
            "det_ssa": ("detailed", "ssa")}
BNG_SEEDS, BNG_SEED0, BNG_BLOCK = 1000, 40000, 100
BNG_MOIS = (1, 2, 3, 4, 5)
BNG_KEEP_T = (40, 80, 120)          # sample indices of t = 20, 40, 60 min (0.5-min grid)


def find_bionetgen():
    """BNG2.pl from $BNG2_PL or $BNGPATH, then BNG2.pl on PATH."""
    import shutil
    for cand in (os.environ.get("BNG2_PL"),
                 os.path.join(os.environ["BNGPATH"], "BNG2.pl") if os.environ.get("BNGPATH")
                 else None,
                 shutil.which("BNG2.pl")):
        if cand and os.path.isfile(cand):
            return os.path.abspath(cand)
    raise FileNotFoundError("BioNetGen not found: set BNG2_PL (or BNGPATH) to a BioNetGen "
                            "2.9.3 BNG2.pl.")


def _bng_block(args):
    """One BNG2.pl run: `arm` at (moi, v0) for `seeds`. Returns per-seed arrays."""
    import shutil
    arm, moi, v0, seeds = args
    model, method = BNG_ARMS[arm]
    src = open(os.path.join(HERE, _MODEL_ALIASES[model])).read()
    acts = [f'setParameter("MOI",{moi})', f'setParameter("v0",{v0})']
    if method == "ssa":
        acts.append("generate_network({overwrite=>1})")
    for i, s in enumerate(seeds):
        if i:
            acts.append("resetConcentrations()")
        sim = "simulate_nf" if method == "nf" else 'simulate'
        extra = "" if method == "nf" else 'method=>"ssa",'
        acts.append(f'{sim}({{{extra}prefix=>"s{s}",t_start=>0,t_end=>{T_END:g},'
                    f'n_steps=>{N_POINTS - 1},seed=>{s}}})')
    text = src[: src.index("begin actions")] + "begin actions\n  " + "\n  ".join(acts) + \
        "\nend actions\n"
    work = tempfile.mkdtemp(prefix=f"bng_{arm}_m{moi}_v{v0}_", dir=HERE)
    try:
        bngl = os.path.join(work, "m.bngl")
        open(bngl, "w").write(text)
        r = subprocess.run(["perl", find_bionetgen(), "--outdir", work, bngl],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
        stacks = {o: [] for o in TRAJ}
        for s in seeds:
            g = os.path.join(work, f"s{s}.gdat")
            with open(g) as fh:
                names = fh.readline().lstrip("#").split()
            Y = np.loadtxt(g, ndmin=2)
            assert Y.shape[0] == N_POINTS and abs(Y[-1, 0] - T_END) < 1e-9, (arm, s, Y.shape)
            assert Y[0, names.index("Obs_DNA")] == moi, (arm, s, "MOI not applied")
            assert Y[0, names.index("Obs_Vol")] == v0 * VSCALE, (arm, s, "v0 not applied")
            for o in TRAJ:
                stacks[o].append(Y[:, names.index(o)])
        tvec = Y[:, 0]
    finally:
        shutil.rmtree(work, ignore_errors=True)
    st = {o: np.vstack(v) for o, v in stacks.items()}
    c = classify(st, tvec)
    cii, q = _conc(st, "Obs_CII"), _conc(st, "Obs_Q")
    fpt = {k: np.array([tvec[np.argmax(row)] if row.any() else np.inf for row in x])
           for k, x in (("t_cii", cii >= CII_T), ("t_q", q >= Q_T))}
    keep = {f"{o}": st[o][:, list(BNG_KEEP_T)].astype(np.int32) for o in TRAJ}
    return arm, moi, v0, np.array(seeds), c, fpt, keep


def bng(n_seeds=BNG_SEEDS, arms=tuple(BNG_ARMS), workers=1):
    """Run every (arm, cell size, MOI) at n_seeds through BNG2.pl and cache
    reference/cortes2017_pergenome_bng.npz."""
    import hashlib
    from multiprocessing import Pool
    blocks = [(a, m, v0, tuple(range(s0, min(s0 + BNG_BLOCK, BNG_SEED0 + n_seeds))))
              for a in arms for v0 in CELL_SIZES.values() for m in BNG_MOIS
              for s0 in range(BNG_SEED0, BNG_SEED0 + n_seeds, BNG_BLOCK)]
    print(f"=== bng: {len(blocks)} BNG2.pl runs of <= {BNG_BLOCK} seeds, arms {list(arms)}, "
          f"{n_seeds} seeds per condition from {BNG_SEED0}, {workers} workers", flush=True)
    acc = {}
    t0 = time.time()
    with Pool(workers) as pool:
        for i, (a, m, v0, seeds, c, fpt, keep) in enumerate(
                pool.imap_unordered(_bng_block, blocks, chunksize=1)):
            acc.setdefault((a, m, v0), []).append((seeds, c, fpt, keep))
            if (i + 1) % 25 == 0:
                print(f"  {i + 1}/{len(blocks)} blocks [{time.time() - t0:.0f}s]", flush=True)
    out = {}
    for (a, m, v0), parts in acc.items():
        parts.sort(key=lambda p: p[0][0])
        seeds = np.concatenate([p[0] for p in parts])
        assert np.array_equal(seeds, np.arange(BNG_SEED0, BNG_SEED0 + n_seeds))
        key = f"{a}__moi{m}__v{v0}"
        for k in ("lyso1", "lysis1", "lyso2", "lysis2"):
            out[f"{key}__{k}"] = np.concatenate([p[1][k] for p in parts])
        for k in ("t_cii", "t_q"):
            out[f"{key}__{k}"] = np.concatenate([p[2][k] for p in parts])
        for o in TRAJ:
            out[f"{key}__{o}"] = np.concatenate([p[3][o] for p in parts])
    print("  %lysogeny, criteria #1 (rows: cell size; columns: MOI 1..5)")
    for a in arms:
        for size, v0 in CELL_SIZES.items():
            vals = []
            for m in BNG_MOIS:
                k = f"{a}__moi{m}__v{v0}"
                n1, n0 = out[f"{k}__lyso1"].sum(), out[f"{k}__lysis1"].sum()
                vals.append(100.0 * n1 / (n1 + n0))
            print(f"    {a:8s} {size:8s}" + "".join(f"{v:7.1f}" for v in vals))
    sha = {m: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
           for m, f in _MODEL_ALIASES.items()}
    path = os.path.join(HERE, "reference", "cortes2017_pergenome_bng.npz")
    np.savez_compressed(path, arms=np.array(arms), mois=np.array(BNG_MOIS),
                        cell_sizes=np.array(list(CELL_SIZES)),
                        v0s=np.array(list(CELL_SIZES.values())), n_seeds=n_seeds,
                        seed0=BNG_SEED0, t_end=T_END, n_points=N_POINTS,
                        keep_t=np.array(BNG_KEEP_T) * T_END / (N_POINTS - 1),
                        bng="BioNetGen 2.9.3 BNG2.pl (run_network SSA; bin/NFsim 1.14.3, -cb)",
                        sha_pergenome=sha["pergenome"], sha_detailed=sha["detailed"], **out)
    print(f"  cached -> reference/{os.path.basename(path)}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sweep"
    if mode == "sweep":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 400
        sweep(n, sys.argv[3] if len(sys.argv) > 3 else "rm")
    elif mode == "traj":
        moi = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        traj(moi, float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
    elif mode == "bng":
        # bng [n_seeds] [workers] [arm,arm,...]
        bng(int(sys.argv[2]) if len(sys.argv) > 2 else BNG_SEEDS,
            tuple(sys.argv[4].split(",")) if len(sys.argv) > 4 else tuple(BNG_ARMS),
            int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    elif mode == "parity":
        parity(int(sys.argv[2]) if len(sys.argv) > 2 else 100)
    else:
        print(f"unknown mode {mode!r}; use sweep | traj | parity | bng")
