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


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sweep"
    if mode == "sweep":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 400
        sweep(n, sys.argv[3] if len(sys.argv) > 3 else "rm")
    elif mode == "traj":
        moi = int(sys.argv[2]) if len(sys.argv) > 2 else 1
        traj(moi, float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
    elif mode == "parity":
        parity(int(sys.argv[2]) if len(sys.argv) > 2 else 100)
    else:
        print(f"unknown mode {mode!r}; use sweep | traj | parity")
