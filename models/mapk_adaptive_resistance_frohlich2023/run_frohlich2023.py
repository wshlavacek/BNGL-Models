#!/usr/bin/env python
"""Check the curated BNGL against the authors' own MARM2.0 model (Fröhlich et al., 2023).

The curated file is a condition-specific BioNetGen translation of MARM2.0, the PySB
energy-pattern model released with Fröhlich F, Gerosa L, Muhlich JL, Sorger PK (2023),
Mol Syst Biol 19:e10988. This driver rebuilds the Fig. 2A simulation from the authors' own
release and compares it with the pinned `.bngl`. The two arms share no model code:

  authors' arm   reference/frohlich2023_marm2_RTKERK__pRAF.py, the authors' flattened PySB
                 model (variant pRAF), at row 0 -- the best-ranked set -- of their calibrated
                 parameter table. PySB generates the network (through BNG2.pl) and compiles
                 the right-hand side from its own rate expressions (ScipyOdeSimulator's RHS
                 builder); the Jacobian is PySB's construction kept sparse; SciPy BDF
                 integrates; PySB's SimulationResult evaluates observables and expressions.
  curated arm    the pinned `.bngl`, expanded by BNG2.pl and integrated by the independent
                 SciPy implementation in independent_frohlich2023.py, with its hand-written
                 analytic Jacobian.

ScipyOdeSimulator.run itself is not used to integrate: VODE (and LSODA) form a dense
1,612 x 1,612 finite-difference Jacobian, and VODE had not finished the 0.1-1 h stretch of
the pre-equilibration after 25 minutes. Analytic Jacobians matter in both arms: once the
state stops moving, BDF with a finite-difference Jacobian is held to small steps.

Two adaptation protocols are run in the authors' arm, because they are not the same:

  authors'  what the authors' pipeline does for Fig. 2A. Rows 0-9 of their processed data
            table (MARM/data/processed_EGF_EGFR_MEKi_PRAFi_RAFi.csv; columns recorded in
            reference/frohlich2023_marm2_upstream.json) set `*_preeq` columns and
            `t_presim` = 0, which AMICI executes as a pre-equilibration to steady state at
            1 uM vemurafenib without EGF; the paper: "model simulations were pre-equilibrated
            to these conditions". (Their finite pre-treatment experiments instead start from
            the drug-free steady state and use t_presim.)
  24 h      what the curated `.bngl` seeds encode: the model's own initial values (free
            monomers at the calibrated abundances, no transcripts, no DUSP, SPRY or EGFR)
            integrated for 24 h at 1 uM vemurafenib without EGF.

Both then re-initialize free EGF at 100 ng/mL (AMICI reinitialization_state_idxs_sim) and
integrate 0-2 h with the authors' solver tolerances (MARM/estimation.py get_solver:
rtol 1e-9, atol 1e-11).

Two steps of the authors' pipeline import AMICI and are re-implemented here rather than run:
MARM.model.get_model_instance, which for the EGF_RAFi instance drops the PRAFi and MEKi
initials; and the parameter mapping of MARM.parameters / MARM.estimation, which reads
`bind_Vemurafenib_RAF_kf` for the model's `bind_RAFi_RAF_kf`. The flat file ships the unused
PRAFi and MEKi binding constants at 0; BNG2.pl 2.9.3 refuses to evaluate their log(0)
activation energies, so they are set to a placeholder of 1, and the driver asserts that no
reaction is generated from any rule that uses them.

Comparisons (tolerances and reasons in TOLERANCES below):

  model     protocol-independent: every calibrated value of row 0 against the `.bngl`; a
            species bijection (canonical labels) with identical reactions per rule; dy/dt
            of the two networks at six states, against each species' gross flux, and their
            two independently derived Jacobians at the same states
  per       the `.bngl` seeds against the authors' state at t = 0; pERK_IF_obs and
  protocol  pMEK_IF_obs, all 30 observables and all 1,612 species over 0-2 h

Verdicts: `model_equivalent` is the model comparisons plus the 24 h protocol; `all_pass`
also requires the authors' protocol. The exit status is 1 unless all_pass.

Usage
-----
    export BNGPATH=/path/to/BioNetGen-2.9.3          # BNG2.pl; run_network is not needed
    uv sync --group marm2
    uv run --group marm2 python models/mapk_adaptive_resistance_frohlich2023/run_frohlich2023.py

Runs offline: it reads only committed files. About 4 minutes on one core (3 min 53 s on a
Xeon E7-8890 v4), half of it BNG2.pl expanding the two networks and PySB compiling the
right-hand side. `--compiler python` avoids the C compiler PySB's Cython backend needs, at
roughly 30x the right-hand-side cost.

Writes into reference/ (overwriting):
    frohlich2023_marm2_fig2a_authors_protocol.csv   authors' model, steady-state adaptation
    frohlich2023_marm2_fig2a_24h_protocol.csv       authors' model, the .bngl's 24 h adaptation
    frohlich2023_marm2_adapted_state.csv            .bngl seeds beside both authors' states
    frohlich2023_marm2_comparison.json              every metric, tolerance and verdict
The outputs are written before the verdict, so a failure can be read.
"""

from __future__ import annotations

import os

# One BLAS thread: nothing here gains from more, and the driver must not saturate a shared
# login node.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
from collections import Counter  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.sparse import csc_matrix  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import independent_frohlich2023 as ind  # noqa: E402

STEM = "mapk_adaptive_resistance_frohlich2023"
BNGL = HERE / f"{STEM}.bngl"
REF = HERE / "reference"
UPSTREAM = REF / "frohlich2023_marm2_upstream.json"
AUTHORS_MODEL = REF / "frohlich2023_marm2_RTKERK__pRAF.py"
AUTHORS_PARAMETERS = REF / "frohlich2023_marm2_RTKERK_pRAF_EGF_EGFR_MEKi_PRAFi_RAFi.csv"
OUT_TRAJECTORY = {"authors": REF / "frohlich2023_marm2_fig2a_authors_protocol.csv",
                  "24h": REF / "frohlich2023_marm2_fig2a_24h_protocol.csv"}
OUT_STATE = REF / "frohlich2023_marm2_adapted_state.csv"
OUT_JSON = REF / "frohlich2023_marm2_comparison.json"

# The authors' solver settings (MARM/estimation.py get_solver), used in both arms.
RTOL, ATOL = 1e-9, 1e-11
TIMES = np.linspace(0.0, 2.0, 201)          # the .bngl protocol: t_end 2 h, 200 steps
DRUG = "Vemurafenib"                          # RAFi in rows 0-9 of the authors' data table
CONDITION_ADAPT = {"RAFi_0": 1.0, "EGF_0": 0.0, "EGFR_crispr": 1.0,
                   "N_Avogadro": 6.02214076e23, "volume": 1e-12, "m_Da_EGF": 6200.0}
CONDITION_SIM = dict(CONDITION_ADAPT, EGF_0=100.0)
# Parameters the calibration never estimated: the experimental condition, and the binding
# constants of the two drugs absent from Fig. 2A (their rules fire no reaction here).
NOT_CALIBRATED = {"N_Avogadro", "volume", "m_Da_EGF", "RAFi_0", "PRAFi_0", "MEKi_0", "EGF_0",
                  "EGFR_crispr", "bind_PRAFi_RAF_kf", "bind_PRAFi_RAF_dG", "bind_PRAFi_RAF_phi",
                  "catalyze_RAFrafiNone_MEKmeki_MEKi_p_kcatr", "bind_MEKi_MEK_kf",
                  "bind_MEKi_MEK_dG", "bind_MEKi_MEK_phi", "ep_pMEK_MEKi_ddG"}
LOG0_PLACEHOLDERS = ("bind_PRAFi_RAF_kf", "bind_MEKi_MEK_kf")
# The curated file renames two constants and carries the Fig. 2A doses as its own values.
BNGL_NAME = {"N_Avogadro": "NA", "volume": "V_ref"}
# Adaptation: checkpoints in hours. The 24 h state is taken at its checkpoint; the steady
# state is declared once a checkpoint interval of at least STEADY_SPAN hours moves no species
# above FLOOR by more than STEADY_DRIFT (relative), and must be reached by the last checkpoint.
ADAPT_24H = 24.0
CHECKPOINTS = [1.0, 3.0, 10.0, ADAPT_24H, 30.0, 100.0, 200.0, 300.0, 400.0, 600.0, 800.0, 1000.0]
STEADY_SPAN, STEADY_DRIFT = 100.0, 1e-8
FLOOR = 1e-9                                  # uM, ~6e-4 molecules per 1 pL cell

TOLERANCES = {
    "parameters_max_rel": (1e-12, "the same IEEE doubles; the .bngl prints them to 17 "
                                  "significant digits, so only last-digit rounding may differ"),
    "rhs_max_rel_to_gross_flux": (1e-10, "identical equations evaluated in a different "
                                         "summation order; the denominator is each species' "
                                         "gross turnover sum_j |S_ij r_j|, so cancellation "
                                         "cannot inflate it"),
    "jacobian_max_rel_to_row_max": (1e-10, "the two analytic Jacobians, PySB-derived and "
                                           "hand-written, at the same six states; each row "
                                           "relative to its largest entry"),
    "seeds_max_rel": (1e-5, "the .bngl prints its seeds to 13 significant digits, and both "
                            "arms integrate at rtol 1e-9; relative to max(|value|, "
                            f"{FLOOR:g} uM)"),
    "pERK_IF_max_rel": (1e-5, "two right-hand-side implementations, each integrated at "
                              "rtol 1e-9 / atol 1e-11, from states that agree to the seed "
                              "tolerance; relative to max(|value|, 1e-6 x peak)"),
    "pMEK_IF_max_rel": (1e-5, "as pERK_IF_obs"),
    "observables_max_rel": (1e-5, "as pERK_IF_obs, for each of the 30 observables"),
    "species_max_rel_to_peak": (1e-5, "per species, max_t |diff| / max(peak_t |x_i|, "
                                      f"{FLOOR:g} uM); below the floor is under a thousandth "
                                      "of a molecule per cell"),
}


# ----------------------------------------------------------------------------- provenance


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_pins():
    """Refuse to run unless the vendored files are byte-identical to the pinned release."""
    pins = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    for rel, info in pins["files"].items():
        digest = sha256(HERE / rel)
        if digest != info["sha256"]:
            sys.exit(f"{rel}: sha256 {digest} != pinned {info['sha256']} "
                     f"({pins['upstream_repository']} @ {pins['pinned_commit']})")
    return pins


# ------------------------------------------------------------------ species canonical form


SPECIES_RE = re.compile(r"[A-Za-z_]\w*\([^)]*\)(?:@\w+)?")


def canonical_species(text):
    """Canonical label of a BioNetGen species string.

    Invariant to molecule order, component order and bond numbering: molecules are grouped
    by their bond-free signature, every ordering within each group is tried with bonds
    renumbered in order of appearance, and the lexicographically smallest string wins. The
    complexes here have at most two copies of any molecule type, so the search is small.
    """
    compartment = None
    match = re.match(r"^@(\w+)::(.*)$", text)
    if match:
        compartment, text = match.group(1), match.group(2)
    fixed = text.startswith("$")
    body = text.lstrip("$")
    tokens = SPECIES_RE.findall(body)
    if ".".join(tokens) != body:
        raise ValueError(f"cannot parse species {text!r}")
    molecules = []
    for token in tokens:
        name, inner, mol_comp = re.match(r"^(\w+)\(([^)]*)\)(?:@(\w+))?$", token).groups()
        components = []
        for comp in inner.split(",") if inner else []:
            comp, _, bond = comp.partition("!")
            comp, _, state = comp.partition("~")
            components.append((comp, state, bond))
        components.sort(key=lambda c: (c[0], c[1]))
        molecules.append((name, mol_comp or compartment, components))

    def signature(mol):
        return mol[0], mol[1], tuple((c, s, bool(b)) for c, s, b in mol[2])

    order = sorted(range(len(molecules)), key=lambda i: signature(molecules[i]))
    groups = [list(g) for _, g in itertools.groupby(order, key=lambda i: signature(molecules[i]))]
    best = None
    for arrangement in itertools.product(*(itertools.permutations(g) for g in groups)):
        relabel, parts = {}, []
        for i in (i for group in arrangement for i in group):
            name, mol_comp, components = molecules[i]
            rendered = []
            for comp, state, bond in components:
                if bond:
                    bond = relabel.setdefault(bond, str(len(relabel) + 1))
                rendered.append(comp + (f"~{state}" if state else "")
                                + (f"!{bond}" if bond else ""))
            parts.append(f"{name}({','.join(rendered)})@{mol_comp}")
        label = ("$" if fixed else "") + f"@{compartment}::" + ".".join(parts)
        if best is None or label < best:
            best = label
    return best


def net_sections(text):
    """Species strings and (reactants, products, rule) triples of a `.net`, in file order."""
    species, reactions, block = [], [], None
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("begin "):
            block = line.split()[1]
            continue
        if line.startswith("end "):
            block = None
            continue
        body, _, comment = line.partition("#")
        fields = body.split()
        if not fields:
            continue
        if block == "species":
            species.append(fields[1])
        elif block == "reactions":
            rule = comment.split()[0]
            reactants = [int(x) - 1 for x in fields[1].split(",") if x != "0"]
            products = [int(x) - 1 for x in fields[2].split(",") if x != "0"]
            reactions.append((reactants, products, rule))
    return species, reactions


def strip_net(text):
    """`.net` content with comments and whitespace removed (house reproducibility rule)."""
    out = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            out.append(" ".join(line.split()))
    return out


# ------------------------------------------------------------------------- curated arm


def curated_network(workdir, bngpath):
    """Expand the pinned `.bngl` with BNG2.pl (generate_network only) and parse it."""
    text = BNGL.read_text(encoding="utf-8")
    start = text.index("begin actions")
    end = text.index("end actions") + len("end actions")
    run_dir = workdir / "curated"
    run_dir.mkdir()
    (run_dir / BNGL.name).write_text(
        text[:start] + "begin actions\n  generate_network({overwrite=>1})\nend actions"
        + text[end:], encoding="utf-8")
    proc = subprocess.run(["perl", str(Path(bngpath) / "BNG2.pl"), BNGL.name], cwd=run_dir,
                          capture_output=True, text=True)
    if proc.returncode:
        sys.exit(f"BNG2.pl failed on the curated model:\n{proc.stdout[-3000:]}\n{proc.stderr}")
    net_path = run_dir / f"{STEM}.net"
    committed = (REF / f"{STEM}.net").read_text(encoding="utf-8")
    if strip_net(net_path.read_text(encoding="utf-8")) != strip_net(committed):
        sys.exit("the regenerated .net does not match reference/ structurally; regenerate the "
                 "committed reference output before running this check")
    version = re.search(r"BioNetGen version ([\d.]+)", proc.stdout)
    return net_path, (version.group(1) if version else "unknown")


# ------------------------------------------------------------------------- authors' arm


def load_authors_model():
    """The authors' PySB model at calibrated row 0, as the Fig. 2A (EGF_RAFi) instance."""
    spec = importlib.util.spec_from_file_location("RTKERK__pRAF", AUTHORS_MODEL)
    module = importlib.util.module_from_spec(spec)
    # No bytecode cache: reference/ must stay flat (skills/bngl/skill.md section 5.4 rule 4).
    previous, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    model = module.model

    with open(AUTHORS_PARAMETERS, newline="") as handle:
        rows = list(csv.DictReader(handle))
    fvals = [float(r["fval"]) for r in rows]
    if fvals.index(min(fvals)) != 0 or rows[0][""] != "0":
        sys.exit("row 0 of the authors' parameter table is not the best-ranked set")
    row0 = {k: float(v) for k, v in rows[0].items() if k not in ("", "fval")}

    def table_name(name):
        # MARM.parameters.specialise_par_name(name, prafi=None, rafi=DRUG, meki=None)
        return name.replace("RAFi", DRUG)

    applied, unmatched = {}, set()
    for parameter in model.parameters:
        key = table_name(parameter.name)
        if key in row0:
            parameter.value = row0[key]
            applied[parameter.name] = row0[key]
        else:
            unmatched.add(parameter.name)
    if unmatched != NOT_CALIBRATED:
        sys.exit(f"unexpected parameters without a calibrated value: "
                 f"{sorted(unmatched ^ NOT_CALIBRATED)}")
    for name, value in CONDITION_ADAPT.items():
        model.parameters[name].value = value
    for name in LOG0_PLACEHOLDERS:
        model.parameters[name].value = 1.0
    # MARM.model.get_model_instance(..., 'EGF_RAFi', ...): keep EGF and RAFi initials only.
    model.initials = [i for i in model.initials if i.value.name not in ("PRAFi_0", "MEKi_0")]
    return model, applied


def authors_network(model, workdir):
    """PySB network generation (BNG2.pl), kept as text so its species can be labelled."""
    import pysb.bng

    text = pysb.bng.generate_network(model)
    path = workdir / "authors_RTKERK__pRAF__EGF_RAFi.net"
    path.write_text(text, encoding="utf-8")
    pysb.bng.load_equations(model, str(path))
    # The placeholders enter only Ea0_bind_PRAFi_RAF and Ea0_bind_MEKi_MEK, i.e. the PRAFi and
    # MEKi binding rules. With no PRAFi or MEKi in the instance those rules must fire nothing.
    fired = {rule for reaction in model.reactions for rule in reaction["rule"]}
    placeholder_rules = {r.name for r in model.rules if "PRAFi" in r.name or "MEKi" in r.name}
    if fired & placeholder_rules:
        sys.exit(f"rules that use the log(0) placeholders fired: "
                 f"{sorted(fired & placeholder_rules)}")
    return text


def authors_ode(model, simulator):
    """PySB's compiled right-hand side for the authors' network, and its analytic Jacobian.

    The right-hand side is ScipyOdeSimulator's own (its RHS builder compiles PySB's rate
    expressions). The Jacobian follows PySB's own construction -- the derivative of the
    kinetics vector with respect to species plus that with respect to observables times the
    observables matrix, left-multiplied by the stoichiometry -- but kept sparse: PySB's dense
    version would be a 17,516 x 1,612 symbolic matrix. EGF_0 enters only the initial value of
    free EGF, so one parameter vector serves both phases.
    """
    import sympy
    from sympy.matrices.expressions.matexpr import MatrixElement

    builder = simulator.rhs_builder
    p_vec = np.array([CONDITION_SIM.get(p.name, p.value) for p in model.parameters_all()],
                     dtype=float)
    e_vec = builder.calc_expressions_constant(p_vec)

    def rhs(t, y):
        return builder.rhs_fn(t, y, p_vec, e_vec)

    entries = {"y": ([], [], []), "o": ([], [], [])}
    for j, rate in enumerate(builder.kinetics):
        for atom in rate.atoms(MatrixElement):
            key = "y" if atom.parent == builder.y else "o" if atom.parent == builder.o else None
            if key:
                rows, cols, exprs = entries[key]
                rows.append(j)
                cols.append(int(atom.i))
                exprs.append(rate.diff(atom))
    n_y = len(entries["y"][2])
    derivatives = sympy.lambdify([builder.y, builder.p, builder.e, builder.o],
                                 entries["y"][2] + entries["o"][2])
    observables = builder.observables_matrix.astype(float)
    stoichiometry = builder.stoichiometry_matrix.astype(float)
    shape_y = (len(model.reactions), len(model.species))
    shape_o = (len(model.reactions), observables.shape[0])

    def jac(t, y):
        values = np.asarray(derivatives(y[:, None], p_vec[:, None], e_vec,
                                        (observables @ y)[:, None]), dtype=float).ravel()
        d_y = csc_matrix((values[:n_y], entries["y"][:2]), shape=shape_y)
        d_o = csc_matrix((values[n_y:], entries["o"][:2]), shape=shape_o)
        return csc_matrix(stoichiometry @ (d_y + d_o @ observables))

    return rhs, jac


def integrate_authors(rhs, jac, t_span, y0, t_eval=None):
    solution = solve_ivp(rhs, t_span, y0, method="BDF", t_eval=t_eval,
                         rtol=RTOL, atol=ATOL, jac=jac)
    if not solution.success:
        sys.exit(f"authors' model integration failed on {t_span}: {solution.message}")
    return solution


def adapt(rhs, jac, y_start, log):
    """Integrate from the model's initial values at 1 uM vemurafenib, no EGF.

    Returns the 24 h state and the steady state. The authors' AMICI pre-equilibration also
    reaches its steady state by simulation (their get_solver sets Newton steps to 0).
    """
    y, t0, at_24h, previous = np.asarray(y_start, dtype=float), 0.0, None, {}
    for t1 in CHECKPOINTS:
        start = time.time()
        y = integrate_authors(rhs, jac, (t0, t1), y).y[:, -1]
        previous[t1] = y
        if t1 == ADAPT_24H:
            at_24h = y.copy()
        span_start = max((t for t in previous if t <= t1 - STEADY_SPAN), default=None)
        drift = (float(np.max(np.abs(y - previous[span_start])
                              / np.maximum(np.abs(y), FLOOR)))
                 if span_start is not None else np.inf)
        log(f"  adaptation t = {t1:g} h ({time.time() - start:.0f} s); max relative change "
            f"since t = {span_start if span_start is not None else '-'} h: {drift:.2e}")
        t0 = t1
        if drift <= STEADY_DRIFT:
            return at_24h, y, {"t_h": t1, "since_h": span_start, "max_rel_change": drift}
    sys.exit(f"no steady state by t = {CHECKPOINTS[-1]:g} h (last change {drift:.2e})")


# ----------------------------------------------------------------------------- metrics


def cython_version():
    try:
        import Cython
    except ImportError:
        return None
    return Cython.__version__


def rel_err(a, b, floor):
    return np.abs(a - b) / np.maximum(np.abs(b), floor)


def verdict(name, value, extra=None):
    tol, why = TOLERANCES[name]
    out = {"value": float(value), "tol": tol, "pass": bool(value <= tol), "why": why}
    if extra:
        out.update(extra)
    return out


def protocol_metrics(net, curated_y, authors_y, authors_y0, authors_result, species_text):
    """Seed and trajectory comparisons for one adaptation protocol (curated species order)."""
    p_c = net["params"]
    seeded = net["y0"] != 0
    seed_err = rel_err(net["y0"], authors_y0, FLOOR)
    perk_c = ind.perk_if(net, curated_y)
    perk_a = np.asarray(authors_result.expressions["pERK_IF_obs"], dtype=float)
    pmek_c = (ind.group_value(net, "pMEK", curated_y) / ind.group_value(net, "tMEK", curated_y)
              * p_c["pMEK_IF_scale"] + p_c["pMEK_IF_offset"])
    pmek_a = np.asarray(authors_result.expressions["pMEK_IF_obs"], dtype=float)
    perk_err = rel_err(perk_c, perk_a, 1e-6 * np.max(np.abs(perk_a)))
    pmek_err = rel_err(pmek_c, pmek_a, 1e-6 * np.max(np.abs(pmek_a)))
    obs_worst, obs_name = 0.0, ""
    for name in net["groups"]:
        ref = np.asarray(authors_result.observables[name], dtype=float)
        mine = ind.group_value(net, name, curated_y)
        peak = float(np.max(np.abs(ref)))
        err = (float(np.max(np.abs(mine))) if peak == 0.0
               else float(np.max(rel_err(mine, ref, 1e-6 * peak))))
        if err > obs_worst:
            obs_worst, obs_name = err, name
    peak_i = np.maximum(np.max(np.abs(authors_y), axis=1), FLOOR)
    species_err = np.max(np.abs(curated_y - authors_y), axis=1) / peak_i
    worst_species = int(np.argmax(species_err))
    worst_seed = int(np.argmax(np.where(seeded, seed_err, 0.0)))
    return {
        "seeds_max_rel": verdict("seeds_max_rel", float(np.max(seed_err[seeded])), {
            "n_seeded": int(np.sum(seeded)), "median": float(np.median(seed_err[seeded])),
            "worst": species_text[worst_seed],
            "unseeded_authors_max_uM": float(np.max(np.abs(authors_y0[~seeded]))),
            "pERK_IF_at_t0": {"bngl": float(perk_c[0]), "authors": float(perk_a[0])}}),
        "pERK_IF_max_rel": verdict("pERK_IF_max_rel", float(np.max(perk_err)), {
            "at_time_h": float(TIMES[int(np.argmax(perk_err))]),
            "max_abs": float(np.max(np.abs(perk_c - perk_a))),
            "authors_peak": float(np.max(perk_a)),
            "authors_peak_time_h": float(TIMES[int(np.argmax(perk_a))])}),
        "pMEK_IF_max_rel": verdict("pMEK_IF_max_rel", float(np.max(pmek_err))),
        "observables_max_rel": verdict("observables_max_rel", obs_worst, {
            "n": len(net["groups"]), "worst": obs_name}),
        "species_max_rel_to_peak": verdict(
            "species_max_rel_to_peak", float(species_err[worst_species]), {
                "n": len(species_err), "worst": species_text[worst_species],
                "median": float(np.median(species_err))}),
    }


def write_trajectory(path, result):
    """Every observable and output expression of the authors' model on the 0-2 h grid."""
    observables = list(result.observables.dtype.names)
    expressions = [n for n in result.expressions.dtype.names if n.endswith("_obs")]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["time_h", *expressions, *observables])
        for k, t in enumerate(TIMES):
            writer.writerow([f"{t:.10g}"]
                            + [f"{result.expressions[n][k]:.15g}" for n in expressions]
                            + [f"{result.observables[n][k]:.15g}" for n in observables])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--compiler", choices=("cython", "python"), default="cython",
                        help="PySB right-hand-side backend (default: cython)")
    parser.add_argument("--bngpath", default=os.environ.get("BNGPATH"),
                        help="BioNetGen directory holding BNG2.pl (default: $BNGPATH)")
    parser.add_argument("--keep", action="store_true", help="keep the scratch directory")
    args = parser.parse_args()
    if not args.bngpath or not (Path(args.bngpath) / "BNG2.pl").exists():
        sys.exit("set BNGPATH (or --bngpath) to a BioNetGen directory containing BNG2.pl")
    os.environ["BNGPATH"] = args.bngpath  # PySB finds BNG2.pl through it

    import pysb
    import scipy
    import sympy
    from pysb.simulator import ScipyOdeSimulator
    from pysb.simulator.base import SimulationResult

    started = time.time()

    def log(message):
        print(f"[{time.time() - started:6.0f} s] {message}", flush=True)

    pins = check_pins()
    log(f"vendored MARM2.0 files match {pins['upstream_repository']} @ "
        f"{pins['pinned_commit'][:7]}")
    workdir = Path(tempfile.mkdtemp(prefix="frohlich2023_marm2_"))
    try:
        # --- curated arm -------------------------------------------------------------------
        net_path, bng_version = curated_network(workdir, args.bngpath)
        net = ind.parse_net(net_path)
        log(f"curated .bngl expanded by BNG2.pl {bng_version}: {len(net['y0'])} species, "
            f"{len(net['reactions'])} reactions (matches reference/{STEM}.net)")
        curated = ind.integrate(net, TIMES, rtol=RTOL, atol=ATOL, analytic_jacobian=True)
        log("curated network integrated 0-2 h from its seeds (independent SciPy BDF)")

        # --- authors' arm: model, network, right-hand side ---------------------------------
        model, applied = load_authors_model()
        authors_text = authors_network(model, workdir)
        log(f"authors' MARM2.0 instance EGF_RAFi expanded by PySB {pysb.__version__}: "
            f"{len(model.species)} species, {len(model.reactions)} reactions")
        simulator = ScipyOdeSimulator(model, compiler=args.compiler)
        rhs_a, jac_a = authors_ode(model, simulator)
        log(f"authors' right-hand side compiled by PySB ({args.compiler}), Jacobian derived")

        # --- species bijection and reactions per rule --------------------------------------
        species_a, reactions_a = net_sections(authors_text)
        species_c, reactions_c = net_sections(net_path.read_text(encoding="utf-8"))
        label_c = [canonical_species(s) for s in species_c]
        label_a = [canonical_species(s) for s in species_a]
        if len(set(label_c)) != len(label_c) or set(label_c) != set(label_a):
            sys.exit("the two networks do not have the same species")
        index_c = {label: i for i, label in enumerate(label_c)}
        to_curated = np.array([index_c[label] for label in label_a])  # authors i -> curated

        def reaction_key(reactants, products, rule, mapping):
            return (tuple(sorted(mapping[i] for i in reactants)),
                    tuple(sorted(mapping[i] for i in products)), rule)

        identity = np.arange(len(label_c))
        if (Counter(reaction_key(*r, identity) for r in reactions_c)
                != Counter(reaction_key(*r, to_curated) for r in reactions_a)):
            sys.exit("the two networks do not have the same reactions")

        def curated_order(x):
            out = np.empty_like(x)
            out[to_curated] = x
            return out

        # --- parameters --------------------------------------------------------------------
        p_c = net["params"]
        worst_param, worst_name = 0.0, ""
        for name, value in applied.items():
            mine = p_c[BNGL_NAME.get(name, name)]
            err = abs(mine - value) / abs(value) if value else abs(mine)
            if err > worst_param:
                worst_param, worst_name = err, name
        for name in ("N_Avogadro", "volume", "m_Da_EGF", "EGFR_crispr"):
            err = abs(p_c[BNGL_NAME.get(name, name)] / CONDITION_SIM[name] - 1.0)
            if err > worst_param:
                worst_param, worst_name = err, name
        doses_ok = (p_c["RAFi_0"] == CONDITION_SIM["RAFi_0"]
                    and p_c["EGF_0"] == CONDITION_SIM["EGF_0"])

        # --- adaptation, both protocols, then the EGF pulse --------------------------------
        y_start = np.asarray(simulator.initials, dtype=float)[0]   # PySB's own initials
        y_24h, y_steady, steady = adapt(rhs_a, jac_a, y_start, log)
        egf_free = species_a.index("@CP::$EGF(rtk)")
        for name, value in CONDITION_SIM.items():
            model.parameters[name].value = value
        egf_value = float(model.expressions["initEGF"].expand_expr().xreplace(
            {p: p.value for p in model.parameters}))
        pulses, t0_states = {}, {}
        for label, adapted in (("authors", y_steady), ("24h", y_24h)):
            y0 = adapted.copy()
            y0[egf_free] = egf_value
            start = time.time()
            solution = integrate_authors(rhs_a, jac_a, (TIMES[0], TIMES[-1]), y0,
                                         t_eval=TIMES)
            pulses[label] = (solution.y, SimulationResult(simulator, [TIMES], [solution.y.T]))
            t0_states[label] = y0
            log(f"authors' model, {label} protocol, integrated 0-2 h "
                f"({time.time() - start:.0f} s)")

        # --- equations: dy/dt at six states ------------------------------------------------
        rhs_c = ind.make_rhs(net)
        jac_c = ind.make_jacobian(net)
        abs_s = abs(net["S"])
        moving = ~net["fixed"]
        egf_c = int(to_curated[egf_free])
        seeds_without_egf = net["y0"].copy()
        seeds_without_egf[egf_c] = 0.0
        rng = np.random.default_rng(20231988)
        states = {"curated seeds": net["y0"],
                  "authors' steady state + EGF": curated_order(t0_states["authors"]),
                  "authors' 24 h state + EGF": curated_order(t0_states["24h"]),
                  "authors' t = 0.1 h (steady-state protocol)":
                      curated_order(pulses["authors"][0][:, 10])}
        for k in range(2):
            states[f"random state {k + 1}"] = 10.0 ** rng.uniform(-6, 0, len(net["y0"]))
        permute = csc_matrix((np.ones(len(to_curated)), (to_curated, np.arange(len(to_curated)))),
                             shape=(len(to_curated), len(to_curated)))
        rhs_worst = jac_worst = 0.0
        for state in states.values():
            f_c = rhs_c(0.0, state)
            f_a = curated_order(np.asarray(rhs_a(0.0, state[to_curated]), dtype=float).ravel())
            gross = abs_s @ np.abs(ind.fluxes(net, state))
            denom = np.maximum(gross[moving], 1e-300)
            rhs_worst = max(rhs_worst, float(np.max(np.abs(f_c - f_a)[moving] / denom)))
            j_c = jac_c(0.0, state).tocsr()
            j_a = (permute @ jac_a(0.0, state[to_curated]) @ permute.T).tocsr()
            row_max = np.maximum(abs(j_c).max(axis=1).toarray().ravel(), 1e-300)
            diff = abs(j_c - j_a).max(axis=1).toarray().ravel()
            jac_worst = max(jac_worst, float(np.max(diff / row_max)))
        # How far from steady the .bngl seeds are, in the measures AMICI uses (get_solver:
        # steady-state rtol 1e-7, atol 1e-9; converged when the WRMS of dy/dt is below 1).
        steady_c = curated_order(y_steady)

        def amici_wrms(state):
            f = rhs_c(0.0, state)[moving]
            return float(np.sqrt(np.mean((f / (1e-7 * np.abs(state[moving]) + 1e-9)) ** 2)))

        f_seeds = rhs_c(0.0, seeds_without_egf)
        above = moving & (seeds_without_egf > FLOOR)

        model_results = {
            "parameters_max_rel": verdict("parameters_max_rel", worst_param, {
                "n_calibrated": len(applied), "worst": worst_name,
                "doses_match": bool(doses_ok)}),
            "network": {"species": len(label_c), "reactions": len(reactions_c),
                        "species_bijection": True, "reactions_identical_per_rule": True,
                        "pass": True},
            "rhs_max_rel_to_gross_flux": verdict("rhs_max_rel_to_gross_flux", rhs_worst, {
                "states": list(states)}),
            "jacobian_max_rel_to_row_max": verdict("jacobian_max_rel_to_row_max", jac_worst, {
                "states": list(states)}),
        }
        protocols = {}
        for label in ("authors", "24h"):
            authors_y, result = pulses[label]
            protocols[label] = protocol_metrics(
                net, curated.y, curated_order(authors_y), curated_order(t0_states[label]),
                result, species_c)
        seeds_steadiness = {
            "note": "dy/dt of the model at the .bngl seeds with free EGF removed, against the "
                    "authors' steady state; AMICI's pre-equilibration stops below WRMS 1",
            "amici_wrms_dydt_at_bngl_seeds": amici_wrms(seeds_without_egf),
            "amici_wrms_dydt_at_authors_steady_state": amici_wrms(steady_c),
            "max_rel_rate_per_h_seeded_above_floor": float(np.max(
                np.abs(f_seeds[above]) / seeds_without_egf[above])),
            "n_unseeded_produced": int(np.sum(moving & (seeds_without_egf == 0)
                                              & (np.abs(f_seeds) > 0))),
        }

        def passed(block):
            return all(v.get("pass", True) for v in block.values() if isinstance(v, dict))

        model_ok = passed(model_results) and doses_ok
        model_equivalent = model_ok and passed(protocols["24h"])
        all_pass = model_equivalent and passed(protocols["authors"])

        # --- outputs -----------------------------------------------------------------------
        for label, path in OUT_TRAJECTORY.items():
            write_trajectory(path, pulses[label][1])
        with open(OUT_STATE, "w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["net_index", "species", "bngl_seed_uM",
                             "authors_steady_state_uM", "authors_24h_uM"])
            steady_c, h24_c = curated_order(t0_states["authors"]), curated_order(t0_states["24h"])
            for i, text in enumerate(species_c):
                writer.writerow([i + 1, text, f"{net['y0'][i]:.15g}", f"{steady_c[i]:.15g}",
                                 f"{h24_c[i]:.15g}"])
        summary = {
            "_comment": "Written by run_frohlich2023.py; do not edit by hand.",
            "bngl_sha256": sha256(BNGL),
            "authors_release": {"repository": pins["upstream_repository"],
                                "commit": pins["pinned_commit"],
                                "files": {k: v["sha256"] for k, v in pins["files"].items()},
                                "parameter_row": 0, "drug": DRUG},
            "protocol": {
                "adaptation_condition": CONDITION_ADAPT, "simulation_condition": CONDITION_SIM,
                "authors": "steady state at the adaptation condition (AMICI pre-equilibration; "
                           "processed data rows 0-9: *_preeq set, t_presim 0)",
                "24h": "24 h at the adaptation condition from the model's initial values",
                "steady_state_reached": steady, "pulse_times_h": [0.0, 2.0, len(TIMES)],
                "rtol": RTOL, "atol": ATOL},
            "software": {"pysb": pysb.__version__, "scipy": scipy.__version__,
                         "numpy": np.__version__, "sympy": sympy.__version__,
                         "bionetgen": bng_version, "pysb_compiler": args.compiler,
                         "cython": cython_version()},
            "model": model_results,
            "protocols": protocols,
            "bngl_seeds_steadiness": seeds_steadiness,
            "model_equivalent": bool(model_equivalent),
            "all_pass": bool(all_pass),
        }
        OUT_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")
        log("wrote " + ", ".join(p.name for p in (*OUT_TRAJECTORY.values(), OUT_STATE, OUT_JSON)))
    finally:
        if args.keep:
            log(f"scratch kept in {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    def show(name, res):
        print(f"  {'PASS' if res['pass'] else 'FAIL'}  {name:28s} {res['value']:.3e}"
              f"  (tol {res['tol']:.0e})")

    print("\nmodel (protocol-independent)")
    for name, res in model_results.items():
        if "value" in res:
            show(name, res)
    print(f"  {'PASS' if doses_ok else 'FAIL'}  {'doses in the .bngl':28s} RAFi_0 "
          f"{p_c['RAFi_0']:g} uM, EGF_0 {p_c['EGF_0']:g} ng/mL")
    print(f"  PASS  {'network':28s} {len(label_c)} species, {len(reactions_c)} "
          f"reactions, identical")
    for label, title in (("24h", "24 h adaptation (what the .bngl seeds encode)"),
                         ("authors", "steady-state pre-equilibration (the authors' protocol)")):
        print(f"\n{title}")
        for name, res in protocols[label].items():
            show(name, res)
    print(f"\nmodel equivalent to the authors' release: {'YES' if model_equivalent else 'NO'}")
    print(f"reproduces the authors' Fig. 2A protocol:  {'YES' if all_pass else 'NO'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
