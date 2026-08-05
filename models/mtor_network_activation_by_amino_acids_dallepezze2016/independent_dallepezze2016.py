"""Independent implementation of the Dalle Pezze et al. (2016) mTOR network model.

Level-1 verification for ``models/mtor_network_activation_by_amino_acids_dallepezze2016``.
Two arms, deliberately different in kind:

``rhs`` / ``integrate``
    A hand-written ODE system: 31 state variables and the 48 mass-action reactions of
    BioModels BIOMD0000000640, transcribed one by one below. It shares no code with
    BioNetGen and never reads the generated ``.net``. This is what tests that the BNGL
    rules expand to the intended network -- a rule that matches the wrong context, or a
    rate constant attached to the wrong context, shows up here and nowhere else.

``sbml_reference``
    A *mechanical* integration of ``BIOMD0000000640.xml`` itself: libSBML parses the
    file, the kinetic laws are evaluated from their MathML, and nothing is transcribed
    by hand. The hand-written arm above and the ``.bngl`` share one failure mode --
    both are somebody reading the SBML -- and this arm is the only thing that can catch
    it. It needs the deposited SBML, which lives under ``dev/papers/`` and is not
    committed; when it is absent the function says where to get it instead of failing.

Reaction numbering in the comments is the SBML's own (``reaction_1`` ... ``reaction_48``),
so any line here can be checked against the deposit in one grep.

Run this file directly for a self-test: it integrates all three published conditions and
prints the peak and t=120 value of every readout.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# --- state vector layout -------------------------------------------------------------
# 31 dynamic species, in the order BIOMD0000000640 declares them.
(
    IR_u, IR_p, IR_r,
    IRS_u, IRS_p, IRS_pS636,
    AMPK_u, AMPK_p,
    Akt_00, Akt_T, Akt_S, Akt_TS,
    TSC_u, TSC_T1462, TSC_S1387,
    mTORC1_u, mTORC1_p,
    mTORC2_u, mTORC2_p,
    S6K_00, S6K_229, S6K_389, S6K_both,
    PRAS_00, PRAS_246, PRAS_183, PRAS_both,
    PI3Kv_u, PI3Kv_p,
    PI3K_u, PI3K_p,
) = range(31)

N_STATES = 31

# --- initial pool sizes (BIOMD0000000640 initial amounts; arbitrary relative units) ---
POOLS = {
    "IR": 50.0, "IRS": 150.0, "AMPK": 50.0, "Akt": 300.0, "TSC": 50.0,
    "mTORC1": 100.0, "mTORC2": 100.0, "S6K": 300.0, "PRAS40": 20.0,
    "PI3K_var": 50.0, "PI3K_PDK1": 50.0,
}

# --- fitted rate constants (BIOMD0000000640 listOfParameters, verbatim) --------------
P = {
    "IRS_phos_by_Amino_Acids": 0.0331672,
    "AMPK_T172_phos_by_Amino_Acids": 17.6284,
    "mTORC2_S2481_phos_by_Amino_Acids": 0.0268658,
    "mTORC1_S2448_activation_by_Amino_Acids": 0.0156992,
    "IR_beta_phos_by_Insulin": 0.0203796,
    "IR_beta_pY1146_dephos": 0.493514,
    "IR_beta_ready": 323.611,
    "IRS_phos_by_IR_beta_pY1146": 2.11894,
    "IRS_p_phos_by_p70_S6K_pT229_pT389": 0.338859859949792,
    "IRS_phos_by_p70_S6K_pT229_pT389": 0.0863775267376444,
    "IRS_pS636_turnover": 25.0,
    "PI3K_PDK1_phos_by_IRS_p": 0.000187226757782201,
    "PI3K_p_PDK1_dephos": 0.18913343080532,
    "PI3K_variant_phos_by_IR_beta_pY1146": 0.000549027801822575,
    "PI3K_variant_p_dephos": 0.108074886441184,
    "AMPK_T172_phos": 0.490602,
    "AMPK_pT172_dephos": 165.704,
    "Akt_T308_phos_by_PI3K_p_PDK1_first": 7.47437,
    "Akt_T308_phos_by_PI3K_p_PDK1_second": 7.47345,
    "Akt_S473_phos_by_mTORC2_pS2481_first": 1.31992e-05,
    "Akt_S473_phos_by_mTORC2_pS2481_second": 0.159093,
    "Akt_pT308_dephos_first": 88.9654,
    "Akt_pT308_dephos_second": 88.9639,
    "Akt_pS473_dephos_first": 0.376999,
    "Akt_pS473_dephos_second": 0.380005,
    "mTORC2_S2481_phos_by_PI3K_variant_p": 0.120736,
    "mTORC2_pS2481_dephos": 1.42511,
    "TSC1_TSC2_S1387_phos_by_AMPK_pT172": 0.00175772,
    "TSC1_TSC2_T1462_phos_by_Akt_pT308": 1.52417,
    "TSC1_TSC2_pS1387_dephos": 0.25319,
    "TSC1_TSC2_pT1462_dephos": 147.239,
    "mTORC1_pS2448_dephos_by_TSC1_TSC2": 0.00869774,
    "p70_S6K_T389_phos_by_mTORC1_pS2448_first": 0.00261303413778722,
    "p70_S6K_T389_phos_by_mTORC1_pS2448_second": 0.110720890919343,
    "p70_S6K_T229_phos_by_PI3K_p_PDK1_first": 0.0133520172873009,
    "p70_S6K_T229_phos_by_PI3K_p_PDK1_second": 1.00000002814509e-06,
    "p70_S6K_pT389_dephos_first": 1.10036057608758,
    "p70_S6K_pT389_dephos_second": 1.10215267954479,
    "p70_S6K_pT229_dephos_first": 1.00000012897033e-06,
    "p70_S6K_pT229_dephos_second": 0.159201353240651,
    "PRAS40_S183_phos_by_mTORC1_pS2448_first": 0.15881,
    "PRAS40_S183_phos_by_mTORC1_pS2448_second": 0.0683009,
    "PRAS40_T246_phos_by_Akt_pT308_first": 0.279344,
    "PRAS40_T246_phos_by_Akt_pT308_second": 0.279401,
    "PRAS40_pS183_dephos_first": 1.8706,
    "PRAS40_pS183_dephos_second": 1.88453,
    "PRAS40_pT246_dephos_first": 11.8759,
    "PRAS40_pT246_dephos_second": 11.876,
}

# The three published conditions. `pi3k` scales both PI3K pools, which is how the
# authors simulated wortmannin ("PI3K inhibition: residual activity 10%", Fig. 1d).
CONDITIONS = {
    "aa": dict(insulin=0.0, amino_acids=1.0, pi3k=1.0),          # Fig. 2c
    "aa_insulin": dict(insulin=1.0, amino_acids=1.0, pi3k=1.0),  # Fig. 2b
    "aa_wortmannin": dict(insulin=0.0, amino_acids=1.0, pi3k=0.1),  # Fig. 2d
}

# Readout -> the .gdat observable it corresponds to, for the notebook's comparison.
OBSERVABLES = (
    "IR_pY1146", "IRS_pS636", "AMPK_pT172", "Akt_pT308", "Akt_pS473",
    "TSC_pS1387", "mTOR_pS2448", "mTOR_pS2481", "S6K_pT229", "S6K_pT389",
    "PRAS40_pT246", "PRAS40_pS183",
)


def initial_state(pi3k: float = 1.0, mtorc1: float = 1.0, mtorc2: float = 1.0,
                  ampk: float = 1.0) -> np.ndarray:
    """Fully dephosphorylated start; the four factors scale a node's total pool."""
    y = np.zeros(N_STATES)
    y[IR_u] = POOLS["IR"]
    y[IRS_u] = POOLS["IRS"]
    y[AMPK_u] = ampk * POOLS["AMPK"]
    y[Akt_00] = POOLS["Akt"]
    y[TSC_u] = POOLS["TSC"]
    y[mTORC1_u] = mtorc1 * POOLS["mTORC1"]
    y[mTORC2_u] = mtorc2 * POOLS["mTORC2"]
    y[S6K_00] = POOLS["S6K"]
    y[PRAS_00] = POOLS["PRAS40"]
    y[PI3Kv_u] = pi3k * POOLS["PI3K_var"]
    y[PI3K_u] = pi3k * POOLS["PI3K_PDK1"]
    return y


def rhs(t: float, y: np.ndarray, insulin: float, amino_acids: float) -> np.ndarray:
    """The 48 reactions of BIOMD0000000640 as a mass-action right-hand side."""
    d = np.zeros(N_STATES)

    def flux(rate, *, consume=(), produce=()):
        for i in consume:
            d[i] -= rate
        for i in produce:
            d[i] += rate

    # Two composite modifiers the SBML builds inside function_16/17/25/26: the total
    # T308-phosphorylated Akt, and the TSC species that inhibit mTORC1 (naive and
    # S1387-phosphorylated -- the Akt-inhibited pT1462 form does not).
    akt_pT308_total = y[Akt_T] + y[Akt_TS]
    tsc_inhibitory = y[TSC_u] + y[TSC_S1387]

    # --- insulin receptor: reaction_4, reaction_5, reaction_6
    flux(P["IR_beta_phos_by_Insulin"] * insulin * y[IR_u], consume=(IR_u,), produce=(IR_p,))
    flux(P["IR_beta_pY1146_dephos"] * y[IR_p], consume=(IR_p,), produce=(IR_r,))
    flux(P["IR_beta_ready"] * y[IR_r], consume=(IR_r,), produce=(IR_u,))

    # --- IRS1: reaction_7, reaction_1, reaction_8, reaction_9, reaction_10
    flux(P["IRS_phos_by_IR_beta_pY1146"] * y[IR_p] * y[IRS_u],
         consume=(IRS_u,), produce=(IRS_p,))
    flux(P["IRS_phos_by_Amino_Acids"] * amino_acids * y[IRS_u],
         consume=(IRS_u,), produce=(IRS_p,))
    flux(P["IRS_p_phos_by_p70_S6K_pT229_pT389"] * y[S6K_both] * y[IRS_p],
         consume=(IRS_p,), produce=(IRS_pS636,))
    flux(P["IRS_phos_by_p70_S6K_pT229_pT389"] * y[S6K_both] * y[IRS_u],
         consume=(IRS_u,), produce=(IRS_pS636,))
    flux(P["IRS_pS636_turnover"] * y[IRS_pS636], consume=(IRS_pS636,), produce=(IRS_u,))

    # --- PI3K activities: reaction_12, reaction_11, reaction_14, reaction_13
    flux(P["PI3K_PDK1_phos_by_IRS_p"] * y[IRS_p] * y[PI3K_u],
         consume=(PI3K_u,), produce=(PI3K_p,))
    flux(P["PI3K_p_PDK1_dephos"] * y[PI3K_p], consume=(PI3K_p,), produce=(PI3K_u,))
    flux(P["PI3K_variant_phos_by_IR_beta_pY1146"] * y[IR_p] * y[PI3Kv_u],
         consume=(PI3Kv_u,), produce=(PI3Kv_p,))
    flux(P["PI3K_variant_p_dephos"] * y[PI3Kv_p], consume=(PI3Kv_p,), produce=(PI3Kv_u,))

    # --- AMPK: reaction_15, reaction_2, reaction_16
    flux(P["AMPK_T172_phos"] * y[IRS_p] * y[AMPK_u], consume=(AMPK_u,), produce=(AMPK_p,))
    flux(P["AMPK_T172_phos_by_Amino_Acids"] * amino_acids * y[AMPK_u],
         consume=(AMPK_u,), produce=(AMPK_p,))
    flux(P["AMPK_pT172_dephos"] * y[AMPK_p], consume=(AMPK_p,), produce=(AMPK_u,))

    # --- Akt: reaction_17..reaction_24. `_first` acts on Akt unmodified at the other
    # site, `_second` on Akt already modified there; the constants are independent.
    flux(P["Akt_T308_phos_by_PI3K_p_PDK1_first"] * y[PI3K_p] * y[Akt_00],
         consume=(Akt_00,), produce=(Akt_T,))
    flux(P["Akt_S473_phos_by_mTORC2_pS2481_first"] * y[mTORC2_p] * y[Akt_00],
         consume=(Akt_00,), produce=(Akt_S,))
    flux(P["Akt_T308_phos_by_PI3K_p_PDK1_second"] * y[PI3K_p] * y[Akt_S],
         consume=(Akt_S,), produce=(Akt_TS,))
    flux(P["Akt_S473_phos_by_mTORC2_pS2481_second"] * y[mTORC2_p] * y[Akt_T],
         consume=(Akt_T,), produce=(Akt_TS,))
    flux(P["Akt_pT308_dephos_first"] * y[Akt_T], consume=(Akt_T,), produce=(Akt_00,))
    flux(P["Akt_pS473_dephos_first"] * y[Akt_S], consume=(Akt_S,), produce=(Akt_00,))
    flux(P["Akt_pT308_dephos_second"] * y[Akt_TS], consume=(Akt_TS,), produce=(Akt_S,))
    flux(P["Akt_pS473_dephos_second"] * y[Akt_TS], consume=(Akt_TS,), produce=(Akt_T,))

    # --- TSC1-TSC2: reaction_25, reaction_26, reaction_27, reaction_28
    flux(P["TSC1_TSC2_S1387_phos_by_AMPK_pT172"] * y[AMPK_p] * y[TSC_u],
         consume=(TSC_u,), produce=(TSC_S1387,))
    flux(P["TSC1_TSC2_T1462_phos_by_Akt_pT308"] * akt_pT308_total * y[TSC_u],
         consume=(TSC_u,), produce=(TSC_T1462,))
    flux(P["TSC1_TSC2_pS1387_dephos"] * y[TSC_S1387], consume=(TSC_S1387,), produce=(TSC_u,))
    flux(P["TSC1_TSC2_pT1462_dephos"] * y[TSC_T1462], consume=(TSC_T1462,), produce=(TSC_u,))

    # --- mTORC1: reaction_30, reaction_29
    flux(P["mTORC1_S2448_activation_by_Amino_Acids"] * amino_acids * y[mTORC1_u],
         consume=(mTORC1_u,), produce=(mTORC1_p,))
    flux(P["mTORC1_pS2448_dephos_by_TSC1_TSC2"] * y[mTORC1_p] * tsc_inhibitory,
         consume=(mTORC1_p,), produce=(mTORC1_u,))

    # --- mTORC2: reaction_31, reaction_3, reaction_32
    flux(P["mTORC2_S2481_phos_by_PI3K_variant_p"] * y[PI3Kv_p] * y[mTORC2_u],
         consume=(mTORC2_u,), produce=(mTORC2_p,))
    flux(P["mTORC2_S2481_phos_by_Amino_Acids"] * amino_acids * y[mTORC2_u],
         consume=(mTORC2_u,), produce=(mTORC2_p,))
    flux(P["mTORC2_pS2481_dephos"] * y[mTORC2_p], consume=(mTORC2_p,), produce=(mTORC2_u,))

    # --- p70-S6K: reaction_33..reaction_40
    flux(P["p70_S6K_T229_phos_by_PI3K_p_PDK1_first"] * y[PI3K_p] * y[S6K_00],
         consume=(S6K_00,), produce=(S6K_229,))
    flux(P["p70_S6K_T389_phos_by_mTORC1_pS2448_first"] * y[mTORC1_p] * y[S6K_00],
         consume=(S6K_00,), produce=(S6K_389,))
    flux(P["p70_S6K_T229_phos_by_PI3K_p_PDK1_second"] * y[PI3K_p] * y[S6K_389],
         consume=(S6K_389,), produce=(S6K_both,))
    flux(P["p70_S6K_T389_phos_by_mTORC1_pS2448_second"] * y[mTORC1_p] * y[S6K_229],
         consume=(S6K_229,), produce=(S6K_both,))
    flux(P["p70_S6K_pT229_dephos_first"] * y[S6K_229], consume=(S6K_229,), produce=(S6K_00,))
    flux(P["p70_S6K_pT389_dephos_first"] * y[S6K_389], consume=(S6K_389,), produce=(S6K_00,))
    flux(P["p70_S6K_pT229_dephos_second"] * y[S6K_both], consume=(S6K_both,), produce=(S6K_389,))
    flux(P["p70_S6K_pT389_dephos_second"] * y[S6K_both], consume=(S6K_both,), produce=(S6K_229,))

    # --- PRAS40: reaction_41..reaction_48
    flux(P["PRAS40_S183_phos_by_mTORC1_pS2448_first"] * y[mTORC1_p] * y[PRAS_00],
         consume=(PRAS_00,), produce=(PRAS_183,))
    flux(P["PRAS40_T246_phos_by_Akt_pT308_first"] * akt_pT308_total * y[PRAS_00],
         consume=(PRAS_00,), produce=(PRAS_246,))
    flux(P["PRAS40_T246_phos_by_Akt_pT308_second"] * akt_pT308_total * y[PRAS_183],
         consume=(PRAS_183,), produce=(PRAS_both,))
    flux(P["PRAS40_S183_phos_by_mTORC1_pS2448_second"] * y[mTORC1_p] * y[PRAS_246],
         consume=(PRAS_246,), produce=(PRAS_both,))
    flux(P["PRAS40_pS183_dephos_first"] * y[PRAS_183], consume=(PRAS_183,), produce=(PRAS_00,))
    flux(P["PRAS40_pT246_dephos_first"] * y[PRAS_246], consume=(PRAS_246,), produce=(PRAS_00,))
    flux(P["PRAS40_pS183_dephos_second"] * y[PRAS_both], consume=(PRAS_both,), produce=(PRAS_246,))
    flux(P["PRAS40_pT246_dephos_second"] * y[PRAS_both], consume=(PRAS_both,), produce=(PRAS_183,))

    return d


def observables(y: np.ndarray) -> dict[str, np.ndarray]:
    """The twelve Fig. 2 readouts.

    For a two-site node the paper's readout is the total phosphorylation of one site --
    "Akt-pS473-obs was associated to the sum of the species Akt-pS473 and
    Akt-pT308-pS473" (Methods) -- so each is a sum over the other site's states.
    """
    return {
        "IR_pY1146": y[IR_p],
        "IRS_pS636": y[IRS_pS636],
        "AMPK_pT172": y[AMPK_p],
        "Akt_pT308": y[Akt_T] + y[Akt_TS],
        "Akt_pS473": y[Akt_S] + y[Akt_TS],
        "TSC_pS1387": y[TSC_S1387],
        "mTOR_pS2448": y[mTORC1_p],
        "mTOR_pS2481": y[mTORC2_p],
        "S6K_pT229": y[S6K_229] + y[S6K_both],
        "S6K_pT389": y[S6K_389] + y[S6K_both],
        "PRAS40_pT246": y[PRAS_246] + y[PRAS_both],
        "PRAS40_pS183": y[PRAS_183] + y[PRAS_both],
    }


def integrate(condition="aa", t_eval=None, t_end=120.0, rtol=1e-10, atol=1e-12, **levels):
    """Integrate one condition and return ``(t, {observable: trajectory})``.

    ``condition`` names an entry of ``CONDITIONS``; ``levels`` overrides the node
    activity factors (``pi3k``, ``mtorc1``, ``mtorc2``, ``ampk``) for the scans.
    """
    cfg = dict(CONDITIONS[condition])
    pi3k = levels.pop("pi3k", cfg["pi3k"])
    y0 = initial_state(pi3k=pi3k, **levels)
    if t_eval is None:
        t_eval = np.linspace(0.0, t_end, 1201)
    sol = solve_ivp(
        rhs, (0.0, float(t_eval[-1])), y0,
        args=(cfg["insulin"], cfg["amino_acids"]),
        method="LSODA", t_eval=t_eval, rtol=rtol, atol=atol,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    return sol.t, observables(sol.y)


def steady_state_residual(condition="aa", **levels) -> float:
    """max |dy/dt| at t = 120, in units of the largest pool.

    A trajectory comparison over a finite window can miss a sign error that only shows
    up at the fixed point, so the notebook checks this too.
    """
    cfg = CONDITIONS[condition]
    t, _ = integrate(condition, **levels)
    pi3k = levels.pop("pi3k", cfg["pi3k"])
    y0 = initial_state(pi3k=pi3k, **levels)
    sol = solve_ivp(rhs, (0.0, 120.0), y0, args=(cfg["insulin"], cfg["amino_acids"]),
                    method="LSODA", rtol=1e-12, atol=1e-14)
    return float(np.max(np.abs(rhs(120.0, sol.y[:, -1], cfg["insulin"], cfg["amino_acids"]))))


# --- the mechanical arm ---------------------------------------------------------------

DEFAULT_SBML = (
    Path(__file__).resolve().parents[2]
    / "dev" / "papers" / "DallePezze2016" / "BIOMD0000000640.xml"
)

SBML_MISSING = (
    "BIOMD0000000640.xml not found at {path}.\n"
    "It is the deposited model and is not committed (dev/ is gitignored). Fetch it with:\n"
    "  curl -sL -o {path} \\\n"
    "    'https://www.biomodels.org/model/download/BIOMD0000000640.3"
    "?filename=BIOMD0000000640_url.xml'"
)


def sbml_reference(condition="aa", t_eval=None, path=None, **levels):
    """Integrate the deposited SBML itself, evaluating its MathML mechanically.

    Nothing here is transcribed by hand, so this is the arm that tests whether the
    hand-written system above -- and the ``.bngl`` beside it -- read the deposit
    correctly. Returns ``(t, {observable: trajectory})`` with the same keys as
    :func:`integrate`.
    """
    import libsbml  # imported lazily: only this arm needs it

    p = Path(path) if path is not None else DEFAULT_SBML
    if not p.exists():
        raise FileNotFoundError(SBML_MISSING.format(path=p))

    model = libsbml.readSBML(str(p)).getModel()

    # Inline every functionDefinition, then read each kinetic law as a Python expression.
    fdefs = {}
    for fd in model.getListOfFunctionDefinitions():
        args = [fd.getArgument(i).getName() for i in range(fd.getNumArguments())]
        fdefs[fd.getId()] = (args, libsbml.formulaToL3String(fd.getBody()))

    species = [s.getId() for s in model.getListOfSpecies()]
    dynamic = [s.getId() for s in model.getListOfSpecies() if not s.getBoundaryCondition()]
    idx = {sid: i for i, sid in enumerate(dynamic)}

    consts = {p_.getId(): p_.getValue() for p_ in model.getListOfParameters()}
    consts["Cell"] = model.getCompartment("Cell").getSize()

    def expand(expr: str) -> str:
        """Substitute functionDefinition calls until none remain."""
        import re
        changed = True
        while changed:
            changed = False
            for fid, (args, body) in fdefs.items():
                m = re.search(rf"\b{fid}\(", expr)
                if not m:
                    continue
                # split the call's arguments at top-level commas
                i = m.end()
                depth, start, actual = 1, i, []
                while depth:
                    if expr[i] == "(":
                        depth += 1
                    elif expr[i] == ")":
                        depth -= 1
                        if depth == 0:
                            actual.append(expr[start:i].strip())
                            break
                    elif expr[i] == "," and depth == 1:
                        actual.append(expr[start:i].strip())
                        start = i + 1
                    i += 1
                sub = body
                for a, v in zip(args, actual, strict=False):
                    sub = re.sub(rf"\b{a}\b", f"({v})", sub)
                expr = expr[: m.start()] + f"({sub})" + expr[i + 1 :]
                changed = True
        return expr

    laws = []
    for r in model.getListOfReactions():
        f = expand(libsbml.formulaToL3String(r.getKineticLaw().getMath()))
        laws.append((
            compile(f, "<sbml>", "eval"),
            [s.getSpecies() for s in r.getListOfReactants()],
            [s.getSpecies() for s in r.getListOfProducts()],
        ))

    cfg = dict(CONDITIONS[condition])
    pi3k = levels.pop("pi3k", cfg["pi3k"])
    scale = {"PI3K_variant": pi3k, "PI3K_PDK1": pi3k,
             "mTORC1": levels.get("mtorc1", 1.0), "mTORC2": levels.get("mtorc2", 1.0),
             "AMPK": levels.get("ampk", 1.0)}
    def initial(sid):
        # The deposit sets initialConcentration, not initialAmount, and Cell has size 1,
        # so the two coincide -- but reading the wrong one silently yields all zeros.
        s = model.getSpecies(sid)
        return s.getInitialConcentration() if s.isSetInitialConcentration() \
            else s.getInitialAmount()

    y0 = np.array([initial(sid) * scale.get(sid, 1.0) for sid in dynamic])
    fixed = {"Insulin": cfg["insulin"], "Amino_Acids": cfg["amino_acids"]}
    fixed.update({s: 0.0 for s in species if s.endswith("_obs")})

    def f(_t, y):
        env = dict(consts)
        env.update(fixed)
        env.update({sid: y[i] for i, sid in enumerate(dynamic)})
        d = np.zeros(len(dynamic))
        for code, reactants, products in laws:
            # the expressions come from the deposit, not from user input
            v = eval(code, {"__builtins__": {}}, env)  # noqa: S307
            for s in reactants:
                if s in idx:
                    d[idx[s]] -= v
            for s in products:
                if s in idx:
                    d[idx[s]] += v
        return d

    if t_eval is None:
        t_eval = np.linspace(0.0, 120.0, 1201)
    sol = solve_ivp(f, (0.0, float(t_eval[-1])), y0, method="LSODA",
                    t_eval=t_eval, rtol=1e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(f"SBML integration failed: {sol.message}")

    g = {sid: sol.y[idx[sid]] for sid in dynamic}
    return sol.t, {
        "IR_pY1146": g["IR_beta_pY1146"],
        "IRS_pS636": g["IRS_pS636"],
        "AMPK_pT172": g["AMPK_pT172"],
        "Akt_pT308": g["Akt_pT308"] + g["Akt_pT308_pS473"],
        "Akt_pS473": g["Akt_pS473"] + g["Akt_pT308_pS473"],
        "TSC_pS1387": g["TSC1_TSC2_pS1387"],
        "mTOR_pS2448": g["mTORC1_pS2448"],
        "mTOR_pS2481": g["mTORC2_pS2481"],
        "S6K_pT229": g["p70_S6K_pT229"] + g["p70_S6K_pT229_pT389"],
        "S6K_pT389": g["p70_S6K_pT389"] + g["p70_S6K_pT229_pT389"],
        "PRAS40_pT246": g["PRAS40_pT246"] + g["PRAS40_pT246_pS183"],
        "PRAS40_pS183": g["PRAS40_pS183"] + g["PRAS40_pT246_pS183"],
    }


if __name__ == "__main__":
    for cond in CONDITIONS:
        t, obs = integrate(cond)
        print(f"--- {cond} ---  steady-state residual "
              f"{steady_state_residual(cond):.3e}")
        for name in OBSERVABLES:
            c = obs[name]
            print(f"   {name:14s} peak={c.max():10.4f} @t={t[c.argmax()]:6.2f}"
                  f"   t120={c[-1]:10.4f}")
