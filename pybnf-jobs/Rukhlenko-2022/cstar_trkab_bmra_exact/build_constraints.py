#!/usr/bin/env python
"""Generate the PAPER-EXACT BMRA-constraint machinery for the JOINT TrkA + TrkB models:
  * each receptor model `cstar_trk{a,b}_bmra_exact.bngl` = the Trk network with the exact Eq.14
    connection-coefficient functions r_<tgt>_<src>() appended (replacing the sign-approx cc_*
    carriers; the conserved totERK observable is kept), and
  * `bmra_{A,B}_rij.prop` = a TWO-SIDED BPSL constraint per constrained edge pinning the model's
    Eq.14 r_ij inside the BMRA confidence interval [mean-std, mean+std] (Table S5, 10-min posterior).

GROUNDING (paper Methods p.22, "Refining parameters of the dynamic model"): the Trk models were
fit to the TRAINING set = the TrkA/TrkB phosphorylation time courses (Western blot) + the **10-min**
RPPA data, UNDER the constraint that "the connection coefficients defined in Eq. 14 must be within
the confidence intervals of the BMRA inferred connections". The **45-min** RPPA is the VALIDATION
set. So the exact r_ij are constrained to the **10-min** BMRA posterior (Table S5). Eq.14 defines
r_ij at STEADY STATE; the Trk ligand response is transient/adaptive (peaks ~10 min, relaxes back),
so the well-defined reference steady state is the BASAL (Lig_on = 0) state -- a true st.st. and the
pre-perturbation reference Eq.25 uses. The constraint experiment therefore evaluates r_ij at the
no-ligand steady state (as the sign-approx slug does).

THE DECOMPOSITION (derived + verified to 1e-10 vs an operational MRA; see VALIDATION.md):
  r_ij = C_i * L_ij,   L_ij = w(g_ij-1)/((1+w)(1+g_ij w)),  w = x_j/K_ij   [x_j = regulator active form]
STANDARD activation/deactivation node i (ERK,AKT,JNK,S6K -- same MM push-pull form as SKMEL):
  C_i = vn_i/((Kn_i+x_i)*sl_i),  sl_i = Jd*Kp_i/(nX(Kp_i+nX)) + vn_i*Kn_i/(Kn_i+x_i)^2,
  Jd = vn_i*x_i/(Kn_i+x_i),  nX = tot_i - x_i.   (x_i = the node active form; ERK uses ppERK.)
RSK node (basal + ERK-linear synthesis, NOT a hyperbolic alpha; only RSK<-ERK): a custom r derived
  from f_RSK = (vpRSKbasal + kpRSK*ppERK/KpRSK)*Kp*nRSK/(Kp+nRSK) - vnRSK*pRSK/(KnRSK+pRSK):
  r_RSK_ERK = [kpRSK*nRSK/(Kp+nRSK)] / [B*Kp^2/(Kp+nRSK)^2 + vnRSK*KnRSK/(KnRSK+pRSK)^2] * ppERK/pRSK,
  B = vpRSKbasal + kpRSK*ppERK/KpRSK, nRSK = RSK_total - pRSK, Kp = KpRSK{A,B} (tag-specific).

NOT constrained (documented in VALIDATION.md): the receptor-dimerization node RTK/AdRTK (dim0/1/2
states, ligand binding, dimer degradation) has no clean closed-form Eq.14 C_i, so the RTK-incoming
edges (RTK<-S6K for A; RTK<-{ERK,RSK} for B) are reported numerically but not constrained -- the same
honest scope used for the SKMEL IRS node. The DPD (STV) row is NOT constrained (SVM-defined, and the
phospho time-course training does not involve Sval; as in the sign-approx slug).

BMRA source: Table S5 10-min posterior (Trk{A,B}_{rm,rs}_10_log_200_5K.csv), embedded below;
proteins order TRK,ERK,AKT,JNK,S6K,RSK,RTK,STV. Values = posterior mean +/- 1 std (k=1).
"""
import re

PROT = ["TRK", "ERK", "AKT", "JNK", "S6K", "RSK", "RTK", "STV"]
STD_NODES = {"ERK", "AKT", "JNK", "S6K"}
MODELS = {"A": "cstar_trka_bmra_exact.bngl", "B": "cstar_trkb_bmra_exact.bngl"}
SRCMODELS = {"A": "../cstar_trkab_bmra/cstar_trka_bmra.bngl", "B": "../cstar_trkab_bmra/cstar_trkb_bmra.bngl"}

# 10-min BMRA posteriors (rows=target, cols=source; order PROT). From the sign-approx slug.
RM = {
"A": [[0,0,0,0,0,0,0,0],
 [1.24351247889868,0,-0.509948152641893,0.93767128409788,-0.461539076696644,0,0.418133198240322,0],
 [-0.511023977984669,1.02223875523467,0,0,-0.23838387617812,0.594807202238467,0,0],
 [0,0,0,0,0,0,0,0],
 [0,0,0.531253008538487,0,0,0,0,0],
 [0,0,0,0,0,0,0,0],
 [0,0,0,0,0.279093531589222,0,0,0],
 [0,0,0.221217487608546,0,0.179367430730783,0,0,0]],
"B": [[0,0,0,0,0,0,0,0],
 [0.0191246120969057,0,0.00103749147808791,0.343008640111969,0,-0.0745063054295733,0.324570294174492,0],
 [0.243570732004696,0.551754058475933,0,-0.0692282927249645,-0.179690967907822,-0.0538342505751875,0.532130672249929,0],
 [0,0.869415924473424,0,0,0,0,-0.0193409085086149,0],
 [0,0,0.489658219194708,0,0,0,0,0],
 [0,0.753117325242822,0,0,0,0,0,0],
 [0.331502512337904,0.616148972485565,0,0,0,0.788269987784397,0,0],
 [0,-0.681908632223304,0.28949814520325,-0.0384204022620735,0.411426250912263,0.467632044422297,0,0]],
}
RS = {
"A": [[0,0,0,0,0,0,0,0],
 [0.53752449251516,0,0.53922737544398,0.542942246327137,0.54036960812057,0,0.549173304557561,0],
 [0.705668582661992,0.698958984231271,0,0,0.718774867110515,0.714307341801992,0,0],
 [0,0,0,0,0,0,0,0],
 [0,0,0.0340740628951638,0,0,0,0,0],
 [0,0,0,0,0,0,0,0],
 [0,0,0,0,0.165099263635929,0,0,0],
 [0,0,0.543270249456481,0,0.552002116815507,0,0,0]],
"B": [[0,0,0,0,0,0,0,0],
 [0.595838025255946,0,0.0664927393800005,0.591991081247234,0,0.607300717335327,0.640289266488317,0],
 [0.530885727173588,0.518234738518658,0,0.528102092502302,0.528678636761685,0.53004511902799,0.524745755977587,0],
 [0,0.283970297888799,0,0,0,0,0.279348439924267,0],
 [0,0,0.0329040075891381,0,0,0,0,0],
 [0,0.0191142364819504,0,0,0,0,0,0],
 [0.539795511059706,0.536473001531127,0,0,0,0.542388751417443,0,0],
 [0,0.441171623145094,0.439902039426225,0.44231215193596,0.446397537695365,0.439246342004745,0,0]],
}

# standard node -> (x_i obs, total, Kp, vn, Kn) -- param names, same in A and B
STD = {
 "ERK": ("ppERK", "ERK_total", "KpERK", "vnERK", "KnERK"),
 "AKT": ("pAKT",  "AKT_total", "KpAKT", "vnAKT", "KnAKT"),
 "JNK": ("pJNKt", "JNK_total", "KpJNK", "vnJNK", "KnJNK"),
 "S6K": ("pS6Kt", "S6K_total", "KpS6K", "vnS6K", "KnS6K"),
}
# source active-form observable the model's alpha_<src><tgt>() uses as x_j
XA = {"TRK": "pTRK", "RTK": "pAdRTK", "ERK": "ppERK", "AKT": "pAKT",
      "JNK": "pJNK", "S6K": "pS6K", "RSK": "pRSK"}
# RSK rule params per model (kpRSK, KpRSK differ by tag; the rest shared)
RSK_KP = {"A": ("kpRSKA", "KpRSKA"), "B": ("kpRSKB", "KpRSKB")}
ORDER = ["ERK", "AKT", "JNK", "S6K", "RSK", "RTK"]  # target grouping for readable output


def split_edge(edge):
    for i in range(2, len(edge) - 1):
        a, b = edge[:i], edge[i:]
        if a in PROT and b in PROT:
            return a, b            # (src, tgt)
    return None


def parse_alpha_edges(tag):
    """(src, tgt, g_id, K_id) for every alpha_<tag>_<src><tgt>() crosstalk multiplier."""
    s = open(SRCMODELS[tag]).read()
    out = []
    for m in re.finditer(rf"alpha_{tag}_(\w+)\(\)\s*=\s*\(1\.0 \+ (g_{tag}_\w+)\*", s):
        edge, g = m.groups()
        sp = split_edge(edge)
        if not sp:
            continue
        src, tgt = sp
        K = g.replace("g_", "K_", 1)
        out.append((src, tgt, g, K))
    return out


def all_edges(tag):
    """Every model signalling edge (alpha-based + the RSK<-ERK rule edge), as (src,tgt,g,K,kind)."""
    edges = [(s, t, g, K, "std" if t in STD_NODES else ("rtk" if t == "RTK" else "other"))
             for (s, t, g, K) in parse_alpha_edges(tag)]
    edges.append(("ERK", "RSK", None, None, "rsk"))     # RSK<-ERK: from the RSK rule, not an alpha
    return edges


def L_expr(g, K, xa):
    w = f"({xa}/{K})"
    return f"{w}*({g}-1.0)/((1.0+{w})*(1.0+{g}*{w}))"


def rsk_expr(tag):
    kp, Kp = RSK_KP[tag]
    nR = "(RSK_total-pRSK)"
    B = f"(vpRSKbasal + {kp}*ppERK/{Kp})"
    num = f"({kp}*{nR}/({Kp}+{nR}))"
    den = f"({B}*{Kp}*{Kp}/(({Kp}+{nR})*({Kp}+{nR})) + vnRSK*KnRSK/((KnRSK+pRSK)*(KnRSK+pRSK)))"
    return f"r_RSK_ERK() = ({num}/{den})*(ppERK/pRSK)"


def emit_rij_funcs(tag):
    out = ["# --- Eq.14 connection coefficients r_ij = C_i*L_ij (build_constraints.py) --------",
           "# r_ij = d log x_i/d log x_j|st.st (Rukhlenko 2022 Eq.13/14), evaluated at the BASAL",
           "# (Lig_on=0) steady state. Standard node C_i = vn/((Kn+x_i)*sl_i); RSK is a custom node.",
           "# bmra_%s_rij.prop constrains each r_<tgt>_<src> two-sided to the Table S5 (10-min) CI." % tag,
           "# Verified analytic == operational MRA to 1e-10 (VALIDATION.md)."]
    # standard-node self-normalizations (only for nodes that are a constrained target)
    edges = all_edges(tag)
    tgt_nodes = {t for (s, t, g, K, kind) in edges if kind == "std" and RM[tag][PROT.index(t)][PROT.index(s)] != 0}
    for tgt in [n for n in ["ERK", "AKT", "JNK", "S6K"] if n in tgt_nodes]:
        x, tot, Kp, vn, Kn = STD[tgt]
        nX = f"({tot}-{x})"
        Jd = f"({vn}*{x}/({Kn}+{x}))"
        sl = f"({Jd}*{Kp}/({nX}*({Kp}+{nX})) + {vn}*{Kn}/(({Kn}+{x})*({Kn}+{x})))"
        out += [f"sl_{tgt}() = {sl}", f"C_{tgt}() = {vn}/(({Kn}+{x})*sl_{tgt}())"]
    # r_ij functions, only for constrained edges (std targets + RSK)
    for tgt in ORDER:
        for (s, t, g, K, kind) in edges:
            if t != tgt:
                continue
            r = RM[tag][PROT.index(t)][PROT.index(s)]
            if kind == "std" and r != 0:
                out.append(f"r_{t}_{s}() = C_{t}()*({L_expr(g, K, XA[s])})")
            elif kind == "rsk" and r != 0:
                out.append(rsk_expr(tag))
    return "\n".join(out) + "\n"


def emit_prop(tag):
    rm, rs = RM[tag], RS[tag]
    p = [f"# PAPER-EXACT BMRA-CI constraints on the Trk{tag} connection coefficients (Table S5, 10-min",
         "# posterior). Each edge: the model's Eq.14 r_ij (function r_<tgt>_<src>) is pinned TWO-SIDED",
         "# inside [mean-std, mean+std] (k=1). Evaluated at the no-ligand (basal) steady state (the",
         "# constraint experiment runs `at 100000`). weight = 1/std. Names are bare (BPSL forbids `()`).",
         ""]
    nconstr = 0
    doc = []
    for tgt in ORDER:
        for (s, t, g, K, kind) in all_edges(tag):
            if t != tgt:
                continue
            r, sd = rm[PROT.index(t)][PROT.index(s)], rs[PROT.index(t)][PROT.index(s)]
            if r == 0 or sd == 0:
                continue
            if kind in ("std", "rsk"):
                lo, hi = r - sd, r + sd
                w = round(1.0 / sd, 3)
                p.append(f"# {t}<-{s}: BMRA r={r:+.3f}+/-{sd:.3f}  ->  r in [{lo:+.3f},{hi:+.3f}]")
                p.append(f"r_{t}_{s} > {lo:.4f} at 100000 weight {w}")
                p.append(f"r_{t}_{s} < {hi:.4f} at 100000 weight {w}")
                nconstr += 1
            elif kind == "rtk":
                doc.append(f"#   RTK<-{s}: BMRA r={r:+.3f}+/-{sd:.3f} (receptor-dimerization node; "
                           f"no closed-form Eq.14 C_i -> reported numerically, NOT constrained; see VALIDATION.md)")
    if doc:
        p += ["", "# DOCUMENTED, NOT CONSTRAINED (the RTK/AdRTK receptor-dimerization node has no clean",
              "# closed-form Eq.14 self-normalization; verified numerically in VALIDATION.md):"] + doc
    return "\n".join(p) + "\n", nconstr


def build_model(tag):
    """cstar_trk{a,b}_bmra_exact.bngl = the sign-approx model with the cc_<tag>_<edge> sign carriers
    replaced by the exact r_<tgt>_<src>() functions (totERK kept)."""
    base = open(SRCMODELS[tag]).read().splitlines()
    out, skip = [], False
    for ln in base:
        if ln.startswith("# --- BMRA-CI connection-sign carriers"):
            skip = True
            out.append("# --- exact Eq.14 connection coefficients (build_constraints.py) ---")
            continue
        if skip and ln.strip() == "end functions":
            skip = False
        if not skip:
            out.append(ln)
    body = "\n".join(out).replace("end functions", emit_rij_funcs(tag) + "end functions")
    body = body.replace("cstar_trk%s_bmra.bngl" % tag.lower(), "cstar_trk%s_bmra_exact.bngl" % tag.lower())
    return body + ("" if body.endswith("\n") else "\n")


def main():
    for tag in ("A", "B"):
        open(MODELS[tag], "w").write(build_model(tag))
        prop, n = emit_prop(tag)
        open(f"bmra_{tag}_rij.prop", "w").write(prop)
        print(f"Trk{tag}: wrote {MODELS[tag]} + bmra_{tag}_rij.prop ({n} two-sided edges = {2*n} constraints)")


if __name__ == "__main__":
    main()
