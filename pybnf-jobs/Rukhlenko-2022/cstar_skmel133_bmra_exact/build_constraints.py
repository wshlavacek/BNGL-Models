#!/usr/bin/env python
"""Generate the PAPER-EXACT BMRA-constraint machinery for the SKMEL-133 job:
  * the model `cstar_skmel133_bmra_exact.bngl` = the SKMEL network with the exact Eq.14
    connection-coefficient functions r_<tgt>_<src>() appended (and the sign-approx cc_S*
    DPD carriers kept), and
  * `bmra_rij.prop` = a TWO-SIDED BPSL constraint per signaling edge pinning the model's
    Eq.14 r_ij inside the BMRA confidence interval [mean-std, mean+std] (Table S10, with_MYC),
    plus a SIGN constraint per DPD driving-force coefficient.

WHAT MAKES THIS "EXACT" (vs. the sibling ../cstar_skmel133_bmra sign approximation):
  The paper (Rukhlenko et al. 2022, Nature; PMC9644236) fit the model with pyBioNetFit under
  the constraint that "the connection coefficients defined in Eq. 14 must be within the
  confidence intervals of the BMRA inferred connections" (Methods p.22). Eq. 13/14 define
      r_ij = d log x_i / d log x_j |st.st = -(df_i/dx_j)/(df_i/dx_i) * (x_j/x_i) ,  r_ii = -1
  the Jacobian-normalized local response of module i to module j. The sign-approx enforces
  only sign(g-1); THIS build encodes the full r_ij as a model function and constrains it
  two-sided to the Table S10 CI -- the real numeric constraint.

THE DECOMPOSITION (derived + triple-verified; see VALIDATION.md): r_ij = C_i * L_ij, where
  L_ij = w(g_ij-1)/((1+w)(1+g_ij w)),  w = x_j/K_ij        [per-edge, x_j the regulator's
          active form -- the same observable the model's alpha_<src><tgt>() uses]
  C_i  = -J_act,i/(x_i * df_i/dx_i)                         [per-target-node self-normalization]
For a standard activation/deactivation node (ERK,AKT,SRC,mTOR):
  J_act = J_deact = vn_i x_i/(Kn_i+x_i) at st.st;  df_i/dx_i = -J_act Kp_i/(nX(Kp_i+nX)) -
  vn_i Kn_i/(Kn_i+x_i)^2,  nX = tot_i - x_i  =>  C_i = J_act/(x_i*[J_act Kp/(nX(Kp+nX)) +
  vn Kn/(Kn+x_i)^2]).
For the CDK module (regulators enter via CycD synthesis; internal CycD relaxes):
  pCDK = CDK_total * fCycD/(KCycCDK+fCycD),  fCycD_ss = vSynCycD*P/vDegCycD  =>  C_CDK =
  KCycCDK/(KCycCDK+fCycD)  (uses the model observable fCycD).
For MYC (constant synthesis - vDegMYC*alpha_CDKMYC*MYCt): C_MYC = -1  => r_MYC_CDK = -L.

NOT constrained here (documented in VALIDATION.md):
  * IRS-as-target (mTOR->IRS, PKC->IRS via inhibitory phosphorylation; ERK->IRS is OFF in the
    model, g_ERKIRS=1): a 4-species module with no clean closed-form Eq.14 r; verified
    numerically (r_IRS,PKC=-0.75 in band; r_IRS,mTOR=-1.28 just outside; both correct sign).
  * DPD row magnitude: Eq. 25 beta_j = r_Sj*(S_ss/x_j_ss). With the model's PHYSICAL S_ss(~10)
    and x_j_ss(~0.4), r_Sj_model ~ 1e-3, vs BMRA r_Sj ~ O(1): the S-vs-x unit ratio differs
    from BMRA's standardized space (a genuine unit-incommensurability, not a fit error), so an
    exact r_Sj magnitude constraint is ill-posed. The DPD is constrained by SIGN (cc_S* carriers).

BMRA source: Table S10 with_MYC (the model matches the withMyc network); values EQUAL the
posterior mean +/- 1 std (k=1) embedded below.
"""
import os

# ---- Table S10 (with_MYC) BMRA r: (mean, std); rows=target, cols=source ----
BMRA = {
 ('ERK','mTOR'):(0.38,0.56),('ERK','SRC'):(0.38,0.56),('ERK','CDK'):(-0.13,0.37),('ERK','IRS'):(0.49,0.56),
 ('AKT','ERK'):(-0.75,0.46),('AKT','CDK'):(0.41,0.46),('AKT','PKC'):(0.97,0.45),('AKT','IRS'):(1.34,0.46),('AKT','MYC'):(0.57,0.46),
 ('mTOR','AKT'):(0.31,0.51),('mTOR','SRC'):(0.18,0.47),('mTOR','CDK'):(0.99,0.51),('mTOR','PKC'):(0.34,0.51),
 ('SRC','ERK'):(0.33,0.72),('SRC','PKC'):(0.40,0.73),('SRC','MYC'):(-0.56,0.70),
 ('CDK','ERK'):(0.29,0.58),('CDK','mTOR'):(0.18,0.58),('CDK','MYC'):(0.67,0.58),
 ('MYC','CDK'):(0.41,0.11),
}
# BMRA DPD (S) row -- for the SIGN constraints on the driving-force coefficients (Eq. 25)
BMRA_DPD = {'mTOR':(1.62,0.46,'+'),'PKC':(0.64,0.46,'+'),'SRC':(-0.27,0.47,'-')}

# regulator active-form observable (== the obs the model's alpha_<src><tgt>() uses)
XA = dict(IRS='pIRS', mTOR='pmTOR', SRC='pSRC', CDK='pCDK', ERK='pERK', PKC='pPKC', MYC='MYCt', AKT='pAKT')
GK = {
 ('ERK','IRS'):('g_IRSERK','K_IRSERK'),('ERK','mTOR'):('g_mTORERK','K_mTORERK'),
 ('ERK','SRC'):('g_SRCERK','K_SRCERK'),('ERK','CDK'):('g_CDKERK','K_CDKERK'),
 ('AKT','ERK'):('g_ERKAKT','K_ERKAKT'),('AKT','PKC'):('g_PKCAKT','K_PKCAKT'),
 ('AKT','CDK'):('g_CDKAKT','K_CDKAKT'),('AKT','IRS'):('g_IRSAKT','K_IRSAKT'),('AKT','MYC'):('g_MYCAKT','K_MYCAKT'),
 ('SRC','ERK'):('g_ERKSRC','K_ERKSRC'),('SRC','PKC'):('g_PKCSRC','K_PKCSRC'),('SRC','MYC'):('g_MYCSRC','K_MYCSRC'),
 ('CDK','ERK'):('g_ERKCDK','K_ERKCDK'),('CDK','mTOR'):('g_mTORCDK','K_mTORCDK'),('CDK','MYC'):('g_MYCCDK','K_MYCCDK'),
 ('mTOR','AKT'):('g_AKTmTOR','K_AKTmTOR'),('mTOR','SRC'):('g_SRCmTOR','K_SRCmTOR'),
 ('mTOR','PKC'):('g_PKCmTOR','K_PKCmTOR'),('mTOR','CDK'):('g_CDKmTOR','K_CDKmTOR'),
 ('MYC','CDK'):('g_CDKMYC','K_CDKMYC'),
}
# standard nodes: target -> (x_i obs, total, Kp, vn, Kn)
STD = {
 'ERK':('pERK','ERK_total','KpERK','vnERK','KnERK'),
 'AKT':('pAKT','AKT_total','KpAKT','vnAKT','KnAKT'),
 'SRC':('pSRC','SRC_total','KpSRC','vnSRC','KnSRC'),
 'mTOR':('pmTOR','mTOR_total','KpmTOR','vnmTOR','KnmTOR'),
}
ORDER = ['ERK','AKT','mTOR','SRC','CDK','MYC']   # target grouping for readable output


def L_expr(g, K, xa):
    w = f"({xa}/{K})"
    return f"{w}*({g}-1.0)/((1.0+{w})*(1.0+{g}*{w}))"


def emit_rij_funcs():
    out = ["# --- Eq.14 connection coefficients r_ij = C_i*L_ij (build_constraints.py) --------",
           "# r_ij = d log x_i/d log x_j|st.st (Rukhlenko 2022 Eq.13/14). L_ij = w(g-1)/((1+w)",
           "# (1+g w)), w=x_j/K_ij; C_i self-normalization per node. bmra_rij.prop constrains each",
           "# r_<tgt>_<src> two-sided to the Table S10 (with_MYC) CI. Verified analytic == numeric",
           "# operational MRA to 5e-9 (VALIDATION.md)."]
    for tgt in STD:
        x, tot, Kp, vn, Kn = STD[tgt]
        nX = f"({tot}-{x})"
        Jd = f"({vn}*{x}/({Kn}+{x}))"
        sl = f"({Jd}*{Kp}/({nX}*({Kp}+{nX})) + {vn}*{Kn}/(({Kn}+{x})*({Kn}+{x})))"
        # C_i = -J_act/(x_i*df_i/dx_i) = Jd/(x_i*sl). Since Jd = vn*x_i/(Kn+x_i), the x_i cancels:
        # C_i = vn/((Kn+x_i)*sl) -- algebraically identical but FINITE at x_i=0 (=> C_i=1), so the
        # constraint experiment (a time course from the all-off seed) has no 0/0 at t=0.
        out += [f"sl_{tgt}() = {sl}",
                f"C_{tgt}() = {vn}/(({Kn}+{x})*sl_{tgt}())"]
    out.append("C_CDK() = KCycCDK/(KCycCDK+fCycD)")
    for tgt in ORDER:
        for (t, s), (g, K) in GK.items():
            if t != tgt:
                continue
            Lij = L_expr(g, K, XA[s])
            if tgt == 'MYC':
                out.append(f"r_{t}_{s}() = -1.0*({Lij})")
            else:
                out.append(f"r_{t}_{s}() = C_{t}()*({Lij})")
    return "\n".join(out) + "\n"


def emit_prop():
    p = ["# PAPER-EXACT BMRA-CI constraints on the SKMEL-133 connection coefficients (Table S10,",
         "# with_MYC). Each signaling edge: the model's Eq.14 r_ij (function r_<tgt>_<src>) is pinned",
         "# TWO-SIDED inside [mean-std, mean+std] (k=1; the table +/- are the BMRA posterior stds).",
         "# Evaluated at the settled no-drug steady state (the constraint experiment runs `at 500000`).",
         "# weight = 1/std (tighter CI -> stronger). Names are bare (BPSL grammar forbids parentheses).",
         ""]
    for tgt in ORDER:
        for (t, s), (g, K) in GK.items():
            if t != tgt:
                continue
            m, sd = BMRA[(t, s)]
            lo, hi = m - sd, m + sd
            w = round(1.0 / sd, 3)
            p.append(f"# {t}<-{s}: BMRA r={m:+.3f}+/-{sd:.3f}  ->  r in [{lo:+.3f},{hi:+.3f}]")
            p.append(f"r_{t}_{s} > {lo:.4f} at 500000 weight {w}")
            p.append(f"r_{t}_{s} < {hi:.4f} at 500000 weight {w}")
    p += ["",
          "# DPD driving-force coefficients beta_j (Eq. 25): the model's PHYSICAL S/x units are not",
          "# commensurate with BMRA's standardized space (r_Sj_model ~ 1e-3 vs BMRA ~ O(1); see",
          "# VALIDATION.md), so the DPD row is constrained by SIGN only. force = +beta_mTOR*pmTOR",
          "# +beta_PKC*pPKC -beta_SRC*pSRC; cc_S<x>()=beta_x*totERK (totERK==1 conserved carrier)."]
    for s, (m, sd, sgn) in BMRA_DPD.items():
        w = round(abs(m) / sd, 3)
        p.append(f"# DPD<-{s}: BMRA r={m:+.3f}+/-{sd:.3f} (force coeff {sgn}beta_{s})")
        p.append(f"cc_S{s} > 0 at 500000 weight {w}")
    p += ["",
          "# PROLIFERATIVE-STATE ANCHOR. The BMRA connection coefficients were inferred in the",
          "# PROLIFERATIVE state (SKMEL-133 is a proliferative BRAF-V600E melanoma line; the cSTAR",
          "# DPD places proliferation at S > 0). The exact Eq.14 r_ij are the steady-state response",
          "# ratios AT that attractor, so the constraint experiment must stay in the proliferative",
          "# basin. Without this anchor a free (g,K) fit can drift the drug-free attractor into the",
          "# DIFFERENTIATED basin (S < 0), satisfying the numeric r_ij bands but in the wrong",
          "# physiological regime (observed: an unanchored bounded fit reached S = -13, pmTOR = 0.055",
          "# vs the published proliferative S = +9.95, pmTOR = 0.409). Sval > 0 keeps the fit on the",
          "# proliferative side of the separatrix; weight 1.0 so, scaled, it dominates any incentive",
          "# to differentiate (the published params satisfy it with wide margin).",
          "Sval > 0 at 500000 weight 1.0"]
    return "\n".join(p) + "\n"


def build_model():
    """cstar_skmel133_bmra_exact.bngl = sign-approx model with the 20 cc_<edge> sign carriers
    replaced by the exact r_<tgt>_<src>() functions (cc_S* DPD carriers + totERK kept)."""
    base = open("../cstar_skmel133_bmra/cstar_skmel133_bmra.bngl").read()
    lines = base.splitlines()
    out, skip = [], False
    for ln in lines:
        if ln.startswith("# --- BMRA-CI connection-sign carriers"):
            skip = True                     # start of the cc_<edge> block
        if skip and ln.startswith("cc_SmTOR()"):
            skip = False                    # keep cc_S* DPD carriers onward
            out.append("# --- DPD driving-force sign carriers (kept; see bmra_rij.prop) ----------------")
        if not skip:
            out.append(ln)
    body = "\n".join(out)
    # insert the r_ij functions just before `end functions`
    body = body.replace("end functions", emit_rij_funcs() + "end functions")
    # retitle the banner
    body = body.replace("cstar_skmel133_bmra.bngl -- cSTAR SKMEL-133 network, BMRA-CI-CONSTRAINED full fit.",
                        "cstar_skmel133_bmra_exact.bngl -- cSTAR SKMEL-133, PAPER-EXACT BMRA r_ij constraints.")
    return body + ("" if body.endswith("\n") else "\n")


def main():
    open("cstar_skmel133_bmra_exact.bngl", "w").write(build_model())
    open("bmra_rij.prop", "w").write(emit_prop())
    ncon = sum(1 for t in ORDER for (a, b) in GK if a == t)
    print(f"wrote cstar_skmel133_bmra_exact.bngl")
    print(f"wrote bmra_rij.prop ({2*ncon} two-sided signaling + {len(BMRA_DPD)} DPD-sign constraints)")


if __name__ == "__main__":
    main()
