"""Shea-Ackers O_R statistical-thermodynamic promoter system for P_R / P_RM.

Regenerates the 40-configuration O_R partition function and the open-complex
activation rates A_PR() / A_PRM() used by lambda_switch_arkin1998_fullcircuit*.bngl,
from the site energies + cooperativities in Table 1/2 of Arkin, Ross & McAdams (1998)
(Shea & Ackers 1985).

Enumeration (matches the committed model exactly):
  three operators OR1,OR2,OR3 each in {O(empty), R(CI2 dimer), C(Cro2 dimer)};
  RNAP may bind P_R (needs OR1,OR2 empty -- occludes them) and/or P_RM (needs OR3
  empty). Loop order OR1 (outer), OR2, OR3 (inner); within one occupancy: base
  (no RNAP), then P_RM-bound (if OR3 empty), then P_R-bound (if OR1,OR2 empty),
  then both (only all-empty). That yields the 40 states in committed order.

  weight energy = sum of [ occupied-site energies (dGR{1,2,3} for CI2, dGC{1,2,3}
  for Cro2) ] + [ RNAP energies dG_Rp (P_R), dG_RMp (P_RM) ] + [ CI2-CI2 adjacency
  cooperativity dG12 (OR1&OR2) , dG23 (OR2&OR3) ]; Cro2 gets no adjacency term.
  weight = exp(-energy/RT) * CI2^(#R) * Cro2^(#C) * RNAP^(#bound).
"""

CI2 = "(Obs_CI2*molec_to_M*dimer_scale)"
CRO2 = "(Obs_Cro2*molec_to_M*dimer_scale)"
RNAP = "(Obs_RNAP*molec_to_M)"


def _configs():
    """Yield the 40 O_R configurations in committed order as dicts."""
    for s1 in "ORC":
        for s2 in "ORC":
            for s3 in "ORC":
                sites = (s1, s2, s3)
                pr_ok = s1 == "O" and s2 == "O"   # RNAP@P_R occludes OR1,OR2
                prm_ok = s3 == "O"                 # RNAP@P_RM occludes OR3
                # within-occupancy order: base, PRM, PR, both
                combos = [(False, False)]
                if prm_ok:
                    combos.append((False, True))
                if pr_ok:
                    combos.append((True, False))
                if pr_ok and prm_ok:
                    combos.append((True, True))
                for pr, prm in combos:
                    yield sites, pr, prm


def _weight_expr(sites, pr, prm):
    s1, s2, s3 = sites
    site_e = {0: ("dGR1", "dGC1"), 1: ("dGR2", "dGC2"), 2: ("dGR3", "dGC3")}
    energy = []
    for i, s in enumerate(sites):
        if s == "R":
            energy.append(site_e[i][0])
        elif s == "C":
            energy.append(site_e[i][1])
    if pr:
        energy.append("dG_Rp")
    if prm:
        energy.append("dG_RMp")
    # Adjacency cooperativity is CI2-CI2 only, and OR2 participates in at most one
    # pair (OR1-OR2 takes precedence): dG12 if OR1,OR2 = CI2; dG23 if OR2,OR3 = CI2
    # AND OR1 != CI2. This matches the committed model (RRR carries dG12 only).
    if s1 == "R" and s2 == "R":
        energy.append("dG12")
    if s2 == "R" and s3 == "R" and s1 != "R":
        energy.append("dG23")
    esum = "+".join(energy) if energy else "0"
    expr = f"exp(-({esum})/RT)"
    nR = sites.count("R")
    nC = sites.count("C")
    nP = int(pr) + int(prm)
    for obs, p in ((CI2, nR), (CRO2, nC), (RNAP, nP)):
        if p == 1:
            expr += f" * {obs}"
        elif p >= 2:
            expr += f" * {obs}^{p}"
    return expr


def generate():
    """Return (lines, A_PR_line, A_PRM_line): the wR_i function lines and the two
    activation-rate function lines, as they appear in the committed .bngl."""
    lines = []
    pr_idx, prm_act_idx, prm_bas_idx, all_idx = [], [], [], []
    for n, (sites, pr, prm) in enumerate(_configs(), start=1):
        all_idx.append(n)
        expr = _weight_expr(sites, pr, prm)
        comment = f"# ({''.join(sites)}, R@PR={int(pr)}, R@PRM={int(prm)})"
        lines.append(f"  wR_{n}() = {expr}  {comment}")
        if pr:
            pr_idx.append(n)
        if prm:
            (prm_act_idx if sites[1] == "R" else prm_bas_idx).append(n)
    Z = " + ".join(f"wR_{i}()" for i in all_idx)
    num_pr = " + ".join(f"wR_{i}()" for i in pr_idx)
    a_pr = f"  A_PR() = kPR * ({num_pr}) / ({Z})"
    act = " + ".join(f"wR_{i}()" for i in prm_act_idx)
    bas = " + ".join(f"wR_{i}()" for i in prm_bas_idx)
    a_prm = f"  A_PRM() = (kPRM1*({act}) + kPRM2*({bas})) / ({Z})"
    return lines, a_pr, a_prm


if __name__ == "__main__":
    lines, a_pr, a_prm = generate()
    for line in lines:
        print(line)
    print(a_pr)
    print(a_prm)
