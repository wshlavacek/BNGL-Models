"""Shea-Ackers activation of P_RE (CII-activated) and repression of P_L (Cro2/CI2)
for lambda_switch_arkin1998_fullcircuit*.bngl, from Table 1 of Arkin, Ross &
McAdams (1998).

These two promoters are small, fixed partition functions (they do not change for
the deferred elongation/translation refinements), so they are emitted verbatim with
the Table-1 parameter names. P_RE: empty / RNAP-only (basal, koc_RE_2) / CII-only /
CII+RNAP (activated, koc_RE_4). P_L: RNAP binds (koc_L_6) and is repressed by Cro2
(OL1,OL2) and CI2, incl. pairwise terms.
"""

CI2 = "(Obs_CI2*molec_to_M*dimer_scale)"
CRO2 = "(Obs_Cro2*molec_to_M*dimer_scale)"
RNAP = "(Obs_RNAP*molec_to_M)"
CII = "(Obs_CII*molec_to_M)"


def generate_pre():
    """Return (wRE lines, A_PRE line)."""
    lines = [
        f"  wRE_2() = exp(-dG_RE_2/RT) * {RNAP}",
        f"  wRE_3() = exp(-dG_RE_3/RT) * {CII}",
        f"  wRE_4() = exp(-dG_RE_4/RT) * {CII}*{RNAP}",
    ]
    a_pre = ("  A_PRE() = (wRE_2()*koc_RE_2 + wRE_4()*koc_RE_4) / "
             "(1 + wRE_2() + wRE_3() + wRE_4())")
    return lines, a_pre


def generate_pl():
    """Return (wL lines, A_PL line)."""
    lines = [
        f"  wL_2() = exp(-dG_L_2/RT) * {CRO2}",
        f"  wL_3() = exp(-dG_L_3/RT) * {CRO2}",
        f"  wL_4() = exp(-dG_L_4/RT) * {CI2}",
        f"  wL_5() = exp(-dG_L_5/RT) * {CI2}",
        f"  wL_6() = exp(-dG_L_6/RT) * {RNAP}",
        f"  wL_7() = exp(-dG_L_7/RT) * {CRO2}^2",
        f"  wL_8() = exp(-dG_L_8/RT) * {CRO2}*{CI2}",
        f"  wL_9() = exp(-dG_L_9/RT) * {CI2}*{CRO2}",
        f"  wL_10() = exp(-dG_L_10/RT) * {CI2}^2",
    ]
    a_pl = ("  A_PL() = wL_6()*koc_L_6 / "
            "(1 + wL_2() + wL_3() + wL_4() + wL_5() + wL_6() + wL_7() + wL_8() + wL_9() + wL_10())")
    return lines, a_pl


if __name__ == "__main__":
    for blk in (generate_pre(), generate_pl()):
        for line in blk[0]:
            print(line)
        print(blk[1])
