#!/usr/bin/env python
"""Generate the BMRA-CI sign-constraint machinery for the SKMEL-133 BMRA job:
  * the `cc_*()` observable-functions appended to the model (one per constrained
    connection), and
  * the `bmra_signs.prop` BPSL file (one inequality per connection).

WHY sign constraints (grounded in Rukhlenko et al. 2022, Nature; PMC9644236):
  - Eq. 24: the crosstalk multiplier alpha_Y^X = (1 + g*Ya/K)/(1 + Ya/K); the paper states
    g (their gamma_Y^X) > 1 => activation, < 1 => inhibition, = 1 => no interaction.
  - Methods p.22/p.24: the model was fit with pyBioNetFit "constrain[ing] the parameters
    using the BMRA inferred connection coefficients within their confidence intervals ...
    the connection coefficients defined in Eq. 14 must be within the confidence intervals
    of the BMRA inferred connections" -> implicit inequality constraints on the params.
  - Eq. 14 r_ij is a Jacobian-normalized ratio (r_ii = -1); PyBNF BPSL constrains model
    OBSERVABLES, not Jacobian elements, so we enforce the robust, CI-supported content of
    each connection: its SIGN (which side of 1 the connection strength g sits on), weighted
    by the BMRA confidence z = |mean|/std of the inferred coefficient. Verified: at the
    published parameters the model's g-signs match all 20 BMRA-inferred signs exactly.

Encoding: BPSL needs an observable, and a parameter-only function is not emitted to the
gdat, so each constraint rides `cc_<edge>() = (g_<edge> - 1) * totERK`, where totERK is a
conserved (constant = 1) observable. sign(cc) = sign(g-1); the constraint pins that sign.
For a connection acting on a DEGRADATION term (CDK->MYC) the effect on the target LEVEL is
inverted, so the g-side is flipped relative to the BMRA r sign.

BMRA source: results_SKMEL_133/SKMEL133_{rm,rs}_log_200_5K_withMyc.csv (the model matches
the withMyc network). rows=target, cols=source, order ERK,AKT,mTOR,SRC,CDK,PKC,IRS,MYC,DPD.
"""

# edge: (g_id, active_form_obs, term 'prod'|'deg', source, target)
EDGES = [
    ("g_IRSERK",  "pIRS",  "prod", "IRS",  "ERK"),
    ("g_mTORERK", "pmTOR", "prod", "mTOR", "ERK"),
    ("g_SRCERK",  "pSRC",  "prod", "SRC",  "ERK"),
    ("g_CDKERK",  "pCDK",  "prod", "CDK",  "ERK"),
    ("g_ERKAKT",  "pERK",  "prod", "ERK",  "AKT"),
    ("g_PKCAKT",  "pPKCt", "prod", "PKC",  "AKT"),
    ("g_CDKAKT",  "pCDK",  "prod", "CDK",  "AKT"),
    ("g_IRSAKT",  "pIRS",  "prod", "IRS",  "AKT"),
    ("g_MYCAKT",  "MYCt",  "prod", "MYC",  "AKT"),
    ("g_ERKSRC",  "pERK",  "prod", "ERK",  "SRC"),
    ("g_PKCSRC",  "pPKCt", "prod", "PKC",  "SRC"),
    ("g_MYCSRC",  "MYCt",  "prod", "MYC",  "SRC"),
    ("g_ERKCDK",  "pERK",  "prod", "ERK",  "CDK"),
    ("g_mTORCDK", "pmTOR", "prod", "mTOR", "CDK"),
    ("g_MYCCDK",  "MYCt",  "prod", "MYC",  "CDK"),
    ("g_AKTmTOR", "pAKT",  "prod", "AKT",  "mTOR"),
    ("g_SRCmTOR", "pSRC",  "prod", "SRC",  "mTOR"),
    ("g_PKCmTOR", "pPKCt", "prod", "PKC",  "mTOR"),
    ("g_CDKmTOR", "pCDK",  "prod", "CDK",  "mTOR"),
    ("g_CDKMYC",  "pCDK",  "deg",  "CDK",  "MYC"),
]

# DPD driving-force coefficients (Eq. 25): signed beta_j; model force = +mTOR +PKC -SRC.
DPD = [  # (beta_id, structural_sign_in_force, source)
    ("beta_mTOR", +1, "mTOR"),
    ("beta_PKC",  +1, "PKC"),
    ("beta_SRC",  -1, "SRC"),
]

PROT = ["ERK", "AKT", "mTOR", "SRC", "CDK", "PKC", "IRS", "MYC", "DPD"]

# BMRA withMyc posterior mean (rm) and std (rs); rows=target, cols=source.
RM = [
 [0,0,0.377367905825401,0.381046697057911,-0.141872762077825,0,0.494427427150514,0,0],
 [-0.738336360505208,0,0,0,0.414682509802654,0.96416992475767,1.32169963486109,0.562585788137663,0],
 [0,0.296546571673435,0,0.186662298048173,0.992655626671292,0.338466710726726,0,0,0],
 [0.340544205986375,0,0,0,0,0.423352432492264,0,-0.536048729948052,0],
 [0.281176867880822,0,0.182421633813682,0,0,0,0,0.672225906342856,0],
 [0,0,0,0,0,0,0,0,0],
 [0.494514729280279,0,-0.55330951824905,0,0,-0.47740759804276,0,0,0],
 [0,0,0,0,0.408461966503752,0,0,0,0],
 [0,0,1.62769665796758,-0.27202980000372,0,0.646820970217049,0,0,0]]
RS = [
 [0,0,0.562638010719323,0.564489642756417,0.391160379542888,0,0.559820216521637,0,0],
 [0.46469341972691,0,0,0,0.457771429214975,0.462304167990833,0.460477062211142,0.461684882184563,0],
 [0,0.51592136413115,0,0.476469196237981,0.506496617414027,0.498901668232723,0,0,0],
 [0.718690942167919,0,0,0,0,0.70951216119684,0,0.712847861663133,0],
 [0.585048089973825,0,0.57877010813476,0,0,0,0,0.583624238737935,0],
 [0,0,0,0,0,0,0,0,0],
 [0.501201036412273,0,0.505566863835427,0,0,0.499053563465224,0,0,0],
 [0,0,0,0,0.107882861882715,0,0,0,0],
 [0,0,0.469095379810841,0.477404634251289,0,0.472662860843123,0,0,0]]


def z_of(source, target):
    i, j = PROT.index(target), PROT.index(source)
    rm, rs = RM[i][j], RS[i][j]
    return rm, rs, (abs(rm) / rs if rs > 0 else float("inf"))


def main():
    funcs, props = [], []
    props.append("# BMRA-CI sign constraints on the SKMEL-133 connection coefficients (withMyc).")
    props.append("# One line per connection: the sign of (g-1) [Eq. 24: g>1 activation, g<1")
    props.append("# inhibition] is pinned to the BMRA-inferred sign, weight = z = |rm|/rs (the")
    props.append("# BMRA confidence). Evaluated at the settled no-drug steady state (at 100000).")
    props.append("# Enforces 'connection coefficients within the BMRA CIs' (Methods p.24) as its")
    props.append("# BPSL-expressible projection. cc_*() = (g_* - 1)*totERK; totERK==1 (conserved).")
    props.append("")
    for gid, act, term, src, tgt in EDGES:
        rm, rs, z = z_of(src, tgt)
        # LEVEL effect sign from BMRA: rm>0 => source raises target.
        # prod term: raise level => g>1 ; deg term: raise level => g<1 (flip).
        level_up = rm > 0
        if term == "deg":
            want_g_gt1 = not level_up
        else:
            want_g_gt1 = level_up
        edge = gid[2:]  # strip 'g_'
        funcs.append(f"cc_{edge}() = ({gid} - 1.0)*totERK")
        op = ">" if want_g_gt1 else "<"
        w = round(z, 3)
        note = f"{src}->{tgt} BMRA r={rm:+.3f}+/-{rs:.3f} (z={z:.2f}); " \
               f"g{'>' if want_g_gt1 else '<'}1 => {'activation' if want_g_gt1 else 'inhibition'}"
        if term == "deg":
            note += " [degradation term: g-side flipped vs r sign]"
        props.append(f"# {note}")
        props.append(f"cc_{edge}() {op} 0 at 100000 weight {w}")
    props.append("")
    props.append("# DPD driving-force coefficients beta_j (Eq. 25): force = +beta_mTOR*pmTOR")
    props.append("# +beta_PKC*pPKC -beta_SRC*pSRC, so the signed force coefficient matches r_Sj.")
    props.append("# betas stay positive; the structural +/- fixes the DPD-row sign. weight = z(r_Sj).")
    for bid, sgn, src in DPD:
        rm, rs, z = z_of(src, "DPD")
        edge = "S" + src
        funcs.append(f"cc_{edge}() = {bid}*totERK")
        w = round(z, 3)
        note = f"DPD<-{src} BMRA r={rm:+.3f}+/-{rs:.3f} (z={z:.2f}); force coeff {'+' if sgn>0 else '-'}{bid}"
        props.append(f"# {note}")
        props.append(f"cc_{edge}() > 0 at 100000 weight {w}")

    with open("_gen_functions.txt", "w") as f:
        f.write("\n".join(funcs) + "\n")
    with open("bmra_signs.prop", "w") as f:
        f.write("\n".join(props) + "\n")
    print(f"wrote _gen_functions.txt ({len(funcs)} cc functions)")
    print(f"wrote bmra_signs.prop ({sum(1 for p in props if p.startswith('cc_'))} constraints)")


if __name__ == "__main__":
    main()
