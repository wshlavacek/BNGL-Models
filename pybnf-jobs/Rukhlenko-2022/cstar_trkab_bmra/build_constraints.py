#!/usr/bin/env python
"""Generate the BMRA-CI sign-constraint machinery for the JOINT TrkA + TrkB BMRA job.

For each receptor model it:
  * parses the alpha_<X>_<edge>() crosstalk multipliers to recover every signalling
    connection (source->target, its g_<X>_<edge> strength, and published value);
  * maps each to the BMRA-inferred connection coefficient r[target][source] (mean rm,
    std rs) from the 10-min posterior (the training timepoint), embedded below from
    BMRA/results/Trk{A,B}_{rm,rs}_10_log_200_5K.csv (proteins TRK,ERK,AKT,JNK,S6K,RSK,RTK,STV);
  * emits a sign constraint g<X>_<edge> ≷ 1 (Eq. 24: g>1 activation, g<1 inhibition)
    ONLY where BMRA determines the sign confidently (z=|rm|/rs >= Z_THRESH) AND the
    published model already carries that sign -- the paper's own inclusion rule ("only
    interactions ... with statistically significant non-zero values are included", p.19).
    Low-confidence edges (CI includes 0) and the rare model/BMRA-mean sign disagreements
    (all low-z, e.g. TrkA TRK->AKT z=0.7) are left UNCONSTRAINED and logged.

The DPD (STV) force coefficients beta_<X>_* are NOT constrained: the DPD is defined by the
SVM state signature, not directly by the BMRA STV row, and the two disagree in sign at
non-trivial confidence (e.g. TrkB beta_B_ERK>0 vs BMRA r_STV,ERK=-0.68, z=1.55). The
signalling connections are the BMRA-constrained content here. (The phospho time-course
training data does not involve Sval, so the betas are held fixed anyway.)

Carrier encoding (as in ../cstar_skmel133_bmra): BPSL needs an observable and a
parameter-only function is not emitted to the gdat, so each constraint rides
cc_<X>_<edge>() = (g_<X>_<edge> - 1)*totERK, with totERK==1 (ERK conserved).
"""
import re

Z_THRESH = 1.0                    # keep edges whose BMRA CI is at least 1 std from zero
PROT = ["TRK", "ERK", "AKT", "JNK", "S6K", "RSK", "RTK", "STV"]
NODES = ["TRK", "ERK", "AKT", "JNK", "S6K", "RSK", "RTK"]

MODELS = {"A": "cstar_trka_bmra.bngl", "B": "cstar_trkb_bmra.bngl"}
SRCMODELS = {"A": "../cstar_trka/cstar_trka.bngl", "B": "../cstar_trkb/cstar_trkb.bngl"}

# 10-min BMRA posteriors (rows=target, cols=source; order PROT).
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


def split(edge):
    for i in range(2, len(edge) - 1):
        a, b = edge[:i], edge[i:]
        if a in NODES and b in NODES:
            return a, b
    return None


def parse(tag):
    s = open(SRCMODELS[tag]).read()
    out = []
    for m in re.finditer(rf"alpha_{tag}_(\w+)\(\)\s*=\s*\(1\.0 \+ (g_{tag}_\w+)\*", s):
        edge, g = m.groups()
        gv = float(re.search(rf"(?m)^\s*{re.escape(g)}\s*=\s*([-0-9.eE]+)", s).group(1))
        out.append((edge, g, gv))
    return out


def main():
    for tag in ("A", "B"):
        rm, rs = RM[tag], RS[tag]
        funcs, props, kept, dropped = [], [], 0, []
        props.append(f"# BMRA-CI sign constraints on the Trk{tag} connection coefficients (10-min posterior).")
        props.append(f"# g>1 activation / g<1 inhibition (Eq. 24). Kept iff the BMRA CI is >= {Z_THRESH:g} std")
        props.append(f"# from zero (statistically significant) AND the published model carries that sign;")
        props.append(f"# weight = z = |rm|/rs. Evaluated at the no-ligand steady state (at 100000).")
        props.append("")
        for edge, g, gv in parse(tag):
            sp = split(edge)
            if not sp:
                continue
            src, tgt = sp
            i, j = PROT.index(tgt), PROT.index(src)
            r, s = rm[i][j], rs[i][j]
            z = abs(r) / s if s > 0 else 0.0
            gsign = 1 if gv > 1 else (-1 if gv < 1 else 0)
            rsign = 1 if r > 0 else (-1 if r < 0 else 0)
            if z < Z_THRESH or rsign == 0 or rsign != gsign:
                reason = ("low-z" if z < Z_THRESH or rsign == 0 else "model/BMRA sign disagree")
                dropped.append(f"{src}->{tgt} g={gv:g} BMRA r={r:+.2f}+/-{s:.2f} z={z:.2f} [{reason}]")
                continue
            e = g[2:]  # strip 'g_' -> A_TRKERK
            funcs.append(f"cc_{e}() = ({g} - 1.0)*totERK")
            op = ">" if gsign > 0 else "<"
            props.append(f"# {src}->{tgt} BMRA r={r:+.3f}+/-{s:.3f} (z={z:.2f}); "
                         f"g{op}1 => {'activation' if gsign > 0 else 'inhibition'}")
            props.append(f"cc_{e} {op} 0 at 100000 weight {round(z, 3)}")
            kept += 1
        with open(f"_gen_functions_{tag}.txt", "w") as f:
            f.write("\n".join(funcs) + "\n")
        with open(f"bmra_{tag}.prop", "w") as f:
            f.write("\n".join(props) + "\n")
        print(f"Trk{tag}: kept {kept} sign constraints; dropped {len(dropped)}:")
        for d in dropped:
            print("    -", d)


if __name__ == "__main__":
    main()
