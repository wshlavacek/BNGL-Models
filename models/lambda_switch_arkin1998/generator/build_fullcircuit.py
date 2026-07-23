"""Assemble lambda_switch_arkin1998_fullcircuit{,_exact}.bngl from its parts.

This is the reconstructed generator for the network-free phage-lambda full circuit
(the original dev/lambda_arkin_fullmodel/ was local-only and lost). It regenerates the
committed models byte-for-byte:

    python build_fullcircuit.py base   > .../lambda_switch_arkin1998_fullcircuit.bngl
    python build_fullcircuit.py exact  > .../lambda_switch_arkin1998_fullcircuit_exact.bngl

Structured, hard-to-hand-edit content is generated from spec/code and spliced into the
verbatim static prose (`_sections.py`):
  - plus_gen.py   : the ~PLUS RNAP molecule type (loc 0..2283) + per-operon elongation
                    and landmark rules  <-- edit OPERONS here to add minor terminators
                    (tR0/tR2/tL2a), promoter occlusion, or per-nt elongation later
  - or_system.py  : the 40-config Shea-Ackers O_R partition function + A_PR/A_PRM
  - pre_pl.py     : the P_RE (CII-activated) and P_L (Cro2/CI2-repressed) functions

Verify a round-trip with:  python build_fullcircuit.py --check
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import or_system  # noqa: E402
import plus_gen  # noqa: E402
import pre_pl  # noqa: E402
from _sections import BASE, EXACT  # noqa: E402

MODEL_DIR = os.path.normpath(os.path.join(HERE, ".."))  # generator/ lives inside the model dir
TARGETS = {
    "base": "lambda_switch_arkin1998_fullcircuit.bngl",
    "exact": "lambda_switch_arkin1998_fullcircuit_exact.bngl",
}


def build(variant):
    """Return the full .bngl text for variant in {'base','exact'}."""
    S = EXACT if variant == "exact" else BASE
    out = []
    out += S["header"]
    out.append("begin model")
    out += S["params"]
    out.append("begin molecule types")
    out.append(plus_gen.molecule_type_line(variant))
    out += S["mt_other"]
    out.append("end molecule types")
    out += S["seed"]
    out += S["obs"]
    out.append("begin functions")
    wr, a_pr, a_prm = or_system.generate()
    out += wr + [a_pr, a_prm]
    wre, a_pre = pre_pl.generate_pre()
    out += wre + [a_pre]
    wl, a_pl = pre_pl.generate_pl()
    out += wl + [a_pl]
    out += S["fn_tail"]
    out.append("end functions")
    out.append("begin reaction rules")
    out += plus_gen.elongation_block(variant)
    out += S["rules_extra"]
    out += S["rules_growth"]
    out.append("end reaction rules")
    out.append("end model")
    out += S["actions"]
    return "\n".join(out) + "\n"


def check():
    """Diff each variant against its committed file; exit non-zero on mismatch."""
    ok = True
    for variant, fname in TARGETS.items():
        path = os.path.join(MODEL_DIR, fname)
        committed = open(path).read()
        generated = build(variant)
        same = committed == generated
        print(f"{variant:6s}: {'IDENTICAL' if same else 'MISMATCH'}  ({fname})")
        if not same:
            ok = False
            cl, gl = committed.splitlines(), generated.splitlines()
            for i in range(max(len(cl), len(gl))):
                c = cl[i] if i < len(cl) else "<none>"
                g = gl[i] if i < len(gl) else "<none>"
                if c != g:
                    print(f"  first diff at line {i + 1}:")
                    print(f"    committed: {c!r}")
                    print(f"    generated: {g!r}")
                    break
    return ok


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "--check"
    if arg == "--check":
        sys.exit(0 if check() else 1)
    elif arg in TARGETS:
        sys.stdout.write(build(arg))
    else:
        sys.stderr.write("usage: build_fullcircuit.py [base|exact|--check]\n")
        sys.exit(2)
