"""~PLUS RNAP-elongation generator for lambda_switch_arkin1998_fullcircuit*.bngl.

Emits the state-increment ("~PLUS collapse") RNAP molecule type and the per-operon
initiation / elongation / landmark reaction rules from a compact operon spec. This is
the part of the model that future refinements touch (minor terminators tR0/tR2/tL2a,
promoter occlusion, per-nt elongation) -- add a landmark to OPERONS or a rule to
operon_rules() rather than editing thousands of enumerated states by hand.

Each operon advances a generic RNAP one nt-equivalent per step (loc~? -> loc~PLUS) at
k_step; landmarks fire location-specific rules: gene-stop emission (rate-separated
k_emit), N-antitermination at the nut site, terminator fall-off, and transcript
completion (k_comp). `at` tags the antitermination state: t (un-antiterminated),
a (antiterminated); the exact variant adds paused nut states n (N-free) and m (N-bound).
"""

MAX_LOC = 2283  # longest operon (P_L); the loc ladder runs 0..MAX_LOC

# operon layout (positions in nt from the transcription start), committed order.
# nut: antitermination landmark; emits: (loc, mRNA) gene stops mid-operon;
# term: (loc, name, k_fo_param) terminator fall-off; done: terminal transcript.
OPERONS = [
    dict(name="PR", L=630, nut=242, emits=[(218, "mRNA_cro")],
         term=(314, "tR1", "k_fo_PR_tR1"), done="mRNA_cII"),
    dict(name="PL", L=2283, nut=48, emits=[(545, "mRNA_N")],
         term=(1022, "tL1", "k_fo_PL_tL1"), done="mRNA_cIII"),
    dict(name="PRM", L=713, nut=None, emits=[], term=None, done="mRNA_cI"),
    dict(name="PRE", L=1116, nut=None, emits=[], term=None, done="mRNA_cI"),
]


def molecule_type_line(variant="base"):
    """The RNAP(...) molecule type declaration line (2 leading spaces)."""
    states = "t~a~n~m" if variant == "exact" else "t~a"
    locs = "~".join([str(i) for i in range(MAX_LOC + 1)] + ["PLUS", "MINUS"])
    return f"  RNAP(op~PR~PL~PRM~PRE,loc~{locs},at~{states})"


def _nut_rules(op, variant):
    """Antitermination rules at the nut landmark."""
    name, nut = op["name"], op["nut"]
    tag = "nutR" if name == "PR" else "nutL"

    def R(loc, at):
        return f"RNAP(op~{name},loc~{loc},at~{at})"

    if variant == "exact":
        # (nutR carries the fuller comment; nutL the short one -- as committed)
        detail = ": k23 escape, k24/k25 load, k26 escape" if name == "PR" else ""
        return [
            f"  # {tag}: slow-step + reversible N loading (Arkin Table 2{detail})",
            f"  Rpause_{name}:   {R(nut, 't')} -> {R(nut, 'n')} k_pause",
            f"  Rescape0_{name}: {R(nut, 'n')} -> {R('PLUS', 't')} k23",
            f"  RNload_{name}:   {R(nut, 'n')} -> {R(nut, 'm')} f_load()",
            f"  RNunload_{name}: {R(nut, 'm')} -> {R(nut, 'n')} k25",
            f"  Rescape1_{name}: {R(nut, 'm')} -> {R('PLUS', 'a')} k26",
        ]
    return [f"  Ranti_{name}: {R(nut, 't')} -> {R(nut, 'a')} f_anti_{name}()"]


def operon_rules(op, variant):
    """All initiation/elongation/landmark rules for one operon, in committed order."""
    name, L = op["name"], op["L"]

    def R(loc, at):
        return f"RNAP(op~{name},loc~{loc},at~{at})"

    out = [f"  # --- operon {name} (L={L} nt, ~PLUS collapse) ---"]
    out.append(f"  Rinit_{name}: DNA() -> DNA() + {R(0, 't')} A_{name}()")
    out.append(f"  RstepT_{name}: {R('?', 't')} -> {R('PLUS', 't')} k_step")
    out.append(f"  RstepA_{name}: {R('?', 'a')} -> {R('PLUS', 'a')} k_step")
    if op["nut"] is not None:
        out += _nut_rules(op, variant)
    for loc, mrna in op["emits"]:
        for at in ("t", "a"):
            out.append(f"  Remit_{name}_{mrna}_{at}: {R(loc, at)} -> "
                       f"{R('PLUS', at)} + {mrna}() k_emit")
    if op["term"] is not None:
        tloc, tname, kfo = op["term"]
        out.append(f"  Rfo_{name}_{tname}: {R(tloc, 't')} -> 0 {kfo} DeleteMolecules")
    for at in ("T", "A"):
        out.append(f"  Rdone{at}_{name}: {R(L, at.lower())} -> "
                   f"{op['done']}() k_comp DeleteMolecules")
    return out


def elongation_block(variant="base"):
    """The full transcription (initiation + elongation + landmarks) rule block."""
    out = []
    for op in OPERONS:
        out += operon_rules(op, variant)
    return out


if __name__ == "__main__":
    import sys
    v = sys.argv[1] if len(sys.argv) > 1 else "base"
    print(molecule_type_line(v)[:120] + "  ...(truncated)")
    for line in elongation_block(v):
        print(line)
