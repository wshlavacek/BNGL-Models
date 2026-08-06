"""Emit the well-mixed counterpart of the Artyomov et al. (2010) lattice model.

The published model is a reaction-diffusion lattice and its SSC listing is in the
paper's `sapp.doc`. This script transcribes the rules of that listing -- the first
one, "SSC code for simulating irreversible dissociation of MHC off the T-cell
surface", which is what Fig. 2A reports -- and converts the rate constants to a
single well-mixed compartment.

The conversion is fixed by the paper's own SI Methods, Eq. S1b: a propensity in
the lattice is `k * n_A * n_B` in per-chamber counts. Averaged over the 10^4
chambers of the 1 um^2 contact, a bimolecular propensity in total counts is
therefore `k/10^4 * N_A * N_B`, while a unimolecular one is unchanged. So every
rule whose reactants are two separate complexes has its constant divided by the
chamber count, and every rule that closes a ring within one complex keeps the
published constant. Nothing else about the rules changes.

Two files are produced:

  primary       the exact mean-field limit just described
  _encounter    the same rules with intermolecular constants further multiplied by
                the paper's own encounter probability, k/(k + m_A + m_B). This
                partially reintroduces the diffusion limitation the mean-field
                limit discards, and is a hybrid rather than a well-mixed model;
                it is here because it is the paper's own correction and it shows
                how much of the gap that correction accounts for.
"""
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
STEM = "lck_delivery_wellmixed_artyomov2010"

# Table 1 of the paper, in units of per second (the starred entries are the
# experimentally derived ones). kon, konCD, konLck and konLck1 are per-chamber
# bimolecular constants; kdiff is a per-molecule hop rate.
PARAMS = [
    ("n_chambers", "1e4", "chambers in the 1 um^2 contact, SI Methods Eq. S3"),
    ("kon_local", "150", "TCR-MHC on rate, per chamber  [Table 1]"),
    ("koffAg", "0.02", "TCR-MHC off, agonist peptide  [Table 1]"),
    ("koffEn", "20", "TCR-MHC off, endogenous peptide  [Table 1]"),
    ("konCD_local", "1000", "MHC-coreceptor on rate, per chamber  [Table 1]"),
    ("koffCD", "20", "MHC-coreceptor off rate; the x axis of Fig. 2A  [Table 1]"),
    ("konLck_local", "1", "Lck(coreceptor)-TCR on rate, per chamber  [Table 1]"),
    ("konLck1_local", "1", "the same, for the case with no MHC involved"),
    ("koffLck", "1", "Lck(coreceptor)-TCR off rate  [Table 1]"),
    ("kdiff", "50", "hop rate of a surface protein  [Table 1]"),
    ("kdestroy", "1e8", "loss of an MHC with both sites free, from the SSC listing"),
    ("count_ag", "10", "agonist pMHC, pre-bound to TCR; not stated in the paper"),
    ("count_en", "0", "endogenous pMHC; not stated in the paper"),
    ("count_tcr", "200", "TCRs, from the SSC listing"),
    ("count_cd4", "100", "coreceptors, from the SSC listing"),
]

MOLECULES = ["MHC(t,c,p~ag~en)", "TCR(m,c,p~p0)", "CD4(m,t,lck~basal)"]

# (name, reactants, products, rate parameter, arity)
#   "bi"  two separate complexes  -> divide the local constant by n_chambers
#   "uni" one complex, ring closure -> keep the local constant
RULES = [
    # Loss of an MHC that is bound to nothing. The SSC listing gives this an
    # enormous rate, which makes dissociation from the surface irreversible.
    ("mhc_lost", "MHC(t,c)", "0", "kdestroy", "uni"),

    # MHC-TCR bond. The four SSC variants enumerate whether each partner already
    # carries a coreceptor, and the fourth is the ring closure through one CD4.
    ("tcr_mhc_on", "MHC(t,c) + TCR(m,c)", "MHC(t!1,c).TCR(m!1,c)", "kon", "bi"),
    ("tcr_mhc_on_mhc_cd", "MHC(t,c!1).CD4(m!1,t) + TCR(m,c)",
     "MHC(t!2,c!1).CD4(m!1,t).TCR(m!2,c)", "kon", "bi"),
    ("tcr_mhc_on_tcr_cd", "MHC(t,c) + TCR(m,c!1).CD4(t!1,m)",
     "MHC(t!2,c).TCR(m!2,c!1).CD4(t!1,m)", "kon", "bi"),
    ("tcr_mhc_close", "MHC(t,c!2).TCR(m,c!1).CD4(t!1,m!2)",
     "MHC(t!3,c!2).TCR(m!3,c!1).CD4(t!1,m!2)", "kon_local", "uni"),

    ("tcr_mhc_off_ag", "TCR(m!1).MHC(t!1,p~ag)", "TCR(m) + MHC(t,p~ag)", "koffAg", "uni"),
    ("tcr_mhc_off_en", "TCR(m!1).MHC(t!1,p~en)", "TCR(m) + MHC(t,p~en)", "koffEn", "uni"),

    # TCR-coreceptor bond, i.e. delivery of Lck to the TCR.
    ("tcr_cd_on", "TCR(c,m) + CD4(t,m)", "TCR(c!1,m).CD4(t!1,m)", "konLck1", "bi"),
    ("tcr_cd_on_via_mhc", "TCR(c,m!1).MHC(t!1,c) + CD4(t,m)",
     "TCR(c!2,m!1).MHC(t!1,c).CD4(t!2,m)", "konLck", "bi"),
    ("tcr_cd_on_cd_mhc", "TCR(c,m) + CD4(t,m!1).MHC(c!1,t)",
     "TCR(c!2,m).CD4(t!2,m!1).MHC(c!1,t)", "konLck", "bi"),
    ("tcr_cd_close", "TCR(m!2,c).CD4(t,m!1).MHC(c!1,t!2)",
     "TCR(m!2,c!3).CD4(t!3,m!1).MHC(c!1,t!2)", "konLck_local", "uni"),

    ("tcr_cd_off_ring", "TCR(c!1,m!2).MHC(t!2).CD4(t!1)",
     "TCR(c,m!2).MHC(t!2) + CD4(t)", "koffLck", "uni"),
    ("tcr_cd_off", "TCR(c!1,m).CD4(t!1)", "TCR(c,m) + CD4(t)", "koffLck", "uni"),

    # MHC-coreceptor bond.
    ("mhc_cd_on", "MHC(c,t) + CD4(m,t)", "MHC(c!1,t).CD4(m!1,t)", "konCD", "bi"),
    ("mhc_cd_on_mhc_tcr", "MHC(c,t!1).TCR(m!1,c) + CD4(m,t)",
     "MHC(c!2,t!1).TCR(m!1,c).CD4(m!2,t)", "konCD", "bi"),
    ("mhc_cd_on_cd_tcr", "MHC(c,t) + CD4(m,t!1).TCR(c!1,m)",
     "MHC(c!2,t).CD4(m!2,t!1).TCR(c!1,m)", "konCD", "bi"),
    ("mhc_cd_close", "MHC(t!2,c).CD4(m,t!1).TCR(c!1,m!2)",
     "MHC(t!2,c!3).CD4(m!3,t!1).TCR(c!1,m!2)", "konCD_local", "uni"),

    ("mhc_cd_off", "MHC(c!1).CD4(m!1)", "MHC(c) + CD4(m)", "koffCD", "uni"),

    # Ring opening, which has to be spelled out. Each of the three off-rules above
    # writes its products with `+`, so BioNetGen applies it only where breaking the
    # bond really does separate the complex. Inside the closed MHC-TCR-coreceptor
    # ring it does not, so without these four rules the ring is generated with no
    # reaction out of it at all -- an absorbing species, and no warning. The two
    # peptide states need separate rules only because they carry different off rates.
    ("ring_open_tcr_mhc_ag",
     "MHC(t!1,c!3,p~ag).TCR(m!1,c!2).CD4(t!2,m!3)",
     "MHC(t,c!3,p~ag).TCR(m,c!2).CD4(t!2,m!3)", "koffAg", "uni"),
    ("ring_open_tcr_mhc_en",
     "MHC(t!1,c!3,p~en).TCR(m!1,c!2).CD4(t!2,m!3)",
     "MHC(t,c!3,p~en).TCR(m,c!2).CD4(t!2,m!3)", "koffEn", "uni"),
    ("ring_open_tcr_cd",
     "MHC(t!1,c!3).TCR(m!1,c!2).CD4(t!2,m!3)",
     "MHC(t!1,c!3).TCR(m!1,c).CD4(t,m!3)", "koffLck", "uni"),
    ("ring_open_mhc_cd",
     "MHC(t!1,c!3).TCR(m!1,c!2).CD4(t!2,m!3)",
     "MHC(t!1,c).TCR(m!1,c!2).CD4(t!2,m)", "koffCD", "uni"),
]

HEADER = """begin model

#@title: Lck delivery by CD4/CD8, well-mixed counterpart (Artyomov et al., 2010){vt}

#@description: |
#  A WELL-MIXED counterpart of the reaction-diffusion lattice model of Artyomov et al.
#  (2010). Agonist and endogenous pMHC on an APC engage TCR on a T cell; the coreceptor
#  CD4 (standing for a constitutively associated CD4-Lck complex) binds MHC and, through
#  the TCR's intracellular domain, delivers Lck. Because MHC, TCR and coreceptor can each
#  bind the other two, the three of them close a ring, and it is that ring that keeps MHC
#  on the T-cell surface after the TCR-pMHC bond breaks.
#
#  How it works: the rules are transcribed from the first SSC listing in the paper's
#  supporting file, the one headed "irreversible dissociation of MHC off the T-cell
#  surface", which is what Fig. 2A reports. An MHC left bound to nothing is destroyed at
#  an enormous rate, so leaving the surface is irreversible and the recorded observable
#  is a survival curve. Its half-time is the paper's "effective half-life".
#
#  What it does: it is NOT expected to reproduce Fig. 2A, and the gap is the point.
#  See the note below.{vd}

#@keyword: |
#  T-cell receptor, coreceptor, CD4, CD8, Lck, pMHC, kinetic proofreading, avidity,
#  ring closure, trimolecular complex, well-mixed limit, mean-field, spatial model,
#  reaction-diffusion, stochastic, Mus musculus, Homo sapiens

#@reference: |
#  Artyomov MN, Lis M, Devadas S, Davis MM, Chakraborty AK (2010). CD4 and CD8 binding to
#  MHC molecules primarily acts to enhance Lck delivery. Proceedings of the National
#  Academy of Sciences 107:16916-16921. doi:10.1073/pnas.1010568107

#@note: |
#  SIMULATION INTENT. Population-based and stochastic: every amount is a molecule count on
#  the 1 um^2 T-cell/APC contact the paper simulates. Stochastic rather than ODE on
#  purpose -- the published model is stochastic too, so running this one stochastically
#  leaves spatial structure as the only difference between them, which is what the folder
#  sets out to measure.

#@note: |
#  THIS IS A COUNTERPART, NOT A REPRODUCTION. The published model is a 100x100 lattice of
#  chambers over the contact area, with per-molecule hopping and bound complexes given
#  diffusion 0. There is no non-spatial variant anywhere in the paper. This file is the
#  well-mixed limit of the same rules, built so that the difference between the two can be
#  measured; it is not an attempt to reproduce Fig. 2A, and it does not.
#
#  The conversion is the paper's own, SI Methods Eq. S1b: a lattice propensity is
#  k * n_A * n_B in PER-CHAMBER counts, so averaged over the 10^4 chambers a bimolecular
#  propensity in total counts is (k/10^4) * N_A * N_B, while a unimolecular one is
#  unchanged. Every rule whose reactants are two separate complexes is therefore divided
#  by n_chambers, and every ring closure keeps the published constant. That is the only
#  change to the rules.

#@note: |
#  WHAT THE COMPARISON ISOLATES. Not "spatial heterogeneity" loosely: this model differs
#  from the published one by finite diffusion and the local depletion and rebinding that
#  come with it. The lattice is diffusion-limited by construction -- kon per chamber is
#  150 and konCD is 1000 against a hop rate of 50 -- and bound complexes do not move at
#  all, so a partner that has just let go is still next door. The well-mixed limit throws
#  all of that away and lets any two molecules meet at the average rate.
#
#  The model carries its own control. With no coreceptor the ring cannot form, the
#  agonist MHC simply falls off at koffAg and is lost, and the survival curve is exactly
#  exp(-koffAg t) whatever the geometry: a half-time of ln2/0.02 = 34.66 s, against the
#  ~35 s the paper reads off Fig. 2A. That arm must agree and does, which is what makes
#  the disagreement in the coreceptor arm attributable to geometry rather than to the
#  translation.

#@note: |
#  TWO COUNTS ARE NOT IN THE PAPER. The SSC listing seeds `count_ag` agonist pMHC
#  pre-bound to TCR and `count_en` endogenous pMHC, and neither number appears in the
#  paper, the SI or the listing. They are parameters here, defaulting to 10 agonist as a
#  tracer against 200 TCR and 100 coreceptors, and no endogenous pMHC, which leaves the
#  coreceptor pool uncontested and the measurement unconfounded. verify_artyomov2010.ipynb
#  scans both and reports how much the half-time moves.

#@figure: Fig. 2 A in Artyomov et al. (2010)
"""


VARIANT_DESCRIPTION = (
    "\n#\n#  This variant multiplies every intermolecular constant by the paper's own\n"
    "#  encounter probability k/(k + m_A + m_B) (its 1/201 in the Discussion), which puts\n"
    "#  back some of the diffusion limitation the mean-field limit discards. It is a\n"
    "#  hybrid, not a well-mixed model: an effective rate carrying spatial information.\n"
    "#  It is here to show how much of the gap that one correction accounts for.")


def build(encounter):
    out = [HEADER.format(
        vt=" -- with the paper's encounter correction" if encounter else "",
        vd=VARIANT_DESCRIPTION if encounter else "")]

    out.append("\nbegin parameters")
    for name, value, comment in PARAMS:
        out.append(f"  {name}  {value}  # {comment}")
    out.append("")
    out.append("  # Mean-field limit of a per-chamber bimolecular constant (SI Methods Eq. S1b).")
    for k in ("kon", "konCD", "konLck", "konLck1"):
        if encounter:
            out.append(f"  # encounter-corrected: {k} = {k}_local/n_chambers * P, "
                       f"P = {k}_local/({k}_local + 2*kdiff)")
            out.append(f"  {k}  {k}_local/n_chambers * {k}_local/({k}_local + 2*kdiff)")
        else:
            out.append(f"  {k}  {k}_local/n_chambers")
    out.append("end parameters")

    out.append("\nbegin molecule types")
    out += [f"  {m}" for m in MOLECULES]
    out.append("end molecule types")

    out.append("\nbegin seed species")
    out.append("  MHC(t!1,c,p~ag).TCR(m!1,c,p~p0)  count_ag  # agonist, seeded already engaged")
    out.append("  MHC(t,c,p~en)                    count_en")
    out.append("  TCR(m,c,p~p0)                    count_tcr")
    out.append("  CD4(m,t,lck~basal)               count_cd4")
    out.append("end seed species")

    out.append("\nbegin observables")
    out.append("  Molecules AgMHC_on_TCR  MHC(t!1,p~ag).TCR(m!1)  # the SSC listing's `record`")
    out.append("  Molecules AgMHC_total   MHC(p~ag)")
    out.append("  Molecules Ring          MHC(t!1,c!2).TCR(m!1,c!3).CD4(t!3,m!2)")
    out.append("  Molecules CD4_free      CD4(m,t)")
    out.append("end observables")

    out.append("\nbegin reaction rules")
    for name, lhs, rhs, rate, arity in RULES:
        out.append(f"  {name}: {lhs} -> {rhs}  {rate}")
    out.append("end reaction rules")

    out.append("\nend model\n")
    out.append("begin actions\n")
    out.append("# The network is small and closes, so it is generated rather than run network-free.")
    out.append("generate_network({overwrite=>1})")
    out.append('simulate({method=>"ssa",t_end=>250,n_steps=>250})\n')
    out.append("end actions")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    drift = False
    for enc, suffix in ((False, ""), (True, "_encounter")):
        path = FOLDER / f"{STEM}{suffix}.bngl"
        text = build(enc)
        if args.check and path.exists() and path.read_text() != text:
            drift = True
            print(f"DRIFT {path.name}")
        path.write_text(text)
        print(f"wrote {path.name}")
    raise SystemExit(1 if drift else 0)


if __name__ == "__main__":
    main()
