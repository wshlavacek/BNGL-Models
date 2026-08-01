"""Emit the house-style BNGL files for models/amyloid_beta_competing_aggregation_pathways_rana2020.

The Rana et al. (2020) EKS network is an explicitly enumerated oligomer chain (A_1..A_12 on
pathway, A'_4..A'_12 off pathway), so every rule has to be written out: BNGL state labels are
opaque tokens and cannot be incremented. This script writes the rule blocks so the 80 reactions
stay consistent across the six protocol files.
"""
from __future__ import annotations

import sys
from pathlib import Path

STEM = "amyloid_beta_competing_aggregation_pathways_rana2020"

# ---------------------------------------------------------------- parameter values
#
# Rate constants: the best fit recovered by refitting this network to the digitized Fig. 1
# traces, because Rana et al.'s Table I is incomplete (five of the thirteen constants
# appear in no column) and its on-pathway pair is incompatible with the reported 25 uM
# peptide concentration. The on-pathway file uses the `on_pathway` job's fit; the five
# competing-pathway files use the `global` job's fit, so one parameterization serves all
# five published experiments. See pybnf-jobs/Rana-2020/ and verify_rana2020.ipynb.
DEFAULTS = dict(
    k_nuon=1.312541129453643,
    k_nuon_=157.19086714816115,
    k_fbon=337.3859751693723,
    k_fbon_=43.18236738549619,
    k_con=39.491291263895995,
    k_con_=0.584679016282277,
    k_nuoff=322.0135857401437,
    k_nuoff_=43.60981468930446,
    k_fboff=0.1910200863130847,
    k_fboff_=0.01540333874969174,
    k_swi=48815.77413936899,
    k_swi_=1.1411755300932513,
    map_on=0.647864332611721, map_off=1.9045703581464974,
    Abeta_0=25.0, Mic_0=0.0, Mic_present=100.0,
)
# The on-pathway-only file is the sub-model Rana et al. fit first (Sec. IV-A) and carries
# that fit's own values; the five competing-pathway files carry the global fit.
# Table I gives the global fit two off-pathway ThT mapping constants, k_off1 and k_off2:
# the addition series (Figs. 1a and 1c) and the removal series (Fig. 1b) are separate
# measurements with their own arbitrary normalization, and Fig. 1b plateaus near 1.5 a.u.
# where Fig. 1c plateaus near 1.0. DEFAULTS["map_off"] is k_off1; the two removal protocols
# carry k_off2 instead.
REMOVAL_MAP_OFF = 2.094255464135117

ON_PATHWAY_FIT = dict(
    k_nuon=3.6676839864546986,
    k_nuon_=82.14226600765491,
    k_fbon=6.724404486241758,
    k_fbon_=9.989681873346637,
    map_on=8.669274980249881,
)

V = dict(DEFAULTS)


def g(name):
    return f"{V[name]:.6g}"


HEADER_KEYWORDS = """#@keyword: |
#  protein aggregation, amyloid-beta, Alzheimer disease, competing pathways,
#  nucleation-elongation, fatty acid micelle, pathway switching, mass action,
#  lag phase, sigmoidal kinetics, Homo sapiens
"""

REFERENCE = """#@reference: |
#  Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P (2020). Global fitting and
#  parameter identifiability for amyloid-beta aggregation with competing
#  pathways. 2020 IEEE 20th International Conference on BioInformatics and
#  BioEngineering (BIBE), pp. 73-78. doi:10.1109/BIBE50027.2020.00020
#
#  Source of the ThT data reproduced here:
#  Ghosh P, Rana P, Rangachari V, Saha J, Steen E, Vaidya A (2020). A game
#  theoretic approach to deciphering the dynamics of amyloid-beta aggregation
#  along competing pathways. R Soc Open Sci 7(8):191814. doi:10.1098/rsos.191814
"""


def parameters_block(*, off_pathway: bool, spikes: bool) -> str:
    L = []
    A = L.append
    A("begin parameters")
    A("  # Constants for conversion to a population-based (molecule-count) unit system.")
    A("  # V_ref is a 100 uL plate well, the format in which the ThT assays of Rana")
    A("  # et al. (2020) are read.")
    A("  NA     6.02214076e23  # molecules/mol")
    A("  V_ref  1e-4           # L/well")
    A("")
    A("  # 25 uM Abeta42 is the peptide concentration stated in Sec. II-A of Rana")
    A("  # et al. (2020).")
    A(f"  Abeta_0      {g('Abeta_0')}  # uM")
    if off_pathway:
        A("")
        A("  # Mic_present is the pseudo-micelle concentration formed by 5 mM C12 fatty")
        A("  # acid, taken as 100 uM (aggregation number 50). The aggregation number is not")
        A("  # reported and is not identifiable from these data: while micelles are in")
        A("  # excess only the product k_con*[L] is constrained, so Mic_present is fixed and")
        A("  # k_con absorbs it. 100 uM is a large excess over the ~6 uM of micelle the off")
        A("  # pathway can consume from 25 uM peptide. Mic_0 is how much of it is present at")
        A("  # t = 0 in this protocol.")
        A(f"  Mic_present  {g('Mic_present')}  # uM")
        A(f"  Mic_0        {g('Mic_0')}  # uM")
    A("")
    A("  # On-pathway nucleation, A_i + A_1 <-> A_i+1; see Eq. 4 and flux H_i in")
    A("  # Rana et al. (2020).")
    A(f"  k_nuon   {g('k_nuon')}  # /uM/h")
    A(f"  k_nuon_  {g('k_nuon_')}  # /h")
    A("")
    A("  # On-pathway fibril binding, A_i + F <-> F; see Eq. 4 and flux I_i in")
    A("  # Rana et al. (2020). F is the 12-mer nucleus, A_12.")
    A(f"  k_fbon   {g('k_fbon')}  # /uM/h")
    A(f"  k_fbon_  {g('k_fbon_')}  # /h")
    if off_pathway:
        A("")
        A("  # Off-pathway condensation on a pseudo-micelle, 4 A_1 + L <-> A'_4; see")
        A("  # Eq. 4-II and flux G'_1 in Rana et al. (2020).")
        A(f"  k_con    {g('k_con')}  # /uM^4/h")
        A(f"  k_con_   {g('k_con_')}  # /h")
        A("")
        A("  # BioNetGen divides the propensity of a rule with n identical reactant")
        A("  # patterns by (n-1)!, so the four A_1 patterns of R_con_f carry a factor")
        A("  # 1/3! that the published flux G'_1 = k_con[A_1]^4[L] does not. k_con_bng")
        A("  # cancels it, keeping the rate law identical to the paper's.")
        A("  k_con_bng  6*k_con   # /uM^4/h")
        A("")
        A("  # Off-pathway nucleation, A'_i + A_1 <-> A'_i+1; see Eq. 4-II and flux H'_i")
        A("  # in Rana et al. (2020).")
        A(f"  k_nuoff   {g('k_nuoff')}  # /uM/h")
        A(f"  k_nuoff_  {g('k_nuoff_')}  # /h")
        A("")
        A("  # Off-pathway oligomer capture, A'_12 + A'_i <-> F'_1; see Eq. 4-II and flux")
        A("  # I'_i in Rana et al. (2020). k_fboff is k_el1f of Eq. 4-II.")
        A(f"  k_fboff   {g('k_fboff')}  # /uM/h")
        A(f"  k_fboff_  {g('k_fboff_')}  # /h")
        A("")
        A("  # Pathway switching at the tetramer, A_4 <-> A'_4; see Eq. 5 and flux J in")
        A("  # Rana et al. (2020).")
        A(f"  k_swi   {g('k_swi')}  # /h")
        A(f"  k_swi_  {g('k_swi_')}  # /h")
    A("")
    A("  # Thioflavin T mapping constants. Rana et al. (2020) write the measured signal")
    A("  # as signal_on + signal_off and fit an off-pathway mapping constant (map' in")
    A("  # Sec. IV-B, k_off / k_off1 / k_off2 in Table I) over 10^0-10^5.")
    A(f"  map_on   {g('map_on')}  # a.u./uM")
    if off_pathway:
        A(f"  map_off  {g('map_off')}  # a.u./uM")
    A("end parameters")
    return "\n".join(L)


def molecule_types(off_pathway: bool) -> str:
    L = ["begin molecule types"]
    if off_pathway:
        L.append("  # Abeta42 oligomer: p is the pathway it belongs to, n its size in monomers.")
        L.append("  # Ab(p~on,n~12) is the on-pathway fibril F; Ab(p~off,n~4) is A'_4.")
        L.append("  Ab(p~on~off,n~1~2~3~4~5~6~7~8~9~10~11~12)")
        L.append("  # Lumped off-pathway 12-23mer, F'_1: kinetically trapped and unable to")
        L.append("  # aggregate further, so it carries no size state.")
        L.append("  Foff()")
        L.append("  # C12 fatty-acid pseudo-micelle, L.")
        L.append("  Mic()")
    else:
        L.append("  # Abeta42 on-pathway oligomer of size n monomers. Ab(p~on,n~12) is the")
        L.append("  # 12-mer nucleus, which Rana et al. (2020) identify with the fibril F.")
        L.append("  # The p component is retained so that this file and the competing-pathway")
        L.append("  # files share one molecule-type declaration.")
        L.append("  Ab(p~on~off,n~1~2~3~4~5~6~7~8~9~10~11~12)")
    L.append("end molecule types")
    return "\n".join(L)


def seed_species(off_pathway: bool) -> str:
    L = ["begin seed species"]
    L.append("  Ab(p~on,n~1)  Abeta_0  # uM")
    if off_pathway:
        L.append("  Mic()         Mic_0    # uM")
    L.append("end seed species")
    return "\n".join(L)


def observables(off_pathway: bool) -> str:
    L = ["begin observables"]
    L.append("  # On-pathway oligomers A_1..A_11.")
    for i in range(1, 12):
        L.append(f"  Molecules  Obs_A{i}  Ab(p~on,n~{i})")
    L.append("")
    L.append("  # On-pathway 12-mer: the nucleus, identified with the fibril F.")
    L.append("  Molecules  Obs_F  Ab(p~on,n~12)")
    if off_pathway:
        L.append("")
        L.append("  # Off-pathway oligomers A'_4..A'_12.")
        for i in range(4, 13):
            L.append(f"  Molecules  Obs_Ap{i}  Ab(p~off,n~{i})")
        L.append("")
        L.append("  # Lumped off-pathway 12-23mer F'_1 and the free pseudo-micelle pool L.")
        L.append("  Molecules  Obs_Fp1  Foff()")
        L.append("  Molecules  Obs_L    Mic()")
    L.append("end observables")
    return "\n".join(L)


def functions(off_pathway: bool) -> str:
    L = ["begin functions"]
    L.append("  # Thioflavin T signal. Rana et al. (2020) report the measured trace as the sum")
    L.append("  # of an on-pathway and an off-pathway contribution; the cross-beta-rich")
    L.append("  # species are the on-pathway fibril F and the off-pathway aggregate F'_1.")
    L.append("  ThT_on()     = map_on*Obs_F")
    if off_pathway:
        L.append("  ThT_off()    = map_off*Obs_Fp1")
        L.append("  ThT_total()  = map_on*Obs_F + map_off*Obs_Fp1")
    L.append("")
    L.append("  # Dimensional conversions (uncomment for molecule-count output)")
    L.append("# F_count()    = Obs_F*1e-6*NA*V_ref  # molecules/well")
    if off_pathway:
        L.append("# Fp1_count()  = Obs_Fp1*1e-6*NA*V_ref  # molecules/well")
    L.append("end functions")
    return "\n".join(L)


def rules(off_pathway: bool) -> str:
    L = ["begin reaction rules"]
    A = L.append
    A("  # --- On pathway, Eq. 4 in Rana et al. (2020) ---")
    A("")
    A("  # Nucleation by monomer addition, A_i + A_1 <-> A_i+1 (flux H_i). Forward and")
    A("  # reverse are written as separate rules because the reverse of a")
    A("  # molecule-deleting rule is a molecule-synthesizing rule.")
    for i in range(1, 12):
        A(f"  R_nuon_{i}_f: Ab(p~on,n~{i}) + Ab(p~on,n~1) -> Ab(p~on,n~{i+1})  k_nuon")
        A(f"  R_nuon_{i}_r: Ab(p~on,n~{i+1}) -> Ab(p~on,n~{i}) + Ab(p~on,n~1)  k_nuon_")
    A("")
    A("  # Fibril binding, A_i + F <-> F (flux I_i). The fibril is catalytic: it takes up")
    A("  # an oligomer without changing its own count, which is how Rana et al. (2020)")
    A("  # write the reaction and the flux.")
    for i in range(1, 12):
        A(f"  R_fbon_{i}_f: Ab(p~on,n~{i}) + Ab(p~on,n~12) -> Ab(p~on,n~12)  k_fbon")
        A(f"  R_fbon_{i}_r: Ab(p~on,n~12) -> Ab(p~on,n~12) + Ab(p~on,n~{i})  k_fbon_")
    if off_pathway:
        A("")
        A("  # --- Off pathway, Eq. 4-II in Rana et al. (2020) ---")
        A("")
        A("  # Condensation of four monomers on a pseudo-micelle, 4 A_1 + L <-> A'_4")
        A("  # (flux G'_1). See the k_con_bng comment in the parameters block.")
        A("  R_con_f: Ab(p~on,n~1) + Ab(p~on,n~1) + Ab(p~on,n~1) + Ab(p~on,n~1) + Mic() \\")
        A("    -> Ab(p~off,n~4)  k_con_bng")
        A("  R_con_r: Ab(p~off,n~4) -> Ab(p~on,n~1) + Ab(p~on,n~1) + Ab(p~on,n~1) \\")
        A("    + Ab(p~on,n~1) + Mic()  k_con_")
        A("")
        A("  # Off-pathway nucleation by monomer addition, A'_i + A_1 <-> A'_i+1 (flux H'_i).")
        for i in range(4, 12):
            A(f"  R_nuoff_{i}_f: Ab(p~off,n~{i}) + Ab(p~on,n~1) -> Ab(p~off,n~{i+1})  k_nuoff")
            A(f"  R_nuoff_{i}_r: Ab(p~off,n~{i+1}) -> Ab(p~off,n~{i}) + Ab(p~on,n~1)  k_nuoff_")
        A("")
        A("  # Capture of a 4-11mer by an off-pathway 12-mer, A'_12 + A'_i <-> F'_1")
        A("  # (flux I'_i). Both oligomers are consumed into the lumped 12-23mer.")
        for i in range(4, 12):
            A(f"  R_fboff_{i}_f: Ab(p~off,n~12) + Ab(p~off,n~{i}) -> Foff()  k_fboff")
            A(f"  R_fboff_{i}_r: Foff() -> Ab(p~off,n~12) + Ab(p~off,n~{i})  k_fboff_")
        A("")
        A("  # --- Pathway switching, Eq. 5 in Rana et al. (2020) ---")
        A("")
        A("  # Rana et al. (2020) allow switching at any oligomer level but their minimal")
        A("  # model realizes it at the tetramer, which is what flux J (Eq. 6) states.")
        A("  R_swi: Ab(p~on,n~4) <-> Ab(p~off,n~4)  k_swi, k_swi_")
    L.append("end reaction rules")
    return "\n".join(L)


def build(title, description, note, *, off_pathway, spikes, actions):
    parts = [
        "begin model",
        "",
        f"#@title: {title}",
        "",
        description.rstrip(),
        "",
        HEADER_KEYWORDS.rstrip(),
        "",
        REFERENCE.rstrip(),
        "",
        note.rstrip(),
        "",
        parameters_block(off_pathway=off_pathway, spikes=spikes),
        "",
        molecule_types(off_pathway),
        "",
        seed_species(off_pathway),
        "",
        observables(off_pathway),
        "",
        functions(off_pathway),
        "",
        rules(off_pathway),
        "",
        "end model",
        "",
        actions.rstrip(),
        "",
    ]
    return "\n".join(parts)


# ---------------------------------------------------------------- per-file specs
COMMON_NOTE = """#@note: |
#  PARAMETER PROVENANCE. Rana et al. (2020) specify the reaction network completely
#  (Eq. 4, Eq. 5 and the flux equations), but Table I tabulates only eight of the
#  thirteen rate constants: k_fbon, k_fbon_, k_con_, k_fboff_ and k_swi_ are never
#  given, and neither is the pseudo-micelle concentration or the size of the
#  addition/removal spikes. The tabulated values also cannot be used as they stand:
#  at the peptide concentration the paper reports (25 uM), the On-Pathway column
#  (k_nuon 22.04 /uM/h, k_nuon_ 12.72 /h) drives the whole nucleation chain to
#  equilibrium within ~0.03 h, which leaves no slow variable and therefore no lag
#  phase, whereas Fig. 1a shows a ~28 h lag. The values here were therefore obtained
#  by refitting this network to the ThT traces digitized from Fig. 1, which is the
#  paper's own procedure; see reference/ and the PyBNF job setups under
#  pybnf-jobs/Rana-2020/. verify_rana2020.ipynb reports both the agreement with the
#  digitized data and the comparison against Table I.
#
#  RULE ENUMERATION. The oligomer size is a BNGL state label, and state labels are
#  opaque tokens that cannot be incremented, so each step of the A_1..A_12 and
#  A'_4..A'_12 chains is written as its own rule. The model is therefore explicit
#  rather than generative; this is a faithful transcription of an ODE model, not a
#  combinatorial rule-based one.
#
#  MASS BALANCE. Total peptide is not conserved, by construction: the published
#  fibril-binding reaction A_i + F <-> F consumes an i-mer without enlarging F, and
#  its reverse emits an i-mer from F. Fibril mass is not a state variable of the
#  published model and is not added here.
#
#  ODE ONLY. 25 uM in a 100 uL well is 1.5e15 molecules, so stochastic simulation of
#  this model is not computationally feasible and no ssa/psa/nfr protocol is given.
#  The assay is a bulk in-vitro measurement for which the deterministic limit is the
#  intended interpretation."""


def actions_simple(t_end, n_steps, figure, protocol, nspecies, nrxn):
    return f"""begin actions
  generate_network({{overwrite=>1}})  # {nspecies} species, {nrxn} reactions

  # {protocol['title']}

{protocol['body']}

  #@figure: {figure}

  simulate({{method=>"ode",suffix=>"ode",t_start=>0,t_end=>{t_end},\\
    n_steps=>{n_steps},print_functions=>1}})
end actions"""


def actions_event(t_event, t_end, n1, n2, figure, protocol, nspecies, nrxn, spikes):
    spike_lines = "\n".join(f"  {s}" for s in spikes)
    return f"""begin actions
  generate_network({{overwrite=>1}})  # {nspecies} species, {nrxn} reactions

  # {protocol['title']}

{protocol['body']}

  #@figure: {figure}

  simulate({{method=>"ode",suffix=>"ode",t_start=>0,t_end=>{t_event},\\
    n_steps=>{n1},print_functions=>1}})
{spike_lines}
  simulate({{method=>"ode",suffix=>"ode",t_start=>{t_event},t_end=>{t_end},\\
    n_steps=>{n2},continue=>1,print_functions=>1}})
end actions"""


DESC_FULL = """#@description: |
#  Ensemble-kinetics model of amyloid-beta 42 aggregation along two competing
#  pathways. On pathway, monomer adds stepwise to a growing oligomer
#  (A_i + A_1 <-> A_i+1, i = 1..11) up to the 12-mer nucleus A_12, which Rana et al.
#  (2020) identify with the fibril F; the fibril then takes up oligomers
#  catalytically (A_i + F <-> F). Off pathway, four monomers condense on a C12
#  fatty-acid pseudo-micelle L (4 A_1 + L <-> A'_4), the off-pathway oligomer grows
#  by monomer addition to A'_12, and A'_12 captures a 4-11mer to give the lumped,
#  kinetically trapped 12-23mer F'_1. The two pathways exchange at the tetramer
#  (A_4 <-> A'_4). All kinetics are mass action. The thioflavin T signal is read as
#  map_on*F + map_off*F'_1.
#  Dynamically the model is a competition for a single monomer pool. Without fatty
#  acid the on pathway alone runs, and the sequential 11-step nucleation chain
#  produces a long lag (~28 h at 25 uM) followed by a steep sigmoidal rise to a
#  fibril plateau. Adding pseudo-micelles opens the off pathway, which has no lag:
#  it consumes monomer immediately, arrests fibril formation and settles at a lower
#  plateau of trapped oligomer. Removing the micelles reverses the diversion and the
#  on pathway resumes."""

DESC_ON = """#@description: |
#  On-pathway-only reduction of the Rana et al. (2020) competing-pathway model: the
#  reactions of Eq. 4 alone, which is the sub-model the paper fits to the fatty
#  acid-free experiment. Monomer adds stepwise to a growing oligomer
#  (A_i + A_1 <-> A_i+1, i = 1..11) up to the 12-mer nucleus A_12, identified with
#  the fibril F, which then takes up oligomers catalytically (A_i + F <-> F). All
#  kinetics are mass action; the thioflavin T signal is map_on*F.
#  With 25 uM monomer and no fatty acid the eleven sequential nucleation steps have
#  to be traversed before any fibril exists, so the fibril trace is flat for ~28 h,
#  rises steeply between 30 and 42 h as the chain empties into A_12, and plateaus
#  once free monomer is exhausted -- the classical nucleation-dependent sigmoid."""


def spec_list(net_full, net_on):
    ns_f, nr_f = net_full
    ns_o, nr_o = net_on
    S = {}
    S[f"{STEM}.bngl"] = dict(
        title="Amyloid-beta competing aggregation pathways (Rana et al., 2020)",
        description=DESC_FULL, off_pathway=True, spikes=False,
        params=dict(Mic_0="MIC", Abeta_0=25.0),
        actions=actions_simple(
            80, 1600, "Fig. 1c (red trace) in Rana et al. (2020)",
            dict(title="Off-pathway reference: 25 uM Abeta42 with 5 mM C12 fatty acid",
                 body="  #@protocol: |\n"
                      "  #  Pseudo-micelles are present from t = 0, so the off pathway is open\n"
                      "  #  throughout. Expect thioflavin T to rise without a lag phase and to\n"
                      "  #  settle near half the fibril plateau of the fatty acid-free\n"
                      "  #  experiment: monomer is diverted into trapped off-pathway oligomer\n"
                      "  #  before the on-pathway nucleation chain can fill."),
            ns_f, nr_f))
    S[f"{STEM}_on_pathway.bngl"] = dict(
        title="Amyloid-beta on-pathway aggregation (Rana et al., 2020)",
        description=DESC_ON, off_pathway=False, spikes=False,
        params=dict(Mic_0=0.0, Abeta_0=25.0),
        actions=actions_simple(
            50, 1000, "Fig. 1a in Rana et al. (2020)",
            dict(title="On pathway: 25 uM Abeta42, no fatty acid",
                 body="  #@protocol: |\n"
                      "  #  Fibril formation from monomer alone. Expect a flat lag phase to\n"
                      "  #  ~28 h, a steep rise between 30 and 42 h, and a plateau by 45 h."),
            ns_o, nr_o))
    for tag, t_ev, t_end, fig, colour in [
            ("micelle_addition_3h", 3, 80, "Fig. 1c (green trace) in Rana et al. (2020)", "green"),
            ("micelle_addition_24h", 24, 80, "Fig. 1c (blue trace) in Rana et al. (2020)", "blue")]:
        S[f"{STEM}_{tag}.bngl"] = dict(
            title=f"Amyloid-beta on-to-off pathway switching at {t_ev} h (Rana et al., 2020)",
            description=DESC_FULL, off_pathway=True, spikes=True,
            params=dict(Mic_0=0.0, Abeta_0=25.0),
            actions=actions_event(
                t_ev, t_end, 20 * t_ev, 20 * (t_end - t_ev), fig,
                dict(title=f"Pseudo-micelle addition event at {t_ev} h (on-to-off switching)",
                     body="  #@protocol: |\n"
                          f"  #  Phase 1: 25 uM Abeta42 alone for {t_ev} h, in the on-pathway lag\n"
                          "  #  phase.\n"
                          f"  #  Phase 2: 5 mM C12 fatty acid is added at {t_ev} h, together with\n"
                          "  #  fresh monomer, and the sample is followed to 80 h. Expect\n"
                          "  #  thioflavin T to climb immediately and without a lag once the\n"
                          "  #  pseudo-micelles open the off pathway, as Rana et al. (2020)\n"
                          "  #  describe for the micelle addition event."),
                ns_f, nr_f,
                ['setConcentration("Mic()","Mic_present")']))
    for tag, t_ev, t_end, fig in [
            ("micelle_removal_5h", 5, 50, "Fig. 1b (green trace) in Rana et al. (2020)"),
            ("micelle_removal_24h", 24, 50, "Fig. 1b (blue trace) in Rana et al. (2020)")]:
        S[f"{STEM}_{tag}.bngl"] = dict(
            title=f"Amyloid-beta off-to-on pathway switching at {t_ev} h (Rana et al., 2020)",
            description=DESC_FULL, off_pathway=True, spikes=True,
            params=dict(Mic_0="MIC", Abeta_0=25.0, map_off="MAP2"),
            actions=actions_event(
                t_ev, t_end, 20 * t_ev, 20 * (t_end - t_ev), fig,
                dict(title=f"Pseudo-micelle removal event at {t_ev} h (off-to-on switching)",
                     body="  #@protocol: |\n"
                          f"  #  Phase 1: 25 uM Abeta42 with 5 mM C12 fatty acid for {t_ev} h, so\n"
                          "  #  the off pathway runs and thioflavin T rises without a lag.\n"
                          f"  #  Phase 2: the sample is diluted at {t_ev} h to bring the fatty\n"
                          "  #  acid below its critical micelle concentration -- the free\n"
                          "  #  pseudo-micelle pool is emptied and fresh monomer is supplied --\n"
                          "  #  and followed to 50 h. Expect a sharp further rise as trapped\n"
                          "  #  off-pathway material and monomer re-enter the on pathway."),
                ns_f, nr_f,
                ['setConcentration("Mic()",0)']))
    return S


SPECS = spec_list((23, 80), (12, 44))


def build_file(name: str) -> str:
    """The exact text of one committed .bngl file."""
    global V
    s = SPECS[name]
    V = dict(DEFAULTS)
    V.update({k: (DEFAULTS["Mic_present"] if v == "MIC" else
                  REMOVAL_MAP_OFF if v == "MAP2" else v)
              for k, v in s["params"].items()})
    if not s["off_pathway"]:
        V.update(ON_PATHWAY_FIT)
    return build(s["title"], s["description"], COMMON_NOTE,
                 off_pathway=s["off_pathway"], spikes=s["spikes"], actions=s["actions"])


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    out.mkdir(parents=True, exist_ok=True)
    for name in SPECS:
        text = build_file(name)
        (out / name).write_text(text)
        bad = [(i + 1, len(ln)) for i, ln in enumerate(text.split("\n")) if len(ln) > 100]
        print(f"{name}: {len(text.splitlines())} lines"
              + (f"  !! lines >100 chars: {bad[:3]}" if bad else ""))


if __name__ == "__main__":
    main()
