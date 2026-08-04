#!/usr/bin/env python3
"""Generate the four BNGL encodings curated from Malleshaiah et al. (2010)."""

from pathlib import Path
from textwrap import wrap


ROOT = Path(__file__).resolve().parent
BASE = "ste5_fus3_ptc1_switch_malleshaiah2010"


def comment_block(text: str) -> str:
    """Wrap prose as structured-comment continuation lines."""
    return "\n".join(f"#  {line}" for line in wrap(text, width=76))


def wrapped_rule(label: str, lhs: str, rhs: str, rate: str) -> list[str]:
    line = f"  {label}: {lhs} -> {rhs} {rate}"
    if len(line) <= 100:
        return [line]
    return [f"  {label}: {lhs} -> \\", f"    {rhs} {rate}"]


def rules(n_sites: int) -> str:
    result: list[str] = []

    result += wrapped_rule(
        "bind_Fus3_dock",
        "Fus3(dock,cat)+Ste5(dock,cat)",
        "Fus3(dock!1,cat).Ste5(dock!1,cat)",
        "f1_K",
    )
    for n in range(n_sites + 1):
        result += wrapped_rule(
            f"release_Fus3_dock_p{n}",
            f"Fus3(dock!1,cat).Ste5(dock!1,cat,n~p{n})",
            f"Fus3(dock,cat)+Ste5(dock,cat,n~p{n})",
            f"b{n + 1}_K",
        )
    for n in range(n_sites):
        result += wrapped_rule(
            f"bind_Fus3_site_p{n}",
            f"Fus3(dock!1,cat).Ste5(dock!1,cat,n~p{n})",
            f"Fus3(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"{n_sites - n}*f2_K",
        )
        result += wrapped_rule(
            f"release_Fus3_site_p{n}",
            f"Fus3(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"Fus3(dock!1,cat).Ste5(dock!1,cat,n~p{n})",
            "b6_K",
        )
        result += wrapped_rule(
            f"release_Fus3_dock_from_site_p{n}",
            f"Fus3(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"Fus3(dock,cat!2).Ste5(dock,cat!2,n~p{n})",
            f"b{n + 1}_K",
        )
        result += wrapped_rule(
            f"rebind_Fus3_dock_p{n}",
            f"Fus3(dock,cat!1).Ste5(dock,cat!1,n~p{n})",
            f"Fus3(dock!2,cat!1).Ste5(dock!2,cat!1,n~p{n})",
            "f3_K",
        )
        result += wrapped_rule(
            f"phosphorylate_docked_Fus3_p{n}",
            f"Fus3(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"Fus3(dock!1,cat).Ste5(dock!1,cat,n~p{n + 1})",
            "k_K",
        )
        result += wrapped_rule(
            f"phosphorylate_site_only_Fus3_p{n}",
            f"Fus3(dock,cat!1).Ste5(dock,cat!1,n~p{n})",
            f"Fus3(dock,cat)+Ste5(dock,cat,n~p{n + 1})",
            "k_K",
        )

    result += wrapped_rule(
        "bind_Ptc1_dock",
        "Ptc1(dock,cat)+Ste5(dock,cat)",
        "Ptc1(dock!1,cat).Ste5(dock!1,cat)",
        "f1_P",
    )
    for n in range(n_sites + 1):
        result += wrapped_rule(
            f"release_Ptc1_dock_p{n}",
            f"Ptc1(dock!1,cat).Ste5(dock!1,cat,n~p{n})",
            f"Ptc1(dock,cat)+Ste5(dock,cat,n~p{n})",
            "b1_P",
        )
    for n in range(1, n_sites + 1):
        result += wrapped_rule(
            f"bind_Ptc1_site_p{n}",
            f"Ptc1(dock!1,cat).Ste5(dock!1,cat,n~p{n})",
            f"Ptc1(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"{n}*f2_P",
        )
        result += wrapped_rule(
            f"release_Ptc1_site_p{n}",
            f"Ptc1(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"Ptc1(dock!1,cat).Ste5(dock!1,cat,n~p{n})",
            "b2_P",
        )
        result += wrapped_rule(
            f"release_Ptc1_dock_from_site_p{n}",
            f"Ptc1(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"Ptc1(dock,cat!2).Ste5(dock,cat!2,n~p{n})",
            "b1_P",
        )
        result += wrapped_rule(
            f"rebind_Ptc1_dock_p{n}",
            f"Ptc1(dock,cat!1).Ste5(dock,cat!1,n~p{n})",
            f"Ptc1(dock!2,cat!1).Ste5(dock!2,cat!1,n~p{n})",
            "f3_P",
        )
        result += wrapped_rule(
            f"dephosphorylate_docked_Ptc1_p{n}",
            f"Ptc1(dock!1,cat!2).Ste5(dock!1,cat!2,n~p{n})",
            f"Ptc1(dock!1,cat).Ste5(dock!1,cat,n~p{n - 1})",
            "k_P",
        )
        result += wrapped_rule(
            f"dephosphorylate_site_only_Ptc1_p{n}",
            f"Ptc1(dock,cat!1).Ste5(dock,cat!1,n~p{n})",
            f"Ptc1(dock,cat)+Ste5(dock,cat,n~p{n - 1})",
            "k_P",
        )

    for n in range(n_sites):
        result += wrapped_rule(
            f"bind_active_Fus3_p{n}",
            f"Fus3a(cat)+Ste5(dock,cat,n~p{n})",
            f"Fus3a(cat!1).Ste5(dock,cat!1,n~p{n})",
            f"{n_sites - n}*f4_K",
        )
        result += wrapped_rule(
            f"release_active_Fus3_p{n}",
            f"Fus3a(cat!1).Ste5(dock,cat!1,n~p{n})",
            f"Fus3a(cat)+Ste5(dock,cat,n~p{n})",
            "b6_K",
        )
        result += wrapped_rule(
            f"phosphorylate_active_Fus3_p{n}",
            f"Fus3a(cat!1).Ste5(dock,cat!1,n~p{n})",
            f"Fus3a(cat)+Ste5(dock,cat,n~p{n + 1})",
            "k_K",
        )
    return "\n".join(result)


def model_text(n_sites: int, molar_units: bool) -> str:
    suffix_bits = []
    if molar_units:
        suffix_bits.append("molar-association-units interpretation")
    if n_sites == 1:
        suffix_bits.append("one-site variant")
    variant = "; ".join(suffix_bits) or "four-site primary"
    states = "~".join(f"p{i}" for i in range(n_sites + 1))
    p_obs = "\n".join(
        f"  Molecules Ste5_p{i} Ste5(n~p{i})" for i in range(n_sites + 1)
    )
    numerator = "+".join(
        ("" if i == 1 else f"{i}*") + f"Ste5_p{i}" for i in range(1, n_sites + 1)
    )
    assoc = {
        "f1_P": 186000 if not molar_units else 0.000186,
        "f1_K": 12000 if not molar_units else 0.000012,
        "f4_K": 109000 if not molar_units else 0.000109,
    }
    discrepancy = (
        "This primary reading preserves Supplementary Table 3 literally: the three "
        "bimolecular association constants are treated as /nM/s. It does not reproduce "
        "the Fig. 3b steady-state switch; see the verification notebook and the "
        "molar-association-units variant. The displayed Fus3_Ste5n_2 loss term "
        "uses k(P); it is encoded with k(K) because both matching gain terms and "
        "Fig. 23b use k(K), while mixing the two would violate conservation."
        if not molar_units
        else
        "This bracketed reading treats the three association-unit labels in Supplementary "
        "Table 3 as /M/s and converts them to /nM/s. It is physically distinct from the "
        "literal table reading and is shipped because the reported units and Fig. 3b cannot "
        "both be recovered from the supplied equations and values. The displayed "
        "Fus3_Ste5n_2 loss term uses k(P); it is encoded with k(K) because both matching "
        "gain terms and Fig. 23b use k(K), while mixing the two would violate conservation."
    )
    variant_comment = comment_block(
        f"This is the {variant}. A Ste5 internal state stores only the number of "
        "phosphorylated identical sites, exactly as in the 38-equation (four-site) "
        "or 14-equation (one-site) state aggregation in the supplement. The model "
        "is concentration-based in nM and is intended for deterministic ODE use. "
        "The source gives no reactor volume; V_ref is a conversion-only nominal "
        "yeast-cell volume and does not enter the ODEs."
    )
    discrepancy_comment = comment_block(discrepancy)
    network_size = "38 species, 80 reactions" if n_sites == 4 else "14 species, 23 reactions"
    return f'''begin model

#@title: Ste5-Fus3-Ptc1 phosphorylation switch (Malleshaiah et al., 2010)

#@description: |
#  Deterministic concentration model of the mating-decision circuit formed by
#  the Ste5 scaffold, Fus3 kinase, and Ptc1 phosphatase in budding yeast. Fus3
#  and Ptc1 compete for a docking site on Ste5 and then bind a catalytic site;
#  either docking bond may break while the enzyme remains catalytic-site-bound.
#  Fus3 docking affinity increases with Ste5 phosphorylation. Alpha-factor
#  recruits Ptc1 and creates a small active cytosolic Fus3 pool through fitted
#  Hill functions. Ste5 carries {n_sites} identical phosphorylation
#  site{'s' if n_sites != 1 else ''}.
#
{variant_comment}

#@keyword: |
#  yeast mating, Ste5, Fus3, Ptc1, scaffold, multisite phosphorylation,
#  zero-order ultrasensitivity, steric hindrance, enzyme competition,
#  two-stage binding, rule-based model, ordinary differential equations

#@reference: |
#  Malleshaiah MK, Shahrezaei V, Swain PS, Michnick SW (2010). The scaffold
#  protein Ste5 directly controls a switch-like mating decision in yeast.
#  Nature 465:101-105. doi:10.1038/nature08946

#@note: |
{discrepancy_comment}

begin parameters

  # Population-conversion constants. V_ref is a nominal 40 fL yeast-cell
  # volume used only to make this concentration model conversion-ready.
  NA 6.02214076e23  # molecules/mol
  V_ref 4e-14       # L/cell

  # Alpha-factor input and protein concentrations (Supplementary Table 3)
  F 0.2                 # dimensionless; PCA calibration, unused in Fig. 3b
  alpha 1             # nM
  Ste5_tot 52         # nM
  Fus3_tot 197        # nM
  Ptc1_max 39         # nM
  EC50_P 240          # nM
  nH_P 2.3            # dimensionless
  Ptc1_0 1.2          # nM
  Fus3active_max 5.8  # nM
  EC50_K 1680         # nM
  nH_K 1.3            # dimensionless

  # Ptc1 two-stage docking, catalytic-site binding, and catalysis
  f1_P {assoc['f1_P']:.12g}  # /nM/s
  f2_P 327             # /s
  f3_P 0.3             # /s
  b1_P 22              # /s
  b2_P 0.12            # /s
  k_P 0.5              # /s

  # Fus3 two-stage docking, catalytic-site binding, and catalysis
  f1_K {assoc['f1_K']:.12g}  # /nM/s
  f2_K 850             # /s
  f3_K 0.1             # /s
  f4_K {assoc['f4_K']:.12g}  # /nM/s
  b1_K 99              # /s; Ste5 p0
  b2_K 42              # /s; Ste5 p1
  b3_K 21              # /s; Ste5 p2
  b4_K 13              # /s; Ste5 p3
  b5_K 10              # /s; Ste5 p4 (unused in one-site variant)
  b6_K 24              # /s
  k_K 1.13             # /s

  # Alpha-factor-dependent inputs (Supplementary Information, Eqs. on p. 17)
  Ptc1_input Ptc1_0+Ptc1_max*alpha^nH_P/(alpha^nH_P+EC50_P^nH_P)  # nM
  Fus3active_input Fus3active_max*alpha^nH_K/(alpha^nH_K+EC50_K^nH_K)  # nM
  Fus3inactive_input Fus3_tot-Fus3active_input  # nM

end parameters

begin molecule types

  Fus3(dock,cat)
  Fus3a(cat)
  Ptc1(dock,cat)
  Ste5(dock,cat,n~{states})

end molecule types

begin seed species

  Fus3(dock,cat) Fus3inactive_input  # nM
  Fus3a(cat) Fus3active_input        # nM
  Ptc1(dock,cat) Ptc1_input          # nM
  Ste5(dock,cat,n~p0) Ste5_tot       # nM

end seed species

begin observables

  # Total Ste5 in each phosphorylation-count state, including complexes
{p_obs}

  # Total Fus3-Ste5 and Ptc1-Ste5 complexes through either binding stage
  Molecules Fus3_docked Fus3(dock!+)
  Molecules Fus3_site_only Fus3(dock,cat!+)
  Molecules Ptc1_docked Ptc1(dock!+)
  Molecules Ptc1_site_only Ptc1(dock,cat!+)

end observables

begin functions

  # Mean number of phosphorylated sites per Ste5 molecule
  Mean_pSites() ({numerator})/Ste5_tot

  # PCA-model complex totals: docked plus catalytic-site-only complexes
  Fus3_Ste5() Fus3_docked+Fus3_site_only
  Ptc1_Ste5() Ptc1_docked+Ptc1_site_only

end functions

begin reaction rules

{rules(n_sites)}

end reaction rules

end model

begin actions

  generate_network({{overwrite=>1}})  # {network_size}

  # Alpha-factor dose-response at steady state

  #@protocol: |
  #  Starting from unphosphorylated, unbound Ste5 at every scan point, vary
  #  alpha-factor from 0.001 to 10 uM (1 to 10000 nM), integrate the
  #  deterministic mass-action equations, and report mean Ste5
  #  phosphorylation together with Fus3-Ste5 and Ptc1-Ste5 complexes.

  #@figure: Fig. 3b in Malleshaiah et al. (2010)

  parameter_scan({{method=>"ode",suffix=>"scan",parameter=>"alpha",\\
    par_min=>1,par_max=>10000,n_scan_pts=>41,log_scale=>1,t_start=>0,\\
    t_end=>100000,n_steps=>200,steady_state=>1,reset_conc=>1,\\
    print_functions=>1,atol=>1e-10,rtol=>1e-8}})

end actions
'''


def filename(n_sites: int, molar_units: bool) -> str:
    suffix = ""
    if molar_units:
        suffix += "_molar_association_units"
    if n_sites == 1:
        suffix += "_1ps"
    return BASE + suffix + ".bngl"


def main() -> None:
    for molar_units in (False, True):
        for n_sites in (4, 1):
            path = ROOT / filename(n_sites, molar_units)
            path.write_text(model_text(n_sites, molar_units), encoding="utf-8")
            print(path.name)


if __name__ == "__main__":
    main()
