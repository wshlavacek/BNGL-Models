"""Emit the Dushek et al. (2011) two-step model as house-style BNGL, for any site count.

The number of phosphorylation sites N is *structural* — it is N components on the substrate
molecule — so a family over N cannot be a `parameter_scan` and would otherwise be twenty
near-identical files. This generator is the `lambda_switch_arkin1998` precedent: the three
committed `.bngl` files are its output at N = 20, and the verification notebook calls it for
N = 1 ... 20 to reproduce all eighty published curves of Fig. 2.

Usage
-----
    python build_dushek2011.py --list
    python build_dushek2011.py --check    # regenerate the committed files, exit 1 on drift
    python build_dushek2011.py --n 6 --refractory --out /tmp/n6.bngl

Representation
--------------
The N sites are written as N *identical* `y~U~P` components, not as one counter site with
N+1 states. BioNetGen canonicalises a species up to permutation of identical components, so
both give the same 7N+3 species and 14N+4 reactions (143 / 284 at N = 20, the count the
paper states), but only the explicit form lets BioNetGen supply the site multiplicity: the
generated network carries `k_on_star, 2*k_on_star, ... 20*k_on_star` where the authors' own
SI file hard-codes `Ekf1 = 20*10/A, Ekf2 = 19*10/A, ...`. That keeps the on-rate a
single-site constant, as `skills/bngl/skill.md` §1.3.1 requires.

One trap this buys: an observable that names every site, `S(y~P,y~P,...)` twenty times over,
is a twenty-component pattern and BioNetGen will not finish matching it. It is never needed —
`Molecules Obs_pY S(y~P!?)` counts phosphorylated sites directly, which is the numerator of
the paper's Eq. 2.

Enzymes are deliberately not molecules. The paper works in the limit where kinase and
phosphatase are in excess of substrate, so E and F enter only through the pseudo-first-order
encounter rates `ET*k_plus` and `FT*k_plus`; the authors' file has no enzyme molecule type
either. Making them explicit would change the model, because a refractory enzyme that
diffuses away is replaced instantly from an infinite pool.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
STEM = "membrane_anchored_multisite_phosphorylation_dushek2011"

CITATION = (
    "Dushek O, van der Merwe PA, Shahrezaei V (2011). Ultrasensitivity in multisite\n"
    "#  phosphorylation of membrane-anchored proteins. Biophysical Journal 100:1189-1197.\n"
    "#  doi:10.1016/j.bpj.2011.01.060"
)

# The two kinetic regimes of Fig. 2. `k_on` for the reaction-limited regime is NOT the
# main text's k+/100 -- see the `#@note:` the generator writes into every file, and the
# verification notebook, which measures it off the published panel.
REGIMES = {
    "diffusion": dict(k_plus=0.1, k_on=10.0),
    "reaction": dict(k_plus=10.0, k_on=1.0),
}


def header(n_sites: int, regime: str, refractory: bool) -> str:
    panels = {
        ("diffusion", True): "Fig. 2 D", ("reaction", True): "Fig. 2 C",
        ("diffusion", False): "Fig. 2 B", ("reaction", False): "Fig. 2 A",
    }[(regime, refractory)]
    species, rxns = 7 * n_sites + 3, 14 * n_sites + 4
    refr_txt = (
        "After catalysis the enzyme is left inactive for a brief refractory period (1/mu)\n"
        "#  before it can act again, which is what makes local rebinding stop short of full\n"
        "#  processivity."
        if refractory else
        "\n#  This file has NO refractory state: after catalysis the enzyme is immediately\n"
        "#  able to act again. It is the paper's own control on the mechanism, not a limit of the\n"
        "#  primary file -- the a~1 state and its two reactivation rules simply do not exist,\n"
        "#  so the network is smaller."
    )
    behaviour = {
        ("diffusion", True): (
            "Ultrasensitivity: the dose-response of total phosphorylation against [E]/[F]\n"
            "#  steepens with the number of sites, and the Hill number grows roughly linearly\n"
            "#  in N (the published inset fits Hill = 0.58 + 0.24 N, giving 5.38 at N = 20).\n"
            "#  This is the paper's headline result."),
        ("reaction", True): (
            "No ultrasensitivity: the Hill number sits at ~1 for every N (published inset fit\n"
            "#  Hill = 0.98 + 0 N). Enzymes are far from the zero-order regime, so adding sites\n"
            "#  does not steepen the response."),
        ("diffusion", False): (
            "Local saturation makes the encounter effectively processive, so a single\n"
            "#  encounter drives many modifications and the dose-response collapses onto the\n"
            "#  single-site curve -- nearly flat at 0.5 across four decades of [E]/[F]."),
        ("reaction", False): (
            "The classical well-mixed result: Hill number ~1 for every N (published inset fit\n"
            "#  Hill = 0.99 + 0 N). The single-site substrate is *subsensitive*, spanning only\n"
            "#  0.47 to 0.53, and multiple sites reduce that subsensitivity."),
    }[(regime, refractory)]
    title = ("Membrane-anchored multisite phosphorylation"
             if refractory else
             "Multisite phosphorylation without enzyme inactivation")
    return f"""begin model

#@title: {title} (Dushek et al., 2011)

#@description: |
#  A membrane-anchored substrate carrying {n_sites} phosphorylation
#  site{"s" if n_sites > 1 else ""}, 
#  modified by a kinase and dephosphorylated by a phosphatase that are both confined to the
#  same membrane and both in excess of the substrate. The paper's two-step scheme: enzyme and
#  substrate first diffuse into an ENCOUNTER COMPLEX at the diffusion-limited rate X*k_plus,
#  and only inside it can they bind (at the local on-rate, lambda*k_on_star) or move apart
#  again (k_minus). {refr_txt}
#
#  How it works: every reaction is a first-order transition of one substrate molecule, so with
#  the enzymes held in excess the whole model is a linear continuous-time Markov chain on
#  {species} states ({rxns} reactions at N = {n_sites}). The site multiplicity lambda is not
#  written into any rate constant -- the {n_sites} sites are {n_sites} identical components
#  and BioNetGen supplies the
#  factors 1..{n_sites} itself.
#
#  What it does: {behaviour}

#@keyword: |
#  signalling, multisite phosphorylation, ultrasensitivity, membrane diffusion,
#  encounter complex, kinetic proofreading, mass action, T-cell receptor,
#  noncatalytic tyrosine-phosphorylated receptor, Homo sapiens

#@reference: |
#  {CITATION}

#@note: |
#  SIMULATION INTENT. Population-based: every amount is a molecule count on a membrane patch
#  of area A_sim = A_ref*f, and at the default A_ref = 1 um^2, f = 1 the counts are
#  numerically equal to the surface densities (um^-2) the paper quotes. The enzymes are NOT
#  molecules -- the paper works in the enzyme-excess limit, so [E] and [F] enter only through
#  the pseudo-first-order encounter rates ET*k_plus and FT*k_plus, exactly as in the authors'
#  own BioNetGen file (Supporting Material pp. 10-12).

#@note: |
#  SITE MULTIPLICITY IS STRUCTURAL, NOT A CONSTANT. The authors' SI file uses one counter site
#  with {n_sites}+1 states and hard-codes the multiplicity into 40 constants (Ekf1 = 20*10/A,
#  Ekf2 = 19*10/A, ...). This file uses {n_sites} identical y components instead. BioNetGen
#  canonicalises species up to permutation of identical components, so the two give the SAME
#  network -- {species} species and {rxns} reactions, the "143 coupled ODEs" the paper reports at
#  N = 20 -- but here the generated rate laws read k_on_star, 2*k_on_star, ... {n_sites}*k_on_star,
#  supplied by BioNetGen rather than by hand, so k_on_star stays a single-site constant
#  (skills/bngl/skill.md §1.3.1).

#@note: |
#  REACTION-LIMITED k_on IS k_plus/10, NOT THE PUBLISHED k_plus/100. The main text defines the
#  reaction-limited regime as "k_on = k+/100", which at k_plus = 10 um^2/s gives k_on = 0.1.
#  That value does not reproduce the published panel: it puts the N = 1 curve of Fig. 2 A at
#  0.313 -> 0.687 across the plotted range, whereas the panel itself measures 0.476 at
#  log10([E]/[F]) = -2. k_on = 1 um^2/s gives 0.471, agreeing to 0.005 of full scale. The
#  Hill-number insets do not discriminate (both give ~1), so the curve shape is the evidence.
#  The digitized measurement is in verify_dushek2011.ipynb. The diffusion-limited regime is
#  unaffected: "k_on = 100*k+" at k_plus = 0.1 gives k_on = 10, which is also the value in the
#  authors' SI file, and it reproduces the Fig. 2 D Hill numbers (5.387 against a published
#  fit of 5.38 at N = 20).

#@figure: {panels} in Dushek et al. (2011)
"""


def build(n_sites: int, regime: str, refractory: bool) -> str:
    p = REGIMES[regime]
    y_types = ",".join(["y~U~P"] * n_sites)
    y_unphos = ",".join(["y~U"] * n_sites)

    def wrap(prefix: str, body: str, width: int = 96) -> str:
        """Wrap a long component list onto continuation lines under the 100-char limit."""
        out, line = [], prefix
        for i, part in enumerate(body.split(",")):
            add = part + ("," if i < len(body.split(",")) - 1 else "")
            if len(line) + len(add) > width:
                out.append(line + "\\")
                line = add
            else:
                line += add
        out.append(line)
        return "\n".join(out)

    refr_params = "" if not refractory else """
  # Reactivation rate of a refractory enzyme; Fig. 2 C,D caption, mu = 1 /s
  mu         1            # /s
"""
    refr_rules = "" if not refractory else """
  # An enzyme is left refractory by catalysis and recovers in place at mu. This is the
  # mechanism the paper isolates: it caps local rebinding short of full processivity.
  R_react_E: S(e~1,a~1) -> S(e~1,a~0)                              mu
  R_react_F: S(f~1,a~1) -> S(f~1,a~0)                              mu
"""
    cat_E = "S(e~1!1,a~0,y~U!1) -> S(e~1,a~1,y~P)" if refractory else \
            "S(e~1!1,a~0,y~U!1) -> S(e~1,a~0,y~P)"
    cat_F = "S(f~1!1,a~0,y~P!1) -> S(f~1,a~1,y~U)" if refractory else \
            "S(f~1!1,a~0,y~P!1) -> S(f~1,a~0,y~U)"
    escape_refr = "" if not refractory else """
  # A refractory enzyme may also simply diffuse out of the encounter complex; the substrate
  # then meets a fresh active enzyme from the excess pool, so `a` resets.
  R_esc_E:   S(e~1,f~0,a~1) -> S(e~0,f~0,a~0)                      k_minus
  R_esc_F:   S(e~0,f~1,a~1) -> S(e~0,f~0,a~0)                      k_minus
"""
    a_states = "a~0~1," if refractory else "a~0,"

    return header(n_sites, regime, refractory) + f"""
begin parameters
  # --- Conversion readiness (skills/bngl/skill.md §1.3). No volumetric species exist; the
  # model lives on a membrane, so the reference quantity is an AREA.
  NA         6.02214076e23   # molecules/mol
  A_ref      1               # um^2/patch
  f          1               # dimensionless
  A_sim      A_ref*f         # um^2/patch
  h_mem      1e-8            # m
  V_ref      A_sim*h_mem*1e-9  # L/patch

  # --- Substrate
  # Total substrate density; SI file `tS 1`, SI Fig. S1 caption S_T = 1 um^-2
  tS         1            # molecules/um^2

  # --- Encounter geometry
  # Binding radius s and encounter area A = pi*s^2; main text p. 1191 and Table 1.
  # The Fig. 2 caption prints A = 0.01^2 um^2; `pdftotext` renders that as "0.012".
  s_bind     0.0056       # um
  A_enc      1e-4         # um^2

  # --- Two-step kinetics ({regime}-limited regime; Fig. 2 caption and Table 1)
  # Diffusion-limited on-rate k+ (macroscopic)
  k_plus     {p['k_plus']}          # um^2/s
  # Bimolecular on-rate k_on (macroscopic), SINGLE SITE
  k_on       {p['k_on']}          # um^2/s
  # Local rates: k- = k+/A and k*on = k_on/A (Table 1)
  k_minus    k_plus/A_enc  # /s
  k_on_star  k_on/A_enc    # /s
  # Unbinding and modification rates; Fig. 2 caption
  koff       1            # /s
  kr         0.1          # /s
{refr_params}
  # --- Enzyme densities. The scanned quantity is the ratio r = [E]/[F] at fixed total.
  # EF_tot is stated only in the SI Fig. S1 caption (E_T + F_T = 1000 um^-2).
  EF_tot     1000         # molecules/um^2
  r          1            # dimensionless
  ET         EF_tot*r/(1+r)  # molecules/um^2
  FT         EF_tot/(1+r)    # molecules/um^2
  # Pseudo-first-order encounter rates X*k+ (main text, Eq. 1)
  kp_E       ET*k_plus    # /s
  kp_F       FT*k_plus    # /s
end parameters

begin molecule types
  # The membrane-anchored substrate. `e`/`f` record whether a kinase / phosphatase is inside
  # the encounter complex, `a` whether that enzyme is active or refractory, and the {n_sites}
  # identical `y` components are the phosphorylation sites. A bond between `e` (or `f`) and a
  # `y` is the bound enzyme-substrate complex, which is why no separate "bound" site is needed.
{wrap("  S(" + a_states + "e~0~1,f~0~1,", y_types + ")")}
end molecule types

begin seed species
  # All substrate starts unphosphorylated with no enzyme in the encounter complex.
{wrap("  S(" + ("a~0," if refractory else "a~0,") + "e~0,f~0,", y_unphos + ")  tS")}   # molecules
end seed species

begin observables
  # Numerator of the paper's Eq. 2: the number of phosphorylated sites, counted over every
  # substrate state, bound or not. `!?` matches a site whether or not an enzyme is on it.
  Molecules Obs_pY        S(y~P!?)
  # Substrate with a kinase / phosphatase inside the encounter complex, active or refractory
  Molecules Obs_enc_E     S(e~1)
  Molecules Obs_enc_F     S(f~1)
  # Substrate actually bound to an enzyme
  Molecules Obs_bound_E   S(e~1!1,y~U!1)
  Molecules Obs_bound_F   S(f~1!1,y~P!1)
end observables

begin reaction rules
  # --- Encounter complex (Eq. 1, left half). X*k+ in, k- out.
  R_enc_E:   S(e~0,f~0,a~0) <-> S(e~1,f~0,a~0)                     kp_E, k_minus
  R_enc_F:   S(e~0,f~0,a~0) <-> S(e~0,f~1,a~0)                     kp_F, k_minus
{escape_refr}
  # --- Binding inside the encounter complex (Eq. 1, right half). The rate constant is the
  # SINGLE-SITE local on-rate; BioNetGen multiplies by the number of matching sites, which is
  # the paper's lambda = N - j for the kinase and lambda = j for the phosphatase.
  R_bind_E:  S(e~1,a~0,y~U) <-> S(e~1!1,a~0,y~U!1)                 k_on_star, koff
  R_bind_F:  S(f~1,a~0,y~P) <-> S(f~1!1,a~0,y~P!1)                 k_on_star, koff

  # --- Catalysis. Distributive: the enzyme releases the site and stays in the encounter
  # complex, so it may rebind locally rather than diffusing away.
  R_cat_E:   {cat_E}   kr
  R_cat_F:   {cat_F}   kr
{refr_rules}end reaction rules

end model
"""


def actions(n_sites: int, regime: str, refractory: bool) -> str:
    species, rxns = 7 * n_sites + 3, 14 * n_sites + 4
    panel = {("diffusion", True): "2 D", ("reaction", True): "2 C",
             ("diffusion", False): "2 B", ("reaction", False): "2 A"}[(regime, refractory)]
    other = "reaction" if regime == "diffusion" else "diffusion"
    q = REGIMES[other]
    other_panel = {("diffusion", True): "2 D", ("reaction", True): "2 C",
                   ("diffusion", False): "2 B", ("reaction", False): "2 A"}[(other, refractory)]
    return f"""
begin actions
  generate_network({{overwrite=>1}})  # {species} species, {rxns} reactions

  # Dose-response against the kinase/phosphatase balance

  #@protocol: |
  #  Steady-state total phosphorylation as [E]/[F] is scanned over four decades at fixed
  #  [E]+[F]. Each scan point is integrated to steady state, which is what the paper plots;
  #  the system is a linear Markov chain so the steady state is unique and independent of the
  #  initial condition.

  #@figure: Fig. {panel} in Dushek et al. (2011)

  parameter_scan({{method=>"ode",parameter=>"r",par_min=>0.01,par_max=>100,n_scan_pts=>41,\\
    log_scale=>1,t_start=>0,t_end=>1e5,n_steps=>10,steady_state=>1,\\
    atol=>1e-12,rtol=>1e-10,suffix=>"scan"}})

  # The other kinetic regime of the same figure row

  #@protocol: |
  #  The {other}-limited regime is the same network with k_plus = {q['k_plus']} um^2/s and
  #  k_on = {q['k_on']} um^2/s. Switching both reproduces Fig. {other_panel}.

  #@figure: Fig. {other_panel} in Dushek et al. (2011)

  setParameter("k_plus",{q['k_plus']})
  setParameter("k_on",{q['k_on']})
  parameter_scan({{method=>"ode",parameter=>"r",par_min=>0.01,par_max=>100,n_scan_pts=>41,\\
    log_scale=>1,t_start=>0,t_end=>1e5,n_steps=>10,steady_state=>1,\\
    atol=>1e-12,rtol=>1e-10,suffix=>"scan_{other}"}})
end actions
"""


COMMITTED = {
    f"{STEM}.bngl": (20, "diffusion", True),
    f"{STEM}_nonrefractory.bngl": (20, "diffusion", False),
}


def emit(n_sites: int, regime: str, refractory: bool, with_actions: bool = True) -> str:
    txt = build(n_sites, regime, refractory)
    return txt + (actions(n_sites, regime, refractory) if with_actions else "")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n", type=int, default=20, help="number of phosphorylation sites")
    ap.add_argument("--regime", choices=sorted(REGIMES), default="diffusion")
    ap.add_argument("--refractory", action="store_true", default=None)
    ap.add_argument("--no-refractory", dest="refractory", action="store_false")
    ap.add_argument("--no-actions", action="store_true")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--list", action="store_true", help="list the committed files")
    ap.add_argument("--check", action="store_true",
                    help="regenerate the committed files in place and report drift")
    a = ap.parse_args()

    if a.list:
        for name, (n, reg, refr) in COMMITTED.items():
            print(f"{name}: N={n}, {reg}-limited, refractory={refr}")
        return 0

    if a.check:
        drift = 0
        for name, (n, reg, refr) in COMMITTED.items():
            want = emit(n, reg, refr)
            p = FOLDER / name
            have = p.read_text() if p.exists() else ""
            if have != want:
                p.write_text(want)
                print(f"REGENERATED {name}")
                drift += 1
            else:
                print(f"ok          {name}")
        return 1 if drift else 0

    refr = True if a.refractory is None else a.refractory
    txt = emit(a.n, a.regime, refr, with_actions=not a.no_actions)
    if a.out:
        a.out.write_text(txt)
        print(f"wrote {a.out}")
    else:
        sys.stdout.write(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
