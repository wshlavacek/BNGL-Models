"""Rukhlenko-2022 slugs checked against the authors' own files stay in step with that check.

Each slug with an `upstream.json` pins one file of github.com/OleksiiR/cSTAR_Nature, which is
not vendored (the repository grants no licence), and commits the comparison its
`compare_upstream.py` wrote, `upstream_comparison.json`. Running that comparison needs
BioNetGen and the authors' file; these tests need neither. They guard that the pin, the slug's
`.bngl` and the committed comparison agree, so an edit to the `.bngl` that is not re-checked
fails here, and they unit-test the pieces of `cstar_upstream.py` the comparison rests on.
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import pytest

PAPER = Path("pybnf-jobs/Rukhlenko-2022")
SLUGS = sorted(p.parent.name for p in PAPER.glob("*/upstream.json"))


@pytest.fixture(scope="module")
def cu():
    sys.path.insert(0, str(PAPER))
    import cstar_upstream

    return cstar_upstream


def load(slug, name):
    return json.loads((PAPER / slug / name).read_text(encoding="utf-8"))


def test_a_slug_is_pinned():
    assert SLUGS


@pytest.mark.parametrize("slug", SLUGS)
def test_pin_is_well_formed_and_not_vendored(cu, slug):
    pin = load(slug, "upstream.json")
    info = pin["file"]
    assert pin["pinned_commit"] == cu.COMMIT
    assert re.fullmatch(r"[0-9a-f]{64}", info["sha256"])
    assert info["raw_url"] == cu.RAW_URL.format(commit=cu.COMMIT, path=info["upstream_path"])
    assert pin["license"]["consequence"].startswith("not vendored")
    assert not (PAPER / slug / Path(info["upstream_path"]).name).exists()


@pytest.mark.parametrize("slug", SLUGS)
def test_comparison_is_of_the_committed_files(slug):
    pin, report = load(slug, "upstream.json"), load(slug, "upstream_comparison.json")
    bngl = PAPER / slug / f"{slug}.bngl"
    assert report["job_file"]["sha256"] == hashlib.sha256(bngl.read_bytes()).hexdigest(), (
        f"{bngl.name} changed since upstream_comparison.json was written; re-run "
        f"{PAPER / slug / 'compare_upstream.py'}"
    )
    assert report["authors_file"]["sha256"] == pin["file"]["sha256"]
    assert report["authors_file"]["commit"] == pin["pinned_commit"]
    assert report["authors_file"]["path"] == pin["file"]["upstream_path"]


@pytest.mark.parametrize("slug", SLUGS)
def test_comparison_passed_every_assertion(slug):
    report = load(slug, "upstream_comparison.json")
    net = report["network"]
    assert net["species"]["bijection"] and net["reactions"]["bijection"]
    checks = [
        v for cond in net["conditions"].values() for v in cond.values() if isinstance(v, dict)
    ]
    for phase in report["trajectories"].values():
        checks += [phase["species_max_rel_to_peak"], phase["observables_max_rel_to_peak"]]
        checks += list(phase.get("gated_species", {}).values())
    assert checks and all(c["pass"] for c in checks)
    assert all(c["value"] <= c["tolerance"] for c in checks if "value" in c)
    assert all(c["pass"] for c in net["documented_differences"].values())
    assert report["verdict"] == "PASS"


def test_canonical_species_ignores_molecule_order_and_bond_numbers(cu):
    a = cu.canonical_species("TRK(Dim!1,Lig!2,Y~P).TRK(Dim!1,Lig,Y~0).TRKl(Rec!2)")
    b = cu.canonical_species("TRKl(Rec!7).TRK(Y~0,Lig,Dim!3).TRK(Lig!7,Y~P,Dim!3)")
    assert a == b
    other = cu.canonical_species("TRK(Dim!1,Lig!2,Y~0).TRK(Dim!1,Lig,Y~P).TRKl(Rec!2)")
    assert other != a
    assert cu.canonical_species("$IERK(ERKBD)") == "$IERK(ERKBD)"


def test_expressions_follow_bionetgen_syntax(cu):
    expr = "if(x<1,2^3,sqrt(4))+if((x<2),0,1)"
    assert cu.evaluate(expr, {"x": 0.5}) == 8.0
    assert cu.evaluate(expr, {"x": 3.0}) == 3.0


NET = """begin parameters
    1 D       1.3122028733013222  # Constant
    2 S0      -4.3636308697015025  # Constant
    3 k       2*D  # ConstantExpression
end parameters
begin species
    1 X() 0
    2 Y() S0
    3 $G() D
end species
begin reactions
    1 1 2 0.5*k #_R1
end reactions
begin groups
    1 Xo                   1
    2 Go                   2*3
end groups
"""


def test_net_reader_and_equations(cu):
    net = cu.parse_net("t.net", text=NET.replace("0.5*k", "0.5*k"))
    assert net.params["k"] == 2 * 1.3122028733013222
    assert list(net.fixed) == [False, False, True]
    state = cu.State(net, [2.0, 1.0, 3.0])
    assert state.names["Go"] == 6.0
    dydt, gross = state.rhs()
    flux = 0.5 * net.params["k"] * 2.0
    assert list(dydt) == [-flux, flux, 0.0] and list(gross) == [flux, flux, 0.0]


def test_translation_keeps_the_authors_expressions(cu):
    authors = cu.parse_net("a.net", text=NET)
    job = cu.parse_net("j.net", text=NET.replace("1.3122028733013222", "0"))
    gates = {"G()": cu.Gate(params=lambda value, text: {"D": text}, keep_species=True)}
    calls = [
        ("setConcentration", '"X()","D"'),
        ("setConcentration", '"Y()","2*S0"'),
        ("setConcentration", '"G()","2*D"'),
        ("writeSBML", ""),
        ("simulate", '{suffix=>"a",t_end=>1}'),
        ("saveConcentrations", ""),
        ("setConcentration", '"G()",0'),
        ("simulate", '{suffix=>"b",t_end=>1}'),
        ("resetConcentrations", ""),
        ("simulate", '{suffix=>"c",t_end=>1}'),
    ]
    out, gated = cu.translate_actions(calls, authors, job, gates)
    assert out == [
        ("setConcentration", '"X()","(1.3122028733013222)"'),  # D differs: written out
        ("setConcentration", '"Y()","2*S0"'),  # S0 the same: kept
        ("setConcentration", '"G()","2*(1.3122028733013222)"'),
        ("setParameter", '"D","2*(1.3122028733013222)"'),
        ("simulate", '{suffix=>"a",t_end=>1}'),
        ("saveConcentrations", ""),
        ("setConcentration", '"G()",0'),
        ("setParameter", '"D","0"'),
        ("simulate", '{suffix=>"b",t_end=>1}'),
        ("resetConcentrations", ""),
        ("setParameter", '"D","2*(1.3122028733013222)"'),
        ("simulate", '{suffix=>"c",t_end=>1}'),
    ]
    assert gated == {
        "a": {"G()": 2 * 1.3122028733013222},
        "b": {"G()": 0.0},
        "c": {"G()": 2 * 1.3122028733013222},
    }
    assert cu.tighten('{suffix=>"a"}', 1e-12) == '{atol=>1e-12,rtol=>1e-12,suffix=>"a"}'
