"""`tools/conformance.py` is the regression guard over `models/`, so it has to hold two
properties at once: it must be green on the collection as it stands, and it must actually fire
when a folder drifts. A validator that has quietly stopped checking looks exactly like a clean
repository. The tests below cover both — the live collection, and synthetic folders built to
break one rule each."""

import shutil
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools import conformance  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
MODELS = REPO_ROOT / "models"

# A folder that passes every rule today, used as the starting point for the drift fixtures.
EXEMPLAR = "genetic_switch_gardner2000"


# --------------------------------------------------------------------------------------------
# The live collection
# --------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def live():
    return conformance.run(REPO_ROOT)


def test_collection_has_no_unbaselined_errors(live) -> None:
    """The guard is green. A new ERROR here means the commit under test broke conformance —
    read the message, which names the house-style section it came from."""
    baseline = conformance.load_baseline()
    reported, _suppressed, _stale = conformance.apply_baseline(live, baseline)
    errors = [f for f in reported if f.severity == conformance.ERROR]
    assert not errors, "\n" + conformance.format_report(errors, show_warnings=False)


def test_baseline_has_no_stale_entries(live) -> None:
    """Every grandfathered entry still corresponds to a live finding. When someone fixes a
    file, deleting its baseline line is part of the fix — this is what stops the baseline
    outliving the debt it records."""
    baseline = conformance.load_baseline()
    _reported, _suppressed, stale = conformance.apply_baseline(live, baseline)
    assert not stale, "baseline entries that no longer fire:\n  " + "\n  ".join(stale)


def test_baseline_holds_only_the_two_known_classes(live) -> None:
    """The baseline is corpus debt in two named classes — `.bngl` files with no structured
    header, and files that never declare their simulation intent. Anything else appearing here
    is a regression that was grandfathered instead of fixed."""
    data = yaml.safe_load(conformance.BASELINE_PATH.read_text()) or {}
    assert set(data) <= {"header-tags", "simulation-intent"}, sorted(data)


def test_every_model_folder_is_checked(live) -> None:
    """Guards the discovery loop: a rule that silently stops visiting folders would leave the
    suite green. Every folder must be reachable, which the README row rule proves per folder."""
    folders = {p.name for p in MODELS.iterdir() if p.is_dir()}
    readme_rows = conformance.parse_readme((REPO_ROOT / "README.md").read_text())
    assert folders <= set(readme_rows), sorted(folders - set(readme_rows))


# --------------------------------------------------------------------------------------------
# Drift fixtures — one broken rule each
# --------------------------------------------------------------------------------------------


@pytest.fixture
def drifted(tmp_path):
    """A one-folder `models/` copied from the exemplar, plus its README row, that a test can
    then damage. Returns (repo_root, folder, edit_metadata)."""
    root = tmp_path / "repo"
    folder = root / "models" / EXEMPLAR
    folder.parent.mkdir(parents=True)
    shutil.copytree(MODELS / EXEMPLAR, folder)

    readme_line = next(
        line
        for line in (REPO_ROOT / "README.md").read_text().splitlines()
        if conformance.README_ROW.match(line)
        and conformance.README_ROW.match(line).group("folder") == EXEMPLAR
    )
    (root / "README.md").write_text(
        f"| File(s) | Scale | Description | Reference(s) |\n{readme_line}\n"
    )

    def edit_metadata(mutate):
        path = folder / "metadata.yaml"
        data = yaml.safe_load(path.read_text())
        mutate(data)
        path.write_text(yaml.safe_dump(data, sort_keys=False))

    return root, folder, edit_metadata


def rules(root: Path) -> set[str]:
    return {f.rule for f in conformance.run(root)}


def test_exemplar_is_clean(drifted) -> None:
    root, _folder, _edit = drifted
    errors = [f for f in conformance.run(root) if f.severity == conformance.ERROR]
    assert not errors, conformance.format_report(errors, show_warnings=False)


def test_catches_missing_metadata(drifted) -> None:
    root, folder, _edit = drifted
    (folder / "metadata.yaml").unlink()
    assert "metadata-present" in rules(root)


def test_catches_id_folder_mismatch(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d.update(id="something_else"))
    assert "metadata-id" in rules(root)


def test_catches_missing_required_key(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d.pop("created"))
    assert "metadata-schema" in rules(root)


def test_catches_point_of_contact_without_email(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d["point_of_contact"].pop("email"))
    assert "point-of-contact" in rules(root)


def test_catches_invented_source_tag(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d["source"]["tags"].append("definitely_not_a_tag"))
    assert "source-tags" in rules(root)


def test_catches_invented_file_role(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d["files"][0].update(role="miscellaneous"))
    assert "file-roles" in rules(root)


def test_catches_two_primaries(drifted) -> None:
    root, _folder, edit = drifted

    def mutate(d):
        for entry in d["files"]:
            if str(entry.get("name", "")).endswith(".bngl"):
                entry["role"] = "primary"

    edit(mutate)
    assert "primary-file" in rules(root)


def test_catches_bngl_without_scale(drifted) -> None:
    root, _folder, edit = drifted

    def mutate(d):
        for entry in d["files"]:
            if str(entry.get("name", "")).endswith(".bngl"):
                entry.pop("scale", None)

    edit(mutate)
    assert "bngl-scale" in rules(root)


def test_catches_resurrected_rating_field(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d["files"][0].update(rating=4))
    assert "no-rating" in rules(root)


def test_catches_manifest_entry_with_no_file(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d["files"].append({"name": "not_here.gdat", "role": "reference"}))
    assert "manifest-paths" in rules(root)


def test_catches_missing_header_tag(drifted) -> None:
    root, folder, _edit = drifted
    bngl = folder / f"{EXEMPLAR}.bngl"
    bngl.write_text(bngl.read_text().replace("#@keyword:", "# keyword:"))
    assert "header-tags" in rules(root)


def test_catches_undeclared_simulation_intent(drifted) -> None:
    root, folder, _edit = drifted
    bngl = folder / f"{EXEMPLAR}.bngl"
    text = bngl.read_text()
    header = conformance.header_text(text)
    scrubbed = conformance.SIMULATION_INTENT.sub("REDACTED", header)
    bngl.write_text(text.replace(header, scrubbed, 1))
    assert "simulation-intent" in rules(root)


def test_catches_missing_verification_figure(drifted) -> None:
    root, folder, _edit = drifted
    for png in folder.glob("verify_*.png"):
        png.unlink()
    assert "verification-figure" in rules(root)


def test_catches_nested_reference_directory(drifted) -> None:
    """The §5.4 rule settled by flattening `ste5_fus3_ptc1_switch_malleshaiah2010`: `reference/`
    holds files and BioNetGen scan directories, and nothing else."""
    root, folder, _edit = drifted
    (folder / "reference" / "per_model_subdir").mkdir()
    assert "reference-layout" in rules(root)


def test_allows_bng_scan_directories(drifted) -> None:
    """The same rule must not fire on the 35 legitimate `*_scan` directories in the collection."""
    root, folder, _edit = drifted
    scan = folder / "reference" / "some_variant_scan_koff"
    scan.mkdir()
    (scan / "some_variant_scan_koff_0001.gdat").write_text("# point\n")
    assert "reference-layout" not in rules(root)


NF_PROTOCOL = 'simulate({method=>"nf",suffix=>"nfr"})'


def _make_network_free(bngl: Path, protocol: str = NF_PROTOCOL) -> None:
    """Rewrite a `.bngl`'s actions block so its committed protocol is network-free."""
    text = bngl.read_text()
    head = text[: text.index("begin actions")]
    bngl.write_text(f"{head}begin actions\n  {protocol}\nend actions\n")


def test_catches_network_free_bngl_without_xml(drifted) -> None:
    """The defect behind issue #41: for a network-free model the XML is what NFsim and
    RuleMonkey actually read, so it is committed the way a generate-first model commits its
    `.net`. Six folders had shipped without one."""
    root, folder, _edit = drifted
    _make_network_free(folder / f"{EXEMPLAR}.bngl")
    assert "network-free-xml" in rules(root)


def test_catches_writexml_only_protocol_without_xml(drifted) -> None:
    """`lambda_switch_arkin1998_fullcircuit`'s shape: the whole actions block is `writeXML()`,
    so the XML is the only output the protocol produces — and was the one thing not committed."""
    root, folder, _edit = drifted
    _make_network_free(folder / f"{EXEMPLAR}.bngl", "writeXML()")
    assert "network-free-xml" in rules(root)


def test_allows_the_protocol_suffixed_xml_name(drifted) -> None:
    """`reference/<stem>.xml` is the default, but eight folders predate it and carry the
    protocol's own suffix (`blbr_dembo1978_nfr.xml`). Both satisfy the rule."""
    root, folder, _edit = drifted
    _make_network_free(folder / f"{EXEMPLAR}.bngl")
    (folder / "reference" / f"{EXEMPLAR}_nfr.xml").write_text("<sbml/>\n")
    assert "network-free-xml" not in rules(root)


def test_network_free_xml_is_not_satisfied_by_a_sibling(drifted) -> None:
    """The trap the suffix allowance opens, live in `blbr_rings_posner1995`: a variant's
    `<stem>_no_rings_nfr.xml` starts with the primary's stem, and must not count for it."""
    root, folder, _edit = drifted
    _make_network_free(folder / f"{EXEMPLAR}.bngl")
    variant = folder / f"{EXEMPLAR}_variant.bngl"
    shutil.copy(folder / f"{EXEMPLAR}.bngl", variant)
    (folder / "reference" / f"{EXEMPLAR}_variant_nfr.xml").write_text("<sbml/>\n")
    findings = [f for f in conformance.run(root) if f.rule == "network-free-xml"]
    assert [f.where for f in findings] == [f"models/{EXEMPLAR}/{EXEMPLAR}.bngl"]


def test_generate_first_bngl_needs_no_xml(drifted) -> None:
    """The rule keys on the protocol, not on the folder: the exemplar generates its network and
    commits `.net`/`.gdat`, and must not be asked for an XML."""
    root, _folder, _edit = drifted
    assert "network-free-xml" not in rules(root)


def test_catches_folder_missing_from_readme(drifted) -> None:
    root, _folder, _edit = drifted
    (root / "README.md").write_text("| File(s) | Scale | Description | Reference(s) |\n")
    assert "readme-row" in rules(root)


def test_catches_readme_scale_drift(drifted) -> None:
    root, _folder, edit = drifted

    def mutate(d):
        for entry in d["files"]:
            if str(entry.get("name", "")).endswith(".bngl"):
                entry["scale"] = "hours"

    edit(mutate)
    assert "readme-scale" in rules(root)


def test_readme_scale_counts_non_bngl_entries(drifted) -> None:
    """A driver script may be the expensive artifact (§5.6), and four folders are `hours` for
    exactly that reason — so the README maximum must include entries that are not `.bngl`."""
    root, _folder, edit = drifted
    edit(
        lambda d: d["files"].append(
            {"name": "metadata.yaml", "role": "tooling", "scale": "cluster"}
        )
    )
    assert "readme-scale" in rules(root)


# --------------------------------------------------------------------------------------------
# The baseline mechanism itself
# --------------------------------------------------------------------------------------------


def test_baseline_suppresses_only_the_named_finding(drifted) -> None:
    root, _folder, edit = drifted
    edit(lambda d: d["source"]["tags"].append("definitely_not_a_tag"))
    findings = conformance.run(root)
    key = next(f.key for f in findings if f.rule == "source-tags")

    reported, suppressed, stale = conformance.apply_baseline(findings, {key})
    assert [f.key for f in suppressed] == [key]
    assert key not in {f.key for f in reported}
    assert not stale


def test_stale_baseline_entry_is_itself_reported(drifted) -> None:
    root, _folder, _edit = drifted
    findings = conformance.run(root)
    _reported, _suppressed, stale = conformance.apply_baseline(
        findings, {"source-tags models/gone"}
    )
    assert stale == ["source-tags models/gone"]


def test_warnings_are_never_suppressed_by_the_baseline(drifted) -> None:
    """WARN findings do not fail CI, so grandfathering one would only hide it."""
    root, _folder, edit = drifted
    edit(lambda d: [e.pop("documentation_target", None) for e in d["files"]])
    findings = conformance.run(root)
    warning = next(f for f in findings if f.rule == "documentation-target")
    reported, suppressed, _stale = conformance.apply_baseline(findings, {warning.key})
    assert warning in reported
    assert not suppressed
