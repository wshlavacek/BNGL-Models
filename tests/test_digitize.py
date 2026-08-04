"""`skills/curate-model/scripts/digitize.py` holds the machinery shared by the committed
figure digitizers, so a change to it can silently move every digitized value in the collection.
Two guards: the module's own self-test over every helper, and a byte-for-byte replay of the two
digitizers whose source PDFs are obtainable — mu2010 (vector, `pdftocairo -svg`) and
malleshaiah2010 (raster, `pdftoppm` plus colour separation), which between them exercise both
routes. The replays skip when the uncommitted PDF or an optional dependency is absent."""

import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path("skills/curate-model/scripts")

MU2010 = Path("models/four_flux_network_isotopomer_labeling_mu2010")
MU2010_PDFS = [
    Path("dev/papers/Mu2010/HandbookChemoinformatics10.pdf"),
    Path("dev/papers/archived_reference_materials/Mu2010/HandbookChemoinformatics10.pdf"),
]
MALLESHAIAH = Path("models/ste5_fus3_ptc1_switch_malleshaiah2010")
MALLESHAIAH_PDF = Path("dev/papers/Malleshaiah2010/nature08946.pdf")


@pytest.fixture(scope="module")
def digitize():
    sys.path.insert(0, str(SCRIPTS))
    import digitize as module

    return module


def test_self_test_covers_every_helper(digitize):
    """The module's `python digitize.py` self-test, run as a test so CI sees a regression."""
    digitize._self_test()


def test_every_public_name_is_importable(digitize):
    missing = [name for name in digitize.__all__ if not hasattr(digitize, name)]
    assert not missing, f"__all__ names nothing defines: {missing}"


def _replay(script: Path, argv: list[str], outputs: list[Path]) -> None:
    """Run a digitizer and assert it rewrote its committed CSVs unchanged.

    The committed bytes are restored either way, so a regression fails the test instead of
    leaving the working tree dirty.
    """
    before = {path: path.read_bytes() for path in outputs}
    try:
        subprocess.run([sys.executable, str(script), *argv], check=True, capture_output=True)
        after = {path: path.read_bytes() for path in outputs}
    finally:
        for path, data in before.items():
            path.write_bytes(data)
    for path in outputs:
        assert after[path] == before[path], (
            f"{script} no longer reproduces {path}; digitize.py changed a committed value. "
            f"Reconcile before committing — the library must not move digitized data."
        )


@pytest.mark.skipif(shutil.which("pdftocairo") is None, reason="poppler pdftocairo not on PATH")
def test_mu2010_digitizer_reproduces_committed_csv():
    pdf = next((p for p in MU2010_PDFS if p.exists()), None)
    if pdf is None:
        pytest.skip("Mu2010 source PDF is not committed and is not present in dev/papers/")
    _replay(
        MU2010 / "digitize_mu2010.py",
        [str(pdf)],
        [MU2010 / "reference" / "mu2010_fig15_4_digitized.csv"],
    )


@pytest.mark.skipif(shutil.which("pdftoppm") is None, reason="poppler pdftoppm not on PATH")
@pytest.mark.skipif(importlib.util.find_spec("PIL") is None,
                    reason="pillow not installed (uv sync --group digitize)")
def test_malleshaiah2010_digitizer_reproduces_committed_csvs():
    if not MALLESHAIAH_PDF.exists():
        pytest.skip("Malleshaiah2010 source PDF is not committed and is not present in dev/papers/")
    _replay(
        MALLESHAIAH / "digitize_malleshaiah2010.py",
        [],
        sorted((MALLESHAIAH / "reference").glob("malleshaiah2010_fig3b_*_digitized.csv")),
    )
