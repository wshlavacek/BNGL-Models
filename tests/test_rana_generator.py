"""The generator in models/amyloid_beta_competing_aggregation_pathways_rana2020/generator/
must regenerate all six committed BNGL files byte-for-byte. The six files are near-identical
transcriptions of one 80-reaction network under six protocols, so they are written by a
generator rather than by hand; this guards against drift between the two."""

import sys
from pathlib import Path

import pytest

MODEL_DIR = Path("models/amyloid_beta_competing_aggregation_pathways_rana2020")
GEN_DIR = MODEL_DIR / "generator"
TARGETS = sorted(p.name for p in MODEL_DIR.glob("*.bngl"))


@pytest.fixture(scope="module")
def build_file():
    sys.path.insert(0, str(GEN_DIR))
    from build_rana2020 import build_file as _build_file

    return _build_file


def test_all_committed_models_are_covered(build_file):
    from build_rana2020 import SPECS

    assert sorted(SPECS) == TARGETS, "generator and models/ disagree on the file set"


@pytest.mark.parametrize("target", TARGETS)
def test_generator_reproduces_committed_model(build_file, target):
    committed = (MODEL_DIR / target).read_text(encoding="utf-8")
    assert build_file(target) == committed, (
        f"generator output for {target} diverged from the committed file; regenerate with "
        f"`python {GEN_DIR}/build_rana2020.py` or reconcile the generator"
    )
