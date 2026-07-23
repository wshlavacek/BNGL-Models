"""The reconstructed generator in dev/lambda_arkin_fullmodel/ must regenerate the two
committed network-free lambda full-circuit models byte-for-byte. This guards against
drift between the generator and the committed .bngl files."""
import sys
from pathlib import Path

import pytest

GEN_DIR = Path("dev/lambda_arkin_fullmodel")
MODEL_DIR = Path("models/lambda_switch_arkin1998")
TARGETS = {
    "base": "lambda_switch_arkin1998_fullcircuit.bngl",
    "exact": "lambda_switch_arkin1998_fullcircuit_exact.bngl",
}


@pytest.fixture(scope="module")
def build():
    sys.path.insert(0, str(GEN_DIR))
    from build_fullcircuit import build as _build

    return _build


@pytest.mark.parametrize("variant", sorted(TARGETS))
def test_generator_reproduces_committed_model(build, variant):
    committed = (MODEL_DIR / TARGETS[variant]).read_text(encoding="utf-8")
    assert build(variant) == committed, (
        f"generator output for {variant!r} diverged from the committed "
        f"{TARGETS[variant]}; regenerate with "
        f"`python {GEN_DIR}/build_fullcircuit.py {variant}` or reconcile the generator"
    )
