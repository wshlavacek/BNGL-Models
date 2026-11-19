"""models/cell_cycle_oscillator_lang2024/generator/ must regenerate the committed `_fitted` file.

`cell_cycle_oscillator_lang2024_fitted.bngl` takes every number from the authors' vendored
PEtab tables and its model body from the pinned `_fused` file, so it is written by a generator
rather than by hand. This guards against drift between the generator, its pinned inputs, and
the committed file. No BioNetGen is needed; the network-level checks live in
verify_lang2024_fitted.ipynb.
"""

import hashlib
import sys
from pathlib import Path

import pytest

MODEL_DIR = Path("models/cell_cycle_oscillator_lang2024")
GEN_DIR = MODEL_DIR / "generator"
TARGET = MODEL_DIR / "cell_cycle_oscillator_lang2024_fitted.bngl"


@pytest.fixture(scope="module")
def gen():
    sys.path.insert(0, str(GEN_DIR))
    import build_lang2024_fitted

    return build_lang2024_fitted


def test_vendored_upstream_files_match_their_pins(gen):
    for name, (_, digest) in gen.PINS.items():
        path = MODEL_DIR / "reference" / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name


def test_generator_reproduces_committed_model(gen):
    assert gen.build() == TARGET.read_text(encoding="utf-8"), (
        f"generator output diverged from {TARGET.name}; regenerate with "
        f"`python {GEN_DIR}/build_lang2024_fitted.py` or reconcile the generator"
    )
