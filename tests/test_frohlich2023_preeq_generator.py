"""models/mapk_adaptive_resistance_frohlich2023/generator/ must regenerate the `_preeq` file.

`mapk_adaptive_resistance_frohlich2023_preeq.bngl` takes its model body from the pinned primary
file and its seed amounts from the authors' steady state that run_frohlich2023.py writes to
reference/, so it is written by a generator rather than by hand. This guards against drift
between the generator, its pinned inputs, and the committed file. Neither BioNetGen nor PySB
is needed; the comparisons with the authors' model live in run_frohlich2023.py.
"""

import hashlib
import sys
from pathlib import Path

import pytest

MODEL_DIR = Path("models/mapk_adaptive_resistance_frohlich2023")
GEN_DIR = MODEL_DIR / "generator"
TARGET = MODEL_DIR / "mapk_adaptive_resistance_frohlich2023_preeq.bngl"


@pytest.fixture(scope="module")
def gen():
    sys.path.insert(0, str(GEN_DIR))
    import build_frohlich2023_preeq

    return build_frohlich2023_preeq


def test_pinned_inputs_match(gen):
    for path, digest in gen.PINS.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, path.name


def test_generator_reproduces_committed_model(gen):
    assert gen.build() == TARGET.read_text(encoding="utf-8"), (
        f"generator output diverged from {TARGET.name}; regenerate with "
        f"`python {GEN_DIR}/build_frohlich2023_preeq.py` or reconcile the generator"
    )


def test_only_header_seeds_and_actions_differ_from_the_primary(gen):
    primary = gen.PRIMARY.read_text(encoding="utf-8")
    variant = TARGET.read_text(encoding="utf-8")
    for name in gen.BODY_BLOCKS + gen.TAIL_BLOCKS:
        p0, p1 = gen.block_span(primary, name)
        v0, v1 = gen.block_span(variant, name)
        assert primary[p0:p1] == variant[v0:v1], name

    def species(text):
        start, end = gen.block_span(text, "seed species")
        return [line.split()[0] for line in text[start:end].splitlines()[1:-1]]

    assert species(primary) == species(variant)
