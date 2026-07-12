#!/usr/bin/env python
"""Gate 2 (model fidelity): compare two BNGL models for structural / functional equivalence.

The fitted rate-constant *values* differ between a curated job and the author's file (they are
being estimated), so a text diff is misleading. This script compares what actually defines the
model's behaviour:

  * PRIMARY (network-generating models): expand each model's reaction network with BNG2.pl and
    compare the SPECIES set and the REACTION TOPOLOGY (each reaction as sorted reactant species
    -> sorted product species), IGNORING rate values. Identical species + reactions => the models
    generate the same network => functionally equivalent (rate *laws* still checked by hand).
  * FALLBACK (network-free / generation times out): compare the reaction-rules, molecule-types,
    seed-species and observables BLOCKS textually, normalized (comments/whitespace stripped), and
    report lines present in only one model.

Each model is expanded using ITS OWN generate_network directive if it has one (so a max_stoich /
max_iter cap is honoured); otherwise a bare generate_network({overwrite=>1}) is appended. Simulate
/ parameter_scan actions are dropped.

    export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
    PY=~/Code/PyBNF/.venv/bin/python
    $PY model_diff.py <our.bngl> <author_reference.bngl> [--timeout 180] [--rules-only]

Exit 0 if the comparison ran (read the verdict); 1 on a read error; 2 if BNGPATH/BNG2.pl missing.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def _bng() -> Path:
    bp = os.environ.get("BNGPATH")
    if not bp or not (Path(bp) / "BNG2.pl").is_file():
        print("FAIL  set BNGPATH to the folder containing BNG2.pl")
        raise SystemExit(2)
    return Path(bp) / "BNG2.pl"


def _block(text: str, name: str) -> list[str]:
    """Return normalized non-comment lines of a `begin <name> ... end <name>` block."""
    m = re.search(rf"begin\s+{name}\b(.*?)end\s+{name}\b", text, re.S | re.I)
    if not m:
        return []
    out = []
    for ln in m.group(1).splitlines():
        s = ln.split("#", 1)[0].strip()
        s = re.sub(r"\s+", " ", s)
        if s:
            out.append(s)
    return out


def generate_net(model: Path, bng: Path, timeout: int) -> str | None:
    """Expand the model to a .net; return its text, or None on timeout/failure."""
    src = model.read_text()
    head = re.split(r"end\s+reaction\s+rules", src, flags=re.I)[0]
    head += "end reaction rules\n"
    if re.search(r"begin\s+model", head, re.I):
        head += "end model\n"
    gn = re.search(r"^\s*generate_network\s*\(.*?\)", src, re.M | re.S)
    gn_line = gn.group(0).strip() if gn else "generate_network({overwrite=>1})"
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "m.bngl"
        p.write_text(head + gn_line + "\n")
        try:
            subprocess.run(["perl", str(bng), str(p)], cwd=d, check=True, timeout=timeout,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
            return None
        nets = list(Path(d).glob("*.net"))
        return nets[0].read_text() if nets else None


def parse_net(net: str):
    species = {}
    for ln in _block(net, "species"):
        parts = ln.split()
        if len(parts) >= 2 and parts[0].isdigit():
            species[parts[0]] = parts[1]
    reactions = set()
    for ln in _block(net, "reactions"):
        parts = ln.split()
        if len(parts) >= 3 and parts[0].isdigit():
            r = tuple(sorted(species.get(i, i) for i in parts[1].split(",") if i != "0"))
            p = tuple(sorted(species.get(i, i) for i in parts[2].split(",") if i != "0"))
            reactions.add((r, p))
    return set(species.values()), reactions


def _report_set(label: str, a: set, b: set) -> bool:
    only_a, only_b = a - b, b - a
    same = not only_a and not only_b
    print(f"  {label}: A={len(a)}  B={len(b)}  common={len(a & b)}  onlyA={len(only_a)}  onlyB={len(only_b)}"
          + ("   IDENTICAL" if same else ""))
    for tag, s in (("onlyA", only_a), ("onlyB", only_b)):
        for ex in list(s)[:4]:
            print(f"      {tag}: {ex}")
        if len(s) > 4:
            print(f"      {tag}: ... (+{len(s) - 4} more)")
    return same


def net_compare(a: Path, b: Path, bng: Path, timeout: int) -> int:
    na = generate_net(a, bng, timeout)
    nb = generate_net(b, bng, timeout)
    if na is None or nb is None:
        which = ", ".join(x.name for x, n in ((a, na), (b, nb)) if n is None)
        print(f"warn  network generation failed/timed out for: {which} -> falling back to rules-level diff")
        return rules_compare(a, b)
    sa, ra = parse_net(na)
    sb, rb = parse_net(nb)
    print("Generated-network comparison (rate values ignored):")
    same_s = _report_set("species ", sa, sb)
    same_r = _report_set("reactions", ra, rb)
    print()
    if same_s and same_r:
        print("PASS (equiv)  identical generated network (species + reaction topology). "
              "Still confirm rate LAWS + units by hand (this compares topology, not rate expressions).")
    else:
        print("REVIEW  networks differ -- inspect the onlyA/onlyB items above. A difference may be "
              "cosmetic (naming), a cap mismatch (max_stoich/iter), or a real structural change. "
              "Record each in VALIDATION.md Gate 2.")
    return 0


def rules_compare(a: Path, b: Path) -> int:
    ta, tb = a.read_text(), b.read_text()
    print("Rules-level comparison (normalized text; network-free or generation skipped):")
    same_all = True
    for blk in ("molecule types", "seed species", "species", "observables", "reaction rules"):
        la, lb = set(_block(ta, blk)), set(_block(tb, blk))
        if not la and not lb:
            continue
        same_all &= _report_set(blk, la, lb)
    print()
    print("PASS (equiv)  identical normalized blocks." if same_all else
          "REVIEW  blocks differ -- inspect onlyA/onlyB; some differences are cosmetic (spacing, "
          "component order). Confirm rule patterns + rate laws by hand. Record in VALIDATION.md.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model_a", type=Path, help="our curated .bngl")
    ap.add_argument("model_b", type=Path, help="the author's reference .bngl (from the parking garage)")
    ap.add_argument("--timeout", type=int, default=180, help="per-model network-generation timeout (s)")
    ap.add_argument("--rules-only", action="store_true", help="skip generation; compare rule blocks textually")
    args = ap.parse_args()
    for p in (args.model_a, args.model_b):
        if not p.is_file():
            print(f"FAIL  no such file: {p}")
            return 1
    if args.rules_only:
        return rules_compare(args.model_a, args.model_b)
    return net_compare(args.model_a, args.model_b, _bng(), args.timeout)


if __name__ == "__main__":
    sys.exit(main())
