"""Regenerate every reference/blasi2016_*.csv and reference/feller2015_*.csv.

The reported quantities this model is verified against are all available as
tabular data, so nothing here is digitized from a plotted figure.

Three sources, none of which is committed to this repository:

1. ``mmc1.pdf`` -- Supplemental Information of Blasi et al. (2016). Table S1
   holds the eight estimated rate constants of the best motif-specific model
   with 95% profile-likelihood confidence intervals, rounded to three decimals.
   The unspecific and site-specific estimates are quoted in the Results of the
   article body (``mmc3.pdf``). Values are transcribed here rather than parsed,
   because Table S1 is laid out as free-floating text runs in the PDF.

2. ``mmc2.xlsx`` -- Table S2 of Blasi et al. (2016): all 32 acetylation rate
   constants and the BIC of the best 100 models. Its rank-1 row is the best
   motif-specific model, and it reports the rate constants to four digits,
   one more than Table S1. Those are the values carried by the BNGL files.
   Parsed straight from the sheet XML so that no spreadsheet engine is needed.

3. The PEtab problem ``Blasi_CellSystems2016`` of the Benchmarking-Initiative
   ``Benchmark-Models-PEtab`` collection, which redistributes the LC-MS motif
   abundances of Feller et al. (2015) that Fig. 2 of Blasi et al. (2016) plots
   as x markers, together with a reference steady-state simulation of the best
   motif-specific model. That the abundance table is the one the paper actually
   fit is confirmed downstream in verify_blasi2016.ipynb, which recovers all
   three BIC values the paper reports from it.

Usage::

    python extract_blasi2016.py --papers ~/Code/BNGL-Models/dev/papers/blasi2016

The PEtab problem is downloaded on demand; pass --petab to point at an already
extracted copy of the ``Blasi_CellSystems2016`` directory instead.
"""

from __future__ import annotations

import argparse
import csv
import io
import itertools
import tarfile
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

BENCHMARK_TARBALL = (
    "https://api.github.com/repos/Benchmarking-Initiative/Benchmark-Models-PEtab/tarball/master"
)
SHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
SITES = ("k5", "k8", "k12", "k16")

# Table S1 in Blasi et al. (2016) and the Results paragraphs of the article
# body: maximum-likelihood estimate and 95% profile-likelihood confidence
# interval for each acetylation rate constant of each scenario. Rate constants
# are relative to the deacetylation rate constant, which is fixed at 1.
REPORTED_RATES = [
    # scenario, parameter, estimate, ci_low, ci_high, source
    ("unspecific", "a_b", 0.182, 0.170, 0.195, "Results, Blasi et al. (2016)"),
    ("site_specific", "a_K5", 0.225, 0.181, 0.280, "Results, Blasi et al. (2016)"),
    ("site_specific", "a_K8", 0.126, 0.104, 0.154, "Results, Blasi et al. (2016)"),
    ("site_specific", "a_K12", 0.417, 0.331, 0.528, "Results, Blasi et al. (2016)"),
    ("site_specific", "a_K16", 0.107, 0.088, 0.131, "Results, Blasi et al. (2016)"),
    ("motif_specific", "a_b", 0.067, 0.063, 0.071, "Table S1, Blasi et al. (2016)"),
    ("motif_specific", "a_0ac_K8", 0.027, 0.022, 0.033, "Table S1, Blasi et al. (2016)"),
    ("motif_specific", "a_K5_K5K12", 2.062, 1.714, 2.500, "Table S1, Blasi et al. (2016)"),
    ("motif_specific", "a_K12_K5K12", 0.552, 0.419, 0.705, "Table S1, Blasi et al. (2016)"),
    ("motif_specific", "a_K16_K12K16", 0.696, 0.595, 0.813, "Table S1, Blasi et al. (2016)"),
    (
        "motif_specific",
        "a_K5K12_K5K8K12",
        0.325,
        0.273,
        0.387,
        "Table S1, Blasi et al. (2016)",
    ),
    (
        "motif_specific",
        "a_K12K16_K8K12K16",
        2.206,
        1.914,
        2.542,
        "Table S1, Blasi et al. (2016)",
    ),
    (
        "motif_specific",
        "a_K8K12K16_4ac",
        3.592,
        3.102,
        4.154,
        "Table S1, Blasi et al. (2016)",
    ),
]

# Bayesian information criterion of the best model in each scenario. The
# motif-specific value is given to four decimals in the rank-1 row of Table S2;
# the other two are quoted to one decimal in the Results.
REPORTED_BIC = [
    ("unspecific", 1, 697.3, "Results, Blasi et al. (2016)"),
    ("site_specific", 4, 644.5, "Results, Blasi et al. (2016)"),
    ("motif_specific", 8, 67.1194, "Table S2 rank 1, Blasi et al. (2016)"),
]


def motif_name(state: tuple[int, ...]) -> str:
    """Motif label of an acetylation state, in the notation of Fig. 1B."""
    marked = [site for site, bit in zip(SITES, state, strict=True) if bit]
    if not marked:
        return "0ac"
    if len(marked) == 4:
        return "4ac"
    return "".join(marked)


def acetylation_edges() -> list[tuple[str, str, str]]:
    """The 32 edges of the 4-cube, as (substrate, product, site acetylated)."""
    out = []
    for state in itertools.product((0, 1), repeat=4):
        for position, site in enumerate(SITES):
            if state[position] == 0:
                product = list(state)
                product[position] = 1
                out.append((motif_name(state), motif_name(tuple(product)), site))
    return out


def read_table_s2_rank1(xlsx: Path) -> dict[str, float]:
    """Rate constants of the rank-1 (best motif-specific) model of Table S2."""
    with zipfile.ZipFile(xlsx) as archive:
        strings = [
            "".join(node.text or "" for node in item.iter(SHEET_NS + "t"))
            for item in ET.fromstring(archive.read("xl/sharedStrings.xml"))
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    def column_index(reference: str) -> int:
        number = 0
        for char in reference:
            if char.isalpha():
                number = number * 26 + ord(char.upper()) - 64
        return number

    rows = []
    for row in sheet.iter(SHEET_NS + "row"):
        cells = {}
        for cell in row.iter(SHEET_NS + "c"):
            value = cell.find(SHEET_NS + "v")
            if value is None:
                continue
            text = strings[int(value.text)] if cell.get("t") == "s" else value.text
            cells[column_index(cell.get("r"))] = text
        if cells:
            rows.append(cells)

    header = rows[0]
    width = max(header)
    names = [header.get(i, "") for i in range(1, width + 1)]
    for row in rows[1:]:
        values = [row.get(i, "") for i in range(1, width + 1)]
        record = dict(zip(names, values, strict=True))
        if record.get("model rank") == "1":
            return record
    raise ValueError(f"no rank-1 row in {xlsx}")


def fetch_petab(destination: Path) -> Path:
    """Download the Blasi PEtab problem from the benchmark collection."""
    with urllib.request.urlopen(BENCHMARK_TARBALL) as response:  # noqa: S310
        payload = response.read()
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        wanted = [
            member
            for member in archive.getmembers()
            if "/Blasi_CellSystems2016/" in member.name and member.isfile()
        ]
        for member in wanted:
            member.name = member.name.split("/Blasi_CellSystems2016/", 1)[1]
            if "/" in member.name:
                continue
            archive.extract(member, destination)
    return destination


def read_petab_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_csv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"wrote {path}")


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--papers",
        type=Path,
        default=here.parents[1] / "dev" / "papers" / "blasi2016",
        help="folder holding mmc1.pdf, mmc2.xlsx and mmc3.pdf",
    )
    parser.add_argument(
        "--petab",
        type=Path,
        default=None,
        help="already extracted Blasi_CellSystems2016 PEtab folder",
    )
    parser.add_argument("--out", type=Path, default=here / "reference")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # --- reported rate constants and BIC ------------------------------------
    write_csv(
        args.out / "blasi2016_reported_rate_constants.csv",
        ["scenario", "parameter", "estimate", "ci_low", "ci_high", "source"],
        REPORTED_RATES,
    )
    write_csv(
        args.out / "blasi2016_reported_bic.csv",
        ["scenario", "n_parameters", "bic", "source"],
        REPORTED_BIC,
    )

    # --- Table S2, rank 1: every edge of the best motif-specific model ------
    rank1 = read_table_s2_rank1(args.papers / "mmc2.xlsx")
    basal = float(rank1["basal rate"])
    edge_rows = []
    for substrate, product, site in acetylation_edges():
        head = "0" if substrate == "0ac" else substrate.replace("k", "K")
        label = f"{head}->{product.replace('k', 'K')}"
        rate = float(rank1[label])
        edge_rows.append(
            [
                substrate,
                product,
                site,
                f"{rate:.6g}",
                "basal" if rate == basal else "motif_specific",
            ]
        )
    write_csv(
        args.out / "blasi2016_tableS2_rank1_edge_rates.csv",
        ["substrate_motif", "product_motif", "site", "rate_constant", "kind"],
        edge_rows,
    )

    # --- LC-MS abundances and the PEtab reference steady state --------------
    with tempfile.TemporaryDirectory(prefix="blasi_petab_") as scratch:
        petab_dir = args.petab or fetch_petab(Path(scratch))
        measurements = read_petab_table(petab_dir / "measurementData_Blasi_CellSystems2016.tsv")
        simulated = read_petab_table(petab_dir / "simulatedData_Blasi_CellSystems2016.tsv")

    counter: dict[str, int] = {}
    measurement_rows = []
    for row in measurements:
        motif = row["observableId"].removeprefix("observable_")
        counter[motif] = counter.get(motif, 0) + 1
        measurement_rows.append([motif, counter[motif], row["measurement"]])
    write_csv(
        args.out / "feller2015_h4_motif_abundances.csv",
        ["motif", "replicate", "relative_abundance"],
        measurement_rows,
    )

    reference = {}
    for row in simulated:
        reference[row["observableId"].removeprefix("observable_")] = row["simulation"]
    write_csv(
        args.out / "blasi2016_petab_reference_abundances.csv",
        ["motif", "relative_abundance"],
        [[motif, value] for motif, value in reference.items()],
    )


if __name__ == "__main__":
    main()
