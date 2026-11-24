#!/usr/bin/env python3
"""Recover Fig. S6 of Chylek et al. (2014) from the publisher's File S2, and Table S1's rows.

Fig. S6 (File S2, page 9) plots, for the 16 pTyr sites used to build and fit the model, the
measured time course (points: mean of three phosphoproteomic experiments, error bars: SD) and
the published simulation (solid line), both "normalized to baseline and log2 transformed".
The figure is MATLAB vector art, so it is read as geometry (digitization.md section 1, PDF
drawings route), not rasterized:

* every panel frame is one stroked 4-segment path; its tick marks are short strokes of
  0.8 pt on the bottom edge (x ticks, t = 0 and 50 s) and on the right edge (y ticks,
  log2 = -2 ... 3). Both axes are fitted to the tick marks (section 2) and the tick labels
  (vector text) are checked against the fit;
* the simulated curve is the one long all-line stroked path inside the frame; its vertices
  are the published simulation;
* data points are filled 4-curve circles; their abscissae are snapped to the stated sampling
  times 0, 5, 15, 30 and 60 s (section 7: an abscissa on a known grid); error bars are the
  vertical strokes plus 1.6-pt caps, and SD = half the cap-to-cap distance.

MATLAB wrote the PDF on a 0.8-pt grid (1/90 in): every coordinate in the panel is a multiple
of 0.8 pt from the frame corner. That quantization is the whole digitization error of this
route, +/-0.4 pt, which is +/-0.029 log2 units on these axes (13.9 pt per log2 unit) and
+/-0.33 s in time (1.20 pt/s). It is recorded per panel in the CSV header.

Table S1 (File S1, TableS1.xlsx) is the authors' source data: log2 ratios at 0/5/15/30/60 s
for each phosphopeptide in each of the three LC-MS/MS experiments. The rows for the 16 sites
are written out with their per-time mean and SD over the experiments that quantified the
site, which is what Fig. S6 says its points are. The figure's markers are then asserted to
match Table S1 within the 0.8-pt grid -- an independent check that the calibration, the panel
order and the site mapping are right. One naming difference is carried, not corrected: the
paper's "LCK Y424" (panel h, and the model's LCK_pY424) is Table S1's LCK position 394
(peptide LIEDNEY(ph)TAR, the activation-loop tyrosine in UniProt P06239 numbering).

Usage:

    python digitize_chylek2014.py [FILE_S2_PDF] [TABLE_S1_XLSX]

with defaults under dev/papers/chylek2014/ (relative to the repository root; dev/ is not
committed). Obtain both from the article's Supporting Information
(doi:10.1371/journal.pone.0104240: s002 = File S2 PDF; s001 = File S1 zip, which holds
TableS1.xlsx). Re-running regenerates the three CSVs in reference/ byte for byte. Needs
PyMuPDF (`uv sync --group digitize`); the xlsx is read with the standard library.
"""

from __future__ import annotations

import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "skills/curate-model/scripts"))
from digitize import Axis, open_pdf, write_csv  # noqa: E402

DEFAULT_PDF = HERE.parents[1] / "dev/papers/chylek2014/journal.pone.0104240.s002.pdf"
DEFAULT_XLSX = HERE.parents[1] / "dev/papers/chylek2014/TableS1.xlsx"
REF = HERE / "reference"
PAGE = 8  # 0-based: File S2 page 9, Fig. S6

# Panel letter -> (label printed under the panel, model observable, Table S1 gene, position).
PANELS = [
    ("a", "CD3D Y149", "TCR_pY149_D", "CD3D", "149"),
    ("b", "CD3G Y171", "TCR_pY171_G", "CD3G", "171"),
    ("c", "TCRZ Y111", "TCR_pY111", "CD247", "111"),
    ("d", "TCRZ Y123", "TCR_pY123", "CD247", "123"),
    ("e", "CD3E Y188", "TCR_pY188_E", "CD3E", "188"),
    ("f", "CD3E Y199", "TCR_pY199_E", "CD3E", "199"),
    ("g", "ZAP70 Y493", "ZAP70_pY493", "ZAP70", "493"),
    ("h", "LCK Y424", "LCK_pY424", "LCK", "394"),
    ("i", "LCK Y192", "LCK_pY192", "LCK", "192"),
    ("j", "ITK Y512", "ITK_pY512", "ITK", "512"),
    ("k", "PTPN6 Y566", "PTPN6_pY566", "PTPN6", "566"),
    ("l", "WAS Y291", "WAS_pY291", "WAS", "291"),
    ("m", "PAG1 Y163", "PAG1_pY163", "PAG1", "163"),
    ("n", "DOK2 Y299", "DOK2_pY299", "DOK2", "299"),
    ("o", "DOK1 Y449", "DOK1_pY449", "DOK1", "449"),
    ("p", "PLCG1 Y783", "PLCG1_pY783", "PLCG1", "783"),
]
TIMES = (0.0, 5.0, 15.0, 30.0, 60.0)
GRID = 0.8  # pt: the PDF's coordinate quantum


def _lines(g):
    return [((it[1].x, it[1].y), (it[2].x, it[2].y)) for it in g["items"] if it[0] == "l"]


def _inside(r, box, pad=1.0):
    return (r.x0 >= box.x0 - pad and r.x1 <= box.x1 + pad
            and r.y0 >= box.y0 - pad and r.y1 <= box.y1 + pad)


def _spans(page):
    out = []
    for b in page.get_text("dict")["blocks"]:
        for ln in b.get("lines", []):
            for s in ln["spans"]:
                if s["text"].strip():
                    out.append((s["text"], s["bbox"]))
    return out


def read_figure(pdf):
    doc = open_pdf(pdf)
    page = doc[PAGE]
    dr = page.get_drawings()
    spans = _spans(page)
    black = [g for g in dr if g["type"] == "s" and g.get("color") == (0.0, 0.0, 0.0)]
    frames = [g for g in black if len(g["items"]) == 4 and 70 < g["rect"].width < 85
              and all(it[0] == "l" for it in g["items"])]
    frames.sort(key=lambda g: (round(g["rect"].y0), g["rect"].x0))
    frames = frames[:16]  # a-p; q and r (extra PAG1 sites, data only) are not model sites
    panels = []
    for (letter, label, obs, gene, pos), fr in zip(PANELS, frames):
        box = fr["rect"]
        # The printed site label sits under the frame.
        under = [t.strip() for t, bb in spans if box.x0 - 5 < bb[0] and bb[2] < box.x1 + 5
                 and box.y1 + 10 < bb[1] < box.y1 + 25]
        if under != [label]:
            raise SystemExit(f"panel {letter}: expected label {label!r}, found {under}")
        inner = [g for g in black if _inside(g["rect"], box) and g is not fr]
        segs = [s for g in inner for s in _lines(g)]
        # Error bars are the strokes of paths smaller than the frame (frame and tick paths
        # span its full width or height).
        bsegs = [s for g in inner if g["rect"].width < box.width - 0.1
                 and g["rect"].height < box.height - 0.1 for s in _lines(g)]
        # Tick stubs: 0.8-pt strokes on the bottom edge (x) and on the right edge (y).
        xt = sorted({round(a[0], 3) for a, b in segs if abs(a[0] - b[0]) < 1e-6
                     and 0.5 < abs(a[1] - b[1]) < 1.2 and max(a[1], b[1]) > box.y1 - 0.01})
        yt = sorted({round(a[1], 3) for a, b in segs if abs(a[1] - b[1]) < 1e-6
                     and 0.5 < abs(a[0] - b[0]) < 1.2 and max(a[0], b[0]) > box.x1 - 0.01})
        if len(xt) != 2 or len(yt) != 6:
            raise SystemExit(f"panel {letter}: {len(xt)} x ticks, {len(yt)} y ticks")
        xax = Axis.from_ticks(xt, [0.0, 50.0])
        yax = Axis.from_ticks(yt, [3.0, 2.0, 1.0, 0.0, -1.0, -2.0])
        # Tick labels (vector text, centred under / left of their ticks) must agree.
        for t, bb in spans:
            cx, cy = (bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2
            if t in ("0", "50") and box.y1 < bb[1] < box.y1 + 8 and box.x0 - 6 < cx < box.x1:
                if abs(xax(cx) - float(t)) > 1.5:
                    raise SystemExit(f"panel {letter}: x label {t} at t={xax(cx):.2f}")
        curves = [g for g in inner if len(g["items"]) > 12
                  and all(it[0] == "l" for it in g["items"])]
        if len(curves) != 1:
            raise SystemExit(f"panel {letter}: {len(curves)} candidate curves")
        cs = _lines(curves[0])
        vx = np.array([cs[0][0][0]] + [b[0] for a, b in cs])
        vy = np.array([cs[0][0][1]] + [b[1] for a, b in cs])
        curve_t, curve_y = xax(vx), yax(vy)
        # Markers: filled 4-curve circles of ~3.2 pt.
        marks = [g for g in dr if g["type"] == "f" and len(g["items"]) == 4
                 and all(it[0] == "c" for it in g["items"]) and 2.5 < g["rect"].width < 4.0
                 and _inside(g["rect"], box, pad=2.5)]
        mc = sorted(((g["rect"].x0 + g["rect"].x1) / 2, (g["rect"].y0 + g["rect"].y1) / 2)
                    for g in marks)
        if len(mc) != len(TIMES):
            raise SystemExit(f"panel {letter}: {len(mc)} markers")
        # Error-bar caps: horizontal 1.6-pt strokes; the bar is centred on its marker.
        caps = [(0.5 * (a[0] + b[0]), a[1]) for a, b in bsegs
                if abs(a[1] - b[1]) < 1e-6 and 1.4 < abs(a[0] - b[0]) < 1.8]
        points = []
        for (mx, my), t_nom in zip(mc, TIMES):
            t = xax(mx)
            if abs(t - t_nom) > 0.5:
                raise SystemExit(f"panel {letter}: marker at t={t:.2f}, expected {t_nom}")
            cy = sorted({round(y, 3) for x, y in caps if abs(x - mx) < 0.05})
            if len(cy) > 2:
                raise SystemExit(f"panel {letter} t={t_nom}: {len(cy)} cap heights")
            # A bar that runs out of the frame keeps one cap; its vertical stroke ends on the
            # frame edge. Then only the in-frame half is measurable: SD = |cap - mean|, flagged.
            bar = [(a, b) for a, b in bsegs if abs(a[0] - b[0]) < 1e-6 and abs(a[0] - mx) < 0.05
                   and abs(a[1] - b[1]) > 0.01]
            clipped = any(min(a[1], b[1]) < box.y0 + 0.05 or max(a[1], b[1]) > box.y1 - 0.05
                          for a, b in bar)
            if len(cy) == 2:
                sd = 0.5 * abs(yax(cy[0]) - yax(cy[1]))
            elif len(cy) == 1 and clipped:
                sd = abs(yax(cy[0]) - yax(my))
            else:  # two caps on one grid row (or none): the bar is shorter than 0.8 pt
                sd = 0.0
            points.append((t_nom, float(yax(my)), float(sd), clipped))
        res_pt = max(xax.residual / abs(xax(1.0) - xax(0.0)), yax.residual / abs(yax(1.0) - yax(0.0)))
        panels.append(dict(letter=letter, label=label, obs=obs, gene=gene, pos=pos,
                           curve_t=curve_t, curve_y=curve_y, points=points,
                           pt_per_log2=abs(1.0 / (yax(1.0) - yax(0.0))),
                           pt_per_s=abs(1.0 / (xax(1.0) - xax(0.0))),
                           y_resid=yax.residual, res_pt=res_pt))
    return panels


def _xlsx_rows(path):
    """Rows of the first sheet of an .xlsx as lists of strings/floats (stdlib only)."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{ns['m']}}}t")))
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = []
    for r in sheet.iter(f"{{{ns['m']}}}row"):
        cells = {}
        for c in r.findall("m:c", ns):
            ref = c.get("r")
            col = 0
            for ch in ref:
                if ch.isalpha():
                    col = col * 26 + (ord(ch) - 64)
            v = c.find("m:v", ns)
            if v is None:
                continue
            val = shared[int(v.text)] if c.get("t") == "s" else v.text
            try:
                val = float(val)
            except (TypeError, ValueError):
                pass
            cells[col - 1] = val
        rows.append([cells.get(i) for i in range(max(cells) + 1)] if cells else [])
    return rows


def read_table_s1(path):
    """Per site: the three experiments' log2 ratios at the five times, from Table S1."""
    rows = _xlsx_rows(path)
    head = rows[0]
    e0 = [i for i, h in enumerate(head) if isinstance(h, str) and h.startswith("Exp")]
    if len(e0) != 3:
        raise SystemExit(f"Table S1: expected three 'ExpN ratios' column groups, found {e0}")
    out = {}
    for row in rows[2:]:
        if len(row) < 5 or row[2] is None:
            continue
        pos = row[0]
        pos = str(int(pos)) if isinstance(pos, float) else str(pos)
        genes = str(row[2]).split(";")
        for letter, label, obs, gene, want in PANELS:
            if gene in genes and pos == want:
                if obs in out:
                    raise SystemExit(f"Table S1: two rows for {gene} {want}")
                exps = []
                for j in e0:
                    vals = []
                    for k in range(5):
                        v = row[j + k] if j + k < len(row) else None
                        vals.append(float(v) if isinstance(v, float) else np.nan)
                    exps.append(vals)
                out[obs] = dict(peptide=row[4], exps=np.array(exps))
    missing = [p[2] for p in PANELS if p[2] not in out]
    if missing:
        raise SystemExit(f"Table S1: no row for {missing}")
    return out


def main(argv):
    pdf = Path(argv[1]) if len(argv) > 1 else DEFAULT_PDF
    xlsx = Path(argv[2]) if len(argv) > 2 else DEFAULT_XLSX
    for p, what in ((pdf, "File S2 (journal.pone.0104240.s002)"),
                    (xlsx, "TableS1.xlsx from File S1 (journal.pone.0104240.s001)")):
        if not p.exists():
            sys.exit(f"{p} not found: download {what} from the Supporting Information of "
                     "doi:10.1371/journal.pone.0104240 (dev/papers/ is not committed).")
    panels = read_figure(pdf)
    table = read_table_s1(xlsx)

    # Cross-check: every Fig. S6 marker against the Table S1 mean and SD (ddof=1, over the
    # experiments that quantified the site) it plots, to one 0.8-pt grid step. Panel h is the
    # exception and is reported, not asserted: its points are close to, but not exactly,
    # Table S1's LCK 394 row (see the docstring).
    worst = {}
    for P in panels:
        ex = table[P["obs"]]["exps"]
        n = np.sum(np.isfinite(ex), axis=0)
        mean = np.nanmean(ex, axis=0)
        sd = np.array([np.nanstd(ex[:, k], ddof=1) if n[k] > 1 else 0.0 for k in range(5)])
        step = GRID / P["pt_per_log2"]
        dm = max(abs(m - mean[k]) for k, (t, m, s, c) in enumerate(P["points"])) / step
        ds = max((s - sd[k]) / step if c else abs(s - sd[k]) / step
                 for k, (t, m, s, c) in enumerate(P["points"]))
        worst[P["letter"]] = (P["label"], dm, ds)
        if P["letter"] != "h" and (dm > 1.0 or ds > 1.0):
            raise SystemExit(f"{P['label']}: markers vs Table S1 off by {dm:.2f} (mean), "
                             f"{ds:.2f} (SD) grid steps")
    rest = [v for k, v in worst.items() if k != "h"]
    print(f"Fig. S6 markers vs Table S1, 15 panels: worst mean {max(v[1] for v in rest):.2f}, "
          f"worst SD {max(v[2] for v in rest):.2f} grid steps (asserted <= 1 step = {GRID} pt)")
    print(f"  panel h ({worst['h'][0]}) vs Table S1 LCK 394: mean {worst['h'][1]:.2f}, "
          f"SD {worst['h'][2]:.2f} grid steps (reported, not asserted)")

    note = ("Chylek LA et al. (2014) PLoS ONE 9:e104240, doi:10.1371/journal.pone.0104240")
    grid = ("MATLAB vector art on a 0.8-pt grid: +/-0.4 pt = +/-"
            f"{0.4 / np.mean([P['pt_per_log2'] for P in panels]):.3f} log2 units, +/-"
            f"{0.4 / np.mean([P['pt_per_s'] for P in panels]):.2f} s")
    write_csv(REF / "chylek2014_figS6_simulated_digitized.csv",
              ["panel", "site", "label", "t_s", "log2_fold_change"],
              [(P["letter"], P["obs"], P["label"], float(t), float(y))
               for P in panels for t, y in zip(P["curve_t"], P["curve_y"])],
              comments=[note,
                        "File S2 Fig. S6 a-p: published simulation (solid line), log2 fold "
                        "change vs baseline; vertices of the vector polyline",
                        grid, "written by digitize_chylek2014.py"], fmt=".6g")
    write_csv(REF / "chylek2014_figS6_data_digitized.csv",
              ["panel", "site", "label", "t_s", "mean_log2", "sd_log2", "sd_clipped"],
              [(P["letter"], P["obs"], P["label"], t, m, s, int(c))
               for P in panels for t, m, s, c in P["points"]],
              comments=[note,
                        "File S2 Fig. S6 a-p: measured points (mean of 3 experiments) and "
                        "error bars (SD); t snapped to the stated 0/5/15/30/60 s",
                        grid, "written by digitize_chylek2014.py"], fmt=".6g")
    write_csv(REF / "chylek2014_tableS1_sites.csv",
              ["site", "label", "gene", "position", "peptide", "experiment",
               "t0_s", "t5_s", "t15_s", "t30_s", "t60_s"],
              [(P["obs"], P["label"], P["gene"], P["pos"], table[P["obs"]]["peptide"], e + 1,
                *[float(v) for v in table[P["obs"]]["exps"][e]])
               for P in panels for e in range(3)],
              comments=[note,
                        "File S1 Table S1: log2 SILAC ratios vs 0 s for the 16 Fig. S6 sites, "
                        "per LC-MS/MS experiment (blank = not quantified)",
                        "Fig. S6's LCK Y424 is Table S1's LCK position 394",
                        "written by digitize_chylek2014.py"], fmt=".6g")


if __name__ == "__main__":
    main(sys.argv)
