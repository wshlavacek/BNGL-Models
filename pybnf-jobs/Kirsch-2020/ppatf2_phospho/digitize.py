#!/usr/bin/env python
"""Digitize the Fig. 7b pp-ATF2(pT69/pT71) panel -> the fit targets for ppatf2_phospho.

The pp-ATF2(T69/T71) PHOSPHORYLATION time course is NOT in the paper's Source Data
workbook (which has sheets Figure_1..6 + S*, but no Figure_7). It is published only as the
Fig. 7b "pp-ATF2 (pT69/pT71)" panel (absolute uM), which the authors plot in two
conditions:
  * CTR      -- control (blue squares); the model's `pT69pT71` at published params.
  * JNK-IN-8 -- JNK inhibitor (red diamonds); the model with k1 = k2 = 0.

Digitization method (validate-pybnf-job / curate-model discipline):
  * Rendered p.10 of the article PDF at 400 dpi (pdftoppm) and cropped the pp-ATF2 panel.
  * Axes are LINEAR. y = pp-ATF2 (uM), ticks 0/0.2/0.4/0.6/0.8/1.0; x = t(min), ticks
    0/10/20/30/40. Data points sit at t = 0, 10, 20, 40 min (the western-blot lanes within
    the 0-40 min window).
  * Points read by eye against the gridline ticks; marker centers to the nearest ~0.01 uM.
  * Tolerance: ~+/-0.03 uM (marker half-width on a clean linear axis). No error bars are
    drawn on this panel, so the fit uses the `sos` objective (unweighted least squares).

t is written in seconds (t[s] = t[min]*60) to match the model's simulation clock.
This script only records the read values and writes the two .exp files; re-run it to
regenerate them. Stdlib only.
"""

# (t_min, pp-ATF2 uM) read from the Fig. 7b pp-ATF2(pT69/pT71) panel:
CTR = [(0, 0.21), (10, 0.55), (20, 0.72), (40, 0.80)]   # blue squares (control)
JNKI = [(0, 0.085), (10, 0.24), (20, 0.42), (40, 0.42)]  # red diamonds (JNK-IN-8: k1=k2=0)


def write_exp(path, points):
    with open(path, "w") as fh:
        fh.write("# time pT69pT71\n")
        for t_min, uM in points:
            fh.write(f"{round(t_min*60)} {uM:.4f}\n")
    print(f"wrote {path} ({len(points)} points)")


if __name__ == "__main__":
    write_exp("ppatf2_phospho_ctr.exp", CTR)
    write_exp("ppatf2_phospho_jnki.exp", JNKI)
