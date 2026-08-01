"""Score this job's fit against the digitized data and draw the reproduction figure.

Runs the committed <slug>.bngl through BioNetGen -- not through the SciPy transcription and
not through PyBNF -- so the figure is evidence that the BNGL file itself, at its nominal
parameter values, reproduces the fit. The protocol is appended as an actions block that
mirrors what the conf synthesizes: a fixed-length unmeasured phase with one micelle
concentration, then setConcentration on Mic() and a measured phase whose clock the script
shifts back to zero to line up with the .exp times.

Usage (with BNGPATH set, from a slug folder):
    python make_reproduction.py                 # nominal values in <slug>.bngl
    python make_reproduction.py params.json     # an alternative parameter set
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
SLUG = HERE.name
BNGPATH = Path(os.environ.get("BNGPATH", Path.home() / "Simulations" / "BioNetGen-2.9.3"))
MIC = 100.0

# stem -> (event time h, micelle before, micelle after, ThT observable, panel, label)
DESIGN = {
    "on_pathway":   (0.0, 0.0, 0.0, "ThT_add", "Fig. 1a", "on pathway, no fatty acid"),
    "fa_control":   (0.0, MIC, MIC, "ThT", "Fig. 1b/c", "5 mM C12 fatty acid throughout"),
    "addition_3h":  (3.0, 0.0, MIC, "ThT_add", "Fig. 1c", "micelles added at 3 h"),
    "addition_24h": (24.0, 0.0, MIC, "ThT_add", "Fig. 1c", "micelles added at 24 h"),
    "removal_5h":   (5.0, MIC, 0.0, "ThT_rem", "Fig. 1b", "micelles removed at 5 h"),
    "removal_24h":  (24.0, MIC, 0.0, "ThT_rem", "Fig. 1b", "micelles removed at 24 h"),
}

def bngl_source(overrides):
    text = (HERE / f"{SLUG}.bngl").read_text()
    for name, value in (overrides or {}).items():
        text, n = re.subn(rf"^(\s*{re.escape(name)}\s+)\S+(\s*#.*)?$",
                          lambda m: f"{m.group(1)}{value:.10g}{m.group(2) or ''}",
                          text, count=1, flags=re.M)
        if not n:
            raise SystemExit(f"parameter {name} not found in {SLUG}.bngl")
    return text


def simulate(text, stem, t_event, mic_before, mic_after, t_end):
    """One protocol through BioNetGen; returns (time shifted to 0, {observable: array})."""
    has_off = "Foff()" in text
    body = text.split("begin actions")[0].rstrip()
    acts = ["begin actions", "  generate_network({overwrite=>1})"]
    if has_off:
        acts.append(f'  setConcentration("Mic()",{mic_before:g})')
    # Both phases write the SAME suffix, exactly as the curated model's actions block does.
    # A `continue=>1` simulate that opens a fresh suffix writes no column header, so a
    # separate "pre" suffix would leave the measured .gdat unlabelled.
    if t_event > 0:
        acts.append(f'  simulate({{method=>"ode",suffix=>"ode",t_start=>0,'
                    f't_end=>{t_event:g},n_steps=>{max(20, int(20 * t_event))}}})')
        if has_off:
            acts.append(f'  setConcentration("Mic()",{mic_after:g})')
    acts.append(f'  simulate({{method=>"ode",suffix=>"ode",t_start=>{t_event:g},'
                f't_end=>{t_event + t_end:g},n_steps=>{max(200, int(20 * t_end))},'
                f'continue=>{1 if t_event > 0 else 0}}})')
    acts.append("end actions")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "m.bngl").write_text(body + "\n\n" + "\n".join(acts) + "\n")
        r = subprocess.run(["perl", str(BNGPATH / "BNG2.pl"), "m.bngl"],
                           cwd=tmp, capture_output=True, text=True)
        gdat = tmp / "m_ode.gdat"
        if not gdat.exists():
            raise SystemExit(f"{stem}: BioNetGen failed\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        cols = open(gdat).readline().lstrip("#").split()
        arr = np.loadtxt(gdat)
    # Drop the sample written AT the event: BioNetGen's first simulate ends there and
    # records the state before setConcentration ran.
    arr = arr[arr[:, 0] > t_event + 1e-9] if t_event else arr
    return arr[:, 0] - t_event, {c: arr[:, i] for i, c in enumerate(cols)}


def main():
    overrides = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}
    text = bngl_source(overrides)
    params = dict(re.findall(r"^\s*(\w+)\s+([0-9eE.+-]+)\s*#", text, flags=re.M))
    exps = sorted(p.stem for p in HERE.glob("*.exp"))

    fig, axes = plt.subplots(1, len(exps), figsize=(4.1 * len(exps), 3.6), squeeze=False)
    total_sse, total_n, rows = 0.0, 0, []
    for ax, stem in zip(axes[0], exps):
        t_event, mb, ma, obs, panel, label = DESIGN[stem]
        data = np.loadtxt(HERE / f"{stem}.exp")
        td, yd = data[:, 0], data[:, 1]
        t, sim = simulate(text, stem, t_event, mb, ma, float(td.max()) * 1.02)

        # The ThT measurement model lives in the conf, not the model file, so it is
        # rebuilt here from the same parameters: map_on*[F] + map_off*[F'_1]. The
        # individual jobs carry one off-pathway constant, the global job two.
        off_key = {"ThT_add": "map_off1", "ThT_rem": "map_off2"}.get(obs, "map_off")
        if off_key not in params:
            off_key = "map_off"
        fp1 = sim.get("Obs_Fp1", np.zeros_like(t))
        pred_curve = (float(params["map_on"]) * sim["Obs_F"]
                      + float(params.get(off_key, 0.0)) * fp1)
        pred = np.interp(td, t, pred_curve)

        sse = float(((pred - yd) ** 2).sum())
        rel = np.abs(pred - yd) / np.maximum(np.abs(yd), 0.05)
        total_sse += sse
        total_n += len(yd)
        rows.append((stem, panel, len(yd), sse, float(np.median(rel))))

        ax.plot(td, yd, ".", ms=3, color="0.35", label="Rana et al. (2020), digitized")
        ax.plot(t, pred_curve, "-", lw=1.6, color="crimson", label="this fit (BioNetGen)")
        ax.set_title(f"{stem}\n{panel} — {label}", fontsize=8)
        ax.set_xlabel("time after event (h)" if t_event else "time (h)", fontsize=8)
        ax.set_ylabel("ThT (a.u.)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="lower right")

    rmse = np.sqrt(total_sse / total_n)
    fig.suptitle(f"{SLUG} — SSE {total_sse:.3f} over {total_n} points (RMSE {rmse:.4f})",
                 fontsize=10)
    fig.tight_layout()
    out = HERE / f"{SLUG}_reproduction.png"
    fig.savefig(out, dpi=140)

    print(f"{'experiment':16s} {'panel':10s} {'n':>4s} {'SSE':>9s} {'median |rel err|':>17s}")
    for stem, panel, n, sse, rel in rows:
        print(f"{stem:16s} {panel:10s} {n:4d} {sse:9.4f} {rel * 100:16.1f}%")
    print(f"{'TOTAL':16s} {'':10s} {total_n:4d} {total_sse:9.4f}   RMSE {rmse:.4f}")
    print(f"\nPyBNF `sos` objective for this fit = SSE/2 = {total_sse / 2:.4f}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
