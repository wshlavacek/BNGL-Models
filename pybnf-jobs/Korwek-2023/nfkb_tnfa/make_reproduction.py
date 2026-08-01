#!/usr/bin/env python3
"""Compare the nfkb_tnfa fit against the published parameterization of Korwek et al. (2023).

Runs the job's own model through BioNetGen twice -- once at the published maximum-likelihood
values (table S1, which are the nominal values in nfkb_tnfa.bngl) and once at the parameters
PyBNF recovered -- under the same protocol the conf synthesizes: equilibrate unstimulated for
1.2e6 s, apply TNF-alpha, follow for 6 h; the A20 KO arm carries h_A20_gene = 0 through both
phases. It then scores both against the .exp data with the job's own chi_sq objective and
writes nfkb_tnfa_reproduction.png.

The published parameterization is the oracle here in the strong sense: nfkb_tnfa.bngl is a
fitting-ready copy of a model that reproduces the authors' own data file S2 bit-for-bit, so
"published" means the authors' actual fit, not a transcription of it.

The best fit is committed as nfkb_tnfa_bestfit.txt (the header and top row of PyBNF's
sorted_params_refine_final.txt), so this runs without re-running the 22-minute fit; point
--results at a fresh output/Results/sorted_params_refine_final.txt to check a new one.

Usage (with BNGPATH set, from this folder):
    python make_reproduction.py [--results <sorted_params file>]
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
BNG2 = Path.home() / "Simulations" / "BioNetGen-2.9.3" / "BNG2.pl"

EQUIL = 1.2e6          # matches the conf's equil_t_end
T_END = 21600.0        # 6 h of TNF-alpha
N_STEPS = 2160

# exp column -> (model observable, scale parameter, arm)
COLUMNS = {
    "NFkBn_nuc": ("NFkB_nuc_total", "s_NFkBn_nuc", "wt"),
    "pIKK_fine": ("IKK_a", "s_pIKK_fine", "wt"),
    "IkBa_fine": ("IkBa_total", "s_IkBa_fine", "wt"),
    "A20_fine":  ("A20", "s_A20_fine", "wt"),
    "pIKK_long": ("IKK_a", "s_pIKK_long", "wt"),
    "IkBa_long": ("IkBa_total", "s_IkBa_long", "wt"),
    "A20_long":  ("A20", "s_A20_long", "wt"),
    "pIKK_ko":   ("IKK_a", "s_pIKK_ko", "ko"),
    "IkBa_ko":   ("IkBa_total", "s_IkBa_ko", "ko"),
}
EXP_FILES = sorted(HERE.glob("nfkb_tnfa_*_r*.exp"))

ACTIONS = """
begin actions
  generate_network({{overwrite=>1}})
  simulate({{method=>"ode",suffix=>"ode",t_end=>{equil},n_steps=>2}})
  setConcentration("TNFa()",1)
  simulate({{method=>"ode",suffix=>"ode",continue=>0,t_start=>0,t_end=>{t_end},\\
    n_steps=>{n_steps}}})
end actions
"""


def read_exp(path):
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, ndmin=2)
    return header, data


def read_gdat(path):
    with open(path) as fh:
        cols = fh.readline().lstrip("#").split()
    return cols, np.loadtxt(path)


def parse_results(path):
    """{param: value} from a PyBNF sorted_params file (header row, best fit first)."""
    lines = [l.split() for l in Path(path).read_text().splitlines() if l.strip()]
    # The header is "#<TAB>Simulation<TAB>Obj<TAB><param>..."; a data row drops the leading
    # empty field, so the row aligns with header[1:] -- ("Simulation", "Obj", params...).
    header, best = lines[0], lines[1]
    if header[0] == "#":
        header = header[1:]
    if len(header) != len(best):
        raise RuntimeError(f"{path}: header has {len(header)} fields, best row has {len(best)}")
    out = {}
    for name, val in zip(header, best):
        if name in ("Simulation", "Obj"):
            continue
        out[name] = float(val)
    return out


def simulate(params, ko, workdir, tag):
    """Run the job model with `params` overridden; returns (times, {observable: trace})."""
    text = (HERE / "nfkb_tnfa.bngl").read_text()
    over = dict(params)
    if ko:
        over["h_A20_gene"] = 0.0
    for name, value in over.items():
        pat = re.compile(rf"^(\s+{re.escape(name)}\s+)\S+(\s*#.*)?$", re.M)
        if pat.search(text) is None:
            continue
        text = pat.sub(lambda m: f"{m.group(1)}{value:.12g}{m.group(2) or ''}", text, count=1)
    text = text.split("begin actions")[0].rstrip() + ACTIONS.format(
        equil=EQUIL, t_end=T_END, n_steps=N_STEPS)
    src = workdir / f"{tag}.bngl"
    src.write_text(text)
    subprocess.run(["perl", str(BNG2), src.name], cwd=workdir, check=True,
                   capture_output=True)
    cols, arr = read_gdat(workdir / f"{tag}_ode.gdat")
    return arr[:, 0], {c: arr[:, i] for i, c in enumerate(cols)}


def predict(params, sims, column, times):
    obs, scale, arm = COLUMNS[column]
    t, traces = sims[arm]
    return params[scale] * np.interp(times, t, traces[obs])


def score(params, sims):
    """chi_sq (the job's objective) plus the per-point relative error."""
    chi2, rel, n = 0.0, [], 0
    for path in EXP_FILES:
        header, data = read_exp(path)
        for j, name in enumerate(header):
            if name == "time" or name.endswith("_SD"):
                continue
            y = data[:, j]
            sd = data[:, header.index(f"{name}_SD")]
            ok = np.isfinite(y)
            if not ok.any():
                continue
            pred = predict(params, sims, name, data[ok, 0])
            chi2 += float(np.sum(((pred - y[ok]) / sd[ok]) ** 2))
            rel += list(np.abs(pred - y[ok]) / y[ok])
            n += int(ok.sum())
    return chi2, n, np.array(rel)


def run_set(params, workdir, tag):
    return {"wt": simulate(params, False, workdir, f"{tag}_wt"),
            "ko": simulate(params, True, workdir, f"{tag}_ko")}


def refit_scales(params, sims):
    """Best per-blot normalization constants at FIXED kinetics.

    With sigma_i proportional to y_i, chi_sq is 16 * sum((s*m_i/y_i - 1)^2) in each column,
    so the optimal scale has the closed form sum(r_i) / sum(r_i^2) with r_i = m_i / y_i. This
    gives the fair intermediate between the published parameterization and the fit: the
    published kinetics seen through the normalization the data itself prefers, which isolates
    how much of the improvement is kinetics rather than blot scaling."""
    out = dict(params)
    for column, (obs, scale, arm) in COLUMNS.items():
        num = den = 0.0
        for path in EXP_FILES:
            header, data = read_exp(path)
            if column not in header:
                continue
            y = data[:, header.index(column)]
            ok = np.isfinite(y)
            if not ok.any():
                continue
            t, traces = sims[arm]
            m = np.interp(data[ok, 0], t, traces[obs])
            r = m / y[ok]
            num += float(np.sum(r))
            den += float(np.sum(r * r))
        if den > 0:
            out[scale] = num / den
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, default=HERE / "nfkb_tnfa_bestfit.txt",
                    help="a PyBNF sorted_params file; defaults to the committed best fit")
    ap.add_argument("--out", type=Path, default=HERE / "nfkb_tnfa_reproduction.png")
    args = ap.parse_args()

    published = {}
    for line in (HERE / "nfkb_tnfa.bngl").read_text().splitlines():
        m = re.match(r"\s+([A-Za-z_]\w*)\s+([-+0-9.eE]+)\s*(#.*)?$", line)
        if m:
            try:
                published[m.group(1)] = float(m.group(2))
            except ValueError:
                pass

    sets = {"published": published}
    if args.results.exists():
        fitted = dict(published)
        fitted.update(parse_results(args.results))
        sets["fitted"] = fitted
    else:
        print(f"note: {args.results} not found -- reporting the published parameterization only")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        results = {}
        for tag, params in list(sets.items()):
            sims = run_set(params, tmp, tag)
            chi2, n, rel = score(params, sims)
            results[tag] = (params, sims, chi2, n, rel)
            if tag == "published":
                # the same kinetics with only the nine blot normalizations refit
                rescaled = refit_scales(params, sims)
                c2, n2, r2 = score(rescaled, sims)
                results["published+scales"] = (rescaled, sims, c2, n2, r2)
        order = ["published", "published+scales", "fitted"]
        for tag in order:
            if tag not in results:
                continue
            _, _, chi2, n, rel = results[tag]
            print(f"{tag:18s} chi_sq = {chi2:9.2f} over {n} points   "
                  f"median |relative error| = {100*np.median(rel):5.1f}%   "
                  f"mean = {100*rel.mean():5.1f}%")

    if "fitted" in results:
        p_chi = results["published"][2]
        s_chi = results["published+scales"][2]
        f_chi = results["fitted"][2]
        print(f"\nchi_sq  published {p_chi:.1f}  ->  published with blot scales refit "
              f"{s_chi:.1f}  ->  fitted {f_chi:.1f}")
        print(f"  the fit improves on the published kinetics by {s_chi / f_chi:.2f}x once the "
              f"blot normalization is taken out of the comparison")
        print("\nrecovered kinetic parameters (fitted / published):")
        for name in ("a_Tak1_by_Tnfa", "a_Ikk", "d_Ikk_1", "d_Ikk_2", "d_Ikk_3",
                     "p_Ikba_by_Ikk", "g_Ikba_u_with_Nfkb", "i_Ikba", "e_Ikba",
                     "tg_Ikba_mrna", "a_Ikba_gene_by_Nfkb__", "s_Ikba",
                     "tg_A20_mrna", "a_A20_gene_by_Nfkb__", "sg_A20"):
            pv, fv = published[name], results["fitted"][0][name]
            print(f"  {name:24s} {fv:12.4g} / {pv:12.4g}  = {fv/pv:7.2f}x")

    # ---- figure -------------------------------------------------------------
    order = ["NFkBn_nuc", "pIKK_fine", "IkBa_fine", "A20_fine",
             "pIKK_long", "IkBa_long", "A20_long", "pIKK_ko", "IkBa_ko"]
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 9.5))
    for ax, column in zip(axes.ravel(), order):
        for tag, style in (("published+scales", dict(color="0.45", ls="--", lw=1.6)),
                           ("fitted", dict(color="tab:blue", ls="-", lw=1.8))):
            if tag not in results:
                continue
            params, sims, *_ = results[tag]
            t, _ = sims[COLUMNS[column][2]]
            grid = np.linspace(0, T_END, 400)
            ax.plot(grid / 3600, predict(params, sims, column, grid), label=tag, **style)
        for path in EXP_FILES:
            header, data = read_exp(path)
            if column not in header:
                continue
            j = header.index(column)
            y = data[:, j]
            sd = data[:, header.index(f"{column}_SD")]
            rep = path.stem.endswith("_r2")
            ax.errorbar(data[:, 0] / 3600, y, yerr=sd, fmt="^" if rep else "s",
                        ms=4, lw=0.9, color="0.6" if rep else "k",
                        label=("replicate blot" if rep else "blot shown"), alpha=.85)
        ax.set_yscale("log")
        ax.set_ylim(0.06, 20)
        ax.set_title(column, fontsize=10)
        ax.set_xlabel("time after TNF-$\\alpha$ (h)", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=.25)
    axes[0, 0].set_ylabel("relative level", fontsize=9)
    axes[0, 0].legend(fontsize=7, loc="lower left")
    chi_txt = "   ".join(f"{k} $\\chi^2$={results[k][2]:.0f}"
                         for k in ("published+scales", "fitted") if k in results)
    fig.suptitle("Korwek et al. (2023) fig. S12 -- NF-$\\kappa$B module refit to the "
                 f"TNF-$\\alpha$ data.  {chi_txt}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.out, dpi=130)
    print(f"\nwrote {args.out.name}")


if __name__ == "__main__":
    main()
