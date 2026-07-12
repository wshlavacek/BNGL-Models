#!/usr/bin/env python
"""Reproduction figure for the igf1r job (Erickson et al. 2019 fit of the Kiselyov 2009 model).

Simulates the IGF1/IGF1R competition assay with the deterministic ODE method at the PyBNF
best-fit parameters, across the Fig-5B cold-ligand titration, and overlays the model on the
Kiselyov data. Hot ("labelled") IGF1 is fixed at 7 pM (= 8852 copies/cell, the seed default);
cold IGF1_cold_conc is scanned over the F5B doses; each dose is integrated to steady state
(t_end = 14400 s) and the readout IGF1_hot_bound is taken at t_end.

This job is NATIVE-ONLY (`normalization = init`, not PEtab-exportable): the simulated
IGF1_hot_bound is normalized to its no-competitor (first-scan-row) value before comparison to
the pre-normalized F5B data -- exactly `Data.normalize_to_init` (pybnf/data.py:389, divide each
column by row 0). The whole titration is run with one BNGL `parameter_scan` (as the classic
model's actions block generated the data). Requires BNGPATH set (uses BNG2.pl) and matplotlib/numpy.
Usage:  BNGPATH=... python make_reproduction.py
"""
import glob
import os
import re
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "igf1r.bngl")

# PyBNF best-fit (this repo's Scatter Search + Simplex run of igf1r.conf; chi_sq recorded in
# README.md). There is no published 3-param (K1/K2/K1prime) best-fit to cite -- the Mitra 2019
# corpus (RuleHub 15-igf1r) fits the fuller 7-rate-constant model to three datasets, and the
# shipped nominals (K1=9.2nM, K2=483nM: Kiselyov 2009 Table-1 KDs) do NOT reproduce the
# normalized F5B curve under this reduced parameterization. These fitted affinities do.
BESTFIT = dict(K1=5.419475967259732e-09, K2=0.00837237152600558, K1prime=1.4377272779785423e-08)
T_END = 14400
HOT_COPIES = 8852        # 7 pM labelled IGF1 (the model's seed default; set explicitly for parity)

# F5B doses (M), in F5B.exp row order -- the scan order, and row 0 is the `init` reference.
DOSES = [2.5822e-12, 6.9016e-12, 1.3065e-11, 1.3178e-12, 2.6498e-11, 6.6104e-11, 1.2953e-10,
         2.5822e-10, 6.5536e-10, 1.2841e-9, 2.6498e-9, 6.7254e-9, 1.3178e-8, 2.5822e-8,
         6.5536e-8, 1.2841e-7, 2.6498e-7, 6.6104e-7]


def scan_hot_bound():
    """One ODE parameter_scan over the F5B cold doses; returns raw IGF1_hot_bound per dose."""
    with open(MODEL) as fh:
        src = fh.read().split("end model")[0] + "end model\n"
    for k, v in BESTFIT.items():   # override the fitted params (nominal Kiselyov KD -> best-fit)
        src = re.sub(rf"(^{k}\s+)[\d.eE+-]+", rf"\g<1>{v:.10g}", src, count=1, flags=re.M)
    vals = ",".join(f"{d:.10g}" for d in DOSES)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "igf1r.bngl")
        with open(path, "w") as fh:
            fh.write(src)
            fh.write(f'\nsetConcentration("IGF1(ds,hs,label~hot)",{HOT_COPIES})\n'
                     f'parameter_scan({{suffix=>"F5B",parameter=>"IGF1_cold_conc",'
                     f'par_scan_vals=>[{vals}],method=>"ode",t_start=>0,t_end=>{T_END},'
                     f'n_steps=>10,steady_state=>0}})\n')
        subprocess.run(["perl", BNG, path], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        scan = np.loadtxt(glob.glob(os.path.join(d, "*.scan"))[0])
    return scan[:, 1]   # IGF1_hot_bound column


def main():
    if BESTFIT["K1"] is None:
        raise SystemExit("Fill BESTFIT with the fitted K1/K2/K1prime (see README.md) first.")
    exp = np.loadtxt(os.path.join(HERE, "F5B.exp"))
    dose, y_data, sd = exp[:, 0], exp[:, 1], exp[:, 2]

    raw = scan_hot_bound()
    y_model = raw / raw[0]          # normalization = init: divide by the first scan row (pybnf/data.py:389)

    chisq = 0.5 * float(np.sum(((y_model - y_data) / sd) ** 2))
    # relative error only where the data are meaningfully above the near-zero tail
    big = y_data > 0.05
    rel = np.abs(y_model[big] - y_data[big]) / y_data[big]
    print(f"igf1r ODE reproduction at the PyBNF best-fit "
          f"(K1={BESTFIT['K1']:.3g} M, K2={BESTFIT['K2']:.3g} M, K1prime={BESTFIT['K1prime']:.3g}):")
    print(f"  chi_sq (fit objective)      = {chisq:.4f}")
    print(f"  median |rel err| (y>0.05)   = {np.median(rel):.3f}")
    print(f"  max    |rel err| (y>0.05)   = {rel.max():.3f}")

    order = np.argsort(dose)   # plot sorted by dose (data rows are not monotonic)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.errorbar(dose[order], y_data[order], yerr=sd[order], fmt="o", color="#222", ms=7,
                capsize=3, zorder=5, label="Kiselyov 2009 Fig. 5B (mean +/- SD)")
    ax.plot(dose[order], y_model[order], "s-", color="#1f77b4", lw=2, ms=6,
            label="ODE @ PyBNF best-fit (chi_sq=%.2f)" % chisq)
    ax.set(xscale="log", xlabel="cold (unlabelled) IGF1  IGF1_cold_conc (M)",
           ylabel="IGF1_hot_bound  (normalized to no-competitor)",
           title="IGF1 / IGF1R competition -- deterministic ODE\n"
                 "Erickson 2019 fit of the Kiselyov 2009 model (PyBNF corpus 15-igf1r)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    out = os.path.join(HERE, "igf1r_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
