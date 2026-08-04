"""Reproduce the six-condition Figure 2B comparison from committed source data."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
CONDITIONS = [
    "il6_1",
    "il6_10",
    "il10_1",
    "il10_10",
    "il6_1_il10_1",
    "il6_10_il10_10",
]
COLORS = {
    "il6_1": "#f48fb1",
    "il6_10": "#d32f2f",
    "il10_1": "#26c6da",
    "il10_10": "#1976d2",
    "il6_1_il10_1": "#ba68c8",
    "il6_10_il10_10": "#512da8",
}

reported = pd.read_csv(HERE / "fig2b_pSTAT_pooled.csv")
model = pd.read_csv(HERE / "author_ensemble1817_six_conditions.csv")
normalizers = {}
il6_10 = model[model.condition == "il6_10"].set_index("time")
for observable in ("Obs_pSTAT3", "Obs_pSTAT1"):
    normalizers[observable] = il6_10.loc[20, observable]

rows = []
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
for axis, (observable, label) in zip(
    axes, (("Obs_pSTAT3", "pSTAT3"), ("Obs_pSTAT1", "pSTAT1"))
):
    for condition in CONDITIONS:
        sim = model[model.condition == condition]
        data = reported[reported.condition == condition]
        prediction = np.interp(data.time, sim.time, sim[observable]) / normalizers[observable]
        sd = data[f"{label}_fit_SD"].to_numpy()
        rows.extend(
            zip(
                [condition] * len(data),
                [label] * len(data),
                data.time,
                data[f"{label}_mean"],
                prediction,
                sd,
            )
        )
        axis.plot(
            sim.time,
            sim[observable] / normalizers[observable],
            color=COLORS[condition],
            lw=1.3,
            label=condition,
        )
        axis.errorbar(
            data.time,
            data[f"{label}_mean"],
            yerr=sd,
            fmt="o",
            mfc="none",
            color=COLORS[condition],
            ms=4,
            capsize=2,
        )
    axis.set(xlabel="time (min)", ylabel=f"normalized {label}", title=label)

comparison = pd.DataFrame(
    rows, columns=["condition", "observable", "time", "data", "model", "SD"]
)
comparison["z"] = (comparison.model - comparison.data) / comparison.SD
comparison["relative_error"] = np.abs(comparison.model - comparison.data) / np.maximum(
    comparison.data, 0.02
)
chi_square = float(np.sum(comparison.z**2))
median_relative_error = float(np.median(comparison.relative_error))

axes[0].legend(fontsize=7, ncol=2)
fig.suptitle(
    "Cheemalavagu et al. (2024), Fig. 2B\n"
    "lines = representative published ensemble member; open markers = pooled data"
)
fig.savefig(HERE / "jak_stat_fig2b_reproduction.png", dpi=180)

print(f"chi-square (84 points): {chi_square:.6g}")
print(f"median relative error: {median_relative_error:.3%}")
assert chi_square < 100
