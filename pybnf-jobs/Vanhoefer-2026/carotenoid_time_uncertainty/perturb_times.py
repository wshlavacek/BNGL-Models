#!/usr/bin/env python
"""Generate the synthetic timing-perturbed datasets for the carotenoid time_error benchmark.

Reproduces the Vanhoefer et al. (bioRxiv 2026.05.09.724053) Fig-6 protocol -- the fully-synthetic
Fig-2 variant applied to the real 13-parameter carotenoid-cleavage model of Bruno et al. 2016.
From the model at theta* (`theta_star_source.txt`, the reference optimum), for each of the 77
reference measurements `(t_k, obs, sigma_k)`:

  * draw a true sampling time   tau_k = clip(t_k + N(0, sigma_t), t0, tmax)
  * value                       ybar_k = x_obs(tau_k, theta*) + N(0, sigma_k)
  * write ybar_k at the REPORTED time t_k

so theta* is the EXACT ground truth and the only corruption is the timing. `sigma_t = 0` gives
tau_k = t_k, i.e. the model's own output at theta* plus measurement noise -- the control. Each
level is deterministic (seeded), so the committed datasets are reproducible.

The trajectory x(tau, theta*) is read off a DENSE PyBNF simulation (`_simulate.conf`, which asks
for a 400-node grid over [0, 200] per condition through the time_error dense-grid action) and
linearly interpolated at tau_k -- the same conditions (init amounts + rate multipliers driven by
theta*) the fit uses, so the synthetic data is self-consistent with the model.

Usage:
    BNGPATH=$HOME/Simulations/BioNetGen-2.9.3 python perturb_times.py               # levels 0,2,5,10
    python perturb_times.py --levels 0 5 --seed 20260509
    python perturb_times.py --inspect                                               # print sim structure
"""
import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REF_DIR = HERE / "reference"
T0, TMAX = 0.0, 200.0


def load_theta_star():
    """{param: value} from the first (best) data row of theta_star_source.txt."""
    lines = [l for l in (HERE / "theta_star_source.txt").read_text().splitlines() if l.strip()]
    names = lines[0].lstrip("#").split()[2:]          # header: Simulation Obj <params...>
    values = lines[1].split()[2:]
    return {n: float(v) for n, v in zip(names, values)}


def dense_trajectories(theta):
    """Simulate _simulate.conf at theta -> {suffix: (times, {obs: values})} on the dense grid."""
    from pybnf.parse import load_config
    from pybnf.pset import PSet
    from pybnf.algorithms import core

    os.chdir(HERE)
    config = load_config("_simulate.conf")
    pset = PSet([v.set_value(theta[v.name]) for v in config.variables])
    with tempfile.TemporaryDirectory() as sim_dir:
        job = core.Job(
            list(config.models.values()), pset, "simstar", sim_dir,
            config.config["wall_time_sim"], None, config.config["normalization"],
            config.postprocessing, True,
            stochastic_seed_policy=config.config["stochastic_seed"])
        result = core.run_job(job, debug=True)
    if getattr(result, "failed", False) or result.simdata is None:
        raise RuntimeError("simulation at theta* failed: "
                           + (getattr(result, "traceback", "") or "no simdata"))
    traj = {}
    for _model, by_suffix in result.simdata.items():
        for suffix, data in by_suffix.items():
            indvar = min(data.cols, key=data.cols.get)
            times = np.asarray(data[indvar], dtype=float)
            obs = {name: np.asarray(data[name], dtype=float)
                   for name in data.cols if name != indvar}
            traj[suffix] = (times, obs)
    return traj


def read_reference(path):
    """(header names, {name: col}, rows) for one reference .exp."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    header = lines[0].lstrip("#").split()
    rows = np.array([[float(x) for x in l.split()] for l in lines[1:]], dtype=float)
    return header, {n: j for j, n in enumerate(header)}, rows


def match_suffix(traj, dataset_token):
    """The trajectory suffix for experiment `experiment____<token>` under its OWN condition.

    PyBNF simulates the model once per (experiment x condition-mutant) pair, so the suffix is the
    experiment suffix followed by the applied mutant's suffix; the scored trajectory is the
    DIAGONAL -- experiment `experiment____model1_dataN` under condition `model1_dataN` -- exactly
    as ``_condition_for_suffix`` matches them at fit time. Off-diagonal pairs (experiment N run
    under condition M) are simulated but never scored, and are not what this dataset's data is."""
    diagonal = f"experiment____{dataset_token}{dataset_token}"
    if diagonal not in traj:
        raise RuntimeError(f"no diagonal trajectory {diagonal!r} in {sorted(traj)}")
    return diagonal


def perturb_dataset(ref_path, times, obs_traj, sigma_t, noise_rng, time_rng):
    """Return the perturbed rows (same shape/columns as the reference) for one dataset.

    ``noise_rng`` is level-INDEPENDENT (the measurement noise, drawn in a fixed order), so every
    sigma_t level shares one measurement-noise realization and the ONLY difference between levels
    is the timing perturbation drawn from ``time_rng`` -- the controlled comparison the standard-vs-
    marginal MSE gate needs. At sigma_t = 0 no timing draw is taken (tau_k = t_k), so sigma_t = 0
    is exactly "the model at theta* plus that same measurement noise"."""
    header, cols, rows = read_reference(ref_path)
    indvar = header[0]
    observables = [n for n in header if n != indvar and not n.endswith("_SD")]
    out = rows.copy()
    for r in range(rows.shape[0]):
        t_k = rows[r, cols[indvar]]
        # One latent sampling time per measurement ROW (a single sample is drawn at one instant,
        # shared across that row's observables); clip to the marginal support [t0, tmax].
        tau_k = t_k if sigma_t == 0 else float(np.clip(t_k + time_rng.normal(0.0, sigma_t), T0, TMAX))
        for obs in observables:
            if obs not in obs_traj:
                raise RuntimeError(f"observable {obs!r} not simulated (have {list(obs_traj)})")
            clean = float(np.interp(tau_k, times, obs_traj[obs]))
            sigma_k = rows[r, cols[obs + "_SD"]]
            out[r, cols[obs]] = clean + noise_rng.normal(0.0, sigma_k)
    return header, out


def write_exp(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("# " + "\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(f"{x:.6g}" for x in row) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=float, nargs="+", default=[0.0, 2.0, 5.0, 10.0],
                    help="sigma_t levels (minutes) to generate")
    ap.add_argument("--seed", type=int, default=20260509, help="base seed (from the bioRxiv date)")
    ap.add_argument("--inspect", action="store_true", help="print the simulated structure and exit")
    args = ap.parse_args()

    theta = load_theta_star()
    print(f"theta* ({len(theta)} params): " + ", ".join(f"{k}={v:.4g}" for k, v in theta.items()))
    traj = dense_trajectories(theta)
    if args.inspect:
        for suffix, (times, obs) in traj.items():
            print(f"  {suffix}: {len(times)} nodes over [{times.min():.1f}, {times.max():.1f}], "
                  f"obs = {sorted(obs)}")
        return 0

    ref_files = sorted(REF_DIR.glob("experiment____*.exp"))
    for level in args.levels:
        # The measurement-noise stream is the SAME across levels (fixed seed, fixed draw order),
        # so the datasets differ ONLY in the timing perturbation (a per-level stream). This isolates
        # the timing effect the standard-vs-marginal comparison measures.
        noise_rng = np.random.default_rng(args.seed)
        time_rng = np.random.default_rng(args.seed + 1 + int(round(level * 1000)))
        outdir = HERE / f"sigmaT{level:g}"
        for ref in ref_files:
            token = ref.stem.replace("experiment____", "")     # model1_dataN
            times, obs_traj = traj[match_suffix(traj, token)]
            header, rows = perturb_dataset(ref, times, obs_traj, level, noise_rng, time_rng)
            write_exp(outdir / ref.name, header, rows)
        print(f"  sigma_t = {level:g}: wrote {len(ref_files)} datasets to {outdir.name}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
