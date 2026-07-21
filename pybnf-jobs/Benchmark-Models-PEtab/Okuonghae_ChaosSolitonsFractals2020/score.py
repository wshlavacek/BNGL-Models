#!/usr/bin/env python
"""Score this PyBNF benchmark job against the Grein et al. 2026 reference J*.

The Grein/Hasenauer optimizer benchmark scores every fit by its optimality gap
    OG = J_pybnf_paperscale - J*,      "solved" iff OG < 1.92   (chi^2, alpha=0.05, 1 dof)
where J* = min over all optimizer runs on Marvin of the Eq. 6 Gaussian negative
log-likelihood (suppl/data/best_fx_marvin.csv), and J_pybnf_paperscale is PyBNF's
best objective put on the paper's Eq. 6 scale.

FIDELITY. Both objectives are the Gaussian NLL with an ESTIMATED sigma (every
subset-I problem estimates its noise as a free `sigma`/`sd_*` parameter). PyBNF
MINIMIZES the reduced objective (it drops the parameter-independent per-point
constants -- 1/2 log(2 pi), and, for a log10-transformed observable, the
change-of-variables Jacobian sum log(y_obs * ln10)). It then reports the FULL
normalized log-likelihood at the best fit in information_criteria.txt, which RESTORES
every dropped constant (it matches scipy.stats.norm.logpdf / lognorm.logpdf). So the
paper's Eq. 6 NLL is exactly
    J_paper == -log_likelihood
for BOTH linear (`observableTransformation = lin` -> gaussian) and log10
(`observableTransformation = log10` -> lognormal, Jacobian included) observables.
This is the primary, transform-agnostic scale used here.

Cross-check (linear observables only): J_paper == J_pybnf + N * 1/2 log(2 pi). For a
log10 observable the identity carries the extra data Jacobian, so it is reported but
not asserted.

Usage:  python score.py                 # shipped provenance (best_fit_params.txt + information_criteria.txt)
        python score.py output          # a fresh run's output/ (reads output/Results/)
"""
import os, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
if len(sys.argv) > 1:
    RES = os.path.join(sys.argv[1], "Results")
    best_path = os.path.join(RES, "sorted_params_final.txt")
    ic_path = os.path.join(RES, "information_criteria.txt")
else:
    best_path = os.path.join(HERE, "best_fit_params.txt")
    ic_path = os.path.join(HERE, "information_criteria.txt")

JSTAR = float(open(os.path.join(HERE, "jstar.txt")).read().split()[0])

# --- PyBNF's best reduced objective (drops per-point constants) ---
best_line = [l for l in open(best_path) if l.strip() and not l.startswith("#")][0]
J_pybnf = float(best_line.split()[1])

# --- PyBNF's full normalized log-likelihood at the best fit (restores every constant) ---
ic = {}
for l in open(ic_path):
    if l.startswith("#") or not l.strip():
        continue
    k_, v_ = l.split()
    ic[k_] = v_
n = int(ic["n"]); k = int(ic["k"])
lnL = float(ic["log_likelihood"])

J_paper = -lnL                          # PyBNF's own paper-scale Eq. 6 NLL (primary)
C = n * 0.5 * math.log(2 * math.pi)     # the 1/2 log 2pi part; = J_paper - J_pybnf only for LINEAR obs
extra = J_paper - J_pybnf               # restored per-point constants: C (linear) or C + Jacobian (log10)

OG = J_paper - JSTAR
solved = OG < 1.92

print(f"N scored points            n = {n}")
print(f"free parameters            k = {k}")
print(f"PyBNF reduced objective       = {J_pybnf:.6f}   (minimized; drops per-point constants)")
print(f"restored constants (-lnL - J) = {extra:.6f}   (= N/2 log2pi = {C:.4f} for linear; + Jacobian for log10)")
print(f"J_paper = -log_likelihood     = {J_paper:.6f}   (PyBNF's full normalized NLL = paper Eq. 6)")
print(f"reference           J* (Grein) = {JSTAR:.6f}")
print(f"OPTIMALITY GAP  OG = J_paper-J* = {OG:.6f}")
print(f"=> {'SOLVED' if solved else 'NOT solved'} (threshold OG < 1.92)")
