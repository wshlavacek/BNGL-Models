# Snippet to add to examples/real-world/_manifest.py, inside EXAMPLES, in the correct
# simulator group (# ODE / # SSA / # NF section comments). See references/real-world-anatomy.md.
#
# The test cross-checks: resolved model.stochastic == stochastic (True iff ssa/nf);
# simulator drives the NF-only guards; heavy=True excludes it from the executable tier-2
# set. `observables` are the .exp column names (functions listed WITHOUT parentheses).

# --- ODE, experimental data (minimal) ---
RealWorldExample(
    folder='<name>', conf='<name>.conf', simulator='ode',
    observables=('<obs1>', '<obs2>'),
    system='<biology> (<First-author Year>, PMC#######); ODE, <protocol>'),

# --- SSA / NFsim variant ---
# RealWorldExample(
#     folder='<name>', conf='<name>.conf', simulator='nf',   # or 'ssa'
#     observables=('<obs1>',),
#     system='<biology> (<First-author Year>, PMC#######); NFsim, <protocol>',
#     stochastic=True, heavy=True),                            # heavy if cluster-scale

# --- Optional parameter-recovery targets (paper's reported best-fit values) ---
# ..., recover={'<param>': <published_value>}, tol=0.5),

# --- BPSL (constraint-bearing) example ---
# A job that attaches a .prop/.con is NATIVE-ONLY (not PEtab-exportable). There is no
# dedicated manifest field for this -- note it in `system` (e.g. "... ; BPSL constraints,
# native-only"), mark it in the README coverage row, and use the export-refused / check
# assertion instead of the PEtab-lint test. See references/bpsl-constraints.md.
