# Verifying a network-free model

How to establish the two `curate-model` verification levels when the reaction network cannot be
generated. The default level-1 check — integrate the same equations independently — is
unavailable, because there are no equations: there is a rule set and a sampler. What replaces
the independent implementation, and how tightly it must agree, is what this document fixes.

The replacement is **RuleMonkey used as the correctness check on NFsim**: an exact network-free
method standing witness for a rejection-based one. Fourteen of the collection's fifty-four model
folders are actively network-free; seven of those document a cross-engine check. The protocol
below is generalized from `tlbr_steric_monine2010` and `erbb_receptor_signaling_creamer2012`,
which run it most completely, together with `egfr_oligomerization_mitra2019`,
`tcr_signaling_chylek2014`, `p53_nhej_dolan2015` and `lambda_switch_cortes2017`.

`skills/nfsim/SKILL.md` is authoritative for engine mechanics — the flags, `+` versus `.`, the
UTL bug, the three constructs NFsim rejects, the ~12% legacy-binary divergence. This document
owns the verification protocol built on them and does not restate them.

## Contents
1. The independence rule, network-free form
2. Choosing the exact witness
3. The three arms, and the reference of record
4. Making the exact arm affordable
5. Ensemble sizing and the metric
6. Determining the flags instead of assuming them
7. The driver-script campaign shape
8. Where the collection stands
9. Worked examples

---

## 1. The independence rule, network-free form

`stochastic-verification.md` §1 states the rule for a finite network: the independent
implementation is transcribed from the paper and never reads the generated network. Network-free
has no generated network, so the rule becomes:

**The independent check must run a different *algorithm*, not a different *build* of the same
one.**

NFsim is rejection-based — it proposes firings from pattern-matched reactant lists and rejects
proposals that violate a constraint. Rejection is exactly where a network-free run can be
silently wrong: the output is a legal-looking trajectory sampled at the wrong rate, with no `.net`
to inspect and usually no closed form to check it against. RuleMonkey computes exact event times
over every particle and rejects nothing. Both sample the same continuous-time Markov chain, so
their ensemble means must agree; where they do not, the disagreement is in the sampler, because
there is nothing else it could be.

Two NFsim builds agreeing proves the build is reproducible. `egfr_oligomerization_mitra2019` is
the demonstration: at 1 nM EGF, BNG2.pl's bundled NFsim v1.14.3 reads 280 phospho-EGFR against
RuleMonkey's 317 — **12.4% low** across the informative doses — while bngsim's NFsim core agrees
with RuleMonkey to **1.8%** on the same rules, same XML, same flags. `-utl 4`, `5` and `6` give
the same means, so traversal is not the cause. Two NFsim builds would have averaged two runs of
the same algorithm and called it agreement.

Hence the house position: **parity of record is bngsim-NFsim ≡ RuleMonkey**, and the legacy
binary is reported as a contrast, not as a target.

## 2. Choosing the exact witness

The same doctrine as `stochastic-verification.md` §2 — prefer an exact calculation to a sampled
one. RuleMonkey is exact in *algorithm* but still sampled in *output*: it carries Monte-Carlo
error just as NFsim does, so a closed form beats it. Use the strongest witness the model admits.

| witness | when | what it proves |
|---|---|---|
| **analytic equilibrium or kinetic theory** | the rules encode a theory with a closed form | rules and rate laws, with no sampling error at all |
| **a generated network at reduced size** | the network is enumerable at some smaller valence, copy number or `max_stoich`, even though the production model's is not | that the rules and the network-free sampler agree, exactly, where both can run |
| **a RuleMonkey ensemble** | no closed form, and no size at which the network can be generated | everything peculiar to the sampler: rejection bookkeeping, molecularity, traversal |
| **a second NFsim build** | never sufficient alone | that the build reproduces |

This ordering is why six of the seven actively network-free folders with no RuleMonkey arm are
nonetheless verified. The blbr/tlbr family is *built* from theory and checked against it —
`tlbr_solution_macken1982` against the paper's ODE kinetics (Eq. 11), its closed-form equilibrium
(Eqs. 13, 17) and its branching-process aggregate-size distribution (Theorem 1, Eq. 21);
`blbr_heterogeneity_goldstein1980` against the Goldstein–Wofsy equilibrium Eqs. 13–14;
`blbr_cooperativity_posner2004` against the Wofsy–Goldstein equilibrium Eqs. 10, 15, 17 solved with
`scipy.optimize`. Do not add a RuleMonkey arm where an exact result already pins the answer.

Two qualifications, both from the collection:

- **A theory arm often covers only part of the parameter range.** `tlbr_yang2008` compares NFsim
  against mean-field ODEs for the ligand bond-count distribution at β = 0.1, and says why that is
  legitimate: below the percolation transition intra-complex rejection is negligible (~0% null
  events, per the paper). Above it, the mean-field check stops applying — and that is precisely
  the regime where an exact sampler would be the only witness left. State the range your theory
  covers, so the uncovered part is visible.
- **Reducing to a generatable network is a real option.** `blbr_rings_posner1995` ships
  `verify_generate_network_ode.ipynb`, a `generate_network` ODE of its own reversible-ring rules
  under `max_stoich` truncation; `camkii_holoenzyme_activation_michalski2012`'s network-free
  dodecamer file is anchored against the exact hexamer ODE of the primary file. This is the
  cheapest exact witness available and it should be tried before RuleMonkey.

The seven folders whose active protocol is ODE or SSA on a finite network, and which merely offer
a commented-out `simulate_nf` alternative — `two_state_gene_expression_noise_munsky2012`,
`three_stage_stochastic_gene_expression_shahrezaei2008`,
`noise_induced_bistable_futile_cycle_samoilov2005`,
`demographic_noise_predator_prey_cycles_mckane2005`,
`bursty_autoregulated_gene_expression_lin2016`,
`prion_nucleated_polymerization_rubenstein2007`, `tnf_tnfr1_sequential_binding_mcmillan2021` —
are not in scope here. Their checks are `stochastic-verification.md`'s, or a plain SciPy
integration, and both are stronger.

## 3. The three arms, and the reference of record

The standard shape is three arms on the same rules and the same protocol:

| arm | role |
|---|---|
| BNG2.pl `simulate({method=>"nf"})` — the bundled legacy NFsim binary | what a reader reproduces from the `.bngl`; a **contrast** |
| bngsim `NfsimSession` | the rejection method under test |
| bngsim `RuleMonkeySession` | the exact witness |

**Say which arm produced the committed reference output, and why.** The `.gdat`/`.scan` in
`reference/` come from BNG2.pl, because that is what the model's own actions block runs — a
reader reproduces the file, not the campaign. The parity claim is a separate artifact and lives
in the driver. Commit the `writeXML()` output too (`reference/tlbr_steric_monine2010.xml`): it is
the single input all engines consume, and it is what makes them comparable.

**Bake parameters into the XML; do not `set_param` after `initialize()`.** bngsim's RuleMonkey
does not re-derive seed-species counts on `set_param`, and NFsim does not refresh a rule's rate
constant post-init (lanl/bngsim#44). The failure is silent and looks exactly like an engine
disagreement: you believe you ran a 10 nM dose, one engine actually ran the default, and the
z-scores blow up on a model that is fine. Every folder here writes the dose, the subvolume and
the ligand switch into the parameter block and calls `writeXML()` before load —
`_dose_xml` in `run_monine2010.py`, `build_xml` in `run_creamer2012.py` and
`run_tcr_chylek2014.py` — so all three arms load byte-identical initial conditions from one file.

## 4. Making the exact arm affordable

RuleMonkey steps every particle, so its cost tracks particle count rather than events at the
reaction center, and it is routinely an order of magnitude above the NFsim run it validates
(`skills/bngl/skill.md` §5.6.1). The published protocol is usually out of reach. Reduce the
*protocol*, apply the identical reduction to all three arms, and leave the rules alone.

| model | axis reduced | published | cross-check | factor |
|---|---|---|---|---|
| `erbb_receptor_signaling_creamer2012` | subvolume `f` | 0.02 | 0.002 | 10× |
| `tcr_signaling_chylek2014` | subvolume `Fx` | 0.02 | 0.005 | 4× |
| `p53_nhej_dolan2015` | time window | 2400 min | 480 min | 5× |
| `egfr_oligomerization_mitra2019` | dose set | 6-point log scan to 100 nM (~10⁶ EGF) | 0.1 / 1 / 10 nM | ~10⁵ in particles |

Also collapse a two-phase protocol to one phase. `creamer2012` and `chylek2014` both publish an
equilibration with ligand binding off followed by a stimulation; the cross-check runs
single-phase, ligand-on from t = 0. The equilibration is a large share of the cost and verifies
nothing about the sampler, because almost nothing is firing.

Three rules:

- **Reduce a parameter, never a rule.** Whatever you change must be a subvolume, a dose or a
  horizon, applied identically to every arm. The moment two arms run different rules the check
  means nothing.
- **Reduce one axis.** Then the reduction is a sentence, and a reader can judge whether the
  reduced regime still exercises the mechanism.
- **Say what you reduced and that the model is unchanged.** `dolan2015`: "A reduced window
  (t_end=480 min) keeps the exact RuleMonkey method and multi-seed ensembles tractable; the model
  is unchanged."

The consequence for `scale:` is that **the cross-check, not the model, is normally the expensive
artifact**. `creamer2012`, `chylek2014` and `dolan2015` all carry `scale: minutes` on the `.bngl`
and `scale: hours` on the driver script. Class the artifact that actually pays
(`skills/bngl/skill.md` §5.6.2).

## 5. Ensemble sizing and the metric

### Sizing

Both arms are sampled, so this is a two-sample comparison of means and the standard error is
`s/√n`. Counts in current use: **24 seeds per engine** (`creamer2012`, `chylek2014`), **20**
(`monine2010`), **12** (`dolan2015`), **6** (`mitra2019`). These are deliberately small. You are
estimating a mean, not a distribution, and the comparison is then repeated across every
observable and every time point — so the *number of comparisons*, not the seed count, does most
of the work.

12–24 seeds per engine is the working range. Go below it only when the observable is nearly
deterministic (`mitra2019`'s 6 seeds against a ~500-cluster count). Go above it when the statistic
is a **binomial fraction**, whose noise is `√(p(1−p)/n)` and does not benefit from averaging over
observables: `lambda_switch_cortes2017`'s parity check runs **100 seeds per engine** because it
compares lysogeny percentages, and states the expected noise scale in the same breath as the
difference.

### The metric

The default is the **pairwise z-score of ensemble means**, computed for every engine pair ×
observable × time point:

```
z = |mean_a − mean_b| / sqrt(se_a² + se_b²)
```

**Do not accept on `max|z| < 3`.** You are taking a maximum over hundreds to thousands of
z-values; if the engines agree perfectly those values are standard normal, and the maximum grows
with the count. Every passing model in the collection would fail that test:

| model | comparisons *N* | expected max\|z\| = Φ⁻¹(1 − 1/2*N*) | observed max\|z\| | fraction \|z\| < 3 |
|---|---|---|---|---|
| `erbb_receptor_signaling_creamer2012` | 1488 (31 obs × 16 t × 3 pairs) | 3.40 | 3.48 | 0.997 |
| `tcr_signaling_chylek2014` | 702 (18 × 13 × 3) | 3.19 | 3.36 | 0.997 |
| `p53_nhej_dolan2015` | 525 (7 × 25 × 3) | 3.10 | 3.29 | 0.994 |

Report **both** numbers and judge on the pair: the fraction below 3 should sit at ≈ 0.997, and
`max|z|` should sit near `Φ⁻¹(1 − 1/2N)` rather than far above it. A real defect does not look
like this — it puts a whole observable systematically off, which collapses the fraction, not just
the maximum.

**Restrict to active observables.** `creamer2012` filters the 55 tracked sites down to the 31
whose ensemble mean ever exceeds one molecule. A z-score between two engines that both report
zero is `0/0`, and carrying those comparisons inflates the fraction-below-3 into meaninglessness.

**When the s.e.m. is tiny, drop z for an absolute difference against a physical criterion.**
`monine2010` is the case, and says so in the driver: equilibrium λ is so tightly determined
(s.e.m. ~0.001–0.005) that a sub-2% difference between independent implementations registers as
`max|z| = 19`. The meaningful statement is the absolute one — **max |Δλ| = 0.011** across
informative doses, against Monine's own acceptance criterion of RMS λ < 0.02 (SI Eq. 11).
`mitra2019` does the same with a relative difference (pEGFR 1.8%, cluster count 0.9%). This is the
`when-the-paper-is-wrong.md` §4 rule applied to agreement rather than disagreement: a number, and
a number normalized by something the paper's authors would recognize. z is the default because it
needs no such anchor; when the model hands you one, prefer it — and state that you swapped, and
why.

## 6. Determining the flags instead of assuming them

`nfsim/SKILL.md` says what `-bscb` and `-utl` do. Whether they matter **for this model** is a
measurement, and the exact arm is the instrument. Run NFsim at default flags, with `-bscb`, and
with `-bscb` plus a raised `-utl`, and compare each against the same RuleMonkey ensemble.

**First check whether the flag can apply at all — a ring in the seed species is not the same as
ring-forming rules.** `-bscb` constrains *bimolecular* reactant patterns, so a model whose rules
never create a bond gives it nothing to act on.
`camkii_holoenzyme_activation_michalski2012_nfsim.bngl` seeds a pre-formed twelve-subunit ring and
thereafter only changes subunit states: both autophosphorylation rules are single connected
two-molecule patterns joined by an explicit bond — unimolecular — so `-bscb` is inert, and the
auto-computed UTL of 2 already covers the largest pattern in the file.
`blbr_rings_posner1995` looks structurally similar and is the opposite case: its rules *build* the
rings, through bimolecular crosslinking and an explicit closure rule, and it requires
`-bscb -utl 5`. Read the rules, not the topology of the seeded species.

Where the flag can apply, three outcomes, all present in the collection:

- **The flags are no-ops — and that is a result worth having.** `creamer2012`'s adaptor layer can
  in principle bridge a receptor dimer through Grb2 and Shc, so the same-complex question is real
  and the driver refuses to assume it away. The measurement finds default, `-bscb` and
  `-bscb -utl 8` producing **bit-identical** ensembles, all at `max|z| = 2.79` against RuleMonkey.
  The conclusion — legacy NFsim at default flags, which is what the authors used, is a valid
  reference for this model — is much stronger than "we passed `-bscb` to be safe", and it is the
  statement the paper actually needs.
- **The flag is load-bearing.** In `monine2010` the acyclic TLBR assumption *is* `-bscb`. Dropping
  it moves equilibrium λ by up to **0.046** (at 0.50 nM) — four times the three-engine agreement
  scale of 0.011, and more than twice Monine's own RMS < 0.02 acceptance. A no-`bscb` run is a
  different model, not a looser run of the same one.
- **It matters in principle and barely in practice.** `mitra2019` carries dedicated ring-closure
  rules, so its crosslink rules must block same-complex binding; at the fitted parameters `-bscb`
  and no-`bscb` differ by a few percent, because higher-order rings are rarely populated. Recorded
  as a measurement, not asserted — and the contrast is committed, so a later reparameterization
  can be re-checked against it.

Give the determination its own driver mode and its own cached artifact (`reference/flags.npz`). It
answers a question a reader will ask of every network-free model in the collection, and it is not
answerable by reading the `.bngl`.

## 7. The driver-script campaign shape

Network-free verification is the canonical driver-script case in `curate-model` §"Verification
Artifact Shape": hours of ensemble runs across three engines will not fit in a notebook. Five
folders have converged on one shape — use it.

```
python run_<author><year>.py reference [N]   # faithful published protocol vs committed reference/
python run_<author><year>.py agreement [N]   # the three-engine cross-check, N seeds/engine
python run_<author><year>.py flags [N]       # -bscb / -utl determination against the exact arm
python run_<author><year>.py plot            # rebuild verify_<author><year>.png from cache alone
python run_<author><year>.py all
```

Four properties make it work:

- **Every expensive mode caches to `reference/*.npz`, and `plot` reads only the cache.** This is
  what makes the figure regenerable in a second on a machine that could not run the campaign at
  all. Cache the **per-seed ensembles**, not the summary statistics — `agreement.npz` holds the
  full `(n_seeds, n_points)` array per engine per observable, so the metric can be recomputed at a
  different threshold, or a new statistic added, without re-running hours of NFsim.
- **The seed count is an argument, not a constant.** `agreement 24`. The committed number is what
  was run; a reviewer with more machine can raise it and get a tighter answer from the same code.
- **The docstring carries the argument.** Each of these drivers opens with a paragraph per mode
  saying what is checked and why that check is the right one here. `mitra2019` and `monine2010`
  both explain *why the legacy binary is not the reference* before any code runs. That paragraph
  is the verification's reasoning and it belongs where a reader hits it first.
- **Name the cross-check mode `agreement`.** `cortes2017` calls the same thing `parity`; either
  word is fine, but one name across the collection is worth more than either.

The driver still owes everything a notebook would: both verification levels, the metrics printed
on the PNG, the verdict stated. The three-engine agreement **is** level 1 — the independent
implementation for a model that has no network to integrate — not a bonus panel.

## 8. Where the collection stands

Twenty-one folders reference `method=>"nf"`. Fourteen are **actively** network-free; the other
seven run ODE or SSA on a finite network and only offer a commented-out `simulate_nf` alternative
(§2). Of the fourteen, seven document a cross-engine check: `tlbr_steric_monine2010`,
`erbb_receptor_signaling_creamer2012`, `egfr_oligomerization_mitra2019`, `tcr_signaling_chylek2014`,
`p53_nhej_dolan2015`, `lambda_switch_cortes2017`, `lambda_switch_arkin1998`.

One of those seven is a claim without an artifact. `lambda_switch_arkin1998`'s `metadata.yaml`
says "NFsim/RuleMonkey parity retained" of the exact full-circuit variant;
`run_fullcircuit.py` takes `method` as an argument so either engine can be driven, but no
comparison is cached, and neither `fullcircuit_fig3.npz` nor `fullcircuit_fig6.npz` records which
engine produced it. **A parity claim with no number is the thing this document exists to
prevent.** When that folder is next touched it needs an `agreement` mode, or the sentence needs to
come out.

Of the seven actively network-free folders with no cross-engine arm, six are the blbr/tlbr theory
family and need none — they are checked against closed-form equilibrium or kinetic theory, which
is the stronger witness (§2). The seventh is the open case:

**`camkii_holoenzyme_activation_michalski2012_nfsim.bngl`.** The file is anchored where an exact
witness exists — its reformulation reproduces the primary file's exact hexamer ODE — but its
reason for existing is the dodecamer-to-100-mer range that `generate_network` cannot reach
(44,368 dodecamer configurations), and there it has no exact witness. It is the cheapest
outstanding cross-check in the collection: the model carries ~2,500 dodecamers over six seconds
of simulated time, so a five-configuration, eight-seed campaign runs in about nine seconds.
It is not a flags problem — see below — it is simply a check that has never been committed.

## 9. Worked examples

- **`tlbr_steric_monine2010`** — the template for an aggregation model with a load-bearing flag.
  Three engines, all `-bscb`, agreeing to max |Δλ| = 0.011 against Monine's own RMS < 0.02
  criterion; dropping `-bscb` moves λ by 0.046. Also the case for abandoning z when the s.e.m. is
  tiny, and for saying so out loud rather than reporting a `max|z|` of 19 as a failure.
- **`erbb_receptor_signaling_creamer2012`** — the template for a large multisite model. 24 seeds
  per engine at a 10×-reduced subvolume, single-phase; `max|z| = 3.48` over 1488 comparisons with
  99.7% below 3; and the flag determination that licenses legacy NFsim at default flags — the
  authors' own configuration — as a valid reference.
- **`egfr_oligomerization_mitra2019`** — the case that proves the independence rule. bngsim NFsim
  ≡ RuleMonkey to 1.8% on pEGFR while BNG2.pl's bundled v1.14.3 sits 12.4% low, with `-utl 4/5/6`
  ruled out as the cause.
- **`tcr_signaling_chylek2014`** — the same shape on a model with no same-complex ambiguity at
  all: the combinatorial complexity is per-protein multisite phosphorylation, so no flags are
  needed, and the driver's opening paragraph says that rather than passing `-bscb` defensively.
- **`p53_nhej_dolan2015`** — reduction on the time axis instead of the copy-number axis, plus a
  mean-field ODE arm that is a genuinely independent third check, with the low-copy species
  expected to sit *above* the mean-field value and that offset explained rather than tolerated.
- **`lambda_switch_cortes2017`** — parity on a **binary** statistic (lysogeny fraction), where the
  right ensemble size is set by `√(p(1−p)/n)` and the answer is 100 seeds per engine, not 24.
