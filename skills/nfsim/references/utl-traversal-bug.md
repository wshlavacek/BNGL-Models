# NFsim Universal Traversal Limit (UTL) Bug

## The bug

NFsim's auto-computed `suggestedTraversalLimit` is set to the maximum template molecule count across all patterns. This value is the BFS depth limit for collecting "product" molecules whose reaction memberships must be updated after a reaction fires. The BFS skips molecules at `currentDepth >= limit`, so a limit of N visits depths 0 through N-1 only.

For models with a multi-molecule unimolecular reactant pattern (e.g., a 4-molecule ring closure), the auto-computed limit can be too shallow. Molecules can silently drop out of reactant lists and never be re-added.

## Source locations

Auto-computation of the limit (`NFinput.cpp`, line ~3250):

```cpp
if((int)tMolecules.size()>suggestedTraversalLimit) {
    suggestedTraversalLimit = (int)tMolecules.size();
}
```

BFS cutoff (`molecule.cpp`, line ~668):

```cpp
if((depth!=ReactionClass::NO_LIMIT) && (currentDepth>=depth)) continue;
```

The limit is applied system-wide to all reactions (`system.cpp`, line ~441):

```cpp
void System::setUniversalTraversalLimit(int utl) {
    this->universalTraversalLimit = utl;
    for(rxnIter = allReactions.begin(); rxnIter != allReactions.end(); rxnIter++ )
        (*rxnIter)->setTraversalLimit(utl);
}
```

## Concrete failure scenario

Model: bivalent ligand (L) / bivalent receptor (R) with reversible ring closure of L-R-L-R chains into cyclic dimers. The ring closure reactant pattern spans 4 molecules. The auto UTL is 4.

1. A 4-molecule chain `L1-R1-L2-R2` exists. `L1` (the endpoint with a free site) is in the ring closure reactant list.
2. Capture extends the chain: free `L3` binds `R2`'s free site, making `L1-R1-L2-R2-L3` (5 molecules). `L1` is correctly removed from ring closure (pattern no longer matches).
3. The extension dissociates: `L3-R2` bond breaks. Dissociation root is `L3`. BFS from `L3`: depth 0 (`L3`) → 1 (`R2`) → 2 (`L2`) → 3 (`R1`) → **4 (`L1`) — SKIPPED**.
4. `L1` is never re-added. The restored 4-molecule chain is invisible to ring closure.

## Quantitative effect

With the reproducer model at RT = 3000 molecules:

| UTL | Cyclic dimers (NFsim) | Cyclic dimers (ODE) | Error |
|-----|----------------------|---------------------|-------|
| 4 (auto) | 1073 ± 12 | 1175 | **-8.7%** |
| 5 | 1170 ± 8 | 1175 | -0.4% |

The error is silent (no warnings), systematic, and grows with system size. At RT = 30 the bug is masked by stochastic noise.

## Workaround

First, check the auto-computed UTL by running NFsim with `-v` (verbose):

```sh
NFsim -xml model.xml -v ... 2>&1 | grep "Traversal Limit"
# Output: Universal Traversal Limit (UTL) set automatically to: 4
```

Then pass `-utl N` with N at least one above the auto value. Re-run and compare outputs — if results change, keep the higher setting:

```sh
NFsim -xml model.xml ... -bscb -utl 5
```

In BNGL:

```bngl
simulate_nf({..., param=>"-bscb -utl 5"})
```

## When to suspect this bug

- Your model has a unimolecular rule with a connected reactant pattern spanning N ≥ 3 molecules.
- Other rules (capture, dissociation, crosslinking) can fire within N bonds of the anchor molecule of that pattern.
- NFsim results show a systematic deficit in the reaction governed by the multi-molecule pattern, compared to the `generate_network` ODE or analytic theory.
- The deficit grows with system size but vanishes (or is masked) at very small molecule counts.

## GitHub issue

Filed at https://github.com/RuleWorld/nfsim — see the bug report in `models/bivalent_ligand_receptor_cyclic_dimers_posner1995/nfsim_utl_bug_report.md` for a full reproducer and suggested fix.
