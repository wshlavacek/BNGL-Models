---
name: nfsim
description: Use when working on BNGL or NFsim models whose behavior depends on NFsim-specific flags, molecularity bookkeeping, species observables, raw `simulate_nf` parameters, or validation of BNGL dynamics against analytic theory or notebooks.
---

# NFsim

Use this skill for BNGL/NFsim work where the result can change because of NFsim execution details rather than the reaction rules alone. Typical triggers are `complex=>1`, `-cb`, `-bscb`, `print_functions=>1`, species observables, cyclic binding models, or questions about whether a BNGL model matches an external theory.

## Workflow

1. Read the BNGL rules and decide which bimolecular rules must forbid intracomplex matches.
2. Check observables and outputs. Species observables and some function outputs depend on complex bookkeeping.
3. Choose the NFsim flags deliberately instead of assuming BNGL defaults are enough.
4. For stochastic validation, compare ensemble means to theory. Do not judge correctness from a single NFsim seed when copy numbers are small.

## Flag Semantics

- `complex=>1` in `simulate_nf` maps to NFsim `-cb`. In BioNetGen 2.9.3 this is already the default for `simulate_nf`, so writing `complex=>1` is usually redundant.
- `complex=>1` or `-cb` turns on complex bookkeeping. It does not, by itself, enforce the full reactant-side molecularity constraint for bimolecular rules.
- `param=>"-bscb"` appends the raw NFsim flag that blocks same-complex binding. Use this when two reactant patterns must come from different complexes, especially in crosslinking or ring-forming models.
- Without `-bscb`, NFsim can fall back to an overlap-only check for bimolecular reactants. That prevents the same reaction center from being reused, but it can still allow two reactant patterns to map into the same connected complex.
- `print_functions=>1` maps to NFsim `-ogf`.
- `simulate({method=>\"nf\", ...})` and `simulate_nf({...})` use the same BioNetGen wrapper, so `param=>"-bscb"` works in either form.

## `+` Versus `.`

- A top-level `+` creates separate BNGL reactant or product patterns. BioNetGen writes one XML `<ReactantPattern>` or `<ProductPattern>` per top-level pattern.
- A top-level `.` keeps the molecules in one pattern. If explicit bonds connect them, that pattern is an ordinary bonded species graph.
- If `.` is used without explicit bonds, as in `A().B()`, NFsim treats it as a disjoint-pattern or `connected-to` pattern. On the left-hand side this means "match these molecules within the same connected complex, not necessarily through a direct bond." NFsim warns that this is dangerous and slower.
- On the left-hand side, `A()+B()` gives two reactants. With `-bscb`, NFsim can require those two reactants to come from different complexes. `A().B()` gives one reactant pattern, so `-bscb` does not split `A` and `B`; they are already part of the same reactant pattern.
- Mixed forms behave exactly the way the pattern count suggests. `A().B()+C()` means one reactant that matches `A` and `B` in the same connected complex, plus a second reactant `C`. With `-bscb`, the `A().B()` complex and `C` must be distinct.
- Multiple plus signs just add more reactant patterns. `A()+B()+C()` yields three reactants. With `-bscb`, NFsim checks those reactant patterns pairwise and rejects a firing if any two map to the same complex.
- On the right-hand side, `+` versus disconnected `.` mostly changes how BioNetGen groups product molecules in XML. NFsim changes state only through explicit operations such as `AddBond`, `DeleteBond`, state changes, compartment moves, molecule addition, or molecule deletion.
- Because of that, `A()+B() -> A().B()` is a no-op unless some other operation is present. BioNetGen warns that the rule has no transformations. Likewise, `A().B() -> A()+B()` is also a no-op if no bond is deleted.
- Unbinding is the subtle case. `A(a!1).B(b!1) -> A(a)+B(b)` and `A(a!1).B(b!1) -> A(a).B(b)` both delete a bond, but BioNetGen adds `ensureConnected="1"` to the second XML form. Current NFsim does not parse that attribute, so do not rely on disconnected product patterns to enforce a post-reaction connectivity constraint.

## Universal Traversal Limit (`-utl`)

NFsim's universal traversal limit (UTL) controls how far the post-reaction BFS walks from the reaction site when updating reaction memberships. The auto-computed value equals the largest template molecule count across all patterns. **This default is too low for models that combine multi-molecule unimolecular patterns (e.g., a 4-molecule ring closure) with simpler rules (e.g., 2-molecule dissociation).**

The problem: when a simple rule fires near a molecule that anchors a large pattern, the BFS may not reach the anchor molecule. That molecule silently drops out of the large rule's reactant list and is never re-added.

- **General guidance:** always check the auto-computed UTL and verify that increasing it does not change results. Run NFsim once with `-v` (verbose) to see the auto value:
  ```
  Universal Traversal Limit (UTL) set automatically to: 4
  ```
  Then re-run with `-utl` set one or two above that value and compare outputs. If results change, keep the higher setting.
- The auto-computed UTL equals the largest template molecule count across all patterns. It can be too low when simpler rules (capture, dissociation) fire near the anchor molecule of a larger pattern.
- For the size-2 ring model (4-molecule ring closure/opening patterns), the auto UTL is 4 and the fix is `-utl 5`.
- This is an NFsim bug (see [references/utl-traversal-bug.md](references/utl-traversal-bug.md)). The `-utl` flag is the workaround.

```bngl
simulate_nf({t_end=>600,n_steps=>240,param=>"-bscb -utl 5"})
```

## Global Molecule Limit (`gml`)

The `gml` parameter caps the total number of molecules NFsim will track. It was originally a guard for RAM-constrained machines. Modern machines have plenty of memory, so routinely set `gml` to the largest 32-bit signed integer (`2147483647`) to avoid silent truncation:

```bngl
simulate({method=>"nf",...,gml=>2147483647,param=>"-bscb -utl 5"})
```

If `gml` is omitted, BioNetGen defaults to 200000, which can silently cap molecule creation in models with large copy numbers.

## Practical Rules

- Do not use reversible shorthand if the reverse reading adds a nonphysical restriction. Write the forward binding rule and dissociation rule explicitly instead.
- If you need species observables, expect slower NFsim runs. NFsim warns because species tracking requires complex bookkeeping.
- When running NFsim directly on XML, `-bscb` is the molecularity-critical flag. It also turns on complex bookkeeping internally, so an extra `-cb` is redundant. Likewise, `complex=>1` is redundant when `-bscb` is already present.
- If your rule meaning depends on molecules being in the same connected complex without an explicit bond, write that intentionally and expect slower matching. Do not use disconnected `.` on the right-hand side as a way to enforce product-side molecularity.
- When running through BNGL actions, prefer:

```bngl
simulate_nf({t_end=>600,n_steps=>240,gml=>2147483647,print_functions=>1,param=>"-bscb -utl 5"})
```

- If you need the generic action form, this is equivalent:

```bngl
simulate({method=>"nf",t_end=>600,n_steps=>240,gml=>2147483647,print_functions=>1,param=>"-bscb -utl 5"})
```

- For direct XML runs, use a command like:

```sh
NFsim -xml model.xml -sim 600 -oSteps 240 -ogf -bscb -utl 5 -gml 2147483647 -o model.gdat
```

## Validation Pattern

- Use explicit observables for the quantities you want to compare to theory.
- For cyclic dimer models, separate chain-like bonds from ring-closing bonds if that makes observables easier to define.
- Validate both kinetics and equilibrium:
  - horizontal equilibrium level from theory,
  - deterministic or analytic time course from theory,
  - NFsim ensemble mean from the BNGL model.
- If the BNGL single-run trajectory looks noisy, increase the number of NFsim seeds before concluding the model is wrong.

## Reference

Read [references/flags-and-molecularity.md](references/flags-and-molecularity.md) when you need source-level confirmation for `-cb`, `-bscb`, species observables, or the size-2 ring validation workflow.

Read [references/plus-vs-dot.md](references/plus-vs-dot.md) when you need source-level details for `+`, `.`, disconnected patterns such as `A().B()`, mixed forms such as `A().B()+C()`, or the unsupported `ensureConnected` product-side case.

Read [references/utl-traversal-bug.md](references/utl-traversal-bug.md) for the NFsim bug where the auto-computed universal traversal limit is too shallow for multi-molecule unimolecular patterns, causing silent reactant list corruption at scale.
