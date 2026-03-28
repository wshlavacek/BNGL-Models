# NFsim Flags And Molecularity

This note is for source-level confirmation of the BNGL/NFsim flag behavior that matters for cyclic crosslinking models.

## BioNetGen 2.9.3 Wrapper

- `Perl2/BNGAction.pm`
  - `simulate_nf` maps `complex => -cb`
  - `simulate_nf` maps `print_functions => -ogf`
  - raw `param` text is appended directly to the NFsim command line
- `simulate({method=>"nf", ...})` forwards to the same NFsim wrapper, so the same mappings and `param` behavior apply there too.

Useful inspection command:

```sh
rg -n "complex|print_functions|param" /Users/wish/Simulations/BioNetGen-2.9.3/Perl2/BNGAction.pm
```

## What `-cb` Actually Does

- NFsim complex bookkeeping is the machinery needed for species observables and complex-aware matching.
- In BioNetGen `simulate_nf`, `complex=>1` only requests NFsim `-cb`.
- `-cb` is not the same thing as enforcing distinct reactant complexes for a bimolecular rule.

## What `-bscb` Adds

- `src/NFsim.cpp`
  - `-bscb` sets `blockSameComplexBinding = true`
  - if either `-cb` or `-bscb` is present, NFsim enables complex bookkeeping for the loaded XML model
- `src/NFreactions/transformations/transformationSet.cpp`
  - with complex bookkeeping enabled for the transformation set, `checkMolecularity()` rejects bimolecular matches when two reactant patterns point to the same complex
  - without that path, NFsim can use the weaker collision-only check, which only rejects overlapping reaction centers

Useful inspection commands:

```sh
rg -n "bscb|cb" /Users/wish/Code/nfsim/src/NFsim.cpp
rg -n "checkMolecularity|complex_bookkeeping|check_collisions" /Users/wish/Code/nfsim/src/NFreactions/transformations/transformationSet.cpp
```

Interpretation:

- Use `-cb` when you need bookkeeping.
- Use `-bscb` when model correctness depends on reactant-side distinct-complex enforcement.
- In direct NFsim XML runs, `-bscb` already implies the bookkeeping needed for that enforcement.

## Species Observables

- `src/NFinput/NFinput.cpp`
  - NFsim warns that species observables require complex bookkeeping
  - if species observables are present and bookkeeping is off, NFsim auto-enables it

Useful inspection command:

```sh
rg -n "Species Observable|Auto-enabling complex bookkeeping" /Users/wish/Code/nfsim/src/NFinput/NFinput.cpp
```

## Size-2 Ring Validation Pattern

For the bivalent ligand / bivalent receptor cyclic dimer model:

- avoid reversible shorthand if the reverse reading would incorrectly restrict dissociation
- add `param=>"-bscb"` to BNGL NFsim actions
- compare NFsim ensemble means against Posner-style kinetic and equilibrium theory, not against one stochastic trajectory
- if needed, keep chain bonds and ring-closing bonds in separate site states so observables can distinguish open chains from closed rings
- add `-utl 5` when the ring closure pattern spans 4 molecules (see [utl-traversal-bug.md](utl-traversal-bug.md))
