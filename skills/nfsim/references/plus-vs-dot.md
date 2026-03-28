# `+` Versus `.` In BNGL/NFsim

This note explains what NFsim actually receives from BioNetGen for top-level `+` and `.` on each side of a rule.

## Core Translation

- BioNetGen stores each top-level reactant pattern separately and writes each one as its own XML `<ReactantPattern>`.
- BioNetGen does the same for products with `<ProductPattern>`.
- Therefore:
  - `A()+B()` means two patterns
  - `A().B()` means one pattern containing two molecules
  - `A().B()+C()` means two patterns
  - `A()+B()+C()` means three patterns

Relevant source:

- `Perl2/RxnRule.pm` writes one XML pattern block per entry in `rr->Reactants` and `rr->Products`.

Useful inspection command:

```sh
sed -n '1768,1800p' /Users/wish/Simulations/BioNetGen-2.9.3/Perl2/RxnRule.pm
```

## Left-Hand Side Semantics

### `A()+B()`

- NFsim receives two reactant patterns.
- It builds one mapping set per reactant pattern.
- Reactant-side molecularity checks operate across those mapping sets.
- With `-bscb`, NFsim rejects a firing if two reactant patterns map to the same complex.

Relevant source:

- `src/NFcore/reactionClass.cpp` checks molecularity before firing.
- `src/NFreactions/transformations/transformationSet.cpp` enforces unique complex ids when complex bookkeeping for molecularity is active.

### `A().B()`

- NFsim receives one reactant pattern containing two molecules with no explicit bond.
- When NFsim reads that pattern, it detects disjoint sets and converts them into `connected-to` constraints.
- At match time, NFsim traverses the bonded neighborhood of the first matched molecule and searches that connected complex for molecules matching the other disconnected pattern pieces.
- So `A().B()` on the left-hand side means "A and B must be in the same connected complex", not "A and B are independent reactants."

Relevant source:

- `src/NFinput/NFinput.cpp` warns about disjoint patterns and adds `connectedTo` links.
- `src/NFcore/reactionClass.cpp` restructures those `connectedTo` links for matching.
- `src/NFcore/templateMolecule.cpp` traverses the bonded neighborhood to satisfy the `connectedTo` part of the match.

Practical consequence:

- `-bscb` does not split apart molecules inside one disconnected reactant pattern.
- If you write `A().B()`, you have already said they belong to the same reactant pattern.

### Mixed Forms

- `A().B()+C()` means:
  - one reactant pattern matching `A` and `B` in the same connected complex
  - one separate reactant pattern matching `C`
- With `-bscb`, NFsim requires the `C` reactant to come from a different complex than the `A().B()` reactant pattern.

### Three Or More `+`

- `A()+B()+C()` creates three reactant patterns.
- With `-bscb`, NFsim checks all reactant patterns for distinct complex ids.
- This is pairwise distinctness across reactant patterns, not a separate "product molecularity" concept.

## Right-Hand Side Semantics

### `A()+B()` Versus `A().B()` With No Bond Operations

- These differ in XML grouping, but not in actual transformations.
- If there is no `AddBond`, `DeleteBond`, state change, move, add, or delete operation, NFsim has no state change to apply.
- BioNetGen warns that such rules have no transformations.

Practical consequence:

- `A()+B() -> A().B()` is not a way to force co-membership in a complex.
- `A().B() -> A()+B()` is not a way to split a complex.

### Explicit Bond Formation

- `A()+B() -> A!1.B!1` and `A().B() -> A!1.B!1` both create an `AddBond`.
- The difference is on the reactant side:
  - in the first rule, the bond forms between two reactant patterns
  - in the second rule, the bond forms within one disconnected reactant pattern that already required both molecules to be in the same connected complex

### Explicit Bond Deletion

- `A!1.B!1 -> A()+B()` produces a `DeleteBond` and two product patterns.
- `A!1.B!1 -> A().B()` produces a `DeleteBond` and one disconnected product pattern.
- BioNetGen marks the second case with `ensureConnected="1"` on the XML `DeleteBond`, meaning the bond deletion should only fire if the endpoints remain in the same product pattern after deletion.

Relevant source:

- `Perl2/RxnRule.pm` adds `ensureConnected="1"` when both post-delete endpoints map into the same product pattern.

## Important Caveat: NFsim Ignores `ensureConnected`

- Current NFsim XML parsing of `DeleteBond` reads only `site1` and `site2`.
- There is no `ensureConnected` handling in the NFsim source tree.
- Therefore, disconnected right-hand-side product grouping is not a reliable way to enforce a product-side connectivity constraint in current NFsim.

Relevant source:

- `src/NFinput/NFinput.cpp` reads `DeleteBond` attributes `site1` and `site2` only.
- searching the NFsim source for `ensureConnected` returns no runtime handling.

Practical rule:

- If correctness depends on a post-reaction connectivity condition, do not rely on disconnected `.` on the product side.
- Prefer explicit bond structure, a separate reverse/forward rule formulation, or external logic.

## Product Patterns In NFsim

- NFsim stores reactant and product templates mainly for reaction-connectivity inference, not for an independent product-side molecularity check during firing.
- The actual fired state change comes from the transformation list: add bond, delete bond, state change, delete molecule, add molecule, move compartment.

Relevant source:

- `src/NFcore/reactionClass.cpp` stores all reactant and product templates for connectivity inference.
- `src/NFreactions/transformations/transformationSet.cpp` uses transformed product templates in `checkConnection()`, which is about reaction-connectivity bookkeeping.
