# BNGL Formal Grammar (ANTLR 3)

## Source

Copied from the BioNetGen repository:

```
~/Code/bionetgen/parsers/BNGParser/src/bngparser/grammars/
```

- **Commit:** `9709520b` — `Fix tfun() parsing to handle complex expressions correctly`
- **Date copied:** 2026-03-19
- **Parser language:** Java (ANTLR 3 target)

These `.g` files are the formal grammar that the BioNetGen ANTLR parser uses.
They are the authoritative specification of BNGL syntax. The Perl parser in
`bng2/Perl2/` is the runtime engine but has no formal grammar file — these
ANTLR grammars are the closest thing to a spec.

## File inventory

| File | Lines | Purpose |
|---|---|---|
| `BNGLexer.g` | 275 | Lexical tokens: operators, keywords, literals |
| `BNGGrammar.g` | 288 | Top-level structure: `begin model`/`end model`, block dispatch |
| `BNGGrammar_Parameters.g` | 43 | `begin parameters` block |
| `BNGGrammar_MoleculeDef.g` | 68 | `begin molecule types` block: sites, state labels |
| `BNGGrammar_SeedSpecies.g` | 344 | `begin seed species` / `begin species`: patterns, bonds, compartment assignment |
| `BNGGrammar_Observables.g` | 156 | `begin observables`: type (Molecules/Species), patterns, relational operators |
| `BNGGrammar_Expression.g` | 204 | Arithmetic expressions, built-in functions, conditionals |
| `BNGGrammar_ReactionRules.g` | 383 | `begin reaction rules`, `begin population maps`: reactants, products, rate laws |
| `BNGGrammar_Actions.g` | 712 | All action directives: simulate, generate_network, parameter_scan, I/O, etc. |
| `BNGTree.g` | 186 | AST tree walker (post-parse evaluation) |

**Total:** ~2,659 lines

## Architecture

```
BNGLexer.g          (tokenizer)
    |
BNGGrammar.g        (top-level parser, imports all sub-grammars)
    |--- BNGGrammar_Parameters.g
    |--- BNGGrammar_MoleculeDef.g
    |--- BNGGrammar_SeedSpecies.g    (also defines species_def used by rules & observables)
    |--- BNGGrammar_Observables.g
    |--- BNGGrammar_ReactionRules.g  (also defines population_maps_block)
    |--- BNGGrammar_Expression.g     (shared by all blocks that accept expressions)
    |--- BNGGrammar_Actions.g
    |
BNGTree.g           (AST walker, partially implemented)
```

## House style subset

The house style (skill.md) uses a subset of the full BNGL grammar. The
relevant grammar rules for each house style block:

| House style block | Grammar file | Key rules |
|---|---|---|
| `begin model` / `end model` | BNGGrammar.g | `prog` |
| `begin parameters` | BNGGrammar_Parameters.g | `parameters_block`, `parameter_def` |
| `begin compartments` | BNGGrammar.g | `compartments_block`, `compartment` |
| `begin molecule types` | BNGGrammar_MoleculeDef.g | `molecule_types_block`, `molecule_def`, `site_def` |
| `begin seed species` | BNGGrammar_SeedSpecies.g | `seed_species_block`, `seed_species_def`, `species_def2`, `species_element`, `DOLLAR` (constant species) |
| `begin observables` | BNGGrammar_Observables.g | `observables_block`, `observable_def_line`, `observable_type` |
| `begin functions` | BNGGrammar.g | `functions_block`, `function_def` |
| `begin energy patterns` | BNGGrammar.g | `energy_patterns_block` |
| ~~`begin population maps`~~ | ~~BNGGrammar_ReactionRules.g~~ | Out of scope |
| `begin reaction rules` | BNGGrammar_ReactionRules.g | `reaction_rules_block`, `reaction_rule_def`, `reaction_def` |
| `begin actions` | BNGGrammar_Actions.g | `actions_block`, `action`, `simulate_method`, `generate_network` |
| Expressions | BNGGrammar_Expression.g | `expression`, `expression2`, arithmetic chain |

### Features NOT used by house style

The following grammar features are valid BNGL but outside house style scope:

- `readFile` action (BNGGrammar.g: `read_file`)
- `generate_hybrid_model` action
- `visualize` action
- `writeMfile`, `writeMexfile`, `writeSSC` actions
- `simulate_pla` (PLA method)
- `population_maps_block` (hybrid particle/population simulation)
- Net grammar mode (`gParent.netGrammar` branches)

### Key syntax details from the grammar

1. **Block order is not enforced** by the parser (`program_block` is a choice,
   not a sequence). House style enforces order via §2.

2. **Comments are discarded** by the lexer (`LINE_COMMENT: '#' ~('\n'|'\r')* {skip();}`).
   Our structured comments (`#@key: value`) are invisible to the parser —
   they're a convention layered on top.

3. **`begin species` and `begin seed species`** are both valid
   (`BEGIN (SEED)? SPECIES`). House style prefers `seed species`.

4. **Line continuations** are handled by the lexer:
   `ULB:('\\'(' ')*'\r'?'\n'(WS)*) {skip();}`

5. **Compartment assignment** uses `@` syntax: `species_element ... AT STRING`.

6. **Bond notation**: `!` followed by `+` (any bond), `?` (wildcard), or
   an integer/string bond name.

7. **State labels**: `~` followed by STRING, INT, or `?`.

8. **Reaction arrows**: `->` (UNI_REACTION_SIGN), `<->` (BI_REACTION_SIGN).

9. **Rate law functions**: Built-in `Sat`, `MM`, `Hill`, `Arrhenius` with
   fixed argument counts. User-defined functions are also supported.

10. **Expression built-ins**: `sin`, `cos`, `tan`, `exp`, `ln`, `log10`,
    `log2`, `sqrt`, `abs`, trig/hyperbolic inverses, `if()`, `min`, `max`,
    `sum`, `avg`, `time()`, `_pi`, `_e`.

## Updating

If the BioNetGen grammar changes upstream, re-copy the `.g` files and update
the commit hash above. The grammar is stable — changes are infrequent.
