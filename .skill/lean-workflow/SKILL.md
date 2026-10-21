---
name: lean-workflow
description: Write, repair, and explain Lean 4 proofs using the project's pinned Lean and Mathlib. Includes finding existing lemmas, choosing tactics, and naming declarations. Use for ordinary Lean coding and formalization, rather than a dedicated performance investigation or source-fidelity audit.
---

# Lean proofs and development

Work in the existing project. Read its instructions, toolchain, Lake configuration,
and the affected module before choosing commands. This skill provides guidance;
it does not install Lean or imply that a Lean server, search service, or terminal
is available. Use the tools actually provided. If compilation is unavailable,
say which claims remain unchecked.

## Work on the immediate goal

1. Identify the exact goal and relevant definitions. For a mathematical source,
   preserve premises, quantifier order, carriers, and constant dependencies.
2. Search the project and pinned Mathlib before constructing a new proof.
   Start with a focused `rg` search or a goal-directed search such as `exact?`.
   Check the candidate's actual type, implicit arguments, and imports locally.
3. Make the smallest coherent proof or API change. Keep intermediate obligations
   in helper lemmas or the proof, not as new assumptions of the claimed theorem.
4. Validate the affected target with the project's command; otherwise use a
   focused command such as `lake build Project.Module` or `lake env lean File.lean`.
   Follow existing formatting and linter conventions. Replace exploratory search
   with its suggested proof where appropriate; ordinary `simp` and automation
   do not need a blanket rewrite into explicit tactics.
5. Report what compiled, what remains incomplete, and any change to the statement
   or trust assumptions. A successful build alone does not establish source fidelity.

Preserve dependency pins and build caches. Do not use `lake clean`, delete
upstream artifacts, or upgrade dependencies as a routine proof repair. Temporary
`sorry` experiments follow the project's draft policy and must not be reported
as proved results. Do not introduce custom axioms to close a goal.

## Read only the reference needed now

These are bundled documents, not separate skills to activate. Resolve links
relative to this skill's directory. Do not read the entire reference tree.

| Current obstacle | Read |
|---|---|
| Finding a lemma, import, or theorem name | [Search](references/search/guide.md) |
| Translating a mathematical step or debugging tactics | [Proof patterns](references/proofs/guide.md) |
| Naming a declaration, choosing a namespace, or style review | [Naming and style](references/style/guide.md) |
| Analysis/PDE/probability APIs, carriers, or recurring proof patterns | [Analysis reference map](references/analysis/guide.md) |

For a dedicated slow-build investigation, choose `lean-elaboration`; for comparing
Lean declarations with a mathematical source, choose `lean-statement-audit`.
For a project plan, module restructuring, coordinated proof campaign, dependency
graph, or public release, choose `lean-project-architecture`. These are separate
tasks, not prerequisites for fixing a proof. If another skill is unavailable,
continue within this skill's scope and report the missing specialized guidance.
