# Write and repair Lean proofs

Inspect the goal and local context before choosing a tactic. Express the
mathematical step directly: a helper lemma, a case split, an induction, an
extensionality argument, or a calculation. Search for existing lemmas before
recreating an API.

Use automation that fits the goal. Keep ordinary `simp`, `ring`, `omega`, and
similar checked tactics when they are clear and stable. Replace exploratory
search with the suggested proof when useful; do not rewrite working proofs solely
to eliminate automation. Narrow expensive automation only after evidence of cost.

Preserve the target statement. An intermediate result belongs inside the proof
or in a separately proved helper, not in a new hypothesis of the advertised
result. Do not use a placeholder definition to disguise a missing construction.
Apply the project's draft and module-visibility conventions.

- Read the [Math-to-Lean phrasebook](references/math-to-lean.md) when translating
  an unfamiliar mathematical step.
- Read the relevant section of [Tactics and automation](references/tactics-and-automation.md)
  when a tactic or normal form is the obstacle.
- For module visibility problems, consult the bundled [Module system](../module-system.md)
  and the pinned toolchain before changing imports or exposure.
