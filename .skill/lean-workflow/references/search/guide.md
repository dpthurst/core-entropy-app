# Find a Lean declaration

Search the project's pinned sources first. Match both the intended conclusion
and its assumptions; a familiar name or a remote search result is not enough.

- For a nearly closed goal, try `exact?` or `apply?`; inspect remaining goals.
- For rewriting, inspect the actual types and use `rw?` or `simp?` to discover
  suitable lemmas.
- For a missing name or import, use a focused `rg` search in the project and
  `.lake/packages/mathlib/Mathlib`, then read nearby declarations.
- If local discovery is insufficient and network tools are available, use a
  type-pattern or natural-language search service. Verify hits against the pinned
  dependency before relying on them. A service mentioned here is not necessarily
  configured in this agent's environment.

Read [Search examples](references/search-and-discovery.md) for naming heuristics,
namespace navigation, and search patterns. Read [Tool index](references/tools.md)
only when selecting an unfamiliar search or inspection tool. Do not exhaust a
fixed number of search methods when the correct lemma is already found.
