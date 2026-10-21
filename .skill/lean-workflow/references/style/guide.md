# Lean naming and style

Follow the repository and pinned Mathlib conventions. Look at adjacent declarations
for namespace, file location, and API shape before renaming anything.

Theorems and proofs generally use `snake_case`; type and proposition constructors
use `UpperCamelCase`; data-valued definitions generally use `lowerCamelCase`.
Use conclusion-oriented theorem names and natural namespaces. Do not rename a
public API or alter a mathematical binder merely to satisfy a cosmetic preference.

Read [Naming and style](references/naming-and-style.md) for detailed conventions.
Treat its Mathlib-specific formatting as a reference, not a mandate to reformat
an unrelated project. Use the project's actual linter command and respect its
supported Lean module conventions. Drafts and proof holes follow local policy;
this reference does not introduce a frozen-anchor process.
