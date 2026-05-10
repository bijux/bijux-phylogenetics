# Portability Audit

This audit enforces portable path semantics for checked-in evidence payloads.

- audited json files: `139`
- report-like files: `186`
- action required: `0`

## Rules

Checked-in evidence paths must stay portable and reviewer-readable.

Allowed locator forms:

- repository-relative paths rooted at governed prefixes such as `evidence-book/`, `packages/`, `docs/`, or `artifacts/`
- explicit external locators such as `external:lund/...`

Forbidden locator forms:

- workstation-local absolute paths tied to one machine layout
- parent traversal such as `../...` that guesses sibling repositories or ambient workspace shape
- ambiguous path-like strings that are neither explicit external locators nor rooted repository-relative paths

Portable checked-in reports and plots must reference their governed sources through those same locator forms.

## Locator Kinds

- `external_locator`: `157`
- `repo_relative`: `509`
- `suspicious_path_like`: `32`

## Issues

- none
