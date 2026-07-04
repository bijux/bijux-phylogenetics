# Changelog

## 0.1.2 - 2026-07-05

### Changed

- Restored the managed badge block on the alias package's single README
  surface so the shorter install-name page keeps release, docs, and package
  links without diverging into a second PyPI-only document.
- Raised the canonical runtime dependency floor to `bijux-phylogenetics
  >=0.1.2,<1.0` so the alias package resolves to the next release line and the
  repository-owned publishability checks enforce the same contract.

## 0.1.1 - 2026-07-05

### Changed

- Raised the alias package metadata bar on PyPI with a clearer summary,
  broader discovery keywords, a richer project-link set, and classifiers for
  Python `3.12`, `3.13`, and `3.14`.
- Reworked the PyPI long description so the alias package publishes a clean
  text-first explanation of what changes and what does not change when users
  choose the shorter distribution name.
- Raised the canonical-runtime dependency floor to `bijux-phylogenetics
  >=0.1.1` so the alias package cannot resolve to the older release line after
  the `0.1.1` publication.

## 0.1.0 - 2026-07-04

### Added

- Added `phylogenetic` as the short compatibility distribution for the
  canonical `bijux-phylogenetics` runtime.
- Added the shorter package and CLI import surface while preserving the same
  benchmark, workflow, and reporting surfaces as the canonical runtime.

### Changed

- Aligned alias-package docs, package metadata, and release navigation with the
  canonical runtime so users can choose the shorter install name without
  losing the public trust and workflow guidance.
