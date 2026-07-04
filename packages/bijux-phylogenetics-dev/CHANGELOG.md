# Changelog

## 0.1.2 - 2026-07-05

### Changed

- Re-centered badge synchronization on the single package `README.md`
  surfaces so package publication stays aligned with the same SSOT pattern
  used in sibling Bijux repositories.
- Added tests that keep package readmes Mermaid-free for PyPI while still
  requiring the generated badge blocks and trusted-publisher-ready release
  workflow defaults.
- Updated publish-readiness policies and guards to expect the
  `bijux-phylogenetics>=0.1.2,<1.0` alias dependency contract for the next
  release line.

## 0.1.1 - 2026-07-05

### Changed

- Tightened the repository truth surface around root and handbook entry pages
  so reader-facing docs describe runtime depth, evidence limits, and release
  stewardship in final public language rather than meta-documentation.
- Continued release-preparation support for `0.1.1` by documenting the PyPI
  publication cleanup and the post-`0.1.0` documentation alignment work.

## 0.1.0 - 2026-07-04

### Added

- Added the repository-owned maintainer tooling package for documentation
  alignment, release readiness, package-bundle checks, publication guards,
  evidence governance, and standards-aware quality gates.
- Added package-install smoke proofs that build the runtime wheel and sdist,
  install them into clean environments, and verify packaged example and
  dataset resources through the public API and CLI.

### Changed

- Centralized repository truth contracts so READMEs, public docs, internal
  docs, release guides, and package metadata are checked as one governed
  maintainer surface.
- Expanded publish-readiness and release-readiness auditing to cover execution
  surfaces, package boundaries, package bundles, evidence-book closure, and
  artifact-output discipline.
