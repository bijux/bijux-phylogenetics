# Changelog

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
