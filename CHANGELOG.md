# Changelog

This file records notable repository-level changes for `bijux-phylogenetics`.

It does not replace package-level release history. Versioning and package-local
release notes belong to each distribution under `packages/`.

Use this changelog for workspace changes that affect multiple packages or
change contributor and maintainer workflows across the repository.

## 0.1.2 - 2026-07-05

### Changed

- Restored managed badge blocks to the two PyPI long-description surfaces so
  the published package pages keep their quick trust and navigation cues
  without falling back to Mermaid-rendered content that PyPI does not support.
- Raised the next release line to `0.1.2` across repository and package
  fallback-version metadata so local builds, tagged releases, and publish
  checks point at one coherent post-`0.1.1` target.
- Advanced the alias-package dependency floor and the repository-owned
  publish-readiness policies to `bijux-phylogenetics>=0.1.2,<1.0` so release
  governance matches the next canonical runtime publication line.
- Added maintainer-side badge and workflow guards that keep PyPI readmes
  synchronized, Mermaid-free, and trusted-publisher-ready ahead of the
  `0.1.2` release.

## 0.1.1 - 2026-07-05

### Changed

- Prepared the `0.1.1` release line across the repository and package
  changelog surfaces so the first post-`0.1.0` maintenance release has clear,
  bounded release notes.
- Strengthened the two published package records on PyPI with more compelling
  summaries, broader discovery metadata, richer project links, and classifier
  declarations through Python `3.14`.
- Reworked the root README and the public handbook landing pages so they speak
  directly to readers about runtime breadth, evidence boundaries, and release
  trust instead of narrating how the documentation ought to be written.
- Prepared the published package descriptions for a Mermaid-free PyPI release
  surface and aligned the alias package so it resolves to the `0.1.1`
  canonical runtime line.

## 0.1.0 - 2026-07-04

### Added

- Established the repository on the shared Bijux Python workspace model with
  governed `configs/`, `makes/`, `tox`, artifact discipline, and synchronized
  standards from `bijux-std`.
- Published the canonical package `bijux-phylogenetics`, the short
  compatibility alias `phylogenetic`, and the repository-owned maintainer
  package `bijux-phylogenetics-dev`.
- Added the first broad runtime surface across tree validation, tree
  comparison, rooting, rendering, alignment diagnostics, trimming, coding
  diagnostics, distance analysis, comparative workflows, ancestral-state
  reconstruction, simulation, and governed report generation.
- Added native finite-state likelihood foundations, native maximum-likelihood
  search surfaces, native benchmark contracts, and public native Bayesian
  inference entry points for owned DNA workflows.
- Added governed adapter workflows for external tools including MAFFT,
  trimAl, IQ-TREE2, FastTree, MrBayes, and BEAST, with explicit provenance and
  boundary reporting.
- Added the evidence-book trust program with five PCM study families, governed
  bundle structure, freshness and integrity ledgers, comparative study guides,
  and public reading surfaces.
- Added packaging, smoke-install, publication-readiness, release-readiness,
  evidence-governance, artifact-governance, package-boundary, and package-bundle
  verification lanes.

### Changed

- Reworked public docs, the root handbook, and package entry surfaces around
  reader questions, trustworthy scope boundaries, and named hyperlinks instead
  of file-path-shaped navigation.
- Clarified the repository contract around native ownership versus
  wrapper-backed orchestration so runtime breadth and evidence closure are
  described honestly rather than as one undifferentiated maturity claim.
- Standardized artifact output discipline so transient local outputs belong
  under `artifacts/`, while governed comparative trust studies stay under
  `evidence-book/`.
- Expanded the comparative analysis program with stronger parity references,
  sensitivity review, multiple-testing summaries, reviewer-facing caveats, and
  packaged example resources that support installation smoke proofs.
- Hardened release and publication governance so docs, package metadata,
  standards sync, build surfaces, and evidence reports are checked together as
  one repository-level publishability story.

### Fixed

- Aligned managed shared make surfaces across sibling repositories so
  repository contracts, full-suite package overrides, and release gates agree
  on one governed build surface.
- Restored public evidence-book and README trust wording after brittle
  summary-style regressions so documentation reflects checked-in governed
  sources instead of fragile snapshot counters.
- Closed artifact-hygiene regressions around root cache directories and
  publish-readiness drift so evidence-governance and release-readiness gates
  pass on the intended repository-owned surfaces.

## Changelog Scope

Use this file for changes such as:

- root governance and contributor policy
- shared automation under `makes/`
- shared configuration under `configs/`
- root handbook and repository navigation
- repository-level CI, publish, and release process changes
- shared evidence-book and trust-surface governance

Do not use this file for changes that only affect one package release stream
unless the repository-level workflow changed too.
