---
title: Release Surface Alignment
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Release Surface Alignment

Release readiness in this repository depends on narrative alignment as much as
on green checks.

## Surfaces That Must Agree

- root `README.md`
- package `README.md` files
- public docs
- maintainer docs where release rules are described
- package `pyproject.toml` metadata
- release-readiness and publish-readiness reports

## Alignment Questions

Ask these before treating a release as serious:

- do the README and public docs name the same product surface
- do package descriptions and classifiers match the public maturity story
- do native inference claims stay bounded by the actual public runtime
  contracts
- do wrapper-backed workflows remain labeled as wrapper-backed where that is
  still the honest implementation boundary
- do benchmark and evidence docs describe differences in review language rather
  than flattening them into one generic claim

## Minimum Validation For Documentation And Metadata Changes

- run `uv run mkdocs build --strict`
- run the narrow README or documentation contract tests when they exist
- run publish-readiness or release-readiness validation when metadata or public
  package identity changes materially

## Failure Mode To Watch

The common failure mode is not missing prose. It is contradiction:

- code says one thing
- package metadata says another
- the public docs say a third

That is the state maintainers should treat as release-blocking.
