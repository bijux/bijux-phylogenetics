---
title: Release and Publication Gates
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Release And Publication Gates

Release readiness is governed, not assumed.

The repository checks package boundaries, shared standards sync, package build
proof, evidence-governance rules, and publication-readiness reports before
claiming a clean release state.

That is especially important here because the repository has both runtime and
evidence responsibilities, and those do not mature at the same speed.

## Release Questions This Repository Has To Answer

- do the documented public packages, READMEs, and handbook pages describe the
  runtime truthfully
- do package metadata, release descriptions, and public docs describe the same
  product surface
- do native inference and benchmark claims stay bounded by explicit runtime
  contracts rather than by implication
- do wrapper-backed workflows remain clearly labeled as wrapper orchestration
  when that is still the honest boundary
- do publication-readiness and evidence-governance reports stay green on the
  repository state being released

## Practical Reading Rule

Release readiness here is not one badge and not one green test lane.

The serious release story has to hold across:

- package metadata
- root and package READMEs
- public and internal docs
- package build verification
- release-readiness and publish-readiness reports
- runtime, parity, benchmark, and evidence-governance surfaces
