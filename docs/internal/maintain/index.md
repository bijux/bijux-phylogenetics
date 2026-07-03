---
title: Maintainer Guide
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Maintainer Guide

Use the maintainer handbook when the question is about repository health rather
than runtime product usage.

The main internal responsibilities covered here are:

- keeping public docs, package READMEs, and metadata aligned with the actual
  runtime surface
- preserving the boundary between native runtime ownership, wrapper-backed
  workflow ownership, and evidence-book closure
- preparing a release story that stays truthful across package metadata, docs,
  and release-readiness reports

## Read By Maintainer Question

- how should public docs and package READMEs map to the runtime surface:
  [documentation ownership](documentation-ownership.md)
- how do native flagship surfaces stay documented across reference, workflow,
  and quality pages:
  [native surface documentation](native-surface-documentation.md)
- how do architecture, evidence, release, and package docs stay aligned as one
  public trust story:
  [public trust surface maintenance](public-trust-surface-maintenance.md)
- how do interface, output, and metadata guides stay aligned as one public
  contract story:
  [public interface and artifact maintenance](public-interface-and-artifact-maintenance.md)
- how do I check that metadata, docs, and release surfaces all say the same
  thing:
  [release surface alignment](release-surface-alignment.md)
