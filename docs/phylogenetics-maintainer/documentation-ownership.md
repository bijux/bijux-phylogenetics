---
title: Documentation Ownership
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-04
---

# Documentation Ownership

This repository has several documentation layers. They do not answer the same
question, and maintainers should not let them blur together.

## Public Surface Responsibilities

- the root `README.md` introduces what the repository publishes and where the
  major boundaries are
- package `README.md` files explain the durable role of each published package
- `docs/phylogenetics/` describes the runtime product surface, public
  interfaces, operations, and quality boundaries
- `docs/phylogenetics-evidence-book/` explains study evidence, parity, and
  scientific trust boundaries
- flagship native surfaces should be visible across reference, workflow, and quality pages rather than only in one inventory page

## Reader-Facing Rule

Public documentation should help an external reader understand the runtime
depth quickly. It should not read like a private repo note, a checklist, or a
thin wrapper disclaimer when the code and tests now support a deeper story.

## Maintainer Surface Responsibilities

- `docs/phylogenetics-maintainer/` is for repository stewardship, release
  operations, and documentation alignment rules
- maintainer docs should explain how to keep the public story truthful, not
  replace the public docs

## When Runtime Ownership Changes

When the runtime gains one meaningful new owned surface, update every affected
surface in the same review pass:

- root `README.md`
- the relevant package `README.md`
- the matching public handbook page
- package metadata if description or maturity language changed
- maintainer handbook guidance if maintainers need a new review rule

That review pass should make the new surface visible in public language, not
just technically present in code or tests.

## Failure Modes To Avoid

- public docs deny a surface that the code already owns
- wrapper-backed workflows are described as if they were native
- the evidence book is used as a substitute for product docs
- package metadata lags far enough behind the docs that the release story
  becomes contradictory
