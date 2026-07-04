---
title: Documentation Ownership
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Documentation Ownership

This repository has several documentation layers, and they do not answer the
same question.

## Public Surface Responsibilities

- the root `README.md` tells a new reader what the repository publishes and
  where the major boundaries are
- package `README.md` files explain the durable role of each published package
- the public docs under `docs/phylogenetics/` describe the runtime
  product surface, public interfaces, operations, and quality boundaries
- the evidence-book docs describe governed scientific trust and study evidence,
  not the entire runtime product surface
- flagship native surfaces should be visible across reference, workflow, and
  quality pages rather than only one API inventory page

## Maintainer Surface Responsibilities

- `docs/phylogenetics-maintainer/` is for maintainers, release operations, and repo-health
  contracts
- maintainer docs should explain how to keep the public story truthful, not
  replace the public docs

## When Runtime Ownership Changes

When the runtime gains one meaningful new owned surface, update all of the
following that are affected:

- root `README.md`
- the relevant package `README.md`
- the matching public handbook page
- package metadata if the package description or maturity claim changed
- maintainer handbook guidance when maintainers need one new review rule

## What To Avoid

- do not leave public docs claiming a surface does not exist after the code now
  owns it
- do not describe wrapper-backed workflows as if they were native
- do not use the evidence-book as one substitute for package or runtime docs
- do not let package metadata lag so far behind the docs that the release story
  becomes internally contradictory
