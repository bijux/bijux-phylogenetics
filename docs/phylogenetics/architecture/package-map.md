---
title: Package Map
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Package Map

The repository intentionally exposes three different package roles, not three
equivalent public products.

## Canonical Runtime Package

`bijux-phylogenetics` is the public center of gravity.

Use it when you want:

- the durable Python API under `bijux_phylogenetics`
- the canonical CLI entry point
- the public runtime documentation surface
- the repository's native inference, workflow, benchmark, and reporting story

This is the package public documentation should describe first, because it
owns the product meaning of the repository.

## Compatibility Alias Package

`phylogenetic` is the shorter compatibility install surface.

Use it when you want:

- the shorter distribution name
- the shorter CLI name
- the same runtime family without a forked implementation story

What it does not mean:

- it is not a second independent runtime
- it is not a separate documentation boundary
- it should not diverge in capability claims from `bijux-phylogenetics`

## Maintainer Package

`bijux-phylogenetics-dev` is repository-owned maintainer tooling.

It exists for:

- documentation contract checks
- release and publish-readiness support
- repository quality and governance automation
- package and docs integrity guardrails

It is not an end-user install surface and should not be described as if it
were part of the public scientific runtime.

## Package Choice Rule

Choose packages by durable responsibility:

- install `bijux-phylogenetics` for normal runtime use
- install `phylogenetic` only when the shorter compatibility name is the real
  requirement
- use `bijux-phylogenetics-dev` only when maintaining this repository

That split keeps the public product story stable instead of letting package
names drift into overlapping meanings.
