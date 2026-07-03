---
title: Artifact Hygiene
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Artifact Hygiene

Repository execution should leave behind governed outputs under `artifacts/`,
not durable clutter at the repository root.

## Canonical Paths

Prefer the canonical artifact-owned paths directly when one maintainer command
needs a concrete executable or output root:

- root environment: `artifacts/root/check-venv/`
- tox work directory: `artifacts/root/tox/`
- hypothesis state: `artifacts/root/hypothesis/`
- benchmark scratch space: `artifacts/root/benchmarks/`

Those paths are the durable contract. Root aliases exist only as disposable
compatibility conveniences.

## Disposable Root Aliases

The root `.venv`, `.tox`, `.hypothesis`, and `.benchmarks` paths should be
treated as disposable alias links into `artifacts/root/`.

- `make setup` may recreate them
- `make clean-root-artifacts` may remove them
- maintainers should not depend on them as the only way to reach the governed
  environment or cache roots

If one root alias becomes a real directory instead of a symlink, that is drift,
not an acceptable steady state.

## Zero-Clutter Rule

The repository root should stay free of ad hoc Python caches, tool caches, and
environment directories.

In practice that means:

- keep `.pytest_cache`, `.ruff_cache`, `.mypy_cache`, coverage outputs, and
  similar scratch state out of the root
- route local command outputs into `artifacts/`
- use `make clean-root-artifacts` when one lane leaves behind disposable alias
  links or stray root caches

## Review Questions

- does one new command write into `artifacts/` rather than the root
- does one new maintainer guide point to the canonical artifact path instead of
  only to `./.venv`
- does root cleanup remove disposable aliases and stray caches without touching
  durable checked-in evidence
