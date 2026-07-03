---
title: Surface Selection Guide
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Surface Selection Guide

Use this guide when you know what kind of work you need to do, but not which
public surface in the repository should own that work.

## Start By Question

- I want one governed end-user workflow or report artifact
  Start with the CLI.
- I want the same governed workflow behavior inside Python code
  Start with `bijux_phylogenetics.api`.
- I want one owned native runtime object or inference contract directly
  Start with the named Python native surfaces.
- I want to judge trust, parity, freshness, or evidence closure
  Start with the evidence-book and quality docs.
- I want a frozen integration contract rather than one live runtime process
  Start with the checked-in API bundle.

## Choose The Main Public Surface

### CLI

Use the CLI when you want:

- reviewer-facing artifacts
- governed workflow defaults
- shell-first execution
- release-facing command surfaces such as `benchmark` or `report release-truth`

Start from [CLI surface](cli-surface.md).

### Workflow Python

Use `bijux_phylogenetics.api` when you want:

- notebook or pipeline control
- typed workflow results
- the same governed report fields the CLI uses
- stable `write_json(...)` or `write_tsv(...)` helpers around workflow outputs

Start from [API surface](api-surface.md) and [Python surface](python-surface.md).

### Native Python

Use the named native Python surfaces when you want:

- one owned in-memory runtime object such as `PhyloTree`
- one native maximum-likelihood result contract
- one native Bayesian run contract
- one benchmark report contract

Start from [native inference and benchmarks](native-inference-and-benchmarks.md)
and the broader [Python surface](python-surface.md).

### Evidence And Trust

Use the evidence-book and quality guides when you want:

- freshness and integrity status
- benchmark interpretation
- release readiness
- bounded claim review instead of only runtime usage

Start from [Evidence Book](../../phylogenetics-evidence-book/index.md),
[native benchmark evidence](../quality/native-benchmark-evidence.md), and
[release review workflow](../quality/release-review-workflow.md).

### Frozen API Bundle

Use the checked-in API bundle when you want:

- one reviewable schema contract
- one pinned OpenAPI surface
- integration planning without assuming a live hosted service

Start from [API surface](api-surface.md).

## Practical Rules

- if the output should be one governed file or report bundle, prefer the CLI
- if the output should be one typed workflow object inside Python, prefer
  `bijux_phylogenetics.api`
- if the output should be one owned native contract rather than one full
  workflow wrapper, prefer the named native Python surfaces
- if the question is whether the repository has proved something, leave the
  runtime handbook and read the evidence or quality surfaces directly

## What Not To Do

- do not assume every internal helper is a public contract
- do not treat wrapper-backed workflows as if they were native implementation
  surfaces
- do not use the evidence-book as a substitute for runtime usage docs
- do not use one runtime page as proof of scientific closure
