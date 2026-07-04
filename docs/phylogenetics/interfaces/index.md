---
title: Public Interfaces
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# Public Interfaces

This section describes the runtime contracts that readers are expected to use
directly.

## Interface Families

- CLI commands for end-to-end workflows and reviewer-facing files
- workflow Python surfaces under `bijux_phylogenetics.api`
- native Python contracts for runtime objects and inference results
- benchmark and evidence-adjacent runtime contracts
- emitted artifact and report families

## Interface Selection Heuristic

| If you need | Start here | Why |
| --- | --- | --- |
| one governed workflow or report | CLI | easiest path to durable artifacts |
| notebook or pipeline composition | workflow Python | typed results and stable writers |
| lower-level runtime control | native Python | direct access to owned contracts |
| trust or schema review | artifact pages | explicit output families and boundaries |

```mermaid
flowchart LR
    A[Need an interface] --> B[CLI]
    A --> C[Workflow Python]
    A --> D[Native Python]
    A --> E[Frozen API bundle]
    B --> F[Files and reports]
    C --> G[Typed workflow results]
    D --> H[Lower-level owned contracts]
    E --> I[Schema review]
```

## Start Here

- If you need help choosing an entry point, use the
  [surface selection guide](surface-selection.md).
- If you already know the entry point, continue to the matching reference
  page below.

## Section Map

- [Surface selection guide](surface-selection.md)
- [CLI surface](cli-surface.md)
- [API surface](api-surface.md)
- [Python surface](python-surface.md)
- [Native inference and benchmarks](native-inference-and-benchmarks.md)
- [Artifact consumption guide](artifact-consumption-guide.md)
- [Artifact contracts](artifact-contracts.md)

The surface selection guide and artifact consumption guide are the two fastest
ways to orient a new reader before they drop into the larger CLI and Python
reference pages.
