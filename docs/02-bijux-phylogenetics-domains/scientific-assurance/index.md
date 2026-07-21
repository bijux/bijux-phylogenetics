---
title: Scientific Assurance
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Scientific Assurance

Scientific assurance connects controlled inputs to bounded claims. Simulation
and curated datasets define what is known or held fixed; validation and
benchmarks test named contracts; reporting preserves the provenance,
diagnostics, and limitations needed for review.

```mermaid
flowchart LR
    corpus["Known-answer · recovery<br/>reference · stress · biological"]
    method["Owned or external<br/>method execution"]
    checks["Invariants · tolerances<br/>diagnostics · metrics"]
    artifact["Structured artifacts"]
    claim["Bounded conclusion"]

    corpus --> method --> checks --> artifact --> claim
```

- [Simulation And Datasets](simulation-and-datasets.md) explains generators,
  known truth, catalog roles, bundle provenance, and circularity controls.
- [Validation, Benchmarks, And Reporting](validation-benchmarks-and-reporting.md)
  explains proof layers, experiment design, implementation comparison,
  artifacts, and promotion into evidence.

## Assurance Rule

Every conclusion names its corpus, method, acceptance rule, result, and
revision. Passing one lane never silently grants structural, numerical,
benchmark, evidence, or release confidence in another lane.
