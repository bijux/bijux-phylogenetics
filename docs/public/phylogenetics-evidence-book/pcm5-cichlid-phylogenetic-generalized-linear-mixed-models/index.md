---
title: PCM5 Cichlid Phylogenetic Generalized Linear Mixed Models
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# PCM5: Cichlid Phylogenetic Generalized Linear Mixed Models

PCM5 is the most demanding lecture family in the current evidence book. It
combines posterior phylogenetic signal, ancestral-state extraction, count
models, binomial-style models, and repeated measures.

```mermaid
flowchart LR
    A[Count models] --> B[Binary-style models]
    B --> C[Repeated measures]
    C --> D[Posterior signal and ancestral summaries]
```

## What PCM5 Answers

PCM5 needs especially careful documentation because it is easy to blur together
distinct model families and distinct scientific questions:

- Gaussian versus non-Gaussian models
- logit-scale versus probability-scale summaries
- species-level versus observation-level analyses
- repository-owned approximations versus the lecture's posterior mixed-model workflow

## Read PCM5 In Order

1. [generalized mixed-model evidence boundary](generalized-mixed-model-evidence-boundary.md)

## What The Repository Already Provides

- a complete checked-in source record for the lecture family
- dataset and provenance registration for the cichlid workflow
- an evidence decomposition that keeps the Poisson, multinomial, ancestral, and
  repeated-measure questions separate
- one governed parity bundle ID for each lecture slice, `evidence-001` through
  `evidence-009`
- nearby binary comparative regression support, though not a full PGLMM
  replacement

## What PCM5 Does Not Prove

- It does not claim a full repository-owned replacement for the lecture's
  posterior PGLMM workflow.
- It does not merge count, binary-style, and repeated-measure results into one
  oversimplified summary.
- It does not treat a nearby binary regression surface as proof of full mixed
  model reproduction.

## Checked-In Evidence Record

The full checked-in record lives under
[`evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/README.md).
