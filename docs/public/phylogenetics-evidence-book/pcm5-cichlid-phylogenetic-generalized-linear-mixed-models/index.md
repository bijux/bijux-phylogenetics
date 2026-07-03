---
title: PCM5 Cichlid Phylogenetic Generalized Linear Mixed Models
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# PCM5: Cichlid Phylogenetic Generalized Linear Mixed Models

PCM5 is the most demanding lecture family in the current evidence book. It
combines posterior phylogenetic signal, ancestral-state extraction, count
models, binomial-style models, and repeated measures.

## Why PCM5 Matters

This study family is where the documentation has to stay especially clear.
Without that clarity, it is easy to blur together:

- Gaussian versus non-Gaussian models
- logit-scale versus probability-scale summaries
- species-level versus observation-level analyses
- repository-owned approximations versus the lecture's posterior mixed-model
  workflow

## What The Repository Already Has

- a complete checked-in study dossier for the lecture family
- dataset and provenance registration for the cichlid workflow
- an evidence decomposition that keeps the Poisson, multinomial, ancestral, and
  repeated-measure questions separate
- nearby binary comparative regression support, though not a full PGLMM
  replacement

## Read This Study Next

1. [generalized mixed-model evidence boundary](generalized-mixed-model-evidence-boundary.md)

## Checked-In Study Dossier

The full dossier lives in the repository under
[`evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/README.md).

## Current Repository State

PCM5 is deployed as a source-backed study family, but the posterior PGLMM core
remains outside the repository's governed parity surface today.
