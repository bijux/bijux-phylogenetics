---
title: PCM4 Cichlid Bayesian Linear and Phylogenetic Mixed Models
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# PCM4: Cichlid Bayesian Linear And Phylogenetic Mixed Models

PCM4 introduces the cichlid parasite dataset and the first Bayesian
mixed-model workflow in the course-derived evidence book.

## Why PCM4 Matters

This lecture family is where the repository stops being able to claim lecture
reproduction simply because it has nearby comparative regression surfaces. The
core analytical story depends on posterior sampling, random-effect uncertainty,
and a phylogenetic inverse covariance fed into `MCMCglmm`.

## What The Repository Already Has

- a complete checked-in study dossier for the cichlid lecture family
- explicit dataset and provenance registration
- a named evidence decomposition showing how the lecture should be reviewed
- a Brownian PGLS comparator surface that matches the lecture's final
  comparison step

## Read This Study Next

1. [mixed-model evidence boundary](mixed-model-evidence-boundary.md)

## Checked-In Study Dossier

The full dossier lives in the repository under
[`evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/README.md).

## Current Repository State

PCM4 is fully documented as a study family. The central Bayesian posterior
workflow is still an external-engine boundary rather than a current
repository-owned parity claim.
