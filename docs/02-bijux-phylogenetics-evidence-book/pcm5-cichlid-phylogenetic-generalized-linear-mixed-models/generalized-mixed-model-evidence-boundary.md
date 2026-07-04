---
title: PCM5 Generalized Mixed-Model Evidence Boundary
audience: public
type: study-guide
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# PCM5 Generalized Mixed-Model Evidence Boundary

PCM5 extends the cichlid workflow into generalized mixed models and posterior
phylogenetic quantities. It is the clearest example of a study family that is
fully documented but not yet fully reproduced inside the repository.

```mermaid
flowchart LR
    A[Study sources and datasets] --> B[Governed evidence bundles]
    B --> C[Documented analytical map]
    C --> D[Posterior PGLMM boundary remains open]
```

## What Is Source-Backed And Deployed

- the lecture scripts and supplementary code references
- the cichlid observation table, phylogeny, and saved workspace
- governed bundle IDs for each checked-in evidence unit
- the evidence map for posterior signal, ancestral values, count models,
  infection models, and repeated-measure random effects

## What Is Not Yet A Repository-Owned Parity Claim

- Poisson phylogenetic mixed-model posteriors
- multinomial mixed-model posteriors
- posterior ancestral-state ledgers
- probability-scale posterior slope transforms
- governed variance partitions matching the lecture's multi-random-effect
  models

## What Readers Can Safely Conclude

- the repository now documents the PCM5 lecture family at a much finer-grained
  evidence-ID level than one generic summary
- the repository preserves the study inputs and analytical decomposition well
  enough for serious review
- the remaining gap is not source intake or documentation quality; it is the
  posterior mixed-model reproduction boundary itself

## Nearby Repository Support

The repository does expose narrower comparative surfaces such as binary-response
comparative regression and reviewer-facing model-selection ledgers. Those are
useful context, but they are not substitutes for the lecture's generalized
mixed-model posterior workflow.

## Reading Rule

Treat PCM5 as a documented study boundary with a strong source map. Do not
read it as proof that the repository already owns a full PGLMM implementation.
