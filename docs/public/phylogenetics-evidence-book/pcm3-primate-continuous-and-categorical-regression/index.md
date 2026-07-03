---
title: PCM3 Primate Continuous and Categorical Regression
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# PCM3: Primate Continuous And Categorical Regression

PCM3 is the multivariable regression study in the evidence book. It stays with
the primate dataset, but the main question changes: how should readers interpret
models once continuous predictors, categorical predictors, interactions, and
phylogenetic correction all appear together?

```mermaid
flowchart LR
    A[Standardized predictors] --> B[Categorical predictors]
    B --> C[Interaction terms]
    C --> D[Phylogenetic regression]
```

## What PCM3 Answers

This study is useful when you need to understand:

- how the intercept changes meaning after standardization or transformation
- how categorical predictors depend on a clear coding choice
- how interactions should be read as joint effects instead of isolated slopes
- how the regression story changes again after phylogenetic correction

## Read PCM3 In Order

1. [model design and phylogenetic regression](model-design-and-phylogenetic-regression.md)

## What The Repository Already Provides

- a checked-in source record for the full lecture family
- repository-owned comparative regression and PGLS surfaces that match much of
  the lecture's analytical shape
- one governed parity bundle ID for each lecture slice, `evidence-001` through
  `evidence-008`
- a study breakdown that keeps each lecture question separate instead of
  collapsing everything into one regression claim

## What PCM3 Does Not Prove

- It does not turn deployed parity-style bundle IDs into a blanket claim of
  full numerical parity.
- It does not mean each PCM3 bundle already closes its lecture-specific
  coefficient or diagnostic verdict.
- It does not replace the need to inspect the checked-in source record when you
  need detailed provenance.

## Checked-In Evidence Record

The full checked-in record lives under
[`evidence-book/studies/primate-continuous-and-categorical-regression/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-continuous-and-categorical-regression/README.md).
