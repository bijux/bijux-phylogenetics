---
title: PCM3 Primate Continuous and Categorical Regression
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# PCM3: Primate Continuous And Categorical Regression

PCM3 is the multivariable regression study in the evidence book. It keeps the
same primate lecture system but moves into categorical predictors, interaction
terms, standardized predictors, and a phylogenetically corrected regression.

## Why PCM3 Matters

This is the lecture family where model interpretation becomes harder:

- the intercept changes meaning when the response is transformed or predictors
  are standardized
- categorical predictors require an explicit encoding scheme
- interaction terms must be interpreted as combinations of coefficients rather
  than one slope
- phylogenetic correction changes the regression story again

## What The Repository Already Has

- a deployed study dossier for the full lecture family
- repository-owned comparative regression and PGLS surfaces that match much of
  the lecture's analytical shape
- explicit evidence decomposition for the lecture sections that still need
  governed parity packaging

## Read This Study Next

1. [model design and phylogenetic regression](model-design-and-phylogenetic-regression.md)

## Checked-In Study Dossier

The full dossier lives in the repository under
[`evidence-book/studies/primate-continuous-and-categorical-regression/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-continuous-and-categorical-regression/README.md).

## Current Repository State

PCM3 is documented and deployed as a study family. What is still missing is a
checked-in governed parity bundle for each lecture slice.
