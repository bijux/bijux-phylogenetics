---
title: PCM3 Model Design and Phylogenetic Regression
audience: public
type: study-guide
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# PCM3 Model Design And Phylogenetic Regression

PCM3 is about model structure, not just model output. The lecture walks through
how one primate response changes interpretation as the regression becomes more
realistic.

## Main Analytical Steps

- Fit a baseline GLS with one continuous and one categorical predictor.
- Transform the response to improve residual behavior.
- Remove the intercept when group-specific estimates are the real target.
- Add an interaction term to ask whether one slope depends on the categorical
  grouping.
- Standardize the continuous predictor to change the intercept meaning.
- Refit the model with a phylogenetic covariance structure through
  `corMartins`.

## What Makes This Lecture Family Valuable

PCM3 is where the repository has to be explicit about:

- how categorical levels are encoded
- which level is treated as the reference group
- how interaction coefficients should be read
- whether taxa were excluded before or after the phylogenetic correction
- how diagnostic improvement is separated from biological interpretation

## What The Repository Can Already Support

Nearby runtime surfaces already exist for:

- categorical contrast reporting
- interaction-coefficient ledgers
- comparative model-selection ledgers
- residual review after phylogenetic regression

Those surfaces make PCM3 more than a documentation-only study, even though the
lecture-specific governed parity bundles are still to be checked in.

## Boundary To Keep In Mind

This page documents the lecture family and the available runtime shape. It does
not yet claim one finished PCM3 parity suite equivalent to the PCM1 and PCM2
governed bundles.
