---
title: Native Benchmark Evidence
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Native Benchmark Evidence

The repository now keeps native benchmark evidence as one public trust surface,
not as one hidden maintainer-only utility.

For the practical command workflow, start with
[native benchmark review](../operations/native-benchmark-review.md). This page
explains how to interpret the benchmark family once those reports exist.

## What The Benchmark Family Covers

The native maximum-likelihood benchmark family covers:

- runtime scaling
- peak memory behavior
- topology and likelihood recovery on governed truth cases
- wrapper comparison against cached external reference lanes
- one unified suite that combines those dimensions into one review surface

The public CLI family is:

- `benchmark native-maximum-likelihood-speed`
- `benchmark native-maximum-likelihood-memory`
- `benchmark native-maximum-likelihood-accuracy`
- `benchmark native-maximum-likelihood-suite`
- `benchmark maximum-likelihood-wrapper-comparison`

## Why This Exists

Parity alone is not enough for a serious native runtime story.

The repository needs to show:

- where native code is stronger
- where native code is slower but still behaving as expected
- where model assumptions explain a difference honestly
- where one difference is a real native bug
- where one case remains unsupported instead of being silently rounded into a
  pass

## How To Read The Unified Suite

`benchmark native-maximum-likelihood-suite` is the main public summary lane.

Its governed status classes are:

- `native-advantage`
- `native-bug`
- `expected-assumption-difference`
- `unsupported-case`

Those categories are intentional. The suite is a review surface, not one
generic performance score.

## What This Does Not Claim

- it does not claim that one green suite closes every scientific question
- it does not claim that wrapper tools are irrelevant once native ownership
  exists
- it does not claim that a slower native result is automatically incorrect
- it does not claim that every model family and every benchmark lane are at
  the same maturity level
