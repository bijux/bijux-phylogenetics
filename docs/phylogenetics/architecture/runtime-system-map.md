---
title: Runtime System Map
audience: public
type: explanation
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Runtime System Map

The public runtime is organized by durable capability families, not by the
sequence in which features were added.

## Core Runtime Families

- trees
  Validation, inspection, rooting, comparison, rendering, tree-set review, and
  the owned `PhyloTree` runtime.
- alignments
  FASTA validation, trimming, coding diagnostics, translation, identity review,
  and durable loaded alignment surfaces such as `DnaBinAlignment`.
- likelihood and inference
  Native finite-state likelihood foundations, native maximum-likelihood
  inference results, supported native Bayesian DNA workflows, and governed
  native benchmark families.
- comparative analysis
  Signal, Brownian and OU model fitting, PGLS, discrete Mk fitting, and
  dataset-backed comparative review.
- ancestral reconstruction and simulation
  Continuous and discrete ancestral-state inference, stochastic mapping, and
  governed trait or history simulation surfaces.
- adapters
  Governed orchestration for external tools such as MAFFT, trimAl, IQ-TREE2,
  FastTree, MrBayes, and BEAST where wrapper ownership remains the honest
  boundary.
- reports and artifacts
  Reviewer-facing bundles, figures, manifests, benchmark outputs, and
  evidence-linked reporting surfaces.

## How To Read This Map

Some families represent owned native implementation surfaces. Others represent
governed orchestration or review surfaces around external tools. The important
architectural rule is to keep those responsibilities explicit instead of
flattening them into one vague runtime claim.

The clearest current examples are:

- native maximum-likelihood and native Bayesian DNA workflows as owned runtime
  surfaces
- external-engine FASTA-to-tree and posterior-engine workflows as wrapper
  surfaces
- the evidence-book as a consumer of runtime outputs rather than one runtime
  family itself

## Read This Alongside

- [package map](package-map.md) when the question is which distribution owns
  which responsibility
- [evidence consumer boundary](evidence-consumer-boundary.md) when the question
  is runtime capability versus trust closure
