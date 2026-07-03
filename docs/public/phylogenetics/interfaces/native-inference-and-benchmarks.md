---
title: Native Inference and Benchmarks
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Native Inference And Benchmarks

This page names the owned native inference and benchmark contracts that now
belong to the public runtime story.

Those surfaces are real, but they are not one blanket claim that every
phylogenetic workflow in the repository is already native-first.

## Native Maximum-Likelihood Surface

The current public native maximum-likelihood result surface is:

- `bijux_phylogenetics.phylo.likelihood.infer_nucleotide_maximum_likelihood_result(...)`
- `bijux_phylogenetics.phylo.likelihood.infer_nucleotide_maximum_likelihood_result_from_alignment(...)`
- `bijux_phylogenetics.phylo.likelihood.NucleotideMaximumLikelihoodResult`

That result contract keeps:

- the selected model and model-selection strategy
- the final tree Newick and topology fingerprint
- the final log likelihood and fitted parameter values
- multi-start run summaries and accepted search trace rows
- optional bootstrap and SH-like support reports
- warning messages and optional wrapper-correspondence benchmark metadata

Today this documented native ML contract is strongest on the owned nucleotide
tree-inference family. The wider likelihood package contains more foundations
than this public page names directly, so callers should treat this result
surface as the stable entry point instead of assuming every lower-level helper
is equally durable.

## Native Bayesian Surface

The current public native Bayesian dispatcher is:

- `bijux_phylogenetics.bayesian.run_bayesian_inference(...)`

That dispatcher accepts one supported native model definition plus one matching
proposal schedule and returns the matching native run report.

Today the documented public boundary is intentionally narrow:

- fixed-topology DNA workflows are supported through the public dispatcher
- joint-topology DNA workflows are supported through the public dispatcher
- the dispatcher rejects unsupported model-definition and proposal-schedule
  pairings explicitly
- the public dispatcher does not route through BEAST, MrBayes, or wrapper
  imports to claim native execution

The wider `bijux_phylogenetics.bayesian` package also exposes owned priors,
proposal builders, checkpoints, adaptive tuning, burn-in policies,
independent-chain diagnostics, posterior summaries, and clock-related model
surfaces. Those are part of the native runtime story, but they should still be
read through their documented families rather than as one undifferentiated
claim that every Bayesian submodule is equally mature.

## Native Benchmark Surface

The current public native benchmark family lives under
`bijux_phylogenetics.benchmark`.

The named entry points include:

- `benchmark_native_maximum_likelihood_speed(...)`
- `benchmark_native_maximum_likelihood_memory(...)`
- `benchmark_native_maximum_likelihood_accuracy(...)`
- `benchmark_native_maximum_likelihood_suite(...)`
- `benchmark_maximum_likelihood_wrapper_comparison(...)`

The unified suite is the main public review surface when one caller wants one
governed summary across:

- native scaling behavior
- wrapper-comparison lanes
- support evidence
- truth-recovery evidence

The suite report keeps four status families explicit:

- `native-advantage`
- `native-bug`
- `expected-assumption-difference`
- `unsupported-case`

That distinction matters because the repository treats benchmark evidence as a
review surface, not as one generic pass or fail badge.

## Public Reading Rule

When you need an owned inference contract, start from these named public native
surfaces first. When you need one broader biological workflow, use the
workflow-level CLI or Python surfaces. When you need trust claims, read those
runtime surfaces together with the evidence-book and quality sections instead
of treating native ownership alone as full scientific closure.
