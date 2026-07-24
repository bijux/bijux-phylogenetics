# PCM4 Runtime Reference Mapping

`lund_pcm4_mixed_model_source.R` is the checked source definition for PCM4.
It names the observation aggregation, Gaussian baseline, Bayesian chains,
scaled tree inverse, phylogenetic mixed model, and Brownian PGLS comparator
that a governed reproduction must preserve.

There is currently no executable reference extraction or repository-owned
mixed-model implementation registered in this directory. The seven bundles
contain contract wrappers and zero governed primary outputs.

## Source Objects And Required Evidence

| Source object | Meaning | Evidence needed for comparison |
| --- | --- | --- |
| `fish` and `fishtree` | repeated observations and pruned tree | input hashes, taxon exclusions, row/tip identity |
| `fishSum` | species-level aggregation | grouping keys, aggregation rules, resulting population |
| `glsFish` | Gaussian ML baseline | formula, coefficients, uncertainty, likelihood, diagnostics |
| `m1`, `m1a`, `m1b`, `m1c` | Gaussian posterior chains | priors, chain state, draws, summaries, convergence diagnostics |
| `inv.fishphylo` | scaled inverse tree covariance | tree/tip order, scaling, matrix identity, random-effect levels |
| `m2` | phylogenetic mixed-model posterior | fixed effects, variance components, draws, diagnostics |
| `pgls6` | Brownian ML comparator | covariance, coefficients, likelihood, diagnostics |

## Ownership Boundary

```mermaid
flowchart LR
    source["PCM4 lecture source"] --> rexec["future governed<br/>MCMCglmm extraction"]
    runtime["future Bijux mixed-model<br/>state and execution"] --> comparison["parameter and diagnostic<br/>comparison"]
    rexec --> comparison
    adjacent["existing native Bayesian,<br/>PGLS and diagnostics"] -. context only .-> runtime
    comparison --> bundles["seven claim-scoped<br/>bundle verdicts"]
```

`MCMCglmm` owns the preserved lecture computation. Bijux owns its current
native Bayesian, comparative, and diagnostic implementations. The evidence
bundle owns correspondence. None of those ownership statements creates the
missing Bijux PCM4 mixed-model estimand.

## Why Adjacent Methods Are Insufficient

A Brownian PGLS fit can qualify the lecture's comparator claim but cannot
produce a random-effect posterior. Generic posterior diagnostics can assess
draws but cannot create the model that generated them. Native Bayesian DNA,
clock, or continuous-trait models have different state and likelihood
contracts. Reusing their sampler components does not make their posteriors
comparable with `m2`.

## Admission Requirements

A future reference path must emit immutable, provenance-bound outputs for all
monitored parameters and chains. A future Bijux path must expose the same
estimands under an explicit state, prior, covariance, and execution contract.
Only after both sides retain paired outputs may the bundle define
normalization, posterior comparison rules, exclusions, and verdicts.

Until that boundary is crossed, contract wrapper execution is useful integrity
evidence but not scientific posterior evidence.

## Compare Posterior Records At Their Own Level

| Comparison | Minimum retained observations |
| --- | --- |
| fixed effects | parameter key, scale, draws or lossless summaries, interval rule and chain denominator |
| variance components | random/residual component identity, prior, boundary behavior and posterior draws |
| chain behavior | per-chain trace identity, retained iterations, ESS, autocorrelation and between-chain diagnostics |
| covariance construction | tree/tip order, inverse entries, scaling and random-effect level mapping |
| PGLS comparator | its own likelihood, coefficients, uncertainty, covariance and optimizer state |

The Gaussian baseline, Bayesian Gaussian model, phylogenetic mixed model, and
PGLS comparator must not share a generic coefficient row. Equal parameter
labels can refer to different conditional or marginal estimands.

A posterior comparison rule must also state whether it evaluates means,
quantiles, intervals, draws, distributional distance, or diagnostic thresholds.
Agreement on one projection does not establish agreement on the others.
