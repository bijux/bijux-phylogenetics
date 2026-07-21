---
title: GEIGER Correspondence
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# GEIGER Correspondence

GEIGER is the live reference for continuous and discrete evolutionary-model
fitting. The governed registry currently contains 48 cases across 17 model
shards. It is designed to expose parameterization and optimizer differences,
not merely compare the name of the selected model.

## Continuous Models

The `fitContinuous` registry covers Brownian motion, OU, early burst, lambda,
kappa, delta, white noise, and explicit candidate-model comparisons. Depending
on the model, the comparison can include log likelihood, information criteria,
ancestral mean, evolutionary variance, attraction strength, optimum, or tree
transformation parameter.

## Discrete Models

The `fitDiscrete` registry covers ER, SYM, and ARD rate structures plus governed
lambda, kappa, delta, and early-burst transformations where reference behavior
is supported. State order and the mapping between rate-vector positions and
transitions are recorded because a numerically similar vector with different
labels is not parity.

| Review artifact | Question it answers |
| --- | --- |
| summary table | Which model shards passed, failed, or skipped? |
| observation table | Which compared value, tolerance, outcome, and mismatch reason belong to each case? |
| optimizer triage | Is disagreement associated with starts, convergence, or boundary estimates? |
| boundary warnings | Which parameters lie on or near constrained boundaries? |
| likelihood policy | Which sign, constant, and comparability rules apply? |
| model confidence | How strong and stable is the declared comparison? |
| parameterization registry | How do Bijux fields correspond to GEIGER parameters? |
| generated report | Which cases and non-claims can a reviewer read without decoding JSON? |

## Run The Live Lane

```bash
RUN_LIVE_PARITY=1 bijux-phylogenetics parity \
  --reference-source geiger-live \
  --geiger-shard geiger-fitcontinuous-model-ou \
  --summary-out artifacts/parity/geiger-summary.tsv \
  --observations-out artifacts/parity/geiger-observations.tsv \
  --optimizer-triage-out artifacts/parity/geiger-optimizer-triage.tsv \
  --parameterization-registry-out artifacts/parity/geiger-parameters.tsv \
  --generated-report-out artifacts/parity/geiger-report.md \
  --json
```

Repeat `--geiger-case` or `--geiger-shard` to select a bounded registry subset.
The live lane requires `Rscript`, GEIGER, and `RUN_LIVE_PARITY=1`.

## Programmatic Surface

| Public locator | Purpose |
| --- | --- |
| `bijux_phylogenetics.parity.geiger:list_geiger_parity_cases` | enumerate 48 registered model cases |
| `bijux_phylogenetics.parity.geiger:list_geiger_parity_shards` | enumerate 17 model-owned shards |
| `bijux_phylogenetics.parity.geiger:run_geiger_parity_cases` | execute the selected live reference and Bijux fits |
| `bijux_phylogenetics.parity.geiger:write_geiger_parity_observation_table` | persist value, tolerance, boundary, and mismatch observations |
| `bijux_phylogenetics.parity.geiger:write_geiger_optimizer_triage_table` | separate optimizer-sensitive disagreements from other outcomes |
| `bijux_phylogenetics.parity.geiger:write_geiger_parameterization_registry_table` | persist the field-to-parameter correspondence |
| `bijux_phylogenetics.parity.geiger:build_generated_geiger_parity_report` | build a reviewer projection from governed records |

The generated report is not a second adjudication layer. Its totals and
categories must be reproducible from the selected cases and observation rows.

## Reading Differences

A likelihood difference is actionable only after tree scale, trait transform,
root treatment, parameter bounds, likelihood convention, and optimizer search
are aligned. Similar AIC values can coexist with materially different parameter
estimates. Conversely, a boundary estimate can be scientifically equivalent
while failing naive relative-error comparison.

The generated GEIGER report keeps these cases separate:

- agreement within the declared value and likelihood policy;
- expected assumption or parameterization difference;
- optimizer-sensitive result requiring triage;
- reference boundary or unsupported case;
- unexplained mismatch requiring investigation.

These are diagnostic and review categories around a `passed`, `failed`, or
`skipped` harness case. They are not Evidence Book verdict names. In particular,
an expected assumption difference must remain visible even when a dedicated
comparison policy accepts the case, and an unexplained harness mismatch needs
claim-scoped review before it can be recorded as `mismatch_unexplained` in a
scientific bundle.

## Denominator And Model Identity

Report selected, executed, passed, failed, and skipped cases by model shard.
For each fit, retain the model family, parameterization registry entry,
likelihood policy, optimizer starts, bounds, and observed boundary state. Two
rows labeled “OU” are not comparable when those identities differ.

The registry does not claim universal GEIGER coverage or that agreement makes
the fitted model appropriate for a biological study.

## Model Parameter Record

Every case retains the GEIGER model name, Bijux model identity, state or trait
ordering, tree and trait fixture, transformation, root treatment, bounds,
starts, convergence state, likelihood convention, free-parameter count, and
reported parameters. Information criteria are compared only when likelihood
and parameter-count conventions are aligned.

If different optima satisfy the declared likelihood rule but parameters differ
materially, retain both facts. Do not collapse the case to “passed” without its
parameter observations and boundary classification.
