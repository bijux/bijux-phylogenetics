---
title: MrBayes And BEAST Correspondence
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-22
---

# MrBayes And BEAST Correspondence

MrBayes and BEAST own their posterior computations. Bijux owns explicit
analysis preparation, governed execution, native artifact parsing, posterior
diagnostics, tree-set review, and evidence packaging. The record always says
whether an analysis was only prepared, parsed from existing outputs, or
executed by the current workflow.

## A Posterior Claim Begins With Its Sample Population

| Required identity | Review requirement |
| --- | --- |
| analytical target | likelihood, priors, clock or tree prior, topology treatment, and monitored parameters |
| execution population | independent chains, resumed segments, seeds, generations or states, and terminal status |
| selected posterior | burn-in, thinning or subsampling, failed-chain policy, and retained draw/tree identities |
| diagnostic population | exact series, estimator, threshold, and chains contributing to ESS, drift, or convergence |
| reported quantity | parameter or clade key, transformation, estimator, interval rule, units, and retained denominator |

A posterior mean without this lineage is not portable evidence. Similar means
can arise from different priors, chains, burn-in windows, transformations, or
failed-sample policies. Compare values only after their sample populations and
estimands align.

```mermaid
flowchart LR
    inputs["Alignment · partitions<br/>tree · dates · calibrations"]
    prepare["NEXUS or XML<br/>validation"]
    engine["MrBayes or BEAST<br/>external posterior"]
    native["Trace · diagnostics<br/>posterior trees"]
    select["Burn-in · subsampling<br/>chain/state identity"]
    assess["ESS · drift · convergence<br/>clade/topology sensitivity"]
    evidence["Methods · tables<br/>checksummed evidence package"]

    inputs --> prepare --> engine --> native --> select --> assess --> evidence
```

## Preparation, Execution, And Review Are Separate

| State | MrBayes | BEAST | Supported statement |
| --- | --- | --- | --- |
| prepared | `mrbayes-prepare` writes NEXUS with model, rates, chains, sampling, and burn-in policy | `beast-prepare` writes XML with substitution model, clock, tree prior, dates, calibrations, chain and logging | configuration is structurally represented; no posterior exists yet |
| inspected | prepared NEXUS and declared inputs are reviewable | `beast-xml`, calibration, and tip-date commands validate the time-tree contract | the supplied configuration passes supported checks |
| executed | `mrbayes-run` invokes `mb` and captures native outputs | `beast-run` invokes `beast` with seed, threads, and overwrite policy | the named executable ran under the retained identity |
| parsed | trace, MCMC, consensus, ESS, parameter and tree commands read supplied artifacts | log, parameter, tree, consensus and diversity commands read supplied artifacts | observations come from those artifacts, not necessarily a current run |
| packaged | report and Bayesian evidence commands join selected records | methods, diagnostics table, and evidence package join selected records | the package is internally inventoried; scientific adequacy remains reviewable |

The report must distinguish a prepared-only analysis, a parsed prior run, and
an engine execution performed by the current workflow. File presence alone is
not execution proof.

## MrBayes Correspondence

Preparation accepts aligned nucleotide or protein FASTA, optional partitions,
model and rate settings, generations, chain count, sampling/printing frequency,
and burn-in fraction. Execution retains the `mb` executable identity, version,
command, process state, timeout, native outputs, streams, manifest, and
checksums.

Review commands expose distinct posterior questions:

- `mrbayes-traces`, `mrbayes-mcmc`, and `mrbayes-ess` parse parameter and MCMC
  records;
- `mrbayes-parameters`, `mrbayes-burnin-sensitivity`, and
  `mrbayes-convergence` evaluate retained trace behavior;
- `mrbayes-trees`, `mrbayes-subsample`, `mrbayes-summarize`, and
  `mrbayes-consensus` preserve generation identity while reviewing topology;
- `mrbayes-report` joins the selected posterior trees and traces into a
  reviewer-facing report.

## BEAST Correspondence

Preparation accepts an alignment plus optional starting tree, fossil
calibrations, and tip dates, with explicit clock model, tree prior, chain
length, and logging interval. XML inspection retains taxon, character,
calibration, tip-date, chain, and logger counts plus validation issues.

Execution records the `beast` executable, version, threads, seed, overwrite
policy, command, timeout, streams, outputs, manifest, and checksums. Review is
split by ownership:

- calibration and tip-date commands verify identifiers and dating inputs;
- `beast-log`, `beast-parameters`, `beast-burnin-sensitivity`, and
  `beast-convergence` retain state-indexed trace behavior;
- `beast-trees`, `beast-subsample`, `beast-consensus`, and `beast-diversity`
  retain state identity while reviewing the posterior tree set;
- Bayesian methods and diagnostics-table commands produce publication-facing
  records from explicitly selected inputs.

## Run Identity And Resume Safety

| Identity component | MrBayes | BEAST |
| --- | --- | --- |
| prepared analysis | NEXUS checksum, taxa, partitions, model/rates, generations, chains, sampling and burn-in | XML checksum, taxa, substitution model, clock, tree prior, dates/calibrations, chain and logging |
| execution | `mb` path/version, command, streams, exit and timeout | `beast` path/version, command, threads, seed, overwrite, streams, exit and timeout |
| posterior | generation-indexed traces and trees, consensus and native sidecars | state-indexed log and trees, consensus/MCC and native sidecars |
| selection | burn-in and subsampling against retained generations | burn-in and subsampling against retained states |
| assessment | parameter ESS, trace drift, chain and clade behavior | parameter ESS, trace drift, topology, calibration and time-tree behavior |

`mrbayes-run --resume` and `beast-run --resume` reuse only a verified completed
manifest whose executable, inputs, configuration, and output hashes still
agree. Timeout markers, killed processes, malformed posterior files, or
changed identity require rejection or an explicit clean restart.

Resume therefore preserves execution lineage; it does not create a new
independent chain and must not increase the denominator in a convergence
claim.

## Run And Review

```bash
bijux-phylogenetics adapter mrbayes-prepare \
  dataset/alignment.fasta \
  --model gtr \
  --rates gamma \
  --ngen 1000000 \
  --out artifacts/bayesian/analysis.nex

bijux-phylogenetics adapter mrbayes-run \
  artifacts/bayesian/analysis.nex \
  --timeout-seconds 7200
```

```bash
bijux-phylogenetics adapter beast-prepare \
  dataset/alignment.fasta \
  --tree dataset/starting-tree.nwk \
  --calibrations dataset/calibrations.tsv \
  --clock-model strict \
  --tree-prior yule \
  --out artifacts/bayesian/analysis.xml

bijux-phylogenetics adapter beast-xml \
  artifacts/bayesian/analysis.xml \
  --json

bijux-phylogenetics adapter beast-run \
  artifacts/bayesian/analysis.xml \
  --threads 4 \
  --seed 11 \
  --timeout-seconds 7200
```

Parser and reporting commands can review governed existing outputs without
claiming that the current machine reran the posterior analysis.

## Diagnostic Interpretation

| Observation | Supports | Does not establish |
| --- | --- | --- |
| adequate ESS under a named estimator and threshold | retained draws contain enough effective information for that parameter rule | stationarity, model adequacy, or independent-chain agreement |
| stable trace after selected burn-in | the named window/drift rule did not detect instability | correct burn-in for every parameter or topology |
| similar independent-chain summaries | selected estimands agree under the comparison rule | convergence in unobserved dimensions or correctness of the target |
| stable clade frequencies | selected topology summaries are insensitive to the named comparison | identical branch lengths, dates, parameters, or biological history |
| parser acceptance | native files satisfy the supported structural grammar | current execution, adequate sampling, or a valid biological model |

A posterior comparison must align priors, likelihood, clock and tree model,
parameterization, chain/state selection, transformations, and estimands before
values are adjudicated. Otherwise the posteriors remain separate governed
observations.

## Preserve Posterior Sample Lineage

Every reported posterior quantity must resolve to the exact retained samples
that contributed to it. A summary detached from chain and state identity is
not repairable from its mean, interval, or consensus tree.

```text
reported posterior quantity
  = engine run
  + chain or state series
  + burn-in rule
  + thinning or subsampling rule
  + transformation
  + estimand
```

| Projection | Required parent record | Completeness condition |
| --- | --- | --- |
| parameter summary | native trace, chain/state column, selected window, transform, and missing/non-finite handling | retained and excluded sample counts reconcile |
| ESS or drift diagnostic | exact series, estimator/rule, threshold, and selected window | diagnostic denominator matches the reported series |
| clade frequency or consensus | posterior tree set, generation/state identity, burn-in and subsampling | accepted trees and discarded trees are both counted |
| independent-chain comparison | distinct run identities, aligned estimand, and per-chain diagnostics | resumed segments are not miscounted as independent chains |
| evidence package result | all parent records plus inventory and checksums | every table or figure resolves to its selected native inputs |

Changing burn-in, subsampling, topology filtering, or parameter transformation
creates a new projection. Preserve it alongside its rule; do not overwrite the
prior summary and keep the old interpretation.

### Admit Samples Before Summarizing Them

Posterior review has an admission boundary before any pooled estimate. For
each requested chain or run, retain whether it started, completed, parsed,
passed the declared diagnostic policy, and contributed samples. Then retain
requested, emitted, excluded, and admitted draw or tree counts for every
contributing series.

An excluded chain remains in the execution population with its exclusion
reason; it is not erased from the study design. A resumed segment remains part
of its parent chain unless the engine contract establishes independent
initialization. The pooled posterior denominator is therefore the admitted
sample population, while coverage and convergence statements still disclose
all requested chains and their terminal states.

## Evidence Package Acceptance

`adapter bayesian-evidence` requires explicit input, configuration, tree, log,
diagnostic, and report paths. The resulting package inventories and checksums
the selected files. Acceptance means the package is complete under that
inventory; it does not upgrade diagnostics into convergence or correspondence
into posterior parity.

## Deliberate Non-Claims

- Native `run_bayesian_inference(...)` has its own supported DNA contracts and
  does not dispatch through MrBayes or BEAST.
- Similar posterior means do not establish parity unless priors, likelihood,
  clocks, proposals, burn-in, and compared estimands are aligned.
- ESS and convergence thresholds qualify sampling behavior; they do not prove
  model adequacy or biological correctness.
- PCM4 and PCM5 retain an `MCMCglmm` external-posterior boundary that these
  adapters do not erase.
