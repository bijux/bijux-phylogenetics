---
title: MrBayes And BEAST Correspondence
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# MrBayes And BEAST Correspondence

MrBayes and BEAST own their posterior computations. Bijux provides explicit
preparation, execution, parsing, diagnostics, tree-set review, and reporting
surfaces. This is external-engine correspondence; it is not a hidden route
through which external samples are presented as native Bayesian inference.

## MrBayes Surface

- prepare a nucleotide or protein NEXUS analysis with a supported model and
  rate configuration;
- execute the prepared analysis and retain executable identity, command,
  standard streams, exit state, and native files;
- parse parameter traces, MCMC diagnostics, posterior trees, sampled
  generations, and consensus posterior probabilities;
- apply burn-in and subsampling while retaining the selection policy;
- compute ESS, trace drift, convergence, parameter summaries, burn-in
  sensitivity, and clade sensitivity;
- render posterior reports and methods text linked to the parsed artifacts.

## BEAST Surface

- prepare and inspect time-tree XML with sequence, clock, tree-prior,
  calibration, and tip-date inputs;
- validate fossil calibrations and tip dates against tree and alignment taxa;
- execute a prepared XML analysis with captured engine provenance;
- parse state-indexed log rows and posterior tree samples;
- summarize parameters, ESS, trace drift, burn-in sensitivity, topology
  diversity, consensus, calibration assumptions, and time-tree readiness;
- bundle configs, logs, trees, diagnostics, tables, and reports for review.

## Shared Review Contract

```mermaid
flowchart LR
    config["NEXUS or XML<br/>plus governed inputs"]
    engine["MrBayes or BEAST"]
    raw["Native traces,<br/>trees and logs"]
    policy["Burn-in · thinning<br/>chain and state identity"]
    diagnostics["ESS · drift · convergence<br/>topology and clade sensitivity"]
    report["Evidence package<br/>and methods text"]

    config --> engine --> raw --> policy --> diagnostics --> report
```

The report must distinguish a prepared-only analysis, a parsed prior run, and
an engine execution performed by the current workflow. File presence alone is
not execution proof.

## Run Identity And Resume Safety

| Identity component | MrBayes | BEAST |
| --- | --- | --- |
| prepared analysis | NEXUS content, partitions, model/rates, generations, chains, sample frequency, burn-in policy | XML content, taxa, substitution model, clock, tree prior, calibrations/tip dates, chain length, logging policy |
| execution | `mb` path/version, command, captured streams, exit and timeout state | `beast` path/version, command, threads, seed, overwrite policy, captured streams, exit and timeout state |
| posterior output | parameter trace, tree samples, generation identity, consensus and native files | state-indexed log, tree samples, state identity, consensus/MCC and native files |
| selection | burn-in and subsampling against retained generations | burn-in and subsampling against retained states |
| diagnostics | parameter ESS, trace drift, chain and clade behavior | parameter ESS, trace drift, topology/calibration/time-tree behavior |

`mrbayes-run --resume` and `beast-run --resume` reuse only a verified completed
manifest whose executable, inputs, configuration, and recorded output hashes
still agree. Timeout markers, killed processes, malformed posterior files, or
changed identity require an explicit clean restart or refusal under the chosen
incomplete-run policy.

Resume therefore preserves execution lineage; it does not create a new
independent chain and must not increase the chain denominator in convergence
claims.

## Use The Adapters

```bash
bijux-phylogenetics adapter mrbayes-prepare \
  dataset/alignment.fasta \
  --out artifacts/bayesian/analysis.nex

bijux-phylogenetics adapter beast-xml \
  artifacts/bayesian/analysis.xml \
  --json

bijux-phylogenetics adapter bayesian-evidence \
  --help
```

Execution commands (`mrbayes-run` and `beast-run`) require the external
executables. Parser and report commands can review governed existing outputs
without claiming that the current machine reran the posterior analysis.

## Diagnostic Interpretation

| Observation | Supports | Does not establish |
| --- | --- | --- |
| adequate ESS under a named estimator | retained draws contain enough effective information for that parameter threshold | stationarity, model adequacy, or agreement between independent chains |
| stable trace after selected burn-in | the named drift or window rule did not detect instability | correct burn-in for every parameter or topology |
| similar independent-chain summaries | selected estimands agree under the comparison rule | convergence of unobserved dimensions or correctness of the posterior target |
| stable clade frequencies | topology summaries are insensitive under the named chain/burn-in comparison | identical branch lengths, dates, parameters, or biological history |
| parser acceptance | native files satisfy the supported structural contract | that the current workflow executed the engine or that sampling is adequate |

A posterior comparison must align priors, likelihood, clock and tree model,
parameterization, state/chain selection, transformations, and estimands before
values are adjudicated. Otherwise the two posteriors remain separate governed
observations.

## Deliberate Non-Claims

- Native `run_bayesian_inference(...)` has its own supported DNA contracts and
  does not dispatch through MrBayes or BEAST.
- Similar posterior means do not establish parity unless priors, likelihood,
  clocks, proposals, burn-in, and compared estimands are aligned.
- ESS and convergence thresholds qualify sampling behavior; they do not prove
  model adequacy or biological correctness.
- PCM4 and PCM5 retain an `MCMCglmm` external posterior boundary that these
  adapters do not erase.
