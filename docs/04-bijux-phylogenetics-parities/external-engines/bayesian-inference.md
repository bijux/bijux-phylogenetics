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

## Deliberate Non-Claims

- Native `run_bayesian_inference(...)` has its own supported DNA contracts and
  does not dispatch through MrBayes or BEAST.
- Similar posterior means do not establish parity unless priors, likelihood,
  clocks, proposals, burn-in, and compared estimands are aligned.
- ESS and convergence thresholds qualify sampling behavior; they do not prove
  model adequacy or biological correctness.
- PCM4 and PCM5 retain an `MCMCglmm` external posterior boundary that these
  adapters do not erase.
