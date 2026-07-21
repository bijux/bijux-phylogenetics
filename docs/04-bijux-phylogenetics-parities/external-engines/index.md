---
title: External Engine Correspondence
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# External Engine Correspondence

MAFFT, trimAl, IQ-TREE2, FastTree, MrBayes, and BEAST remain external
computational owners. Bijux adapters make their runs inspectable by validating
inputs, recording executable identity and configuration, capturing failures,
checking required outputs, parsing results, and writing normalized artifacts.

```mermaid
flowchart LR
    input["Validated Bijux input"]
    prepare["Engine-specific<br/>configuration"]
    engine["External executable"]
    capture["Native output and<br/>failure capture"]
    parse["Typed Bijux result"]
    bundle["Manifest and<br/>review artifacts"]

    input --> prepare --> engine --> capture --> parse --> bundle
```

| Stage | External owner | Correspondence surface |
| --- | --- | --- |
| [multiple-sequence alignment](alignment-and-trimming.md) | MAFFT | input preparation, command capture, valid non-empty alignment |
| [alignment trimming](alignment-and-trimming.md) | trimAl | configured trimming, output validation, before/after diagnostics |
| [likelihood inference and support](tree-inference.md) | IQ-TREE2 | model, topology, likelihood, support and run-artifact parsing |
| [approximate large-tree inference](tree-inference.md) | FastTree | topology and SH-like local-support parsing with explicit scale |
| [Bayesian tree inference](bayesian-inference.md) | MrBayes | NEXUS execution, trace/tree parsing, burn-in and consensus review |
| [Bayesian time-tree inference](bayesian-inference.md) | BEAST | XML execution, log/tree parsing, calibration and posterior review |

Cross-engine comparison is valid only for fields with aligned meanings.
IQ-TREE UFBoot, FastTree SH-like local support, and posterior clade probability
can share a normalized display scale while remaining distinct measures.
