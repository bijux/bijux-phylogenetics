---
title: Simulation And Datasets
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Simulation And Datasets

Simulation supplies known generative conditions. Dataset bundles supply stable,
inspectable inputs with biological or methodological context. Together they
support examples, regression tests, recovery experiments, sensitivity studies,
and benchmarks without pretending that one corpus can answer every question.

```mermaid
flowchart LR
    question["Scientific or<br/>validation question"]
    generator["Tree · alignment<br/>trait generator"]
    truth["Parameters · histories<br/>true topology"]
    dataset["Files · metadata<br/>taxon identities"]
    workflow["Analysis workflow"]
    comparison["Recovery · robustness<br/>reference comparison"]

    question --> generator
    generator --> truth
    generator --> dataset
    dataset --> workflow
    truth --> comparison
    workflow --> comparison
```

## Native Simulation Families

| Family | Controlled conditions | Outputs to retain |
| --- | --- | --- |
| trees | birth–death-like tree generation and multispecies-coalescent scenarios | topology, times or lengths, lineage/species mapping, parameters, seed |
| alignments | characters evolved on a supplied tree under a declared state model | alignment, tree identity, model parameters, site count, seed |
| continuous traits | Brownian, correlated, and related continuous processes | tip values, covariance or rate parameters, starting state, seed |
| discrete traits | finite-state trajectories and stochastic histories | tip states, branch histories, rate matrix, root state, seed |

A seed is necessary for replay but is not a full simulation contract. Generator
version, random-number algorithm, tree and model conventions, parameter values,
and output schema must also be recorded. The realized dataset should be retained
when exact byte-for-byte replay matters.

## Dataset Roles

Choose a dataset by the question it can adjudicate:

| Role | What it is good for | What it cannot establish alone |
| --- | --- | --- |
| known-answer fixture | exact parsing, scoring, reconstruction, and small numerical contracts | realistic performance or broad biological behavior |
| recovery suite | whether an estimator recovers declared generative modes or parameters | behavior outside the simulated model family |
| reference corpus | comparison with a named implementation or published result | universal equivalence across versions and settings |
| stress corpus | missingness, identifier errors, scale, degeneracy, and failure behavior | ordinary-case scientific accuracy |
| biological study bundle | an interpretable workflow on a named organism and question | known truth or general performance |

The package catalog spans known-answer references, continuous and discrete mode
recovery, macroevolution recovery, data-quality stress, primate comparative
traits, avian reproductive traits, seashore flora, pathogen and rabies studies,
Influenza A HA, Pleistocene bear cytochrome-b fragments, and a gnathostome
ortholog protein benchmark. These are distinct assets with distinct licenses,
provenance, schemas, and permitted claims—not a generic pool of example files.

## A Dataset Bundle Is More Than Data Files

A durable bundle identifies:

- stable dataset and revision identifiers;
- upstream source, citation, license, and retrieval or derivation history;
- checksums for governed input files;
- taxon naming and reconciliation policy;
- columns, units, state coding, missing-value semantics, and transformations;
- tree rootedness, branch-length meaning, partitions, and alignments;
- known truth or expected properties, when those genuinely exist;
- exclusions, quality warnings, and intended analytical roles.

Exports, panels, and demonstration workflows are projections of this contract.
They should not become alternative sources of truth for the same dataset.

## Designing Recovery Experiments

Recovery means comparing an estimate with a condition known because it was
generated or independently fixed. Define the estimand and acceptance rule
before running the analysis. Use replicated seeds, cover interior and boundary
parameters, and separate estimator bias from optimizer or search failure.

For topology, compare clades or distances as well as an exact-match indicator.
For continuous parameters, retain signed error, absolute or relative error, and
interval coverage where applicable. For discrete histories, distinguish tip
state recovery, internal-state uncertainty, and branch-transition recovery.

## Preventing Simulation Circularity

Agreement is weakest when generator and estimator share the same defect or
hidden assumption. Stronger review combines several sources:

1. hand-checkable fixtures for exact contracts;
2. simulation recovery under multiple parameter regimes;
3. independent reference implementations where assumptions align;
4. biological datasets for realism and failure discovery;
5. stress corpora for explicit limits.

No simulated success licenses a broader claim about empirical adequacy.
