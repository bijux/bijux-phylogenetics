---
title: PCM1 Primate Longevity Signal
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# PCM1: Primate Longevity Signal

PCM1 establishes the data contract reused by the later primate studies. Raw
trait records and a primate phylogeny are normalized, reconciled by species,
and reduced to governed derivatives before phylogenetic signal or ancestral
values are interpreted.

```mermaid
flowchart LR
    raw["Workbook and<br/>primate tree"] --> repair["Type repair, missingness,<br/>duplicate aggregation"]
    repair --> reconcile["Tree–trait taxon<br/>reconciliation"]
    reconcile --> signal["Phylogenetic signal"]
    signal --> ancestral["Ancestral reconstruction"]
```

## Scientific Question

Do primate longevity traits retain phylogenetic structure after the raw table
and tree are made analytically compatible, and how are reconstructed ancestral
values tied to that same prepared state?

The preparation is scientifically material. Type coercion, missing-data
handling, duplicate-species aggregation, pruning, and label reconciliation can
change both the observations and the covariance structure used downstream.

## Study Inputs

| Input | Role |
| --- | --- |
| `primate_raw.xlsx`, `primate.RData`, `primate.csv` | source trait records and lecture workspace |
| `primatetree.nex`, `trimmed_primatetree.nex` | source and lecture-prepared phylogenies |
| `reference_primate.csv` | governed normalized trait derivative |
| `reference_trimmed_primatetree.nwk` | governed reconciled tree derivative |

Input identity and provenance are recorded in the
[dataset catalog](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-longevity-signal/datasets/README.md)
and [provenance record](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-longevity-signal/provenance/README.md).

## Evidence Decomposition

| Bundle | Evaluated surface |
| --- | --- |
| `evidence-001` | end-to-end study summary |
| `evidence-002` | workbook loading |
| `evidence-003` | type repair |
| `evidence-004` | missing-data accounting |
| `evidence-005` | duplicate-species aggregation |
| `evidence-006` | tree import |
| `evidence-007` | tree diagnostics |
| `evidence-008` | tree–data correspondence |
| `evidence-009` | processed exports |

The dossier reports `analytical_alignment_complete` with the plotting boundary
kept separate. Bundle deployment does not turn presentation parity into a
numerical claim.

## Interpreting The Analysis

Phylogenetic signal asks whether observed trait similarity is structured by
the tree under the chosen model; it is not a measure of causal inheritance.
Ancestral reconstruction estimates internal-node states conditional on the
prepared tips, tree, branch lengths, and model. Neither analysis repairs a bad
taxon mapping retroactively.

## Review Priorities

- Confirm that excluded and aggregated observations are enumerated.
- Confirm identical taxon sets and label normalization for traits and tree.
- Separate source-tree pruning from model-driven transformations.
- Inspect signal and ancestral outputs against the same governed derivatives.
- Treat figures as projections of recorded results, not independent evidence.

## Scope

PCM1 supports the registered preparation, correspondence, export, signal, and
ancestral workflow claims for this primate study family. It does not establish
general performance across datasets, validate every comparative method, or
close plotting equivalence beyond the declared boundary.

Continue with [data and tree preparation](data-and-tree-preparation.md),
[signal and ancestral results](signal-and-ancestral-results.md), or
[PCM2](../pcm2-primate-pgls-and-signal/index.md) for regression and evolutionary
mode comparisons on the primate data.
