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

## Analytical Population Dependency Graph

```mermaid
flowchart LR
    workbook["evidence-002<br/>workbook load"]
    types["evidence-003<br/>type repair"]
    missing["evidence-004<br/>missingness"]
    duplicates["evidence-005<br/>species aggregation"]
    tree["evidence-006<br/>tree import"]
    diagnostics["evidence-007<br/>tree diagnostics"]
    alignment["evidence-008<br/>tree–trait alignment"]
    exports["evidence-009<br/>governed derivatives"]
    aggregate["evidence-001<br/>signal and ancestral result"]

    workbook --> types --> missing --> duplicates --> alignment
    tree --> diagnostics --> alignment
    alignment --> exports --> aggregate
```

The graph is the analysis population contract. A change in coercion,
missingness, duplicate aggregation, tree identity, pruning, or label mapping
changes the downstream population and reopens signal and ancestral evidence.
The aggregate result cannot be detached from the preparation chain merely
because its final taxa count remains 75.

## Governed Observations

The aggregate result for the prepared 75-taxon analysis records these values:

| Observation | Bijux result | Interpretation boundary |
| --- | ---: | --- |
| longevity range | 103.0–657.6 | values use the prepared source trait and its recorded units |
| mean longevity | 317.362133 | descriptive value for the reconciled taxa |
| estimated Pagel λ | 0.803 | fitted signal at the precision retained by the PCM1 aggregate record |
| λ=0 likelihood-ratio statistic | 41.030408 | comparison with the registered zero-signal model |
| λ=0 p-value | 1.4988 × 10⁻¹⁰ | model-based evidence against that registered null, not a causal result |
| named MRCA ancestral estimate | 520.078423 | conditional estimate for the five-taxon descendant signature retained in the result |
| branches with positive reconstructed change | 72 | count under the declared ancestral reconstruction |
| branches with reconstructed increase greater than 12 | 55 | thresholded model output, not 55 observed evolutionary events |

The source values live in
[`bijux_reference_results.json`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-longevity-signal/evidence-001/results/bijux_reference_results.json).
The aggregate scalar table contains 31 comparisons: 16 `matched` and 15
`matched_with_tolerance`. The focused preparation bundles separate workbook
loading, type repair, missingness, duplicate aggregation, tree import,
diagnostics, correspondence, and export so a passing headline cannot hide
which transformation produced the analytical population.

## Citation Packet For PCM1

A citation of the signal or ancestral result should retain the PCM1 claim and
bundle identifier, source and prepared input identities, 75-taxon analytical
population, tree and trait derivative checksums, model and software identity,
the exact result row, comparison rule and tolerance, bundle verdict, freshness,
and the plotting limitation. For a branch count, also retain the eligible
branch population and threshold definition; for an ancestral value, retain the
descendant-taxon node key.

Do not cite the 0.803 lambda estimate from this guide alone. The guide is the
interpretation route; the governed JSON observation and its manifests are the
auditable record.

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
