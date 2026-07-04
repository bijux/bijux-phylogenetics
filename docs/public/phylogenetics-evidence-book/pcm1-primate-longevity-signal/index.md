---
title: PCM1 Primate Longevity Signal
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# PCM1: Primate Longevity Signal

PCM1 is the opening comparative case study in the evidence book. It follows
the primate dataset from raw workbook and tree inputs through data repair,
species alignment, phylogenetic signal, and ancestral reconstruction.

## Why PCM1 Matters

Later studies only make sense if these earlier steps are already clear. PCM1
shows where the primate data came from, how the tree was trimmed, and how the
first comparative claims were checked.

## What You Will Learn Here

- how the repository prepares one lecture dataset into a durable study input
- how phylogenetic signal is assessed on that prepared primate dataset
- how ancestral-state outputs are kept tied to the same checked inputs

## Read PCM1 In Order

1. [data and tree preparation](data-and-tree-preparation.md)
2. [signal and ancestral results](signal-and-ancestral-results.md)

## Checked-In Study Dossier

The full study dossier lives in the repository under
[`evidence-book/studies/primate-longevity-signal/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-longevity-signal/README.md).
That dossier includes the shared datasets, source provenance, and the smaller
evidence slices that support this public guide.

## Checked-In Inputs

- lecture source copies:
  `primate_raw.xlsx`, `primate.RData`, `primate.csv`, `primatetree.nex`, and
  `trimmed_primatetree.nex`
- governed derivatives:
  `reference_primate.csv` and `reference_trimmed_primatetree.nwk`
- study-local catalogs:
  [dataset catalog](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-longevity-signal/datasets/README.md)
  and
  [provenance guide](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-longevity-signal/provenance/README.md)

## Evidence Bundle Atlas

- `evidence-001`: end-to-end governed parity summary
- `evidence-002` to `evidence-005`: workbook repair and processed-table trust
- `evidence-006` to `evidence-008`: tree import, diagnostics, and tree-data
  correspondence
- `evidence-009`: processed export contract for downstream reuse

## What PCM1 Does Not Prove

- it does not close every comparative method in the repository
- it does not turn one primate teaching family into one general claim about all
  datasets
- it does not replace the broader evidence-book freshness, integrity, or
  coverage reports

## Current Repository State

- Governed parity bundles are deployed for `evidence-001` through
  `evidence-009`.
- The checked-in evidence record includes the summary bundle plus the
  component-level parity bundles for preprocessing, tree handling, signal, and
  ancestral reporting.
- The repository does not claim plot-for-plot identity with the R plotting
  stack.
- PCM1 remains an introductory study; it does not stand in for the broader
  regression and mixed-model workflows covered later.

## Read This Study Next

- [PCM2: primate PGLS and signal](../pcm2-primate-pgls-and-signal/index.md)
