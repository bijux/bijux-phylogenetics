---
title: PCM2 Primate PGLS and Signal
audience: public
type: study-guide
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-04
---

# PCM2: Primate PGLS And Signal

PCM2 is the deeper primate comparative study. It covers transformed trees,
comparative model fitting, likelihood-ratio comparisons, baseline GLS,
Pagel-lambda regression, and signal diagnostics on the same lecture family.

## Why PCM2 Matters

PCM1 establishes the prepared study inputs. PCM2 asks whether more demanding
comparative methods still behave coherently once phylogenetic correction and
model comparison enter the picture.

## Read PCM2 In Order

1. [Model comparisons and PGLS](model-comparisons-and-pgls.md)
2. [Signal diagnostics and open boundaries](signal-diagnostics-and-open-boundaries.md)

## Checked-In Study Dossier

The full checked-in record lives under
[`evidence-book/studies/primate-pgls-and-signal/`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-pgls-and-signal/README.md).

## Inputs And Evidence Map

- lecture source copies:
  `primate_raw.xlsx`, `primate.RData`, and `primatetree.nex`
- governed derivatives:
  `reference_primate.csv` and `reference_trimmed_primatetree.nwk`
- catalogs:
  [dataset catalog](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-pgls-and-signal/datasets/README.md)
  and
  [provenance guide](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/primate-pgls-and-signal/provenance/README.md)
- evidence bundles:
  `evidence-001` through `evidence-010`

## What To Take From PCM2

- baseline GLS and phylogenetically corrected PGLS can be read as one
  comparative sequence
- transformed trees and model comparisons belong in the same analytical story,
  not in separate islands
- open boundaries can improve credibility when they remain visible

## What PCM2 Does Not Prove

- PCM2 does not close every comparative model in the repository.
- It does not erase the remaining open boundaries named inside the study.
- It does not substitute for release review or broader evidence-book coverage
  reports.

## Read Next

- [PCM3: primate continuous and categorical regression](../pcm3-primate-continuous-and-categorical-regression/index.md)
