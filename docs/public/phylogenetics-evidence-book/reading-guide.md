---
title: Evidence Book Reading Guide
audience: public
type: evidence-navigation
status: canonical
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Evidence Book Reading Guide

Use this section when you want to understand the comparative case studies as
studies, not just as API surfaces or repository status reports.

## Choose A Starting Point

- If you want the introductory primate workflow, start with PCM1.
- If you want comparative model fitting and PGLS, continue to PCM2.
- If you want multivariable regression with interactions and categorical
  predictors, read PCM3.
- If you want Bayesian mixed-model lecture families, read PCM4 and PCM5
  together because they share the cichlid dataset and the same `MCMCglmm`
  dependency.

## Core Reading Order

For a serious first pass through the course-derived studies, read:

1. `overview`
2. `PCM1`
3. `PCM2`
4. `PCM3`
5. `PCM4`
6. `PCM5`

That order mirrors the course progression. It starts with tree and data
preparation, moves through comparative regression, and ends at the Bayesian
mixed-model lectures where the repository boundary is currently sharpest.

## When To Leave These Pages

- Leave for the main phylogenetics handbook when you need installation,
  commands, or public runtime entry points.
- Leave for the checked-in study dossier in `evidence-book/studies/` when you
  need dataset registries, source provenance, or the full evidence map.
- Leave for `evidence-book/index/` when you need repository-wide freshness,
  integrity, or coverage reports rather than one study family.

## What These Pages Deliberately Do

- They explain the scientific workflow in reader-facing language.
- They separate a documented source lecture from a reproduced repository
  workflow.
- They keep open boundaries visible, especially for PCM4 and PCM5, where the
  main posterior workflows still depend on external tooling.
