---
title: Evidence Book
audience: public
type: evidence-navigation
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-04
---

# Evidence Book

This section explains the comparative study record behind the public
phylogenetics runtime.

The five study families come from course material that the repository has
turned into checked-in, reviewable evidence. All five now have governed
parity bundle IDs. Some studies already carry stronger numerical
parity records, while others keep parity verdicts open because the
scientific workflow still depends on an external engine or an unclosed
lecture-specific boundary.

```mermaid
flowchart LR
    A[Source lecture] --> B[Public study guide]
    B --> C[Checked-in evidence dossier]
    C --> D[Freshness and integrity reports]
    C --> E[Study-specific evidence bundles]
```

## What This Guide Is For

- to explain what each PCM study is trying to show in reader-facing language
- to point directly to the checked-in dossier, datasets, and governed evidence
  bundles behind that study
- to keep open boundaries visible instead of flattening every study into the
  same level of trust
- to show that even this first evidence cluster already covers preparation,
  signal, comparative regression, and bounded mixed-model study families

## Start Here

- checked-in study dossiers:
  [`evidence-book/README.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/README.md)
- repository-wide review reports:
  [`evidence-book/index/`](https://github.com/bijux/bijux-phylogenetics/tree/main/evidence-book/index)

Treat those checked-in reports as the current source of truth when you need
live repository status rather than the stable narrative on this page.

## At A Glance

- studies: `5`
- evidence bundles: `43`
- teaching studies: `2`
- migration studies: `2`
- freshness statuses: `current=43`
- coverage gaps still tracked: `28`
- foundational numerical trust: `bounded`
- reviewer readiness: `bounded`
- maturity tier: `reviewable_but_incomplete`

## What The Checked-In Record Already Contains

- five study families with local public study guides
- governed evidence IDs across PCM1, PCM2, PCM3, PCM4, and PCM5
- checked-in datasets, dossier reports, freshness reports, and integrity
  reports
- explicit coverage-gap tracking where the repository has not yet closed a
  lecture or workflow boundary
- reviewer-facing documentation that distinguishes strong numerical parity from
  still-open scientific scope

## Why This Evidence Group Matters

This evidence-book is not trying to prove every scientific surface in the
repository at once. It does already prove that the repository can preserve one
serious comparative-study family with:

- durable source provenance
- named datasets and derivative inputs
- decomposed evidence IDs instead of one vague study label
- bounded trust language where the lecture workflow still exceeds repository
  closure

## Core Reading Order

1. `overview`
2. `PCM1: primate longevity signal`
3. `PCM2: primate PGLS and signal`
4. `PCM3: primate continuous and categorical regression`
5. `PCM4: cichlid Bayesian linear and phylogenetic mixed models`
6. `PCM5: cichlid phylogenetic generalized linear mixed models`

That sequence mirrors the course progression. It starts with data and tree
preparation, moves through comparative regression, and ends with the
cichlid mixed-model studies where external posterior tooling remains the
largest open boundary.

## Study Atlas

| Study | Main question | Evidence IDs | Current repository state |
| --- | --- | --- | --- |
| [Primate Longevity Signal](./pcm1-primate-longevity-signal/index.md) | How do primate traits and tree preparation become a usable signal and ancestral-state workflow? | `evidence-001` to `evidence-009` (9 bundles) | Governed parity bundles deployed; plotting boundary remains explicit |
| [Primate PGLS And Signal](./pcm2-primate-pgls-and-signal/index.md) | How do baseline GLS, PGLS, transformed trees, and signal tests behave on the same primate dataset? | `evidence-001` to `evidence-010` (10 bundles) | Governed parity bundles deployed; one coverage boundary remains explicit |
| [Primate Continuous And Categorical Regression](./pcm3-primate-continuous-and-categorical-regression/index.md) | How should continuous predictors, categorical predictors, and interactions be interpreted together in comparative regression? | `evidence-001` to `evidence-008` (8 bundles) | Governed parity bundles deployed; lecture-specific verdicts remain open |
| [Cichlid Bayesian Linear And Phylogenetic Mixed Models](./pcm4-cichlid-bayesian-linear-and-phylogenetic-mixed-models/index.md) | Where does the cichlid workflow move from standard regression into Bayesian mixed models with a phylogenetic effect? | `evidence-001` to `evidence-007` (7 bundles) | Governed parity bundles deployed; posterior workflow still external |
| [Cichlid Phylogenetic Generalized Linear Mixed Models](./pcm5-cichlid-phylogenetic-generalized-linear-mixed-models/index.md) | How do count models, binary-style models, repeated measures, and posterior signal estimates fit together in one study family? | `evidence-001` to `evidence-009` (9 bundles) | Governed parity bundles deployed; posterior workflow still external |

## How To Read A Study Dossier

1. Start with PCM1 if you need the full preparation path from raw primate
   traits and tree files to signal estimation and ancestral reconstruction.
2. Continue to PCM2 and PCM3 when the question shifts from preparation into
   comparative regression and model interpretation.
3. Read PCM4 and PCM5 together when the question is the cichlid mixed-model
   sequence, because they share the same dataset family and the same external
   posterior boundary.
4. Follow the linked dossier and repository-wide reports whenever you need the
   governed bundle record rather than the narrative summary.

## When To Use These Pages

- Use the public study guides when you want the scientific story in plain language.
- Use the checked-in evidence-book when you need provenance, registries, or bundle-level proof.
- Use the index reports when you need freshness, integrity, or repository-wide gap tracking.

## Reading Boundaries

- These pages do not claim that every lecture has already been reproduced end to end inside the repository.
- They do not hide external-tool dependencies just because nearby Python surfaces exist.
- They do keep those limits visible, especially for PCM4 and PCM5.

## Review Reports

- teaching and migration index: [`evidence-book/index/teaching-and-migration.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/teaching-and-migration.md)
- parity dashboard: [`evidence-book/index/parity-dashboard.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/parity-dashboard.md)
- freshness report: [`evidence-book/index/freshness-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/freshness-report.md)
- integrity report: [`evidence-book/index/integrity-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/integrity-report.md)
- coverage gaps: [`evidence-book/index/coverage-gaps.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/coverage-gaps.md)
- claim re-audit: [`evidence-book/index/claim-reaudit.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/claim-reaudit.md)
- closure criteria: [`evidence-book/index/closure-criteria.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/closure-criteria.md)

The study record is useful because it keeps strong results visible without
pretending that every workflow is fully reproduced inside the repository.
