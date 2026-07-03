---
title: Evidence Book
audience: public
type: evidence-navigation
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Evidence Book

This section documents the five lecture-derived comparative study families
tracked in `bijux-phylogenetics`.

It is the best place to understand how the repository moves from course
examples into checked-in evidence: what has already been reproduced, what has
been cataloged but not yet turned into governed parity, and which analytical
questions each study is meant to answer.

## Start Here

- [evidence book reading guide](reading-guide.md)
- repository evidence root:
  [`evidence-book/README.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/README.md)
- repository-wide review reports:
  [`evidence-book/index/`](https://github.com/bijux/bijux-phylogenetics/tree/main/evidence-book/index)

## Study Atlas

| Study | Main question | Current repository state |
| --- | --- | --- |
| [PCM1: primate longevity signal](pcm1-primate-longevity-signal/index.md) | How do raw primate traits and trees become a usable comparative dataset, signal estimate, and ancestral reconstruction? | Governed parity deployed |
| [PCM2: primate PGLS and signal](pcm2-primate-pgls-and-signal/index.md) | How do transformed trees, comparative model fitting, PGLS, and diagnostics behave on the same primate system? | Governed parity deployed |
| [PCM3: primate continuous and categorical regression](pcm3-primate-continuous-and-categorical-regression/index.md) | How should multivariable regression, categorical predictors, interactions, and phylogenetic correction be interpreted together? | Study dossier deployed; lecture parity not yet packaged as governed bundles |
| [PCM4: cichlid Bayesian linear and phylogenetic mixed models](pcm4-cichlid-bayesian-linear-and-phylogenetic-mixed-models/index.md) | How does the workflow move from Gaussian regression to Bayesian mixed models with a phylogenetic random effect? | Study dossier deployed; core posterior workflow still external |
| [PCM5: cichlid phylogenetic generalized linear mixed models](pcm5-cichlid-phylogenetic-generalized-linear-mixed-models/index.md) | How are count models, binomial-style models, posterior signal estimates, and repeated measures handled together? | Study dossier deployed; core posterior workflow still external |

## Studies

Each study page does three jobs:

- explain the analytical story in plain language
- point to the checked-in evidence dossier under `evidence-book/studies/`
- state clearly whether the repository already reproduces the lecture or is
  still documenting an external boundary

## Practical Boundary

Use these study pages when the question is:

- what a given lecture family is trying to show
- which repository surfaces have already been checked against that lecture
- where a study remains documented but not yet fully reproduced

Use the main phylogenetics handbook when the question is:

- how to install or run the package
- which API or CLI surfaces are public
- how native and wrapper-backed runtime families are organized
