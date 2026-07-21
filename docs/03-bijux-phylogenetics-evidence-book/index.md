---
title: Evidence Book
audience: public
type: evidence-navigation
status: canonical
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-21
---

# Evidence Book

The Evidence Book is the review layer between scientific source material and
public claims. It records what was run, with which inputs and versions, which
checks were applied, and whether each bounded claim passes, fails, remains
open, or is blocked.

```mermaid
flowchart LR
    source["Course or scientific<br/>source"] --> study["Study dossier"]
    study --> bundle["Evidence bundle"]
    bundle --> checks["Declared checks<br/>and tolerances"]
    checks --> verdict["Claim verdict"]
    verdict --> indexes["Freshness, integrity,<br/>coverage, parity indexes"]
    indexes --> docs["Bounded public claim"]
```

## Two Reading Layers

Use these public pages to understand study questions and limitations. Use the
checked-in [`evidence-book/`](https://github.com/bijux/bijux-phylogenetics/tree/main/evidence-book)
for primary bundle records, registries, sources, and repository-wide review
reports. The narrative explains; the checked-in record adjudicates.

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

## Review A Claim In This Order

1. Find the study and claim identifier.
2. Open the bundle manifest and input manifest.
3. Inspect source provenance and runtime or engine versions.
4. Read the declared checks and acceptance criteria.
5. Inspect machine-readable results before the reviewer summary.
6. Confirm the verdict in the current indexes and freshness report.
7. Carry exclusions and open boundaries into any downstream statement.

## Study Atlas

| Study family | Scientific focus | Governed bundles | Material boundary |
| --- | --- | ---: | --- |
| [Cichlid Bayesian Linear And Phylogenetic Mixed Models](./pcm4-cichlid-bayesian-linear-and-phylogenetic-mixed-models/index.md) | Bayesian linear and phylogenetic mixed models | 7 | posterior computation remains external |
| [Cichlid Phylogenetic Generalized Linear Mixed Models](./pcm5-cichlid-phylogenetic-generalized-linear-mixed-models/index.md) | count, binary, repeated-measure models | 9 | posterior computation remains external |
| [Primate Continuous And Categorical Regression](./pcm3-primate-continuous-and-categorical-regression/index.md) | predictors, interactions, comparative regression | 8 | lecture-specific verdicts remain open |
| [Primate Longevity Signal](./pcm1-primate-longevity-signal/index.md) | preparation, signal, ancestral reconstruction | 9 | plotting boundary remains explicit |
| [Primate PGLS And Signal](./pcm2-primate-pgls-and-signal/index.md) | GLS/PGLS, transformed trees, signal tests | 10 | one coverage boundary remains explicit |

The bundle count describes inventory, not number of passing claims. Consult
the live indexes for verdicts.

## Repository-Wide Control Reports

| Report | Question it answers |
| --- | --- |
| [evidence index](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/evidence-index.json) | which bundles and claims exist? |
| [integrity report](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/integrity-report.md) | are declared records structurally coherent? |
| [freshness report](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/freshness-report.md) | do records still match their material dependencies? |
| [parity dashboard](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/parity-dashboard.md) | which declared parity comparisons pass or remain open? |
| [coverage gaps](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/coverage-gaps.md) | where are claim or workflow proofs incomplete? |
| [false-confidence audit](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/false-confidence-audit.md) | where could presentation overstate evidence strength? |
| [scientific debt register](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/scientific-debt-register.md) | which unresolved boundaries require durable follow-through? |

## What A Bundle Proves

A bundle proves only that its declared checks and acceptance rules produce its
recorded verdict for identified inputs, sources, versions, and configuration.
It may support a larger claim only when the claim mapping says so. It does not
certify untested datasets, models, external-engine versions, or biological
interpretations.

## Open Evidence Is First-Class

PCM4 and PCM5 remain valuable precisely because their external posterior
boundary is visible. An open or blocked verdict identifies the missing work
more honestly than a polished narrative that implies completion. Use
[evidence honesty and limits](../01-bijux-phylogenetics/quality/evidence-honesty-and-limits.md)
when translating these records into public language.
