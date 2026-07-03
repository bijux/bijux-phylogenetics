---
title: Evidence Book
audience: public
type: evidence-navigation
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-03
---

# Evidence Book

The evidence-book is the repository trust surface for governed parity, study
evidence, freshness review, and explicit coverage boundaries.

It should be read as a bounded review system around the runtime, not as a
replacement for the runtime handbook.

## Start Here

- [evidence book reading guide](reading-guide.md)
- repository evidence root:
  [`evidence-book/README.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/README.md)
- claim re-audit:
  [`evidence-book/index/claim-reaudit.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/claim-reaudit.md)
- analytical surface coverage:
  [`evidence-book/index/analytical-surface-coverage.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/analytical-surface-coverage.md)
- completion gates:
  [`evidence-book/index/completion-gates.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/completion-gates.md)

## Live Status Reports

The evidence-book has volatile status surfaces that change as studies are
re-audited. Treat those checked-in reports as the current source of truth
rather than relying on hard-coded summary counts in this page.

- freshness report:
  [`evidence-book/index/freshness-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/freshness-report.md)
- integrity report:
  [`evidence-book/index/integrity-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/integrity-report.md)
- coverage gaps:
  [`evidence-book/index/coverage-gaps.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/coverage-gaps.md)
- maturity scorecard:
  [`evidence-book/index/evidence-maturity-scorecard.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/evidence-maturity-scorecard.md)
- closure criteria:
  [`evidence-book/index/closure-criteria.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/closure-criteria.md)

## Studies

The public study guides currently documented here are:

- [PCM1: primate longevity signal](pcm1-primate-longevity-signal/index.md)
- [PCM2: primate PGLS and signal](pcm2-primate-pgls-and-signal/index.md)

Those pages explain the analytical story of the studies. The repository-root
`evidence-book/` directory and linked index reports carry the broader governed
freshness, integrity, and gap surfaces.

## Practical Boundary

Use the evidence-book when the question is:

- which runtime claims are reviewed or still open
- where freshness or coverage debt remains
- what one study demonstrates and what it does not demonstrate

Use the main phylogenetics handbook when the question is:

- how to install or run the package
- which API or CLI surfaces are public
- how native and wrapper-backed runtime families are organized
