# Evidence Book

This directory is the checked-in study atlas for `bijux-phylogenetics`.
It records which comparative case studies are represented in the repository,
which evidence slices have already been deployed, and where the remaining
runtime boundaries still sit.

Use it for study-grounded review. Use `artifacts/` for scratch outputs and
local runs that should not become repository trust surfaces on their own.

## Read This Directory In Two Layers

- Reader-facing study overviews live under
  [`../docs/public/phylogenetics-evidence-book/`](../docs/public/phylogenetics-evidence-book/)
- Checked-in study dossiers, datasets, and provenance live under
  [`studies/`](./studies/)

The public docs explain the analytical story. The `evidence-book/studies/`
directories hold the durable repository record for each study family.

## Study Atlas

| Study | Lecture family | Repository state |
| --- | --- | --- |
| [`primate-longevity-signal`](./studies/primate-longevity-signal/) | PCM1 plots, signal, and ancestral reconstruction | Governed parity bundles deployed |
| [`primate-pgls-and-signal`](./studies/primate-pgls-and-signal/) | PCM2 evolutionary modes, PGLS, and signal testing | Governed parity bundles deployed |
| [`primate-continuous-and-categorical-regression`](./studies/primate-continuous-and-categorical-regression/) | PCM3 multivariable regression with categorical predictors and interactions | Study dossier deployed, runtime parity surface partially available |
| [`cichlid-bayesian-linear-and-phylogenetic-mixed-models`](./studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/) | PCM4 Bayesian linear models and phylogenetic mixed models | Study dossier deployed, runtime parity boundary still open |
| [`cichlid-phylogenetic-generalized-linear-mixed-models`](./studies/cichlid-phylogenetic-generalized-linear-mixed-models/) | PCM5 phylogenetic generalized linear mixed models | Study dossier deployed, runtime parity boundary still open |

## Repository Reports

These index pages remain the repository-wide review surfaces:

- Teaching and migration landing page:
  [`index/teaching-and-migration.md`](./index/teaching-and-migration.md)
- Freshness tracking:
  [`index/freshness-report.md`](./index/freshness-report.md)
- Integrity tracking:
  [`index/integrity-report.md`](./index/integrity-report.md)
- Coverage gaps:
  [`index/coverage-gaps.md`](./index/coverage-gaps.md)
- Claim re-audit:
  [`index/claim-reaudit.md`](./index/claim-reaudit.md)
- Analytical surface coverage:
  [`index/analytical-surface-coverage.md`](./index/analytical-surface-coverage.md)
- Closure criteria:
  [`index/closure-criteria.md`](./index/closure-criteria.md)
- Maturity scorecard:
  [`index/evidence-maturity-scorecard.md`](./index/evidence-maturity-scorecard.md)
- Completion gates:
  [`index/completion-gates.md`](./index/completion-gates.md)

## Reading Order

Start with the study README that matches the lecture family you care about.
Then follow its dataset registry, source provenance, and evidence map.

When you need repository-wide status rather than one study story, move from
the study dossier into the index reports above.
