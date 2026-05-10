# Primate PGLS and signal evidence study Teaching Guide

- study id: `primate-pgls-and-signal`
- categories: `teaching-study, migration-study`
- course source: `external:lund/pcm2-modes-pgls/script`
- concept tags: `ancestral-gaps, baseline-model, course-context, coverage-boundary, diagnostics, eb, estimated-lambda, fixed-lambda, gls, heteroscedasticity, lambda-zero, likelihood-ratio, ou, pagel-lambda, pgls, phylogenetic-signal, qq, regression, reproducibility, residuals, workspace-reload`

## Workflow contracts

- family id: `workflow-contracts`
- verdict: `matched`
- coverage: `covered`
- concept tags: `workspace-reload, course-context, reproducibility`
- evidence ids: `evidence-001`
- fragment ids: `workspace-reload-contract`
- teaching narrative: The lecture's one-line reload assumption is made explicit so students know which objects and paths the rest of the workflow depends on.

## Baseline regression

- family id: `baseline-regression`
- verdict: `matched`
- coverage: `covered`
- concept tags: `gls, baseline-model, regression`
- evidence ids: `evidence-002`
- fragment ids: `baseline-gls-fit`
- teaching narrative: Baseline GLS is isolated before phylogenetic covariance is introduced, which makes the teaching sequence and the trust sequence line up.

## Phylogenetic regression

- family id: `phylogenetic-regression`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `pgls, pagel-lambda, fixed-lambda, estimated-lambda`
- evidence ids: `evidence-003`
- fragment ids: `pagel-lambda-regression`
- teaching narrative: Fixed-lambda and estimated-lambda regression are kept in one evidence family so the lecture transition from GLS to PGLS remains reviewable.

## Phylogenetic signal

- family id: `phylogenetic-signal`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `phylogenetic-signal, lambda-zero, likelihood-ratio`
- evidence ids: `evidence-004`
- fragment ids: `phylogenetic-signal-test`
- teaching narrative: Signal testing is treated as its own teaching family so the lecture's intercept-only workflow does not get mixed into broader model-comparison claims.

## Diagnostics

- family id: `diagnostics`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `residuals, qq, heteroscedasticity, diagnostics`
- evidence ids: `evidence-005`
- fragment ids: `baseline-gls-diagnostics, estimated-lambda-diagnostics`
- teaching narrative: Residual and QQ diagnostics are converted into machine-recorded summary values so the teaching surface is auditable even without figure matching.

## Coverage boundaries

- family id: `coverage-boundaries`
- verdict: `not_comparable`
- coverage: `coverage-gap`
- concept tags: `eb, ou, ancestral-gaps, coverage-boundary`
- evidence ids: `evidence-006`
- fragment ids: `transformed-tree-workflows, continuous-model-comparison, ancestral-mode-comparison, mode-linked-intercept-models`
- teaching narrative: Open EB and ancestral-mode gaps stay explicit so students are not taught a false sense of completeness.
