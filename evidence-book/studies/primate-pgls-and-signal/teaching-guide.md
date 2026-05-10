# Primate PGLS and signal evidence study Teaching Guide

- study id: `primate-pgls-and-signal`
- categories: `teaching-study, migration-study`
- course source: `external:lund/pcm2-modes-pgls/script`
- concept tags: `ancestral, ancestral-gaps, ancestral-reconstruction, baseline-model, continuous, continuous-model-fitting, course-context, coverage-boundary, diagnostics, eb, estimated-lambda, fitting, fixed-lambda, gls, heteroscedasticity, lambda-zero, likelihood, likelihood-ratio, likelihood-ratio-tests, model, ou, pagel-lambda, pgls, phylogenetic-signal, qq, ratio, reconstruction, regression, reproducibility, residuals, tests, transformed, transformed-tree-workflows, tree, workflows, workspace-reload`

## Workflow contracts

- family id: `workflow-contracts`
- verdict: `matched`
- coverage: `covered`
- concept tags: `workspace-reload, course-context, reproducibility`
- evidence ids: `evidence-001`
- fragment ids: `workspace-reload-contract`
- teaching narrative: The lecture's one-line reload assumption is made explicit so students know which objects and paths the rest of the workflow depends on.

## Transformed tree workflows

- family id: `transformed-tree-workflows`
- verdict: `matched`
- coverage: `covered`
- concept tags: `transformed, transformed-tree-workflows, tree, workflows`
- evidence ids: `evidence-006`
- fragment ids: `transformed-tree-workflows`
- teaching narrative: This evidence family is indexed for teaching review, but its dedicated teaching narrative has not been curated yet. Review the bundle-level reports and source fragments directly until the family guide is strengthened.

## Continuous model fitting

- family id: `continuous-model-fitting`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `continuous, continuous-model-fitting, fitting, model`
- evidence ids: `evidence-007`
- fragment ids: `continuous-model-comparison`
- teaching narrative: This evidence family is indexed for teaching review, but its dedicated teaching narrative has not been curated yet. Review the bundle-level reports and source fragments directly until the family guide is strengthened.

## Likelihood-ratio tests

- family id: `likelihood-ratio-tests`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `likelihood, likelihood-ratio-tests, ratio, tests`
- evidence ids: `evidence-008`
- fragment ids: `evolutionary-mode-likelihood-ratios`
- teaching narrative: This evidence family is indexed for teaching review, but its dedicated teaching narrative has not been curated yet. Review the bundle-level reports and source fragments directly until the family guide is strengthened.

## Ancestral reconstruction

- family id: `ancestral-reconstruction`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `ancestral, ancestral-reconstruction, reconstruction`
- evidence ids: `evidence-009`
- fragment ids: `ancestral-mode-comparison`
- teaching narrative: This evidence family is indexed for teaching review, but its dedicated teaching narrative has not been curated yet. Review the bundle-level reports and source fragments directly until the family guide is strengthened.

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
- evidence ids: `evidence-010`
- fragment ids: `mode-linked-intercept-models`
- teaching narrative: Open EB and ancestral-mode gaps stay explicit so students are not taught a false sense of completeness.
