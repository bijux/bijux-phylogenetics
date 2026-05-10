# Primate PGLS and signal evidence study Migration Guide

- study id: `primate-pgls-and-signal`
- course source: `external:lund/pcm2-modes-pgls/script`

## Comparison Modes

- `direct_parity`: `6`

## Side-By-Side Examples

### Non-phylogenetic GLS fit

- example id: `pcm2-baseline-gls`
- family id: `baseline-regression`
- fragment id: `baseline-gls-fit`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm2-modes-pgls/script#L122-L136`
- Bijux locators: `bijux_phylogenetics.comparative.pgls:run_pgls, bijux_phylogenetics.evidence.studies.primate_pgls_and_signal:build_primate_pgls_signal_bundles`
- Bijux summary: The baseline regression surface stays visible as a governed comparison instead of disappearing once lambda-based models start.
- Why migrate: R users can see the exact non-phylogenetic baseline they already know, then compare how Bijux layers explicit evidence on top.

### Fixed-lambda and estimated-lambda PGLS fits

- example id: `pcm2-pgls`
- family id: `phylogenetic-regression`
- fragment id: `pagel-lambda-regression`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm2-modes-pgls/script#L138-L179`
- Bijux locators: `bijux_phylogenetics.comparative.pgls:run_pgls, bijux_phylogenetics.comparative.signal:estimate_pagels_lambda`
- Bijux summary: Bijux records the regression and lambda-fitting outputs with explicit tolerance rules and verdicts.
- Why migrate: The migration value is that regression parity and tolerance logic become explicit review artifacts instead of silent assumptions.

### Intercept-only PGLS and lambda-zero likelihood-ratio testing

- example id: `pcm2-signal-test`
- family id: `phylogenetic-signal`
- fragment id: `phylogenetic-signal-test`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm2-modes-pgls/script#L181-L192`
- Bijux locators: `bijux_phylogenetics.comparative.signal:compute_phylogenetic_signal_test, bijux_phylogenetics.comparative.signal:estimate_pagels_lambda`
- Bijux summary: The intercept-only signal workflow is carried into a governed Python surface with explicit verdicts and tolerances.
- Why migrate: Students can compare the same signal decision path while seeing which numerical differences remain bounded and reviewable.

### Estimated-lambda residual and fitted diagnostics

- example id: `pcm2-diagnostics`
- family id: `diagnostics`
- fragment id: `estimated-lambda-diagnostics`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm2-modes-pgls/script#L168-L179`
- Bijux locators: `bijux_phylogenetics.evidence.studies.primate_pgls_and_signal:build_primate_pgls_signal_scalar_parity_table, bijux_phylogenetics.comparative.pgls:run_pgls`
- Bijux summary: Bijux exposes diagnostics as recorded scalar checks rather than requiring manual plot comparison to understand parity.
- Why migrate: The migration benefit is reviewer-visible diagnostics that can be versioned and compared, not just re-plotted.
