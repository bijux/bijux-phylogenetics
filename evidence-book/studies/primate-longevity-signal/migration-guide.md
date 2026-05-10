# Primate Longevity Signal Migration Guide

- study id: `primate-longevity-signal`
- course source: `external:lund/pcm1-plots-signal/script`

## Comparison Modes

- `direct_parity`: `9`

## Side-By-Side Examples

### Raw primate preprocessing to checked-in analysis table

- example id: `pcm1-data-preparation`
- family id: `data-preparation`
- fragment id: `primate-data-preprocessing`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm1-plots-signal/script#L23-L79`
- Bijux locators: `bijux_phylogenetics.evidence.studies.primate_pcm1_component_bundles:build_primate_pcm1_component_bundles, bijux_phylogenetics.evidence.studies.primate_longevity_signal:build_primate_source_fragment_map`
- Bijux summary: Bijux records the cleaned primate table as a governed evidence surface instead of burying preprocessing behind later model output.
- Why migrate: The Python-side evidence makes type repair, missing-data accounting, and grouped species decisions machine-checkable.

### Tree import, checking, node labels, and pruning

- example id: `pcm1-tree-operations`
- family id: `tree-operations`
- fragment id: `tree-import-and-pruning`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm1-plots-signal/script#L80-L120`
- Bijux locators: `bijux_phylogenetics.core.pruning:prune_tree_to_requested_taxa, bijux_phylogenetics.evidence.studies.primate_pcm1_component_bundles:build_primate_pcm1_component_bundles`
- Bijux summary: The tree-trimming and alignment surface is governed by structural parity outputs rather than visual inspection.
- Why migrate: Students can move from interactive R tree surgery to a scriptable tree-review workflow with explicit structural checks.

### Primate longevity lambda fit

- example id: `pcm1-signal-fitting`
- family id: `comparative-signal`
- fragment id: `primate-lambda-fit`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm1-plots-signal/script#L347-L354`
- Bijux locators: `bijux_phylogenetics.comparative.signal:estimate_pagels_lambda, bijux_phylogenetics.validation_corpus:build_scientific_validation_report`
- Bijux summary: Bijux turns the lecture's lambda-fitting path into a tracked parity bundle with explicit tolerance rules.
- Why migrate: The migration value is not shorter syntax; it is a governed record of where the fitted values align and where boundaries remain.

### Continuous ancestral point estimates

- example id: `pcm1-ancestral-reconstruction`
- family id: `ancestral-reconstruction`
- fragment id: `continuous-ancestral-point-estimates`
- comparison mode: `direct_parity`
- R locators: `external:lund/pcm1-plots-signal/script#L395-L399`
- Bijux locators: `bijux_phylogenetics.validation_corpus:build_scientific_validation_report, bijux_phylogenetics.evidence.studies.primate_longevity_signal:build_primate_scalar_parity_table`
- Bijux summary: Ancestral-state checks are rendered as parity rows and reviewable evidence rather than isolated notebook output.
- Why migrate: The Python-side evidence helps students compare node-level claims without manually reading long console output.
