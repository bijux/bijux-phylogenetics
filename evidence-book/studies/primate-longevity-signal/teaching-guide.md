# Primate Longevity Signal Teaching Guide

- study id: `primate-longevity-signal`
- categories: `teaching-study, migration-study`
- course source: `external:lund/pcm1-plots-signal/script`
- concept tags: `ancestral-states, confidence-intervals, course-context, data-import, figure-boundaries, handoff, likelihood-ratio, missing-data, node-estimates, package-loading, pagel-lambda, plotting, processed-exports, pruning, random-seeds, reproducibility, signal-inputs, signal-testing, simulation, species-aggregation, teaching-visuals, topology, tree-data-alignment, tree-import, type-repair, workspace-artifacts`

## Workflow contracts

- family id: `workflow-contracts`
- verdict: `not_comparable`
- coverage: `covered`
- concept tags: `reproducibility, package-loading, course-context`
- evidence ids: `evidence-001`
- fragment ids: `environment-and-package-contract`
- teaching narrative: Students see which package and environment assumptions belong to the lecture setup before any numerical claim is made.

## Data preparation

- family id: `data-preparation`
- verdict: `matched`
- coverage: `covered`
- concept tags: `data-import, type-repair, missing-data, species-aggregation`
- evidence ids: `evidence-001, evidence-002, evidence-003, evidence-004, evidence-005, evidence-008`
- fragment ids: `primate-data-preprocessing, primate-longevity-vector-assembly, tip-order-alignment`
- teaching narrative: The lecture data-cleaning path is broken into explicit steps so students can inspect what must match before comparative inference is trusted.

## Tree operations

- family id: `tree-operations`
- verdict: `matched`
- coverage: `covered`
- concept tags: `tree-import, pruning, topology, tree-data-alignment`
- evidence ids: `evidence-001, evidence-006, evidence-007, evidence-008`
- fragment ids: `extract-clade-node-77, rotate-nodes-behavior, tree-import-and-pruning, treeio-node-mapping-and-join, unrooted-tree-demo`
- teaching narrative: Tree loading, pruning, and topology operations are separated from plotting so students can review structure before interpretation.

## Visual surfaces

- family id: `visual-surfaces`
- verdict: `not_comparable`
- coverage: `covered`
- concept tags: `plotting, figure-boundaries, teaching-visuals`
- evidence ids: `evidence-001`
- fragment ids: `ancestral-colored-tree-plot, ape-alternate-layouts, ape-longevity-overlay, ape-plotting-basics, ggtree-tree-visualization, joined-ggtree-trait-plotting, lambda-zero-visual-comparison, phytools-tree-plotting, primate-longevity-visual-inspection, random-simulation-plotting`
- teaching narrative: Plotting-oriented blocks remain visible as course material while the repo stays honest that rendered-figure equivalence is not yet claimed.

## Simulation inputs

- family id: `simulation-inputs`
- verdict: `matched`
- coverage: `covered`
- concept tags: `simulation, random-seeds, signal-inputs`
- evidence ids: `evidence-001`
- fragment ids: `random-simulation-inputs`
- teaching narrative: The lecture's seeded random inputs are frozen so comparisons focus on signal-fitting behavior rather than hidden simulation drift.

## Comparative signal

- family id: `comparative-signal`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `pagel-lambda, likelihood-ratio, signal-testing`
- evidence ids: `evidence-001`
- fragment ids: `lambda-zero-covariance-and-lrt, primate-lambda-fit, random-signal-lambda-fits`
- teaching narrative: Students can follow how lambda estimation and lambda-zero testing are checked numerically instead of being accepted from screenshots or prose.

## Ancestral reconstruction

- family id: `ancestral-reconstruction`
- verdict: `matched_with_tolerance`
- coverage: `covered`
- concept tags: `ancestral-states, node-estimates, confidence-intervals`
- evidence ids: `evidence-001`
- fragment ids: `ancestral-table-assembly, bonobo-gibbon-mrca-estimate, continuous-ancestral-intervals, continuous-ancestral-point-estimates, lifespan-increase-counts`
- teaching narrative: Ancestral-state outputs stay tied to explicit node estimates and interval summaries so reviewers can inspect the exact numerical claim.

## Artifact provenance

- family id: `artifact-provenance`
- verdict: `not_comparable`
- coverage: `covered`
- concept tags: `processed-exports, workspace-artifacts, handoff`
- evidence ids: `evidence-001, evidence-009`
- fragment ids: `final-workspace-artifact, processed-analysis-artifacts`
- teaching narrative: Saved tables and trees are treated as governed handoff artifacts so downstream lessons do not depend on undocumented files.
