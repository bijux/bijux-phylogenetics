# Fragile Example Audit

This audit lists evidence bundles or bundle fragments that still depend on
narrow assumptions, seeded inputs, plot-only interpretations, or explicit
coverage boundaries.

- fragile examples: `20`

## Fragility Kinds

- `artifact_only`: `2`
- `coverage_gap`: `4`
- `instability`: `1`
- `model-boundary`: `1`
- `plot_only`: `10`
- `seeded_input_only`: `1`
- `workflow_only`: `1`

## Examples

- `ancestral-colored-tree-plot` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: The ancestral-state branch-color figure is tracked as a visual rendering surface only.
- `ape-alternate-layouts` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: These layout variants are tracked as visual surfaces only.
- `ape-longevity-overlay` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: This is a rendered trait-overlay surface and is tracked separately from ordering correctness.
- `ape-plotting-basics` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: This is a visual exploration block; the current report does not claim figure-equivalence for base `ape` plots.
- `environment-and-package-contract` — `workflow_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: This setup block is documented for reproducibility, but it is not an analysis-equivalence target.
- `final-workspace-artifact` — `artifact_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: The lecture script saves an `.RData` workspace; this report saves explicit machine-readable evidence artifacts instead.
- `ggtree-tree-visualization` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: These `ggtree` examples are tracked as visual surfaces only.
- `joined-ggtree-trait-plotting` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: Joined-data `ggtree` figures are tracked here as visual outputs only.
- `lambda-zero-visual-comparison` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: The side-by-side real-tree versus lambda=0 plots are tracked as visual outputs only.
- `ou-identifiability-boundary-cases` — `model-boundary`
  Path: `studies/comparative-trust-boundaries/evidence-003/ou-identifiability-audit.json`
  Detail: OU fits on the governed reference cases still sit on identifiability boundaries, so the repository must keep these warnings reviewer-visible instead of treating them as ordinary successful fits.
- `phytools-tree-plotting` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: The `phytools::plotTree()` surface is tracked, but no rendered-figure claim is made here.
- `primate-longevity-signal-evidence-001-pcm1-artifact-provenance-tracking` — `coverage_gap`
  Path: `studies/primate-longevity-signal/evidence-001`
  Detail: Processed files and saved workspaces are indexed as provenance outputs rather than overstated as analytical matches.
- `primate-longevity-signal-evidence-001-pcm1-reproducibility-contract-tracked` — `coverage_gap`
  Path: `studies/primate-longevity-signal/evidence-001`
  Detail: The lecture setup and package context remain visible for reviewers without being misreported as a numerical parity claim.
- `primate-longevity-signal-evidence-001-pcm1-visual-surface-tracking` — `coverage_gap`
  Path: `studies/primate-longevity-signal/evidence-001`
  Detail: Plotting examples remain indexed and reviewer-visible while figure-equivalence claims stay intentionally out of scope.
- `primate-longevity-visual-inspection` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: These are visual inspection surfaces only.
- `primate-pgls-and-signal-evidence-006-pcm2-coverage-boundary-explicit` — `coverage_gap`
  Path: `studies/primate-pgls-and-signal/evidence-006`
  Detail: Unimplemented EB, transformed-tree, and ancestral parity surfaces are indexed as open trust boundaries rather than silently implied.
- `processed-analysis-artifacts` — `artifact_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: This block writes the processed CSV and trimmed tree that the Python evidence pass later consumes.
- `random-simulation-inputs` — `seeded_input_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: For credibility, the report freezes these R-generated simulation inputs with `set.seed(1)` and reuses the resulting artifacts on both sides.
- `random-simulation-plotting` — `plot_only`
  Path: `studies/primate-longevity-signal/evidence-001/scientific_debt_register.json`
  Detail: These are visual teaching plots and are tracked separately from the simulation inputs and signal fits.
- `weak-signal-threshold-instability` — `instability`
  Path: `studies/comparative-trust-boundaries/evidence-002/result-instability-audit.json`
  Detail: The weak-signal case crosses the 0.05 decision boundary across governed reruns, so any single-run significance claim would be overconfident.
