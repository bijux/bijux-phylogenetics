# Continuous Mode Recovery Panel

This packaged simulation panel benchmarks Bijux continuous `fitContinuous`-style recovery against deterministic simulation truth and stored local `geiger` recovery references.

The panel contains:

- `trees/reference-tree-twelve-taxa.nwk`: the governed twelve-taxon rooted review tree for Brownian, OU, early-burst, and weak-OU identifiability cases.
- `trees/reference-tree-twenty-four-taxa.nwk`: the governed twenty-four-taxon ultrametric review tree for Pagel-lambda, Pagel-kappa, and Pagel-delta transformed-branch cases.
- `simulation-cases.tsv`: deterministic case definitions, including seeds, true parameters, declared recovery tolerances, expected model-choice outcomes, and expected warning kinds.
- `expected/`: governed reviewer-facing outputs regenerated from the packaged trees and case table, including Bijux-versus-`geiger` parameter comparisons and the stored `geiger` summary ledger used by the benchmark.

The governed `expected/` bundle contains these reviewer-facing artifacts:

- `workflow-summary.tsv`: one package-level metrics row for the panel.
- `recovery-summary.tsv`: one recovery summary row per governed case.
- `parameter-recovery.tsv`: one truth-versus-fit parameter row per engine and parameter, including OU optimum/root-state recovery.
- `parameter-comparison.tsv`: one paired Bijux-versus-`geiger` comparison row per governed parameter, including the OU optimum/root-state surface.
- `model-choice.tsv`: candidate-model ranking rows for both recovery engines.
- `execution-review.tsv`: fit and model-comparison execution status rows.
- `warning-review.tsv`: governed identifiability and boundary-warning rows.
- `geiger-reference.tsv`: stored local `geiger` reference summaries keyed by case.
- `simulated-traits/`: one deterministic simulated trait table per governed case.

The cases cover:

- Brownian sigma-squared recovery.
- OU alpha, sigma-squared, and optimum recovery.
- Early-burst rate-change recovery.
- Weak-OU identifiability review, where warning transparency and Brownian-like model support are the expected outcome.
- Pagel-lambda transformed-branch review.
- Pagel-kappa transformed-branch review.
- Pagel-delta transformed-branch review.
