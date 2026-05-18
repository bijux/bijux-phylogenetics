# Continuous Mode Recovery Panel

This packaged simulation panel benchmarks Bijux continuous `fitContinuous`-style recovery against deterministic simulation truth and stored local `geiger` recovery references.

The panel contains:

- `trees/reference-tree-twelve-taxa.nwk`: the governed twelve-taxon rooted review tree for Brownian, OU, early-burst, and weak-OU identifiability cases.
- `trees/reference-tree-twenty-four-taxa.nwk`: the governed twenty-four-taxon ultrametric review tree for Pagel-lambda, Pagel-kappa, and Pagel-delta transformed-branch cases.
- `simulation-cases.tsv`: deterministic case definitions, including seeds, true parameters, declared recovery tolerances, expected model-choice outcomes, and expected warning kinds.
- `expected/`: governed reviewer-facing outputs regenerated from the packaged trees and case table, including Bijux-versus-`geiger` parameter comparisons and the stored `geiger` summary ledger used by the benchmark.

The cases cover:

- Brownian sigma-squared recovery.
- OU alpha, sigma-squared, and optimum recovery.
- Early-burst rate-change recovery.
- Weak-OU identifiability review, where warning transparency and Brownian-like model support are the expected outcome.
- Pagel-lambda transformed-branch review.
- Pagel-kappa transformed-branch review.
- Pagel-delta transformed-branch review.
