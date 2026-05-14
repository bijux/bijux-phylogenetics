# Continuous Mode Recovery Panel

This packaged simulation panel validates Brownian, Ornstein-Uhlenbeck, and early-burst continuous trait models against deterministic simulation truth on one shared rooted tree.

The panel contains:

- `reference-tree.nwk`: the rooted tree used for every simulation-recovery case.
- `simulation-cases.tsv`: deterministic case definitions, including seeds, true parameters, declared recovery tolerances, expected model-choice outcomes, and expected warning kinds.
- `expected/`: governed reviewer-facing outputs regenerated from the packaged tree and case table.

The cases cover:

- Brownian sigma-squared recovery.
- OU alpha, sigma-squared, and optimum recovery.
- Early-burst rate-change recovery.
- Weak-OU identifiability review, where warning transparency and Brownian-like model support are the expected outcome.
