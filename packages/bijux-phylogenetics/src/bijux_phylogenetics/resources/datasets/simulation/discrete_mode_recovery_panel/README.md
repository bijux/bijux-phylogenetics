# Discrete Mode Recovery Panel

This packaged simulation panel benchmarks Bijux discrete Mk recovery against known transition-rate truth and stored local `geiger::fitDiscrete` recovery references.

The packaged Python surface lives under
`bijux_phylogenetics.datasets.discrete_mode_recovery` and is intentionally
split by ownership: `panel` exposes dataset metadata, `scenarios` owns the
case-table parser, `workflow` reruns the benchmark, `bundle` writes the
reviewer-facing ledgers, and `demo` materializes the full review bundle.

The panel contains:

- `trees/reference-tree-twelve-taxa.nwk`: the governed twelve-taxon rooted review tree for the overparameterized five-state ARD failure surface.
- `trees/reference-tree-twenty-four-taxa.nwk`: the governed twenty-four-taxon rooted review tree for recoverable ER and SYM cases plus a weak-identification ARD review surface.
- `simulation-cases.tsv`: deterministic case definitions, including state vocabularies, explicit transition-rate matrices, optional branch-length transform truth, tolerated recovery error, transform-parameter tolerances, expected model-choice outcomes, and overparameterized-review expectations.
- `expected/`: governed reviewer-facing outputs regenerated from the packaged trees and case table, including recovery summaries, rate and transform-parameter ledgers, and the stored `geiger` summary ledger used by the benchmark.

The cases cover:

- Three-state equal-rates recovery.
- Three-state symmetric recovery.
- Three-state all-rates-different weak-identification review, where the benchmark must report rate error and model-choice failure honestly instead of claiming successful recovery.
- Five-state overparameterized ARD review, where the benchmark must report overparameterized failure evidence instead of claiming successful rate recovery.
