# Test Suite Layout

The runtime test suite is organized by durable product domains rather than a
single flat namespace.

## Directory Boundaries

- `alignment/`: alignment transforms, partitioning, and sequence layout review
- `ancestral/`: ancestral reconstruction, summaries, and sensitivity analysis
- `bayesian/`: Bayesian inference models, priors, samplers, and diagnostics
- `comparative/`: comparative biology, phylogenetic regression, and trait review
- `datasets/`: checked-in study datasets, reference bundles, and benchmark inputs
- `dating/`: time trees, calibrations, tip dates, and dating workflows
- `distance/`: distance matrices, clustering, and minimum-evolution methods
- `evidence/`: evidence-book contracts, teaching surfaces, and release truth
- `likelihood/`: likelihood search, optimization, and branch-scoring workflows
- `parity/`: live and governed reference-parity validation surfaces
- `parsimony/`: parsimony search, reconstructions, and character scoring
- `publication/`: publication packages, supplementary tables, and golden outputs
- `state_models/`: finite-state kernels, substitution models, and transition math
- `topology/`: topology neighborhoods, rooted move families, and pruning refresh
- `trees/`: tree structures, tree-set workflows, clade/support review, and taxonomy
- `workflows/`: workflow orchestration, CLI/runtime integration, and engine lanes

## Shared Fixtures

Each domain directory exposes a local `fixtures` link back to the shared
`tests/fixtures/` tree. Tests should treat that link as the stable source of
checked-in fixture material instead of assuming the suite stays flat forever.
