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
- `real_local/`: opt-in local engine lanes that depend on workstation-installed tools
- `smoke/`: repository gateway and command-surface smoke coverage
- `state_models/`: finite-state kernels, substitution models, and transition math
- `topology/`: topology neighborhoods, rooted move families, and pruning refresh
- `trees/`: tree structures, tree-set workflows, clade/support review, and taxonomy
- `workflows/`: workflow orchestration, CLI/runtime integration, and engine lanes

## Nested Boundaries

- `parity/ape`, `parity/geiger`, and `parity/phytools` own toolkit-specific parity cases; the `parity/` root is reserved for cross-tool access policy.
- `bayesian/sampling`, `bayesian/priors`, `bayesian/proposals`, `bayesian/models`, `bayesian/posterior`, and `bayesian/engines` separate sampling mechanics, prior families, proposal kernels, native model definitions, posterior analysis, and wrapper-facing execution surfaces.
- `publication/integrity` owns release-surface equivalence and schema checks, `publication/packages` owns publication-bundle validation, and `publication/supplementary` owns table-generation surfaces.
- `likelihood/search` owns topology-search heuristics, convergence, replay, and multi-start orchestration, `likelihood/starting_trees` owns starting-tree construction and validation surfaces, and `likelihood/inference` owns end-to-end inference and bootstrap workflows.
- `state_models/ctmc` owns continuous-time Markov chain normalization, stationary-distribution, and transition-matrix tests so the broader state-model surface can separate likelihood and mixture behavior from core matrix math.
- `state_models/compiled` owns compiled candidate-matrix and compiled finite-state analysis surfaces so generated kernels and proof-style equivalence checks do not stay mixed with native model behavior.
- `state_models/finite_state` owns generic finite-state kernel behavior such as protocol, posterior, cache, and branch or topology update invariants, leaving family-specific likelihoods out of the shared mechanics surface.
- `state_models/discrete_traits` owns discrete-trait evolution, Mk-model, ascertainment, and history-simulation workflows so those trait-focused lanes do not remain interleaved with nucleotide or protein family tests.
- `state_models/families` owns native nucleotide, protein, and codon family likelihood surfaces together with substitution-model selection and site-report exports.
- `state_models/mixtures` owns invariant-site and discrete-gamma family behavior across nucleotide, protein, codon, and generic finite-state models so those cross-cutting mixture surfaces do not stay interleaved with base likelihood tests.
- `state_models/optimization` owns branch-length, exchangeability, kappa, and state-frequency search behavior, while `state_models/parameters` owns simplex, bounded-parameter, and parameter-layout invariants that support those optimizers.
- `trees/reporting` owns tree-review reports, publication packages, and reader-facing CLI surfaces so the `trees/` root can focus on structural tree behavior, fixtures, and single-tree algorithms.
- `trees/pruning` owns shared-taxon pruning workflows and minimal pruning proofs, `trees/reconciliation` owns duplication-loss-transfer study surfaces, and `trees/taxonomy` owns taxon table and taxonomy iteration review.
- `trees/formats` owns Newick, Nexus, and support-label parsing or writing behavior, while `trees/references` owns shared fixture and reference-governance checks for tree-adjacent study inputs.
- `workflows/engines` owns external-engine contracts, `workflows/pipelines` owns end-to-end orchestration and replay surfaces, and `workflows/benchmarks` owns governed scale and performance lanes.
- `trees/tree_sets` owns posterior, consensus, quartet, and uncertainty coverage for multi-tree analyses so single-tree structure tests can stay focused at the `trees/` root.

## Shared Fixtures

Each domain directory exposes a local `fixtures` link back to the shared
`tests/fixtures/` tree. Tests should treat that link as the stable source of
checked-in fixture material instead of assuming the suite stays flat forever.
