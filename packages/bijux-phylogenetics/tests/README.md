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
- `comparative/assessment` owns covariance audits, reference validation, uncertainty review, and maturity-assessment surfaces so analytical review stays separate from comparative model execution and publication contracts.
- `comparative/biogeography` owns geographic state, host-association, and phylogeography workflows so comparative regression and trait-evolution tests do not remain mixed with geography-specific surfaces.
- `comparative/contracts` owns comparative CLI, runtime-export, and artifact-schema surfaces so public interfaces are separated from internal analyses and audit helpers.
- `comparative/foundations` owns dataset readiness, shared comparative math, and model-matrix construction surfaces so low-level preparation and encoding behavior is separated from regression, optimization, and publication review.
- `comparative/macroevolution` owns diversification, disparity, ecological-niche, and multispecies-coalescent review surfaces so those study-style analyses are separated from generic comparative reporting.
- `comparative/optimization` owns shared bounded-search routing and scalar-fit diagnostics so comparative optimizers are separated from domain analyses that consume them.
- `comparative/publication` owns comparative report bundles, figure packages, methods-summary exports, and trait-tree publication surfaces so public-output contracts are separated from analytical model checks.
- `comparative/regression` owns PGLS, phylogenetic regression, clade-stability, residual, and signal-review surfaces so comparative statistics are separated from trait-evolution and geography-specific workflows.
- `comparative/trait_evolution` owns trait-dynamics models, mode-recovery checks, and simulation surfaces so comparative regression and reporting tests do not remain mixed with evolution-specific workflows.
- `bayesian/sampling`, `bayesian/priors`, `bayesian/proposals`, `bayesian/models`, `bayesian/posterior`, and `bayesian/engines` separate sampling mechanics, prior families, proposal kernels, native model definitions, posterior analysis, and wrapper-facing execution surfaces.
- `publication/integrity` owns release-surface equivalence and schema checks, `publication/packages` owns publication-bundle validation, and `publication/supplementary` owns table-generation surfaces.
- `parsimony/reconstruction` owns Fitch, Sankoff, ACCTRAN, DELTRAN, Dollo, Wagner, Camin-Sokal, and constrained discrete reconstruction surfaces so ancestral-state criteria stay separated from search heuristics and score summaries.
- `parsimony/placement` owns parsimony phylogenetic-placement workflows, ranked alternative ledgers, and CLI/reporting surfaces so query placement does not remain mixed with search or reconstruction lanes.
- `parsimony/comparison` owns deep-coalescence counting and reporting surfaces so gene-tree versus species-tree comparison contracts stay separated from character-parsimony workflows.
- `likelihood/search` owns topology-search heuristics, convergence, replay, and multi-start orchestration, `likelihood/starting_trees` owns starting-tree construction and validation surfaces, and `likelihood/inference` owns end-to-end inference and bootstrap workflows.
- `likelihood/optimization` owns bounded search, simplex parameterization, and fixed-topology optimization surfaces so lower-level optimizer behavior does not stay mixed with placement, support, and fit-assessment workflows.
- `likelihood/assessment` owns fit diagnostics, information-criteria review, likelihood-ratio reporting, and result-validation surfaces so evaluation contracts are separated from optimization and search mechanics.
- `likelihood/placement` owns query-placement workflows, `likelihood/support` owns approximate topology and SH-like support review, and `likelihood/wrappers` owns wrapper-correspondence and wrapper-benchmark comparison surfaces.
- `topology/comparison` owns agreement-subtree workflows, topology-distance review, and governed topology reference CLIs so cross-tree comparison surfaces are separated from move mechanics and construction helpers.
- `topology/rearrangements` owns rooted NNI, SPR, and TBR move application, neighborhood enumeration, validity review, and neighborhood-summary surfaces so topology edit mechanics stay grouped by shared move semantics.
- `topology/construction` owns deterministic rooting and stepwise-addition construction surfaces so topology-building helpers are separated from topology comparison and rearrangement analysis.
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
- `trees/structures` owns core tree-model, clade-signature, and structural-parity behavior, `trees/metrics` owns node-depth, tip-distance, and ultrametricity analysis, and `trees/simulation` owns random-tree generation and simulation-runtime surfaces.
- `workflows/engines` owns external-engine contracts, `workflows/pipelines` owns end-to-end orchestration and replay surfaces, and `workflows/benchmarks` owns governed scale and performance lanes.
- `trees/tree_sets` owns posterior, consensus, quartet, and uncertainty coverage for multi-tree analyses so single-tree structure tests can stay focused at the `trees/` root.

## Shared Fixtures

Each domain directory exposes a local `fixtures` link back to the shared
`tests/fixtures/` tree. Tests should treat that link as the stable source of
checked-in fixture material instead of assuming the suite stays flat forever.
