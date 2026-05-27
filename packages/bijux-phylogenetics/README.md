# bijux-phylogenetics

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/bijux-phylogenetics/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml/badge.svg?branch=main)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://img.shields.io/badge/release-pypi%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://img.shields.io/badge/release-ghcr%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://img.shields.io/badge/release-github%20workflow-2563EB?logo=githubactions&logoColor=white)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)
[![phylogenetic](https://img.shields.io/pypi/v/phylogenetic?label=phylogenetic&logo=pypi)](https://pypi.org/project/phylogenetic/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)
[![phylogenetic](https://img.shields.io/badge/phylogenetic-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetic)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
[![phylogenetic docs](https://img.shields.io/badge/docs-phylogenetic-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

Runtime package for the bijux-phylogenetics repository.

This package provides the Python API and CLI for tree validation, inspection,
comparison, metadata linkage, comparative trait analysis, ancestral-state
reconstruction, discrete-state evolution analysis, external engine
orchestration, Bayesian posterior summarization, diversification and
macroevolution analysis, evidence bundle creation, explicit parsimony scoring,
and HTML report generation.

## Install

`bijux-phylogenetics` supports Python 3.11 and newer.

```bash
python3.11 -m pip install bijux-phylogenetics
bijux-phylogenetics --help
```

The installed runtime also ships small packaged example inputs. You can copy
them into one writable directory without relying on a source checkout:

```python
from pathlib import Path

from bijux_phylogenetics.core import copy_example_inputs

copy_example_inputs(Path("artifacts/example-inputs"))
```

## Python Workflow Surface

The stable notebook-and-pipeline surface now lives under
`bijux_phylogenetics.api`.

It exposes the same serious runtime objects the CLI already uses for:

- FASTA validation
- multiple-sequence alignment
- alignment trimming
- full FASTA-to-tree execution
- maximum-likelihood tree inference
- branch-support estimation
- topology comparison
- PGLS comparative modeling
- discrete ancestral reconstruction
- reviewer-facing report generation
- one-command config-driven workflow execution

```python
from pathlib import Path

from bijux_phylogenetics.api import (
    render_report_workflow,
    run_comparative_model_workflow,
    run_sequence_to_tree_workflow,
)

workflow = run_sequence_to_tree_workflow(
    Path("dataset/sequences.fasta"),
    out_dir=Path("artifacts/sequence-to-tree"),
    sequence_type="dna",
)

comparative = run_comparative_model_workflow(
    Path("dataset/tree.nwk"),
    Path("dataset/traits.tsv"),
    response="response",
    predictors=["predictor_one"],
    lambda_value=1.0,
)

report = render_report_workflow(
    tree_path=workflow.output_paths["tree"],
    alignment_path=workflow.output_paths["trimmed_alignment"],
    traits_path=Path("dataset/traits.tsv"),
    metadata_path=Path("dataset/metadata.tsv"),
    out_path=Path("artifacts/sequence-to-tree/report.html"),
)

workflow.write_json(Path("artifacts/sequence-to-tree/workflow.json"))
workflow.write_tsv(Path("artifacts/sequence-to-tree/workflow.tsv"))
comparative.write_tsv(Path("artifacts/comparative-model.tsv"))
```

The Python workflow helpers return typed workflow result objects. Each wrapper
delegates the underlying CLI-grade report fields directly and adds stable
`write_json(...)` serialization plus `write_tsv(...)` where a tabular summary
is meaningful.

## Native Tree Model

The owned in-memory tree surface now lives on `bijux_phylogenetics.PhyloTree`
and `bijux_phylogenetics.TreeNode`.

That native model carries:

- stable node IDs
- parent-child edges
- branch lengths
- rooted or unrooted state
- internal and terminal labels
- node and edge metadata
- deterministic preorder and postorder traversal
- deep-copy support without shared mutation
- native canonical Newick parsing and writing through `PhyloTree.to_newick()` and `PhyloTree.from_newick(...)`, including branch lengths, internal labels, quoted labels, multi-tree loading, and location-aware parse errors

Tree import, tree transformation, native distance review, comparative
covariance, ancestral reconstruction, and simulation surfaces now operate on
that same owned structure instead of inventing per-method node identity rules.
Posterior and bootstrap tree-set summaries now do the same for Newick tree
sets: consensus building, clade-frequency mapping, clade-support mapping,
topology clustering, and posterior tree-set comparison all read one
`PhyloTree` record per Newick statement, including `.trees` files that contain
plain Newick records. Strict consensus and support workflows fail explicitly
when a tree set does not share one exact taxon set, while tolerant inspection
surfaces still skip malformed records and report the skipped-record count.

## Method Tiers

Serious workflow and report surfaces now expose an explicit runtime method
tier so the package does not overclaim what it has done.

The governed tiers are:

- `supported`: validated by reference parity or real-engine validation
- `experimental`: usable for exploratory work, but emitted with an explicit warning
- `advisory`: review or audit output that should not be mistaken for inference
- `parser-only`: summaries of external-engine artifacts that do not claim Bijux performed the inference

Current examples include:

- `adapter fasta-to-tree` as `supported`
- `comparative logistic` as `experimental`, with its approximation method and explicit `ape::compar.gee` non-claim reported
- `report tree-package` as `advisory`
- Bayesian report builders such as `adapter mrbayes-report` and `adapter beast-calibration-report` as `parser-only`

JSON CLI payloads expose `method_tier`, `method_inference_mode`,
`method_validation_basis`, and `method_approximation` when applicable. Python
workflow result objects expose the same information through `result.report.method_tier`
or the underlying report object directly.

## Stable Artifact Schemas

The canonical reviewer-facing TSV and JSON artifacts are schema-governed and
tested so downstream scripts do not have to infer headers or top-level keys
from incidental examples.

The governed stable families now include:

- FASTA-to-tree `prefix.model.tsv`
- FASTA-to-tree `prefix.support.tsv`
- clade review tables such as `clade-table.tsv`
- branch review tables such as `host-switch-branches.tsv`
- derived comparative trait tables such as `comparative-traits.tsv`
- comparative model-output tables such as `comparative-summary.tsv`
- event ledgers such as `biogeography/event-table.tsv`
- workflow and report manifests such as `prefix.manifest.json` and `comparative-report.manifest.json`

When one of those schemas changes, the repository test suite now fails on the
exact header or key drift instead of only failing later through broader fixture
differences.

## Current Scope

- parse Newick trees and FASTA alignments
- inspect tree shape and branch-length health
- inspect internal child counts, singleton nodes, missing internal versus terminal branch lengths, branch-length outlier nodes, support normalization, rootedness confidence, and tree-assumption compatibility
- classify tree validity, biological safety, unsafe external labels, node-label conflicts, and downstream forensic readiness for topology, time-tree, comparative, visualization, and publication use
- normalize unsafe taxon labels and audit normalization collisions
- prune trees from explicit taxa, exclusions, traits, or metadata tables
- quantify pruning information loss across taxa, clades, branch length, and metadata retention
- classify internal node labels as support-like or name-like and detect suspicious or mixed support scales
- compare shared clades, clade changes, and shared-split branch lengths between trees
- validate trait and metadata linkage against tree tips
- check comparative readiness for rooted trees and numeric traits
- compute phylogenetic independent contrasts, Blomberg's K, Pagel's lambda, and permutation-based signal tests with explicit missing-value pruning, seeded reproducibility, constant-trait rejection, and rooted non-ultrametric acceptance reporting
- fit standalone Brownian or OU continuous-trait models with confidence intervals, residual diagnostics, model comparison, and leave-one-taxon-out sensitivity
- fit phylogenetic generalized least-squares models with numeric, categorical, and interaction predictors through explicit formula auditing
- adjust repeated comparative hypothesis tests with Benjamini-Hochberg correction and emit integrated comparative audit, influence, tree-comparison, pruning-comparison, and reviewer-facing report outputs
- score unordered discrete character matrices on one tree with Fitch parsimony, emit per-character step ledgers plus per-node candidate-state sets, accept explicit per-character weight tables, and fail explicitly on missing taxa, unknown states, empty matrices, or invalid weights
- score ordered discrete character matrices on one tree with Wagner parsimony, emit weighted per-character costs plus internal-node ordered cost vectors, accept explicit per-character weight tables, and require either ordinal labels or one explicit state order
- score discrete character matrices on one tree with Sankoff parsimony under one user-supplied transition-cost matrix, validate missing states, negative costs, symmetry-by-default, nonzero diagonals, and unused states before scoring, emit minimum-cost ledgers plus per-node cost vectors and optimal-state tie sets, accept explicit per-character weight tables, and allow asymmetric costs only when requested explicitly
- score binary character matrices on one tree with Dollo parsimony, emit one-gain and many-loss branch ledgers, accept explicit per-character weight tables, and fail explicitly when multistate traits have not been binarized first
- score binary character matrices on one rooted tree with irreversible Camin-Sokal parsimony, assume ancestral absence, allow repeated `0→1` gains, forbid `1→0` reversals, emit exact gain-branch ledgers, accept explicit per-character weight tables, and fail explicitly when multistate traits have not been binarized first
- reconstruct unordered parsimony ambiguities with ACCTRAN, emit resolved node states plus governed `ancestral_states.tsv` and branch-mapped change ledgers with stable `branch_id` and ambiguity flags, accept explicit per-character weight tables, and prefer earlier branch placements over later ones while preserving the same minimum tree length
- reconstruct unordered parsimony ambiguities with DELTRAN, emit resolved node states plus governed `ancestral_states.tsv` and branch-mapped change ledgers with stable `branch_id` and ambiguity flags, accept explicit per-character weight tables, and delay ambiguous changes toward terminal branches while preserving the same minimum tree length
- bootstrap small-taxon parsimony inference by character resampling with replacement, infer one exact-search replicate tree per resample, emit governed replicate-tree, replicate-draw, replicate-score, clade-frequency, consensus-tree, and clade-support artifacts, and map support onto the inferred reference tree by descendant clade identity rather than node index
- jackknife small-taxon parsimony inference by character subsampling without replacement, record retained-character counts and retained-character ledgers per replicate, infer one exact-search replicate tree per subsample, emit governed replicate-tree, score, clade-frequency, consensus-tree, and clade-support artifacts, and keep the without-replacement retention policy explicit through `retain_probability`
- hill-climb one rooted binary starting tree by rooted NNI, score every neighbor through the same governed parsimony methods, accept only score-improving moves, and emit one deterministic `search_trace.tsv` plus `start_tree.nwk`, `final_tree.nwk`, and `run.json` with the start score, every accepted move, the final score and topology, and the stopping reason
- hill-climb one rooted binary starting tree by rooted subtree-prune-regraft, reject self-regrafts by generating targets only from the pruned remainder tree, score every legal regraft through the same governed parsimony methods, accept only score-improving moves, and emit one deterministic `search_trace.tsv` plus `start_tree.nwk`, `final_tree.nwk`, and `run.json` with the start score, every accepted move, the final score and topology, and the stopping reason
- run one deterministic parsimony ratchet over rooted SPR search by temporarily upweighting a seeded random character subset each cycle, restoring the normal weights after the perturbed search, and emitting governed `cycle_history.tsv`, `best_tree_history.tsv`, `start_tree.nwk`, `final_tree.nwk`, `best_tree.nwk`, and `run.json` so the exact perturbation cycles, best scores, and best-tree history are reproducible
- summarize generic parsimony tree length with per-character raw scores, explicit character weights, and one governed total score across Fitch, Wagner, Sankoff, Dollo, Camin-Sokal, ACCTRAN, or DELTRAN methods
- compute per-character and aggregate consistency index across supported parsimony methods, exclude constant characters by explicit `0/0` policy, keep nonconstant uninformative characters in the score, and leave arbitrary Sankoff step-matrix minima outside the current owned surface
- compute per-character and aggregate retention index for unordered Fitch-style methods, emit `null` plus explicit zero-range reasons when `max = min`, and avoid claiming ordered, irreversible, or arbitrary step-matrix maxima that the package does not yet own
- compute per-character and aggregate rescaled consistency index on the common supported CI and RI method set, emit rows with `character_id`, `ci`, `ri`, `rc`, and `undefined_reason`, and derive RC directly from the tested CI and RI outputs instead of recomputing separate logic
- run one generic finite-state Felsenstein pruning kernel over rooted trees, keep the conditional-likelihood recursion independent from any one Mk trait surface, and use the same postorder owner as the native foundation for upcoming DNA and protein maximum-likelihood models
- compress identical aligned site columns into stable integer-weighted site patterns and sum weighted site log likelihoods through the same owned likelihood foundation, so repeated-column alignments do not change total JC69-, HKY-, or GTR-style likelihoods
- evaluate fixed-topology JC69 DNA likelihood natively from explicit branch lengths with the shared pruning kernel and site-pattern compression, and improve one supplied topology by bounded coordinate-wise branch-length optimization instead of reusing corrected-distance JC69
- evaluate fixed-topology K80 DNA likelihood natively with explicit transition versus transversion probabilities, fit one fixed-topology `kappa` value by bounded likelihood search on transition-biased data, and keep the equal-frequency K80 ownership boundary separate from upcoming F81 and HKY85 surfaces
- evaluate fixed-topology F81 DNA likelihood natively with supplied or empirically estimated unequal base frequencies, expose the resulting likelihood and AIC surface on one supplied topology, and keep the unequal-frequency ownership boundary separate from equal-frequency JC69 and K80
- evaluate fixed-topology HKY85 DNA likelihood natively with both unequal base frequencies and transition bias active, fit one fixed-topology `kappa` value while carrying supplied or empirically estimated base frequencies through the report, and expose the resulting likelihood and AIC surface without reducing HKY85 to K80 or F81
- evaluate fixed-topology GTR DNA likelihood natively with six identifiable exchangeabilities plus unequal base frequencies, anchor the exchangeability report at `AC=1`, fit one fixed-topology relative exchangeability surface by bounded coordinate search, and expose the resulting likelihood and AIC surface without reducing GTR to HKY85
- optimize fixed-topology nucleotide substitution parameters through one selected JC69, K80, F81, HKY85, or GTR surface that records model-specific starts, search bounds, convergence state, boundary warnings, and governed recovery behavior instead of hardcoding one parameter set per model
- evaluate fixed-topology nucleotide nested likelihood-ratio tests only for one declared set of nested model pairs across JC69, K80, F81, HKY85, and GTR surfaces, emitting null and alternative fits, log-likelihoods, statistic, degrees of freedom, p-value, and explicit boundary caveat instead of comparing nonnested pairs
- compare fixed-topology nucleotide substitution candidates across JC69, K80, F81, HKY85, GTR, and selected `+G` or `+I` variants in one ranked model-selection table that keeps failed fits visible and emits per-row log-likelihood, parameter count, AIC, AICc, BIC, delta-AIC, Akaike weight, and warning ledgers
- hill-climb one rooted binary nucleotide starting tree by likelihood-improving rooted NNI moves, reoptimize branch lengths on every candidate topology under one explicit coordinate-search branch policy, and emit governed input-tree, start-tree, final-tree, trace, and run payload artifacts that make the local optimum and every accepted move auditable
- export fixed-topology nucleotide site log likelihoods as explicit expanded site rows with stable pattern identifiers, pattern weights, and one governed `site_log_likelihoods.tsv` surface whose row sum matches the tested total likelihood
- compute fixed-topology nucleotide marginal ancestral sequence probabilities as explicit internal-node × site × state posterior rows across selected JC69, K80, F81, HKY85, and GTR surfaces, preserving stable node identities, expanded site positions, and full posterior mass instead of collapsing directly to one most-likely state per site
- export thresholded marginal ancestral nucleotide FASTA sequences for internal nodes together with one governed uncertainty table that preserves per-site posterior support, marks low-confidence sites by explicit probability policy, and avoids emitting FASTA alone without the paired probability artifact
- reconstruct fixed-topology nucleotide joint ancestral sequences as one globally optimal internal-node assignment per site across selected JC69, K80, F81, HKY85, and GTR surfaces, preserving explicit node-site joint-state rows and keeping the joint optimum distinct from independent marginal maxima
- evaluate fixed-topology 20-state protein Poisson likelihood natively with the shared pruning kernel and site-pattern compression, keep amino-acid observation policy separate from DNA likelihood surfaces, and carry explicit gap and missing-state policy through the protein likelihood report instead of reusing a nucleotide alphabet model
- evaluate fixed-topology empirical 20-state protein likelihoods from one validated 20×20 rate matrix plus one optional root prior, keep matrix validation explicit at the owned boundary, and expose a matrix-sensitive protein likelihood surface where changing the empirical matrix changes the likelihood on the same alignment
- evaluate fixed-topology empirical protein likelihoods under discrete-gamma rate heterogeneity by emitting explicit category rates, equal weights, alpha, and per-site mixture likelihood rows instead of accepting alpha without changing the likelihood surface
- evaluate fixed-topology empirical protein likelihoods under one invariant-site mixture by fitting a nonzero invariant proportion on invariant-rich alignments, emitting per-site invariant versus variable component rows, and keeping invariant-site mass inside the likelihood surface instead of only reporting invariant-site diagnostics
- evaluate fixed-topology empirical protein likelihoods under one combined discrete-gamma plus invariant mixture by keeping both gamma categories and invariant mass active in the same site likelihood, reporting alpha, invariant proportion, boundary warnings, and total likelihood from one owned `+G+I` surface instead of merging standalone `+G` and `+I` outputs
- optimize branch lengths on one fixed empirical protein topology under one selected fixed-rate, `+G`, `+I`, or `+G+I` likelihood surface by bounded coordinate search, keeping every branch nonnegative and reporting the exact optimized tree plus per-branch start versus optimized lengths
- reconstruct continuous ancestral states under Brownian or OU-style trait models
- reconstruct discrete ancestral states under Fitch parsimony or likelihood-style ER, SYM, and ARD models with explicit ambiguity and low-confidence reporting, root-prior controls, and fitted transition-rate ledgers
- validate discrete ancestral likelihood surfaces against governed `ape::ace` references plus root-prior, ordered, irreversible, and ambiguity policy checks
- compare continuous ancestral reconstructions across two supported models, summarize ancestral sensitivity across model, tree, pruning, or coding choices, and package publication-ready ancestral figures
- validate discrete geographic state coding, detect incomplete ordered vocabularies, estimate ancestral node states under ordered or unordered assumptions, compare equal-rates, symmetric, and all-rates-different models, export node and transition tables, highlight model-sensitive ancestral regions, simulate seeded stochastic character maps from a fitted discrete-state CTMC with branch-segment and time-in-state ledgers, and render discrete-state HTML reports
- audit alignment inference readiness, validate model-selection outputs against engine artifacts, verify inferred-tree taxa against the alignment, inspect metadata-group clustering, classify inference failures, and validate bootstrap tree sets before interpreting engine outputs
- estimate lineage-through-time curves, Pybus-Harvey gamma statistics, simple Yule or birth-death diversification rates, sampling-aware corrections, clade outlier summaries, trait-linked diversification tables, and publication-oriented diversification figure packages for rooted ultrametric trees, while explicitly excluding `geiger::medusa` parity because the current diversification surface does not yet implement stepwise rate-shift search or shift-count model growth, and explicitly excluding `geiger::bd.ms` birth-death parity because the current Yule and birth-death surfaces are heuristic summaries rather than one geiger-matched diversification likelihood contract
- run governed MAFFT-, trimAl-, IQ-TREE-, and FastTree-style external workflows with captured commands, versions, logs, and warning summaries
- prepare BEAST XML analyses, validate prepared XML assumptions, parse BEAST posterior logs or annotation-rich tree sets, compute burn-in-aware ESS and posterior summaries, and distinguish prepared versus parsed versus recorded prior inference evidence in reviewer-facing diagnostics
- prepare and run deterministic MrBayes analyses, summarize posterior trees after burn-in filtering, parse parameter traces, and compute per-parameter ESS values
- compare fast approximate and maximum-likelihood trees through the same deterministic tree-comparison report surface
- resume inference only when saved manifests, inputs, and outputs still match, and render standalone HTML inference workflow reports from those manifests
- export joined metadata rows and missing trait-value diagnostics
- inspect alignment alphabets, composition, GC content, duplicates, raw-sequence length outliers, sliding-window quality, suspicious alignment regions, coding stop codons, and frameshift-like sequence lengths
- classify FASTA inputs as aligned, raw-sequence, or equal-length-but-shape-ambiguous and report method-specific alignment readiness
- detect mixed coding versus noncoding behavior inside the same nucleotide dataset
- define named alignment-filtering profiles, generate cleaned alignments, compare original versus cleaned versions, and warn when filtering removes signal or biases taxon groups
- score alignment quality with transparent components and emit one-shot alignment forensic reports
- detect low-information alignments, ambiguity-heavy columns, duplicate-handling policy needs, and per-sequence quality rankings before inference
- quantify per-taxon locus occupancy, per-locus taxon occupancy, site-coverage fractions, low-coverage taxa, and low-coverage loci from concatenated multi-locus alignments with retained-matrix exports
- audit tree, metadata, traits, alignment, tip dates, and calibrations together through one-shot dataset readiness decisions
- render dedicated reviewer-facing alignment, dataset, phylo-input, and taxonomy HTML reports with machine-readable sidecars
- validate checked-in Level 1 reference fixtures, aggregate workflow coverage, document known failure cases, classify workflow maturity, and render reviewer-facing workflow-validation or release-gate reports
- aggregate pytest JUnit outputs, method tiers, flagship dataset inventory, parity coverage, and stress coverage into one reviewer-facing release-truth report
- generate taxon crosswalk tables, completeness matrices, exclusion tables, ordering-drift audits, pruning-step retention summaries, and named readiness levels for reviewer-facing dataset inspection
- trim all-gap or all-missing columns and remove high-missingness sequences
- translate coding nucleotide alignments to amino-acid alignments and export pairwise identity matrices
- compute raw or p-distance, Jukes-Cantor, Kimura 2-parameter, Felsenstein 81, Tamura-Nei 93, or amino-acid p-distance matrices with explicit gap-handling, ambiguity policies, and model-parameter reporting
- compute rooted or unrooted pairwise tip-distance matrices from branch-length trees with explicit taxon order and owned wide or long-form exports
- audit saturated pairs, unusually divergent pairs, and low-information pairs before distance-based tree building
- build Neighbor-Joining, BIONJ, single-linkage, complete-linkage, or owned rooted-ultrametric UPGMA and WPGMA trees from computed distance matrices, bootstrap site-resampled trees, summarize clade support, and write reproducibility bundles
- validate imported long-form distance matrices, detect nonmetric violations, and build trees from imported distances
- audit NJ and UPGMA method assumptions, including explicit UPGMA ultrametric-clock violations for computed or imported distance matrices
- load posterior tree sets, compute consensus trees, and export clade-frequency or pairwise tree-distance summaries
- cluster identical rooted topologies, detect unstable taxa or clades, and compare two posterior tree sets
- simulate birth-death or coalescent trees, Brownian, OU, or early-burst continuous traits, discrete traits, and DNA or protein alignments
- ship governed internal recovery panels that check whether Brownian, OU, and early-burst trait models recover known simulation truth with explicit parameter-error and warning ledgers
- benchmark validation, tree comparison, large-tree rendering/report generation, large-alignment trimming/distance/readiness, and alignment diagnostics across increasing problem sizes
- root trees on explicit outgroups or reroot them by midpoint
- audit rooting, ordering, clade extraction, and pruning transforms with before/after summaries and retained-versus-removed taxon reasoning
- validate tree roundtrips across Newick, Nexus, and phyloXML formats with topology-preservation checks, support-label audits, and semantic-loss warnings
- audit ambiguous taxon identities, synonym candidates, namespace mixing, workflow taxon loss, and cross-run taxon stability before downstream comparison or linkage
- produce HTML reports and file-level evidence manifests

The dataset audit runtime is also available as a typed Python package under
`bijux_phylogenetics.core.dataset`. Use `summarize_dataset_readiness()` for
the comparative preflight, `build_dataset_crosswalk()`,
`build_dataset_completeness_matrix()`, and
`build_dataset_mismatch_report()` for reviewer-facing taxon-surface ledgers,
`audit_dataset_taxon_ordering()` for silent order-drift detection, and
`audit_dataset_inputs()` for the integrated blocker, risk, exclusion, and
minimal-fix workflow. Inside that package, shared dataclasses live in
`models`, dataset loading and surface introspection live in `context`,
cross-surface ledgers live in `crosswalk`, order-drift logic lives in
`ordering`, comparative preflight lives in `readiness`, audit policy lives in
`audit_policy`, and the integrated orchestration lives in `workflow`, so the
package root stays a curated API gateway instead of remaining the real
implementation file.

## Example CLI Runs

```bash
bijux-phylogenetics alignment classify sequences.fasta --json
bijux-phylogenetics alignment profiles --json
bijux-phylogenetics alignment windows alignment.fasta --window-size 50 --step-size 10 --json
bijux-phylogenetics alignment readiness alignment.fasta --json
bijux-phylogenetics alignment quality alignment.fasta --json
bijux-phylogenetics alignment low-information alignment.fasta --json
bijux-phylogenetics alignment duplicate-policy alignment.fasta --identity-threshold 0.99 --json
bijux-phylogenetics alignment ambiguous-columns alignment.fasta --threshold 0.5 --json
bijux-phylogenetics alignment sequence-ranking alignment.fasta --json
bijux-phylogenetics alignment concatenate loci/gene-alpha.fasta loci/gene-beta.fasta loci/gene-gamma.fasta --out artifacts/supermatrix.aln.fasta --partitions-out artifacts/supermatrix.partitions.txt --matrix-out artifacts/supermatrix.matrix.tsv --json
bijux-phylogenetics alignment occupancy supermatrix.fasta partitions.txt --taxon-coverage-threshold 0.6 --locus-coverage-threshold 0.6 --minimum-locus-occupancy 0.75 --taxa-out artifacts/occupancy/taxa.tsv --loci-out artifacts/occupancy/loci.tsv --matrix-out artifacts/occupancy/matrix.tsv --filtered-alignment-out artifacts/occupancy/filtered.fasta --filtered-partitions-out artifacts/occupancy/filtered-partitions.txt --json
bijux-phylogenetics alignment length-outliers sequences.fasta --json
bijux-phylogenetics alignment forensic alignment.fasta --json
bijux-phylogenetics alignment filter alignment.fasta --profile moderate --out cleaned.fasta --json
bijux-phylogenetics alignment compare alignment.fasta cleaned.fasta --json
bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4
bijux-phylogenetics alignment distance-matrix alignment.fasta --model raw --gap-handling complete-deletion --out distances.tsv
bijux-phylogenetics alignment distance-matrix alignment.fasta --model f81 --components-out distance-components.tsv --parameters-out distance-parameters.tsv --out distances.tsv
bijux-phylogenetics alignment distance-quality alignment.fasta --model jukes-cantor --json
bijux-phylogenetics alignment distance-saturation alignment.fasta --model jukes-cantor --json
bijux-phylogenetics alignment distance-additivity alignment.fasta --model p-distance --out-dir artifacts/distance-additivity --json
bijux-phylogenetics alignment distance-ultrametricity alignment.fasta --model p-distance --tolerance 1e-6 --json
bijux-phylogenetics alignment distance-suitability alignment.fasta --model jukes-cantor --json
bijux-phylogenetics alignment distance-assumptions alignment.fasta --model p-distance --json
bijux-phylogenetics alignment build-tree alignment.fasta --method bionj --out bionj-tree.nwk
bijux-phylogenetics alignment build-tree alignment.fasta --method neighbor-joining --missing-distance-policy nearest-valid --out alignment-tree.nwk
bijux-phylogenetics alignment build-tree alignment.fasta --method upgma --out upgma-tree.nwk
bijux-phylogenetics alignment build-tree alignment.fasta --method complete-linkage --out complete-linkage-tree.nwk
bijux-phylogenetics distance additivity distances.tsv --out-dir artifacts/distance-additivity --json
bijux-phylogenetics distance ultrametricity distances.tsv --tolerance 1e-6 --json
bijux-phylogenetics distance build-tree distances.tsv --method neighbor-joining --missing-distance-policy triangle-bound --out distance-tree.nwk
bijux-phylogenetics distance minimum-evolution matrix.tsv fixed-topology.nwk --out fitted-tree.nwk
bijux-phylogenetics distance fitch-margoliash matrix.tsv fixed-topology.nwk --out fitted-tree.nwk --json
bijux-phylogenetics distance ordinary-least-squares matrix.tsv fixed-topology.nwk --out fitted-tree.nwk --json
bijux-phylogenetics distance nonnegative-least-squares matrix.tsv fixed-topology.nwk --out fitted-tree.nwk --json
bijux-phylogenetics distance patristic-residuals matrix.tsv fixed-topology.nwk --out-dir artifacts/patristic-residuals --json
bijux-phylogenetics distance taxon-influence matrix.tsv reference-tree.nwk --method neighbor-joining --missing-distance-policy mean-impute --out-dir artifacts/distance-taxon-influence --json
bijux-phylogenetics distance taxon-jackknife matrix.tsv --method neighbor-joining --missing-distance-policy mean-impute --out-dir artifacts/distance-taxon-jackknife --json
bijux-phylogenetics distance method-comparison matrix.tsv --out-dir artifacts/distance-method-comparison --json
bijux-phylogenetics distance bme-nni-search matrix.tsv --start-method bionj --out-dir artifacts/distance-bme-nni --json
bijux-phylogenetics alignment compare-distance-to-tree alignment.fasta inferred-tree.nwk --method neighbor-joining --json
bijux-phylogenetics alignment bootstrap-tree alignment.fasta --method neighbor-joining --replicates 200 --support-out artifacts/distance-support.tsv --tree-set-out artifacts/distance-bootstrap.trees --draws-out artifacts/distance-bootstrap-draws.tsv --json
bijux-phylogenetics alignment distance-support-summary alignment.fasta --method neighbor-joining --replicates 50 --json
bijux-phylogenetics alignment distance-models alignment.fasta --json
bijux-phylogenetics alignment distance-gap-sensitivity alignment.fasta --model p-distance --json
bijux-phylogenetics alignment distance-method-report alignment.fasta --method neighbor-joining --replicates 50 --json
bijux-phylogenetics alignment distance-maturity alignment.fasta --method neighbor-joining --replicates 50 --json
bijux-phylogenetics alignment distance-bundle alignment.fasta --method neighbor-joining --replicates 200 --out-dir artifacts/distance-bundle --json
bijux-phylogenetics distance validate distances.tsv --json
bijux-phylogenetics distance quality distances.tsv --json
bijux-phylogenetics distance assumptions distances.tsv --json
bijux-phylogenetics distance reference --json
bijux-phylogenetics report alignment --alignment alignment.fasta --out artifacts/alignment-report.html --json
bijux-phylogenetics report workflow-validation --out artifacts/workflow-validation-report.html --json
bijux-phylogenetics report release-gate --out artifacts/level-1-release-gate.html --json
bijux-phylogenetics report release-truth --test-report artifacts/pytest/full-suite.xml --real-engine-test-report artifacts/pytest/real-engine.xml --out artifacts/release-truth-report.html --json
bijux-phylogenetics report supplementary-alignment-table --alignment alignment.fasta --filtered-alignment filtered-alignment.fasta --out artifacts/supplementary-alignment.tsv --json
bijux-phylogenetics report supplementary-ancestral-state-table --tree tree.nwk --traits traits.tsv --trait habitat --reconstruction-kind discrete --model equal-rates --out artifacts/supplementary-ancestral-states.tsv --json
bijux-phylogenetics phylo parsimony fitch tree.nwk discrete-character-matrix.tsv --character-weights character-weights.tsv --out-dir artifacts/parsimony-fitch --json
bijux-phylogenetics phylo parsimony wagner tree.nwk ordered-character-matrix.tsv --state-order low,medium,high,very_high --character-weights character-weights.tsv --out-dir artifacts/parsimony-wagner --json
bijux-phylogenetics phylo parsimony sankoff tree.nwk discrete-character-matrix.tsv transition-costs.tsv --character-weights character-weights.tsv --out-dir artifacts/parsimony-sankoff --json
bijux-phylogenetics phylo parsimony sankoff tree.nwk discrete-character-matrix.tsv transition-costs.tsv --allow-asymmetric-costs --out-dir artifacts/parsimony-sankoff-asymmetric --json
bijux-phylogenetics phylo parsimony dollo tree.nwk binary-character-matrix.tsv --character-weights character-weights.tsv --out-dir artifacts/parsimony-dollo --json
bijux-phylogenetics phylo parsimony camin-sokal tree.nwk binary-character-matrix.tsv --character-weights character-weights.tsv --out-dir artifacts/parsimony-camin-sokal --json
bijux-phylogenetics phylo parsimony acctran tree.nwk discrete-character-matrix.tsv --character-weights character-weights.tsv --out-dir artifacts/parsimony-acctran --json
bijux-phylogenetics phylo parsimony deltran tree.nwk discrete-character-matrix.tsv --character-weights character-weights.tsv --out-dir artifacts/parsimony-deltran --json
bijux-phylogenetics phylo parsimony bootstrap discrete-character-matrix.tsv --method fitch --replicate-count 200 --seed 7 --out-dir artifacts/parsimony-bootstrap --json
bijux-phylogenetics phylo parsimony jackknife discrete-character-matrix.tsv --method fitch --replicate-count 200 --seed 7 --retain-probability 0.75 --out-dir artifacts/parsimony-jackknife --json
bijux-phylogenetics phylo parsimony nni-search start-tree.nwk discrete-character-matrix.tsv --method fitch --out-dir artifacts/parsimony-nni-search --json
bijux-phylogenetics phylo parsimony spr-search start-tree.nwk discrete-character-matrix.tsv --method fitch --out-dir artifacts/parsimony-spr-search --json
bijux-phylogenetics phylo parsimony ratchet start-tree.nwk discrete-character-matrix.tsv --method fitch --cycle-count 20 --seed 7 --perturbed-character-count 4 --perturbation-factor 2.0 --out-dir artifacts/parsimony-ratchet --json
bijux-phylogenetics phylo parsimony tree-length tree.nwk discrete-character-matrix.tsv --method fitch --character-weights character-weights.tsv --out-dir artifacts/parsimony-tree-length --json
bijux-phylogenetics phylo parsimony consistency-index tree.nwk discrete-character-matrix.tsv --method fitch --out-dir artifacts/parsimony-consistency --json
bijux-phylogenetics phylo parsimony retention-index tree.nwk discrete-character-matrix.tsv --method fitch --out-dir artifacts/parsimony-retention --json
bijux-phylogenetics phylo parsimony rescaled-consistency-index tree.nwk discrete-character-matrix.tsv --method fitch --out-dir artifacts/parsimony-rescaled-consistency --json
bijux-phylogenetics report supplementary-batch-summary-table --workflow-bundle-root workflow/ --out artifacts/supplementary-batch-summary.tsv --json
bijux-phylogenetics report supplementary-clade-support-table --tree tree.nwk --comparison-tree-set posterior-trees.nwk --out artifacts/supplementary-clade-support.tsv --json
bijux-phylogenetics report supplementary-comparative-model-table --tree tree.nwk --traits traits.tsv --formula 'response ~ predictor_one' --formula 'response ~ predictor_one + predictor_two' --lambda-value 0.0 --out artifacts/supplementary-comparative-model.tsv --json
bijux-phylogenetics report supplementary-diversification-table --tree tree.nwk --metadata metadata.tsv --clade-model birth-death --out artifacts/supplementary-diversification.tsv --json
bijux-phylogenetics report supplementary-model-selection-table --iqtree-report run.iqtree --model-sidecar run.model --out artifacts/supplementary-model-selection.tsv --json
bijux-phylogenetics report supplementary-tree-table --tree tree.nwk --out artifacts/supplementary-tree.tsv --json
bijux-phylogenetics report package-comparison artifacts/package-a/rabies-cross-host-geography-package.manifest.json artifacts/package-b/rabies-cross-host-geography-package.manifest.json --out-dir artifacts/package-comparison --json
bijux-phylogenetics report package-revalidation artifacts/rabies-cross-host-geography-package.manifest.json --out-dir artifacts/package-revalidation --json
bijux-phylogenetics report tree-inference-methods-summary workflow.manifest.json --out artifacts/tree-inference-methods-summary.md --json
bijux-phylogenetics report alignment-filtering-methods-summary alignment.fasta --profile moderate --group-table metadata.tsv --group-column region --out artifacts/alignment-filtering-methods-summary.md --json
bijux-phylogenetics report ancestral-methods-summary tree.nwk traits.tsv --trait habitat --kind discrete --model equal-rates --out artifacts/ancestral-methods-summary.md --json
bijux-phylogenetics report reviewer-audit-checklist artifacts/tree-report.manifest.json --out artifacts/reviewer-audit-checklist.tsv --json
bijux-phylogenetics report tree-validation-methods-summary tree.nwk --out artifacts/tree-validation-methods-summary.md --json
bijux-phylogenetics report taxonomy --tree tree.nwk --synonym-table taxonomy.tsv --metadata metadata.tsv --traits traits.tsv --alignment alignment.fasta --reported-taxa reviewer-table.tsv --out artifacts/taxonomy-report.html --json
bijux-phylogenetics report supplementary-taxon-table --tree tree.nwk --metadata metadata.tsv --traits traits.tsv --alignment alignment.fasta --filtered-alignment filtered-alignment.fasta --inference-tree inferred.nwk --reported-taxa reviewer-table.tsv --out artifacts/supplementary-taxa.tsv --json
bijux-phylogenetics taxonomy synonyms tree.nwk --synonym-table synonyms.tsv --json
bijux-phylogenetics taxonomy resolve-synonyms tree.nwk --synonym-table synonyms.tsv --out normalized-tree.nwk --mapping-out synonym-map.tsv --json
bijux-phylogenetics taxonomy namespaces tree.nwk --json
bijux-phylogenetics taxonomy loss tree.nwk --metadata metadata.csv --traits traits.csv --alignment alignment.fasta --filtered-alignment filtered.fasta --inference-tree inferred.nwk --reported-taxa reported.csv --json
bijux-phylogenetics taxonomy stability --run tree=tree.nwk --run alignment=alignment.fasta --run filtered=filtered.fasta --json
bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative discrete-mk tree.nwk traits.tsv --trait habitat --taxon-column species --model equal-rates --summary-out artifacts/discrete-mk-summary.tsv --rates-out artifacts/discrete-mk-rates.tsv --json
bijux-phylogenetics comparative discrete-mk tree.nwk traits.tsv --trait habitat --taxon-column species --model equal-rates --transform lambda --summary-out artifacts/discrete-mk-lambda-summary.tsv --rates-out artifacts/discrete-mk-lambda-rates.tsv --json
bijux-phylogenetics comparative brownian tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative compare-models tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative model-comparison-package tree.nwk traits.tsv --trait height_cm --out-dir artifacts/model-comparison-package --json
bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json
bijux-phylogenetics comparative pgls tree.nwk traits.tsv --formula "height_cm ~ body_mass * habitat" --json
bijux-phylogenetics comparative phylogenetic-residuals tree.nwk traits.tsv --response brain_mass --predictor body_mass --summary-out artifacts/phylogenetic-residual-summary.tsv --residuals-out artifacts/phylogenetic-residuals.tsv --coefficients-out artifacts/phylogenetic-residual-coefficients.tsv --excluded-taxa-out artifacts/phylogenetic-residual-excluded.tsv --json
bijux-phylogenetics comparative phylogenetic-anova tree.nwk traits.tsv --response trait_value --group habitat --simulations 199 --seed 17 --summary-out artifacts/phylogenetic-anova-summary.tsv --groups-out artifacts/phylogenetic-anova-groups.tsv --pairwise-out artifacts/phylogenetic-anova-pairwise.tsv --simulations-out artifacts/phylogenetic-anova-simulations.tsv --excluded-taxa-out artifacts/phylogenetic-anova-excluded.tsv --json
bijux-phylogenetics comparative covariance-audit tree.nwk traits.tsv --analysis pgls --formula "height_cm ~ body_mass + habitat" --summary-out artifacts/covariance-audit-summary.tsv --candidates-out artifacts/covariance-audit-candidates.tsv --excluded-taxa-out artifacts/covariance-audit-excluded.tsv --json
bijux-phylogenetics comparative multiple-testing tree.nwk traits.tsv --responses height_cm range_km --predictors body_mass log_range --json
bijux-phylogenetics comparative report tree.nwk traits.tsv --formula "height_cm ~ body_mass + habitat" --out artifacts/comparative-report.html --json
bijux-phylogenetics comparative report tree.nwk traits.tsv --formula "height_cm ~ body_mass + habitat" --methods-summary-out artifacts/comparative-methods-summary.md --json
bijux-phylogenetics comparative compare-trees tree-a.nwk tree-b.nwk traits.tsv --response height_cm --predictors body_mass log_range --json
bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --json
bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --estimator anc-ml --json
bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --estimator fast-anc --json
bijux-phylogenetics ancestral discrete tree.nwk traits.tsv --trait habitat --model equal-rates --root-prior-mode empirical --summary-out artifacts/ancestral-discrete-summary.tsv --probabilities-out artifacts/ancestral-discrete-probabilities.tsv --transitions-out artifacts/ancestral-discrete-transitions.tsv --json
bijux-phylogenetics ancestral sensitivity tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --json
bijux-phylogenetics ancestral report tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --out artifacts/ancestral-report.html
bijux-phylogenetics ancestral package tree.nwk traits.tsv --trait habitat --kind discrete --model symmetric --state-ordering ordered --ordered-states low,medium,high --out-dir artifacts/ancestral-package --json
bijux-phylogenetics discrete-evolution model tree.nwk geography.tsv --trait region --model symmetric --state-ordering ordered --ordered-states north,south,island --node-table-out artifacts/node-states.tsv --transitions-out artifacts/transitions.tsv --json
bijux-phylogenetics discrete-evolution stochastic-map tree.nwk geography.tsv --trait region --model symmetric --replicates 200 --collection-out artifacts/geography-maps.json --summary-out artifacts/geography-stochastic-summary.tsv --state-times-out artifacts/geography-stochastic-state-times.tsv --branch-occupancy-out artifacts/geography-stochastic-branch-occupancy.tsv --segments-out artifacts/geography-stochastic-segments.tsv --count-matrix-out artifacts/geography-stochastic-counts.tsv --aggregate-matrix-out artifacts/geography-stochastic-aggregate.tsv --branch-transition-out artifacts/geography-stochastic-branch-transitions.tsv --events-out artifacts/geography-stochastic-events.tsv --json
bijux-phylogenetics discrete-evolution summarize-maps artifacts/geography-maps.json --summary-out artifacts/geography-stochastic-summary.tsv --state-times-out artifacts/geography-stochastic-state-times.tsv --branch-occupancy-out artifacts/geography-stochastic-branch-occupancy.tsv --json
bijux-phylogenetics discrete-evolution count-maps artifacts/geography-maps.json --count-matrix-out artifacts/geography-stochastic-counts.tsv --aggregate-matrix-out artifacts/geography-stochastic-aggregate.tsv --branch-transition-out artifacts/geography-stochastic-branch-transitions.tsv --events-out artifacts/geography-stochastic-events.tsv --json
bijux-phylogenetics discrete-evolution density-maps artifacts/geography-maps.json --focal-state island --branch-probabilities-out artifacts/geography-density-branches.tsv --density-branches-out artifacts/geography-density-envelope.tsv --density-slices-out artifacts/geography-density-slices.tsv --out artifacts/geography-density.html --json
bijux-phylogenetics discrete-evolution report tree.nwk geography.tsv --trait region --compare-model symmetric --out artifacts/geography-report.html
bijux-phylogenetics diversification gamma-stat tree.nwk --metadata sampling.tsv --out artifacts/diversification-gamma-statistic.tsv --json
bijux-phylogenetics diversification estimate tree.nwk --metadata sampling.tsv --model birth-death --json
bijux-phylogenetics diversification bd-ms tree.nwk --metadata sampling.tsv --json
bijux-phylogenetics diversification medusa tree.nwk --metadata sampling.tsv --json
bijux-phylogenetics diversification methods-summary tree.nwk --metadata sampling.tsv --traits traits.tsv --trait habitat --out artifacts/diversification-methods-summary.md --json
bijux-phylogenetics diversification package tree.nwk --metadata sampling.tsv --out-dir artifacts/diversification-figures --json
bijux-phylogenetics diversification report tree.nwk --metadata sampling.tsv --traits traits.tsv --trait habitat --out artifacts/diversification-report.html --methods-summary-out artifacts/diversification-methods-summary.md --json
```

The same diversification surface is available as a typed Python package under
`bijux_phylogenetics.comparative.diversification`. Use
`validate_time_tree_for_diversification()` and
`detect_incomplete_taxon_sampling_metadata()` for the rooted-ultrametric and
sampling contracts, `compute_lineage_through_time_curve()`,
`compute_diversification_gamma_statistic()`, and
`estimate_diversification_rate()` for the core analytical lanes,
`build_diversification_method_report()` and
`render_diversification_report()` for reviewer-facing summaries, and
`build_diversification_figure_package()` for the publication-style SVG bundle.
The package is intentionally split by ownership into `trees`, `sampling`,
`lineage`, `rates`, `clades`, `traits`, `reporting`, and `figure_package`
modules instead of one flat diversification grab-bag.

```bash
bijux-phylogenetics adapter align unaligned.fasta --out aligned.fasta --mode linsi --json
bijux-phylogenetics adapter align coding-cds.fasta --out coding.aligned.fasta --mode linsi --codon-aware --json
bijux-phylogenetics adapter trim aligned.fasta --out trimmed.fasta --mode automated1 --json
bijux-phylogenetics alignment sequence-type raw-sequences.fasta --json
bijux-phylogenetics alignment validate-input raw-sequences.fasta --json
bijux-phylogenetics alignment repair-input raw-sequences.fasta --out artifacts/raw-sequences.repaired.fasta --normalize-identifiers --remove-invalid-records --json
bijux-phylogenetics adapter fasta-to-tree raw-sequences.fasta --out-dir artifacts/fasta-to-tree --prefix mammals --alignment-mode einsi --trimming-mode strictplus --normalize-identifiers --remove-invalid-records --bootstrap-replicates 1000 --json
bijux-phylogenetics adapter model-select alignment.fasta --out-dir artifacts/model-select --prefix mammals --json
bijux-phylogenetics adapter infer-ml alignment.fasta --out-dir artifacts/ml --model GTR+G --prefix mammals --json
bijux-phylogenetics adapter bootstrap alignment.fasta --out-dir artifacts/bootstrap --model GTR+G --replicates 1000 --prefix mammals --json
bijux-phylogenetics adapter consensus artifacts/bootstrap/mammals.ufboot --out-dir artifacts/consensus --prefix mammals --json
bijux-phylogenetics adapter infer-fast alignment.fasta --out artifacts/fasttree.nwk --json
bijux-phylogenetics adapter compare --fast-tree artifacts/fasttree.nwk --ml-tree artifacts/ml/mammals.treefile --out artifacts/engine-comparison.html --split-table-out artifacts/engine-comparison-splits.tsv --json
bijux-phylogenetics adapter beast-prepare alignment.fasta --out artifacts/beast/analysis.xml --tree tree.nwk --calibrations calibrations.tsv --tip-dates tip-dates.tsv --clock-model strict --tree-prior yule --json
bijux-phylogenetics adapter beast-xml artifacts/beast/analysis.xml --json
bijux-phylogenetics adapter beast-run artifacts/beast/analysis.xml --resume --json
bijux-phylogenetics adapter beast-log artifacts/beast/analysis.1.log --burnin-fraction 0.1 --summary-out artifacts/beast/log-summary.tsv --json
bijux-phylogenetics adapter beast-trees artifacts/beast/analysis.1.trees --burnin-fraction 0.1 --tree-set-out artifacts/beast/postburnin.nwk --json
bijux-phylogenetics adapter beast-consensus artifacts/beast/analysis.1.trees --burnin-fraction 0.1 --out artifacts/beast/consensus.nwk --clade-table-out artifacts/beast/clades.tsv --json
bijux-phylogenetics adapter bayesian-methods artifacts/beast/analysis.1.trees --log artifacts/beast/analysis.1.log --analysis-xml artifacts/beast/analysis.xml --out artifacts/beast/methods.md --json
bijux-phylogenetics adapter bayesian-evidence --out-dir artifacts/beast/evidence --inputs alignment.fasta calibrations.tsv tip-dates.tsv --configs artifacts/beast/analysis.xml --trees tree.nwk artifacts/beast/consensus.nwk --logs artifacts/beast/analysis.1.log --diagnostics artifacts/beast/log-summary.tsv --reports artifacts/beast/methods.md --json
bijux-phylogenetics adapter mrbayes-prepare alignment.fasta --out artifacts/mrbayes/analysis.nex --ngen 20000 --samplefreq 100 --json
bijux-phylogenetics adapter mrbayes-run artifacts/mrbayes/analysis.nex --resume --json
bijux-phylogenetics adapter mrbayes-summarize artifacts/mrbayes/analysis.run1.t --burnin-fraction 0.25 --json
bijux-phylogenetics adapter mrbayes-traces artifacts/mrbayes/analysis.run1.p --json
bijux-phylogenetics adapter mrbayes-ess artifacts/mrbayes/analysis.run1.p --json
bijux-phylogenetics adapter report artifacts/mrbayes/analysis.manifest.json --out artifacts/mrbayes/inference-report.html --json
bijux-phylogenetics tree-set inspect posterior.trees --json
bijux-phylogenetics tree-set diversity posterior.trees --out artifacts/posterior.rf-distribution.tsv --json
bijux-phylogenetics tree-set consensus posterior.trees --out consensus.nwk --method majority-rule --clade-frequencies-out artifacts/posterior.clade-frequencies.tsv
bijux-phylogenetics tree-set support-map reference-tree.nwk posterior.trees --out artifacts/reference-tree-support.tsv --json
bijux-phylogenetics tree-set methods-summary posterior.trees --out artifacts/tree-set-uncertainty-methods-summary.md --json
bijux-phylogenetics tree-set report posterior.trees --out artifacts/tree-uncertainty-report.html --max-tree-count 5000 --max-report-table-rows 50 --memory-warning-threshold-bytes 134217728
bijux-phylogenetics demo gnathostome-ortholog-protein-benchmark --out artifacts/gnathostome-ortholog-protein-benchmark --json
bijux-phylogenetics demo rabies-method-sensitivity-panel --out artifacts/rabies-method-sensitivity-panel --json
bijux-phylogenetics demo rabies-cross-host-geography-panel --out artifacts/rabies-cross-host-geography-panel --config src/bijux_phylogenetics/resources/datasets/pathogens/rabies_cross_host_geography_panel/workflow-config.json --json
bijux-phylogenetics demo catarrhine-data-quality-stress-panel --out artifacts/catarrhine-data-quality-stress-panel --json
bijux-phylogenetics demo known-answer-reference-panel --out artifacts/known-answer-reference-panel --json
bijux-phylogenetics demo real-dataset-macroevolution --out artifacts/real-dataset-macroevolution --json
bijux-phylogenetics simulate tree-birth-death --tree-count 5 --tip-count 16 --out simulated.trees
bijux-phylogenetics simulate tree-random --tree-count 64 --tip-count 12 --out simulated-random.trees --record-table-out artifacts/random-tree-records.tsv --envelope-table-out artifacts/random-tree-envelope.tsv --json
bijux-phylogenetics simulate tree-coalescent --tree-count 64 --tip-count 12 --out simulated-coalescent.trees --record-table-out artifacts/coalescent-tree-records.tsv --envelope-table-out artifacts/coalescent-tree-envelope.tsv --json
bijux-phylogenetics simulate traits-early-burst tree.nwk --root-state 1.0 --sigma 0.5 --rate-change 4.0 --out simulated-early-burst.tsv --json
bijux-phylogenetics demo continuous-mode-recovery-panel --out artifacts/continuous-mode-recovery-panel --json
bijux-phylogenetics simulate alignment-dna tree.nwk --sequence-length 500 --out simulated-alignment.fasta
bijux-phylogenetics benchmark tree-comparison --replicates 3 --json
bijux-phylogenetics benchmark large-tree-scaling --replicates 1 --tip-count 512 --tip-count 1024 --json
bijux-phylogenetics benchmark large-alignment-scaling --replicates 1 --sequence-count 256 --alignment-length 512 --sequence-count 512 --alignment-length 1024 --json
bijux-phylogenetics benchmark large-tree-set-scaling --replicates 1 --tree-count 128 --tip-count 48 --tree-count 256 --tip-count 64 --json
bijux-phylogenetics benchmark large-tree-model-fitting --tier small --json
bijux-phylogenetics benchmark real-dataset-macroevolution --json
bijux-phylogenetics benchmark workflow-practical-limits --replicates 1 --stress-tier heavy --json
bijux-phylogenetics report production-scale-readiness --replicates 1 --tree-tip-count 512 --tree-tip-count 1024 --sequence-count 256 --alignment-length 512 --sequence-count 512 --alignment-length 1024 --posterior-tree-count 128 --tree-set-tip-count 48 --posterior-tree-count 256 --tree-set-tip-count 64 --stress-tier heavy --out artifacts/production-scale-readiness.html --json
bijux-phylogenetics report alignment-package alignment.fasta --out-dir artifacts/alignment-quality-package --json
bijux-phylogenetics render tree.nwk --metadata metadata.tsv --label-column species --metadata-strip-columns location --traits traits.tsv --categorical-column habitat --continuous-column height_cm --heatmap-columns height_cm,status --layout phylogram --support-labels --package-dir artifacts/tree-publication-package --out artifacts/tree-publication-package/figure.svg --json
bijux-phylogenetics report trait-tree-package tree.nwk --metadata metadata.tsv --traits traits.tsv --label-column species --categorical-column habitat --continuous-column height_cm --metadata-strip-columns location --heatmap-columns height_cm --support-labels --out-dir artifacts/trait-tree-package --json
bijux-phylogenetics diagnose assumptions tree.nwk --metadata metadata.tsv --json
bijux-phylogenetics alignment translate coding.fasta --out translated.fasta --codon-validation-out artifacts/codon-validation.tsv --excluded-sequences-out artifacts/translation-exclusions.tsv
bijux-phylogenetics report dataset tree.nwk metadata.tsv traits.tsv --alignment alignment.fasta --tip-dates tip-dates.tsv --calibrations calibrations.tsv --out artifacts/dataset-report.html --json
bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk
bijux-phylogenetics phylo preflight --workflow fasta-to-tree --json
bijux-phylogenetics phylo run workflow-config.yaml --json
```

For reviewer-facing scale checks on owned tree workflows, use
`benchmark large-tree-scaling`. It benchmarks `tree-validation`,
`tree-comparison`, `tree-rendering`, and `tree-reporting` on governed
high-taxon synthetic trees so validation, SVG rendering, and HTML report costs
stay visible before larger production claims are made.

For reviewer-facing scale checks on owned alignment workflows, use
`benchmark large-alignment-scaling`. It benchmarks
`alignment-diagnostics`, `alignment-trimming`, `distance-analysis`, and
`alignment-readiness` on governed large aligned FASTA classes so diagnostics,
trimAl-backed workflow cost, distance-report cost, and inference-readiness
review stay visible before larger production claims are made.

For governed continuous `fitContinuous` recovery review, use
`demo continuous-mode-recovery-panel`. It reruns a deterministic seven-case
panel across two rooted trees, writes Bijux-versus-stored-`geiger`
parameter-recovery ledgers, including the OU optimum/root-state recovery
surface, records execution and warning review rows, and keeps
weak-identifiability transformed-model cases explicit instead of forcing them
into binary pass-or-fail claims.

The same packaged panel is available as a typed Python surface under
`bijux_phylogenetics.datasets.continuous_mode_recovery`. Use
`load_continuous_mode_recovery_panel_dataset()` to inspect the packaged trees
and case table, `run_continuous_mode_recovery_panel_workflow()` to rerun the
governed recovery report in memory,
`write_continuous_mode_recovery_panel_workflow_bundle(...)` to emit the
reviewer-facing ledgers, or `run_continuous_mode_recovery_panel_demo(...)` to
materialize the dataset copy plus the workflow outputs in one call.

For governed discrete `fitDiscrete` recovery review, use
`demo discrete-mode-recovery-panel`. It reruns the packaged discrete Mk panel,
writes recovery summaries plus rate-recovery, rate-comparison,
transform-parameter recovery, and transform-parameter comparison ledgers,
records execution and warning review rows, and keeps review-only ARD cases
explicit instead of flattening weak identification into success claims.

The same packaged panel is available as a typed Python surface under
`bijux_phylogenetics.datasets.discrete_mode_recovery`. Use
`load_discrete_mode_recovery_panel_dataset()` to inspect the packaged trees
and case table, `run_discrete_mode_recovery_panel_workflow()` to rerun the
governed recovery report in memory,
`write_discrete_mode_recovery_panel_workflow_bundle(...)` to emit the
reviewer-facing ledgers, or `run_discrete_mode_recovery_panel_demo(...)` to
materialize the dataset copy plus the workflow outputs in one call.

For governed owned known-answer recovery review, use
`demo known-answer-reference-panel`. It materializes the packaged simulation
dataset, rebuilds the distance tree, reruns the continuous and discrete
recovery surfaces on stored truth artifacts, and writes explicit recovery
ledgers for parameters, internal nodes, branch events, and declared threshold
checks. The same panel is exposed as a typed Python surface under
`bijux_phylogenetics.datasets.known_answer_reference` through
`load_known_answer_reference_dataset()`,
`run_known_answer_reference_workflow()`,
`write_known_answer_reference_workflow_bundle(...)`, and
`run_known_answer_reference_demo(...)`.

For governed large-tree macroevolution fitting review, use
`benchmark large-tree-model-fitting`. The small tier runs one 100-taxon
Pagel-lambda fit and one 100-taxon binary discrete ER fit, records owned
runtime, peak memory, and optimizer-step budgets, and compares the resulting
fit summaries against stored local `geiger` references. The heavy tier adds
one 512-taxon Brownian continuous-fit review case so 500-plus-taxon fitting is
still exercised on a real owned surface even when transformed-branch optimizers
are materially slower, and it keeps slowdown explicit through threshold-review
rows instead of silently treating the larger surface as proven.

For governed real-data macroevolution review, use
`benchmark real-dataset-macroevolution` or the release-style bundle command
`demo real-dataset-macroevolution`. This benchmark runs on the published
Central European seashore flora subset from Dryad (`10.5061/dryad.0st06f0`),
compares a continuous `seed_mass` model table across Brownian, white-noise,
Pagel-lambda, OU, and early-burst against stored local `geiger` references,
compares a discrete `lifeform` model table across ER, SYM, and ARD against
stored local `geiger` references, records data provenance, and writes both a
model-table parity ledger and a governed missing-and-mismatched-taxon review
ledger. The benchmark is intentionally cautious: the selected OU surface
remains near the lower alpha boundary, so continuous interpretation is
reported as a weak covariance-shape preference rather than strong
adaptive-process proof, and the sparse five-state discrete surface is reported
with selection-level agreement to `geiger` without pretending that SYM and ARD
row scores are clean owned-reference parity on this real dataset.

For reviewer-facing scale checks on posterior and bootstrap tree-set review,
use `benchmark large-tree-set-scaling`. It benchmarks
`tree-set-consensus`, `pairwise-rf-diversity`, `topology-clustering`, and
`uncertainty-summaries` on governed large tree-set classes so consensus cost,
pairwise RF aggregation, topology-mode collapse, and uncertainty-summary review
stay visible before larger production claims are made.

For one governed summary of the largest workflow classes the repository
currently exercises, use `benchmark workflow-practical-limits`. It aggregates
the tested maxima from the large-tree, large-alignment, large-tree-set, and
stress-suite lanes so reviewer-facing limits for taxa, aligned sites, tree
count, and posterior size stay explicit instead of being inferred from
individual benchmark tables.

For one reviewer-facing classification of which owned workflows are currently
safe at small, medium, large, and HPC scale, use
`report production-scale-readiness`. It turns the governed practical-limit
evidence into one HTML and JSON report that declares each workflow's applicable
scale dimensions, highest proven readiness tier, and the workflows currently
ready at each scale threshold.

For a journal-oriented single-tree figure bundle, use `render` with
`--package-dir`. The package writes one vector `figure.svg`, a structured
`figure-caption.md` draft, a machine-readable `figure-legend.tsv`, the tip
annotation ledger, and a manifest that records support-label validation,
scale-bar audit, legend completeness, and label-legibility review before the
figure is treated as publication-ready.

Every publication-oriented figure package now also writes one dedicated
reproducibility manifest. That artifact records the exact input files, any
explicit filters, the model or analysis mode, the reviewer-facing settings,
the generated figure surfaces, and the generated tables so downstream review
does not depend on reverse-engineering package-specific JSON.

For one reviewer-facing supplementary alignment ledger, use
`report supplementary-alignment-table`. It writes one TSV with one row per
sequence identifier, keeping original missingness, gaps, ambiguity burden,
alignment-wide variable and parsimony-informative site counts, low-information
status, and optional filtered-alignment retention status together on one
reviewable surface.

For one reviewer-facing supplementary tree diagnostics ledger, use
`report supplementary-tree-table`. It writes one TSV with one row per tree
source, keeping topology shape, branch-length health, support-class counts,
rootedness and root-confidence classification, ultrametricity, safety flags,
and the merged warning ledger together on one reviewable surface.

For one reviewer-facing tree-validation methods summary, use
`report tree-validation-methods-summary`. It writes one Markdown summary over
the tree validation, inspection, and forensic safety surfaces, naming the
checks that ran, the governing numeric thresholds, the repair-required labels
or exclusion conditions, and the downstream analyses that remain allowed or
blocked.

For one reviewer-facing alignment-filtering methods summary, use
`report alignment-filtering-methods-summary`. It writes one Markdown summary
over one profile-driven filtering pass, naming the selected cleaning profile,
removed site and sequence reasons, retained alignment dimensions, signal-loss
warnings, and optional metadata or trait group-retention skew on one
reviewable surface.

For one reviewer-facing tree-inference methods summary, use
`report tree-inference-methods-summary`. It writes one Markdown summary over
one governed `fasta-to-tree` workflow manifest, naming the alignment and
trimming steps, the selected substitution model, the maximum-likelihood and
bootstrap-support engine steps, the final tree handoff, and the workflow
warnings that constrain interpretation.

For one reviewer-facing ancestral-reconstruction methods summary, use
`report ancestral-methods-summary`. It writes one Markdown summary over one
continuous or discrete ancestral reconstruction, naming the analyzed trait or
states, taxon pruning, reconstruction model and policy surface, node-level
uncertainty contract, and interpretation risks that remain attached to the
internal-node estimates.

For one reviewer-facing tree-set uncertainty methods summary, use
`tree-set methods-summary`. It writes one Markdown summary over one posterior,
bootstrap, or other governed tree set, naming the retained tree count, shared
taxon contract, consensus rule, clade-support dispersion, topology
multimodality, and instability warnings that remain attached to the tree
collection instead of collapsing the uncertainty review into a consensus tree
alone.

Reviewer-facing analysis HTML reports now keep one explicit limitations
section in the durable report body and machine manifest. That contract applies
to the governed Bayesian inference review surfaces, `ancestral report`,
`comparative report`, `diversification report`, and `tree-set report`, so
warnings, invalid interpretations, exclusion boundaries, and unresolved model
or topology uncertainty remain visible instead of being implied only by side
tables or caveat prose elsewhere in the package.

For one reviewer-facing audit checklist over a supported package manifest, use
`report reviewer-audit-checklist`. It writes one TSV with pass, risk, and
blocked items plus explicit evidence and artifact pointers so tree,
comparative, ancestral, and alignment review packages expose one durable
checklist for reproducibility and validity review rather than scattering that
burden across HTML sections, methods text, and raw manifests alone.

For one reviewer-facing supplementary clade-support ledger, use
`report supplementary-clade-support-table`. It writes one TSV with one row per
reference-tree clade, keeping descendant taxa, direct tree-label support,
support class, source-tree ownership, and optional tree-set clade frequency
mapping together on one reviewable surface.

For one reviewer-facing supplementary model-selection ledger, use
`report supplementary-model-selection-table`. It writes one TSV with one row
per candidate model from parsed IQ-TREE artifacts, keeping likelihood,
parameter count, AIC, AICc, BIC, winner flags, and the selected model decision
together on one reviewable surface.

For one reviewer-facing supplementary comparative-model ledger, use
`report supplementary-comparative-model-table`. It writes one TSV with one row
per fitted coefficient across the declared comparative candidate models,
keeping the shared exclusion ledger, coefficient uncertainty, model ranking,
and fitted diagnostics together on one reviewable surface.

For one reviewer-facing supplementary ancestral-state ledger, use
`report supplementary-ancestral-state-table`. It writes one TSV with one row
per internal node across one continuous or discrete ancestral reconstruction,
keeping descendant taxa, model settings, uncertainty, instability flags, and
the shared warning and exclusion context together on one reviewable surface.

For one reviewer-facing supplementary batch ledger, use
`report supplementary-batch-summary-table`. It reads one written workflow
bundle and writes one TSV with one dataset row plus one row per variant,
keeping task status, output freshness, recovery action, merge state,
reproducibility outcome, output references, and batch warnings together on
one reviewable surface.

For one reviewer-facing supplementary diversification ledger, use
`report supplementary-diversification-table`. It writes one TSV with one row
per scanned clade, repeating the tree-wide Yule and birth-death estimates and
the actual sampling-correction evidence so clade-rate outliers, model ranking,
and sampling warnings stay on one reviewable surface.

For one reviewer-facing supplementary taxon ledger, use
`report supplementary-taxon-table`. It writes one TSV that keeps taxon IDs
across the tree, alignment, metadata, traits, tip dates, and calibrations
together with prefixed metadata and trait values, explicit dataset inclusion
status, dataset exclusion causes, workflow reporting loss, and the exact
reason later workflow surfaces dropped a taxon from the final reviewer-facing
set.

For a journal-oriented annotated trait tree package, use
`report trait-tree-package`. It builds on the figure package and adds one
review HTML surface, one annotation-coverage ledger, one annotation-surface
summary ledger, and one dedicated package manifest. Publication readiness stays
blocked until every requested trait, metadata strip, and heatmap surface covers
the tree taxa completely and the underlying caption, legend, and legibility
audits all pass.

For a journal-oriented ancestral-state figure package with explicit uncertainty,
use `ancestral package`. It writes one SVG, PNG, and figure HTML surface plus
one reviewer HTML audit, one node-uncertainty review ledger, and one manifest
that blocks publication readiness until internal node states and their
uncertainty remain visible on the figure itself. Continuous packages render
estimate `+/-` standard-error labels on internal nodes, while discrete packages
pair internal probability pies with top-state confidence labels so ambiguity is
interpretable without leaving the figure surface.

For a journal-oriented discrete biogeography figure bundle, use
`biogeography report`. The package now keeps the ancestral-region tree, the
geographic transition map, a machine-readable `figure-legend.tsv`, and a
structured `figure-caption.md` together under one manifest-backed audit.
Publication readiness stays blocked until the tree renders one probability pie
and one probability label per internal node, the map retains visible
transition-support lines, and every inferred region remains represented through
the shared tree-and-map state palette.

For a journal-oriented dated phylogeny package with visible uncertainty, use
`report time-tree-package`. It materializes a retained posterior tree set,
builds one maximum clade credibility time tree, renders one `time-tree.svg`
with median node ages and 95% HPD whiskers on every internal node, writes one
`node-age-intervals.tsv` ledger and one reviewer HTML summary, and keeps
publication readiness blocked whenever intervals are missing, the dated tree is
not ultrametric, or supplied tip-date or calibration evidence fails the
time-tree readiness audit.

For a journal-oriented Brownian-versus-OU comparative review, use
`comparative model-comparison-package`. It builds one explicit
`model-comparison-criteria.svg`, `model-comparison-likelihood.svg`,
`model-comparison-parameters.svg`, and `model-comparison-fit-summary.svg`
surface, pairs them with machine-readable criteria, likelihood, parameter, and
fit ledgers, and writes one reviewer HTML summary plus one manifest-backed
audit. Publication readiness stays blocked whenever one candidate model loses a
finite AICc value or the selected model remains within an AICc delta of `2.0`
from the runner-up, but the full package still renders so ambiguous support is
reviewed honestly instead of being collapsed into a forced winner.

For a journal-oriented tree-set uncertainty figure bundle, use
`tree-set package`. It builds one consensus tree plus one explicit
`consensus-tree.svg`, `clade-support-plot.svg`, `unstable-taxa-plot.svg`, and
`topology-clusters-plot.svg` surface, then pairs those figures with
`figure-legend.tsv`, `figure-caption.md`,
`tree-set-uncertainty-methods-summary.md`, `uncertainty-review.html`, and a
manifest-backed publication audit. Publication readiness stays blocked until
the consensus support labels are rendered through the support-scale audit, the
instability panel remains explicit even when no unstable taxa are detected, and
the topology-cluster panel keeps alternative rooted modes visible instead of
burying them in tables alone.

For a journal-oriented alignment quality figure bundle, use
`report alignment-package`. It builds one explicit
`alignment-missingness-heatmap.svg`, `alignment-site-quality-summary.svg`, and
`alignment-sequence-quality-panel.svg` surface, pairs them with machine-readable
heatmap, window, and ranking ledgers, and writes one reviewer HTML summary plus
one manifest-backed audit. Publication readiness stays blocked when the
alignment remains suspicious, contains invalid characters, or falls below the
review quality threshold, but the full figure package still renders so missing
or ambiguity-heavy alignments can be reviewed honestly instead of being hidden
behind a boolean failure.

`demo rabies-cross-host-geography-panel` is the repository's flagship public
biological workflow surface. In addition to the dataset and workflow
subdirectories, the demo now writes one reviewer-facing package layer at the
output root:

- `dataset/source-accessions.tsv` for machine-readable accession provenance
- `rabies-cross-host-geography-overview.html` for one direct public handoff
- `rabies-cross-host-geography-artifacts.tsv` for one checksum-backed inventory
  over the exported study inputs, workflow outputs, and overview surfaces
- `rabies-cross-host-geography-reproducibility-checklist.tsv` for one
  machine-produced pass/risk/blocked reviewer checklist over the study package
- `rabies-cross-host-geography-package.manifest.json` for the biological
  question, short answer, config provenance, package-control artifacts, output
  checksums, and key metrics

Its packaged `workflow-config.json` now also carries explicit resource-budget
controls alongside the scientific settings: `iqtree_threads`,
`timeout_seconds`, `max_bootstrap_tree_count`, `max_report_table_rows`, and
`memory_warning_threshold_bytes`. The runtime records observed workflow and
bootstrap-review runtime or memory where available and emits structured budget
failures or warnings instead of silently overrunning those limits.

Its JSON metrics also surface the same `biological_question` and
`short_answer` directly so reviewers do not need to open the nested HTML report
to understand the intended scientific claim.

For a stored flagship package, use `report package-revalidation` on
`rabies-cross-host-geography-package.manifest.json`. It rereads the stored
manifest, verifies that the package artifact inventory and reproducibility
checklist still match their recorded checksums and row counts, rechecks every
inventory-listed dataset or workflow file against its stored checksum and size,
and emits:

- `publication-package-revalidation-artifacts.tsv` for one per-artifact
  existence and checksum ledger
- `publication-package-revalidation-checks.tsv` for one reviewer-facing
  `pass` / `risk` / `blocked` decision table
- `publication-package-revalidation-summary.json` for one machine-readable
  verdict surface
- `publication-package-revalidation-report.html` for one direct reviewer handoff

`all_original_artifacts_match=true` means the original study inputs and outputs
still match the stored package contract. An overall `risk` result means the
original package still matches but the package root now contains undeclared
extra files. An overall `blocked` result means some declared artifact,
inventory, checklist, or manifest-declared checksum no longer matches.

For two stored flagship package versions, use `report package-comparison` on
their two `rabies-cross-host-geography-package.manifest.json` files. It
compares artifact inventory parity, governed workflow settings, input accession
and sequence drift, alignment surfaces, rooted-tree structure, model-selection
surfaces, reviewer-facing figures and reports, and the stored scientific
findings or short-answer summary, then emits:

- `publication-package-comparison-artifacts.tsv` for one per-artifact
  `same` / `changed` / `left_only` / `right_only` ledger
- `publication-package-comparison-checks.tsv` for one reviewer-facing drift
  decision table across inputs, config, alignments, trees, models, and
  conclusions
- `publication-package-comparison-summary.json` for one machine-readable
  comparison verdict surface
- `publication-package-comparison-report.html` for one direct reviewer handoff

An overall `pass` result means the stored study package versions remain
scientifically and operationally aligned across those governed surfaces. An
overall `risk` result means at least one governed study surface drifted between
the two versions. An overall `blocked` result means the two manifests no longer
describe the same governed study package contract at all.

The workflow bundle now also writes one `workflow/conclusion-stability/`
surface that scores:

- key clade stability across rooted bootstrap trees and method variants
- support-value stability across rooted bootstrap trees and method variants
- ancestral-state stability across rooted bootstrap trees and method variants
- comparative coefficient stability across topology bootstrap trees and method variants

The reviewer-facing HTML report separates `stable`, `weak`, and `unstable`
conclusions directly, and the CLI JSON metrics report
`conclusion_stable_count`, `conclusion_weak_count`, and
`conclusion_unstable_count`.

`demo rabies-method-sensitivity-panel` is the repository's governed public
workflow for method-sensitivity review on the same compact rabies panel. It
reruns four declared alignment-and-trimming variants, compares IQ-TREE and
FastTree on each variant, roots both engine trees on the packaged outgroup,
and now executes those isolated variant roots through one declared
`parallel_workers` budget from the packaged config unless the CLI overrides it.
The resulting bundle keeps one workflow manifest, one parallel-execution
summary, and one task log per variant so reviewers can see which isolated
outputs were executed together without mixing their logs. The workflow summary
records the variant count, stable-clade count, changed-clade count, rooted
preprocessing change count, rooted engine change count, and the number of
variants that still show serious unrooted engine conflicts before rooting.
Concurrent reuse of the same workflow output root is rejected explicitly while
one run is active, and the raw workflow root now writes one
`rabies-method-sensitivity-panel.run.json` execution record so parallel worker
count, execution mode, successful variants, failed variants, and per-variant
task logs remain auditable even when one isolated task fails.
The bundle now also writes one reproducibility audit over the current rabies
inputs, resolved settings, workflow manifest, report manifest, task logs, and
per-variant output directories. That audit ships as
`reproducibility-checks.tsv`, `reproducibility-variants.tsv`, and
`reproducibility-audit.json` so reviewers can verify that the current batch
outputs still correspond to the declared inputs and settings.
It now also writes one Slurm-ready planning surface over the same declared
variants: `slurm-job-plan.tsv` records one estimated job per variant with
suggested `sbatch` options, `slurm-estimation-assumptions.tsv` makes the
resource-sizing rules explicit, and `slurm-planning-summary.json` carries the
same planning contract for machine review. Those estimates are grounded in the
current workflow footprint instead of placeholder job sizes, so the bundle
surfaces the planned job count, total estimated core-hours, maximum estimated
memory, maximum wallclock, and total scratch/output estimates directly.
On top of that plan, the bundle now writes one real job-array partitioning
strategy: `slurm-array-partitions.tsv` groups compatible jobs by dataset-size
class, method group, and resource class, `slurm-array-members.tsv` maps each
variant to its array index, `slurm-array-strategy.json` carries the structured
partition contract, and `workflow/slurm-arrays/*.sbatch` contains executable
per-partition scripts that call the repository CLI with one selected variant
per array task.
The bundle also writes one resumable workflow-status surface:
`slurm-job-status.tsv` classifies each planned job as completed, failed,
pending, or stale from the execution record, task logs, running marker, and
durable variant outputs; `slurm-partition-status.tsv` rolls those counts up by
array partition; and `slurm-workflow-status.json` carries the same summary for
machine review.
It now also writes one retained-storage planning surface before scaling the
same workflow up further: `slurm-storage-categories.tsv` breaks the current
bundle into workflow outputs, canonical logs, tree artifacts, posterior
samples, and reviewer-facing reports; `slurm-storage-variants.tsv` records the
same retained footprint per variant; `slurm-storage-report.json` carries the
machine-readable totals; and `slurm-storage-report.html` gives reviewers one
compact summary of the expected retained storage burden. The posterior-sample
category remains explicit even when it is zero, so the workflow does not hide
the fact that this governed rabies surface currently emits no Bayesian chain
or posterior tree outputs.
The same bundle now also writes one output-explosion warning surface before the
workflow grows further: `slurm-output-explosion-checks.tsv` records the global
and per-variant checks behind the warning contract,
`slurm-output-explosion-variants.tsv` records each variant's risk class and
dominant retained categories, `slurm-output-explosion-report.json` carries the
machine-readable summary, and `slurm-output-explosion-report.html` gives
reviewers one compact explanation of whether current retained outputs, tree
artifacts, posterior samples, or report files are starting to scale badly.
Right after that, the same workflow now writes one tree-retention policy
surface for deciding when multi-tree artifacts should be thinned or
compressed safely: `slurm-tree-retention-checks.tsv` records the consistency
checks behind the policy, `slurm-tree-retention-files.tsv` records one
per-file retention decision, `slurm-tree-retention-policy.json` carries the
machine-readable summary, and `slurm-tree-retention-policy.html` gives
reviewers one compact explanation of whether any retained tree sets need
interval thinning or gzip compression. The current governed rabies bundle has
only single-tree outputs, so the policy explicitly records that thinning is
not currently applicable.
It now also writes one per-job evidence surface for independent debugging:
`slurm-job-evidence.tsv` indexes one provenance package per planned job,
`slurm-job-evidence-summary.json` carries the workflow-wide summary, and
`workflow/slurm-job-evidence/<variant_id>/` stores one local task log copy,
step manifests, a machine-readable evidence JSON, and a reviewer-facing HTML
summary for that specific job. That package is meant to be sufficient to debug
one failed or suspicious job without reopening the whole batch history first.
On top of that, the same bundle now writes one stale-output invalidation
surface: `slurm-output-freshness.tsv` records whether each planned job still
matches the current packaged inputs and output-affecting workflow settings,
`slurm-output-freshness-checks.tsv` exposes the exact checksum and setting
checks behind that verdict, and `slurm-output-freshness.json` carries the same
summary for machine review. That means changed input files or changed declared
settings no longer look like clean completed outputs: the freshness ledger and
the resumable status ledger both mark those results as stale until the batch is
rerun.
The same governed bundle now also writes one failure-recovery surface over
those status and freshness ledgers: `slurm-failure-recovery-jobs.tsv` records
one rerun decision per planned job, `slurm-failure-recovery-partitions.tsv`
rolls that up by array partition, `slurm-failure-recovery-report.json` carries
the machine-readable recovery summary, and
`slurm-failure-recovery-report.html` gives reviewers one compact explanation of
which jobs are rerunnable, which ones are merely blocked by a still-live
workflow, and which likely failure cause was inferred from the governed task
logs and status evidence.
Once those per-job ledgers exist, the bundle now also writes one global merge
surface over the complete distributed run:
`slurm-merge-checks.tsv` records the batch-level merge checks,
`slurm-merge-variants.tsv` records each variant's merge eligibility and linked
evidence package, `slurm-merge-report.json` carries the machine-readable merge
decision, and `slurm-merge-report.html` gives reviewers one compact summary of
whether the current completed outputs are actually ready to be merged into one
coherent workflow result.
Its reviewer-facing HTML report is now intentionally compact: it surfaces one
summary card set plus explicit links to the governed TSV and JSON ledgers
instead of embedding those tables directly, and the bundle now includes
`workflow/report-artifacts/rabies-method-sensitivity-report.manifest.json`
with linked-artifact checksums and byte counts. The JSON metrics for the demo
also report `report_linked_artifact_count`, `report_html_size_bytes`,
`report_linked_artifact_bytes`, `report_total_output_bytes`,
`slurm_job_count`, `slurm_total_estimated_core_hours`,
`slurm_maximum_estimated_memory_mib`,
`slurm_maximum_estimated_wallclock_minutes`,
`slurm_total_estimated_scratch_mib`, `slurm_total_estimated_output_mib`,
`slurm_array_partition_count`, `slurm_array_script_count`, and
`slurm_array_largest_partition_size`,
`slurm_job_evidence_file_count`,
`slurm_job_evidence_total_runtime_seconds`, and
`slurm_job_evidence_total_output_byte_count`,
`slurm_output_explosion_status`,
`slurm_output_explosion_global_issue_count`,
`slurm_output_explosion_warning_variant_count`, and
`slurm_output_explosion_high_risk_variant_count`,
`slurm_tree_retention_status`,
`slurm_tree_set_file_count`,
`slurm_tree_posterior_sample_file_count`,
`slurm_tree_thinning_recommended_file_count`,
`slurm_tree_thinning_required_file_count`,
`slurm_tree_compression_recommended_file_count`, and
`slurm_tree_compression_required_file_count`,
`slurm_output_freshness_check_count`,
`slurm_output_freshness_failed_check_count`,
`slurm_fresh_output_job_count`, and `slurm_stale_output_job_count`,
`slurm_failure_recovery_status`,
`slurm_failure_recovery_rerunnable_job_count`,
`slurm_failure_recovery_blocked_job_count`, and
`slurm_failure_recovery_partition_count`,
`reproducibility_passed`, `reproducibility_check_count`,
`reproducibility_failed_check_count`, and
`reproducibility_failed_variant_count`.

`tree-set report` now follows the same scaling contract. The HTML report keeps
top-level uncertainty summaries in the page body, writes large reviewer tables
to one sibling `*.artifacts/` directory as TSV, JSON, or Markdown, links those
artifacts explicitly, keeps one durable
`tree-set-uncertainty-methods-summary.md` beside the reviewer tables, records
report mode as either `full-review` or `scaled-summary`, and reports
`methods_summary_warning_count`, `html_size_bytes`, `linked_artifact_bytes`,
and `total_output_bytes`. Tree sets with `1,000+` trees stay reviewable by
switching the highest-cost supplemental sensitivity passes to linked note
artifacts instead of rerunning them inline.

`demo gnathostome-ortholog-protein-benchmark` is the repository's governed
public amino-acid workflow benchmark. It packages one small real ortholog
protein FASTA panel, reruns MAFFT, trimAl, and IQ-TREE over those amino-acid
inputs, writes the aligned matrix, trimmed matrix, tree, model table, support
table, log, manifest, and overview, and adds one explicit
`workflow/molecular-assumptions.tsv` ledger so reviewers can see that this
surface uses `-st AA`, searches protein models only, does not translate coding
DNA, and does not rely on nucleotide-specific assumptions such as codon
position or GC interpretation.

`demo known-answer-reference-panel` is the repository's governed internal
simulation truth suite. It exports one deterministic tree, alignment,
Brownian-trait, OU-trait, discrete-trait, host-association, and geographic
panel with stored node truths, branch-transition truths, and explicit recovery
thresholds, then reruns the owned recovery workflow and writes parameter,
internal-node, branch-event, and threshold-evaluation ledgers. Its JSON
metrics now surface not only topology and ancestral recovery but also
`parameter_row_count`, `threshold_pass_count`, `threshold_row_count`,
`host_event_accuracy`, and `geographic_event_accuracy`, so the trust contract
is grounded in known simulation answers instead of only one inferred summary.

`demo catarrhine-data-quality-stress-panel` is the repository's governed
dirty-data gauntlet. It combines one already aligned catarrhine matrix, one
raw FASTA validation surface with duplicate identifiers, illegal characters,
empty sequences, and length outliers, one raw coding FASTA surface with a
frame error and a premature stop codon, one raw trait-linkage mismatch
surface, and one rooted tree with zero, negative, and extreme branch lengths.
The demo writes explicit review ledgers for every defect class plus one cleaned
comparative subset, and its JSON metrics expose the same counts directly
through `duplicate_sequence_identifier_count`,
`raw_sequence_length_outlier_count`, `coding_frame_error_count`,
`coding_internal_stop_count`, `raw_trait_missing_from_traits_count`,
`raw_trait_extra_taxon_count`, and `tree_negative_branch_count`.

The same packaged panel is available as a typed Python surface under
`bijux_phylogenetics.datasets.data_quality_stress`. Use
`load_catarrhine_data_quality_stress_panel_dataset()` to inspect the packaged
dirty-input panel, `export_catarrhine_data_quality_stress_panel_dataset(...)`
to materialize that panel without rerunning the workflow,
`run_catarrhine_data_quality_stress_panel_workflow()` to rerun the governed
cleanup report in memory,
`write_catarrhine_data_quality_stress_panel_workflow_bundle(...)` to emit the
reviewer-facing ledgers, or
`run_catarrhine_data_quality_stress_panel_demo(...)` to materialize the
dataset copy plus the workflow outputs in one call.

`comparative signal` keeps its input policy explicit in JSON output: rooted
trees are accepted whether or not they are ultrametric, overlapping missing
trait values are pruned and reported, permutation rows are reproducible from
the supplied seed, and constant post-pruning trait vectors fail with a typed
comparative-method error instead of returning misleading scalar summaries.

`comparative phylogenetic-residuals` is the owned review surface for one
continuous response and one continuous predictor when reviewers need
tree-aware fitted values and residuals instead of one full multivariate PGLS
report. The command can either fix Brownian covariance or estimate Pagel's
lambda before fitting, writes summary, taxon, coefficient, and exclusion
ledgers, preserves one explicit missing-data audit, and reports large
standardized residuals without reducing the review to one opaque outlier flag
alone. `comparative phylogenetic-anova` is the matching categorical review
surface for one continuous response and one group column; it runs one seeded
Brownian null simulation, emits group, pairwise, simulation, and exclusion
ledgers, keeps unequal group sizes explicit, and warns when one or more
groups have fewer than three taxa after pruning.

`comparative multivariate` is the owned shared-taxon review surface for
multiple continuous responses under one predictor design. It writes
response-model, coefficient, residual-covariance, residual-correlation,
residual-association, and exclusion ledgers, keeps the shared complete-case
missing-value policy explicit, and warns when residual degrees of freedom are
weak or the residual covariance matrix is singular. When estimated
response-level Pagel lambda fits diverge materially, the report now says so
explicitly, because the shared residual covariance and correlation then
compare residuals fit under different phylogenetic error assumptions rather
than one common lambda surface.

`comparative covariance-audit` is the pre-fit review surface for Brownian
trait models, OU trait models, and PGLS covariance choices when reviewers
need to know whether the tree and trait overlap can support comparative
fitting before they trust any coefficient or likelihood. It reports matched
and missing taxa, duplicate taxon keys, zero-length and negative branch
counts, matrix dimension, matrix rank, condition number, and the governed
solve path. `fit_strategy=exact` means the raw covariance is already positive
definite and well-conditioned enough for direct inversion. `fit_strategy=`
`regularization` means the raw covariance is singular or ill-conditioned and
the governed diagonal stabilization path would be required before inversion.
`fit_strategy=failure` means blockers or unrecoverable covariance invalidity
stop fitting before any model coefficient should be interpreted.

`comparative discrete-mk` is the governed standalone discrete trait-evolution
fit surface for one rooted tree and one categorical tip trait when reviewers
need the fitted ER, SYM, or ARD likelihood surface directly rather than one
ancestral reconstruction report. For the ER baseline, the owned runtime now
matches the governed live `phytools::fitMk(model='ER')` lane over binary,
multistate, and missing-value-pruned fixtures, and the SYM surface now also
matches the governed live `phytools::fitMk(model='SYM')` lane over multistate
clean and missing-value-pruned fixtures. The ARD surface now also matches the
governed live `phytools::fitMk(model='ARD')` lane over binary and
missing-value-pruned binary fixtures at full rate-row parity, and over clean
and missing-value-pruned multistate fixtures at summary parity when the
optimizer flags weakly identified boundary rates. The same owned surface now
also supports `--transform lambda` for unordered ER, SYM, and ARD fits on
ultrametric rooted trees, records the fitted Pagel-lambda value separately
from the transition-rate parameters, and warns when the transformed
branch-length support stays close to the zero-signal or untransformed
boundary. The command reports log-likelihood, parameter count, AIC, AICc, one
explicit missing-value pruning audit, one directed transition-rate ledger,
overparameterization status, weak-identifiability warnings, optimizer
diagnostics, and ER baseline comparison metrics instead of reducing the fit
to one scalar alone. The same owned runtime is also reusable directly from
Python through
`fit_discrete_mk_model_from_dataset(...)` once one `AncestralDiscreteDataset`
has already been loaded. Inside the comparative package, that discrete Mk
surface now lives as a domain package with separate `models`, `transforms`,
`fitting`, `comparison`, and `tables` modules so transform search policy,
single-fit execution, model-ranking workflow, and ledger writing do not
collapse back into one flat implementation file.

The PGLS surface is likewise available as a typed Python package under
`bijux_phylogenetics.comparative.pgls`. Use `inspect_pgls_inputs()` and
`build_pgls_model_matrix()` for formula audit and encoded design-matrix
inspection, `run_pgls()` for one comparative regression fit, and
`run_pgls_multiple_testing()` for family-wise coefficient review across many
response traits. Inside that package, formula parsing, input/design audit,
numerical fitting, and multiple-testing policy now live in separate
`formula`, `design`, `fitting`, and `multiple_testing` modules, with shared
dataclasses in `models`, so the package root stays a curated API surface
instead of regressing into one flat implementation file.

The BEAST adapter surface now makes its evidence state explicit. `adapter
beast-prepare` only prepares XML, `adapter beast-log`, `beast-trees`, and
`beast-consensus` parse existing posterior outputs, and reviewer-facing
diagnostics such as `adapter bayesian-methods` state when they are only
summarizing a prepared XML versus parsing existing log/tree outputs. When a
matching `analysis.manifest.json` from `adapter beast-run` is present, those
diagnostics identify the posterior log and tree set as outputs from a recorded
prior BEAST inference rather than implying the report executed BEAST itself.

The MrBayes surface now follows that same owned package shape under
`bijux_phylogenetics.bayesian.mrbayes`. Preparation policy lives in
`preparation`, engine orchestration in `execution`, trace and MCMC table
parsing in `tabular`, posterior-tree parsing and consensus summarization in
`posterior_trees`, and ESS, convergence, parameter-summary, and burn-in
sensitivity reporting in `diagnostics`, with shared dataclasses in `models`
and structured artifact errors in `artifacts`. The package root is now only a
curated re-export surface instead of another long mixed-responsibility module.

The Bayesian runtime controls are intentionally strict and now aligned across
BEAST and MrBayes. `adapter beast-run` and `adapter mrbayes-run` leave an
explicit `.incomplete.json` marker not only for timeouts and nonzero exits but
also when the engine exits yet the emitted posterior files fail validation.
That marker records the failure reason plus the observed posterior output state
so reviewers can see whether the run left missing files, empty files, or other
partial artifacts before deciding to clean or rerun.
That shared safety path now lives in
`bijux_phylogenetics.bayesian.posterior_execution`, so both posterior engines
follow the same resume, reject, clean, and output-validation sequence instead
of maintaining separate orchestration code.
`--resume` reuses only one verified completed manifest from the same command,
same checked inputs, and same recorded engine version, `--incomplete-run-policy
clean` is the governed way to discard that partial state, and a missing
executable stops before any incomplete-run marker is written because no engine
run started.

That strictness now applies across every governed external-engine workflow.
MAFFT and trimAl runs succeed only when they emit non-empty valid alignments,
IQ-TREE runs succeed only when the required `.iqtree`, `.log`, tree, model,
and support artifacts are present and parseable for the selected workflow,
FastTree succeeds only when the tree and local-support annotations can be
parsed, and BEAST or MrBayes succeed only when the expected posterior artifact
set exists. Missing required files, empty required files, missing model
results, and missing required support annotations now surface stable structured
errors such as `engine_required_output_missing`, `engine_output_empty`,
`engine_model_result_missing`, and `engine_support_values_missing` before any
workflow manifest or reviewer-facing report is written.

The same reviewer-facing engine evidence now has one product-owned matrix
surface under `bijux_phylogenetics.engines.validation`. Use
`run_alignment_engine_validation_matrix()` for MAFFT, trimAl, IQ-TREE, and
FastTree, `run_bayesian_engine_validation_matrix()` for MrBayes plus live or
governed-fallback BEAST validation, and `run_external_engine_validation_matrix()`
to merge both families into one JSON report written by
`write_external_engine_validation_matrix()`.

The external-engine and release-install trust surface now has three distinct
verification lanes. Fast `engine_contract` tests keep fake-executable and
parser behavior stable in routine verification. `tests/real_local` carries the
governed `engine_real` lane for installed MAFFT, trimAl, IQ-TREE, FastTree,
and MrBayes executables plus the checked-in real BEAST XML/log/tree corpus.
That same `real_local` suite also carries one installability smoke lane that
builds wheel and sdist artifacts, installs each into a clean virtual
environment, copies the packaged example inputs through the installed runtime
API, runs core CLI commands against those writable copies, and verifies that
packaged resources are present in the built distributions. The engine-real
tests still write one combined external validation matrix JSON artifact over
all governed engines, while preserving focused alignment and Bayesian matrix
artifacts for narrower debugging. Each matrix records reviewer-facing engine
names, validation modes, executable paths, version text, commands, exit codes,
runtime, output paths, and output hashes.

Use `phylo preflight` before any external-engine workflow when you need to know
whether the local environment is actually runnable. The command inspects MAFFT,
trimAl, IQ-TREE, FastTree, MrBayes, and BEAST, reports the resolved executable
path and detected version for each engine, classifies each engine as
`tested`, `untested`, `unsupported`, or `missing`, and then summarizes which
governed workflows are `ready`, `caution`, or `blocked`.

When you pass `--workflow`, the command becomes a real gate instead of a passive
inventory. A blocked selected workflow exits early, while a runnable one returns
both `selected_workflow_status` and `overall_status` in JSON output. That split
is intentional: `selected_workflow_status` answers whether the chosen workflow
can run now, while `overall_status` still reflects the health of the broader
external-engine environment.

That diagnostic surface is not optional glue anymore. The governed direct run
surfaces such as `adapter align`, `adapter trim`, `adapter model-select`,
`adapter infer-ml`, `adapter bootstrap`, `adapter infer-fast`,
`adapter fasta-to-tree`, `adapter beast-run`, `adapter mrbayes-run`, and the
Python workflow APIs now run the same compatibility gate before they start
writing workflow outputs. Missing executables still raise the engine-specific
availability error on single-engine surfaces, while unsupported versions or
multi-engine blockers fail as one explicit preflight workflow block instead of
surfacing halfway through a run.

Every governed external-engine workflow now writes a durable manifest that
captures the workflow identifier, input checksums, structured config, resolved
engine commands, detected engine versions, seeds, runtime, and output
checksums. Use that manifest as the review anchor for provenance, reruns, and
downstream evidence bundles instead of reconstructing those details from logs.

Use `phylo run` when you want one serious workflow to start from a single YAML
or JSON config instead of stitching together CLI flags manually.

```bash
bijux-phylogenetics phylo run workflow-config.yaml --json
```

The current config-driven surface targets the governed `fasta-to-tree`
workflow. One config file can declare the input FASTA, optional metadata and
traits tables, engine executables, alignment and trimming settings, inference
seed and threads, output directory, optional result-bundle directory, and
timeout or incomplete-run controls. Invalid config files fail before engine
preflight or alignment starts. A valid run executes the same governed
tree-building workflow as `adapter fasta-to-tree`, then exports and validates
one complete result bundle carrying the resolved workflow config plus copied
config source, metadata, and traits files alongside the inference outputs and
engine artifacts.

When one of those governed workflows fails under `--json`, the error payload
now stays scientific instead of stopping at a bare Python exception.
Structured failures expose `failure_reason`, `scientific_explanation`,
`likely_causes`, `actionable_fixes`, and `evidence`. Invalid FASTA inputs name
duplicate or empty records and illegal characters, trimming failures
distinguish empty retained alignments from missing artifacts, tree-inference
failures distinguish missing trees from unparsable tree outputs, comparative
taxon-linkage failures list the missing or extra taxa explicitly, and BEAST or
MrBayes parser failures identify the missing file, header, sampled row, or
tree block section that prevented scientific review.

Use `phylo replay` when you need to rerun one governed manifest and verify that
the new outputs are still scientifically equivalent.

```bash
bijux-phylogenetics phylo replay \
  artifacts/fasta-to-tree/example.manifest.json \
  --out-dir artifacts/fasta-to-tree-replay \
  --json
```

Replay refuses to run when any recorded input checksum has changed, records
engine-version drift when the local executable no longer matches the manifest,
and compares replayed outputs with tolerant workflow-specific checks instead of
requiring byte-identical files. Tree workflows compare topology and support
semantically, while alignment and model-selection workflows compare the durable
scientific result for that workflow surface. When the manifest already records
one concrete engine path, replay reuses that path by default instead of
guessing a generic `mafft`, `iqtree2`, `FastTree`, `mb`, or `beast` binary.
Pass explicit executable overrides only when you intentionally want to test one
different local engine install against the recorded workflow evidence.

Use `phylo bundle` when you need one portable handoff directory for review or
bundle-local reruns.

```bash
bijux-phylogenetics phylo bundle \
  artifacts/fasta-to-tree/example.manifest.json \
  --out-dir artifacts/fasta-to-tree-bundle \
  --json
bijux-phylogenetics phylo validate-bundle \
  artifacts/fasta-to-tree-bundle \
  --json
```

That bundle copies the workflow manifest, extracted config, bundle-local rerun
ledger, reviewer-facing HTML report, copied inputs when they still exist,
reviewer-facing outputs, and step-level engine artifacts plus step manifests.
`phylo validate-bundle` then checks both checksum integrity and workflow
completeness, including the presence of the required report, outputs, and
declared step manifests.

The rabies demonstration bundle now publishes one governed reproducibility and
review layer alongside the biological outputs: `workflow-config-audit.tsv`,
`workflow-config.resolved.json`, one rooted-ML-versus-bootstrap-consensus
comparison under `bootstrap-review/`, and `scientific-findings.tsv`. The demo
JSON payload surfaces the config-check count, rooted consensus RF distance, and
scientific-finding count so reviewers can quickly tell whether the run stayed
within the expected contract.

For aligned multi-locus datasets, `alignment concatenate` is now the canonical
supermatrix assembly surface. It preserves taxon identifiers across loci,
inserts `?` blocks for absent taxa, writes a remapped partition file, and can
materialize the taxon-by-locus occupancy matrix in the same run. When a locus
contains residues that are alphabet-ambiguous across DNA and protein codes, use
repeated `--data-type` flags so the partition file declares the intended
datatype honestly instead of guessing from overlapping symbols alone.

Run `alignment occupancy` before inference when the supermatrix contains partial
fragments, short recovered loci, or uneven locus completeness. By default, a
locus counts as covered when at least one non-missing site is present. Raise
`--minimum-locus-occupancy` when thin fragments should count as absent for
taxon/locus thresholding instead. The occupancy TSV outputs include both
binary-coverage summaries and `site_coverage_fraction` columns so reviewers can
see the difference between locus presence and overall retained signal.

The `adapter fasta-to-tree` workflow is the canonical one-command bridge from
raw FASTA to a supported inference bundle. It accepts DNA and protein FASTA
inputs, runs alignment, trimming, model selection, maximum-likelihood
inference, and bootstrap support estimation, then writes:

- `prefix.aln`
- `prefix.trimmed.aln`
- `prefix.tree`
- `prefix.log`
- `prefix.model.tsv`
- `prefix.support.tsv`
- `prefix.manifest.json`

Engine-specific intermediate manifests and working files stay under
`out-dir/engine-artifacts/prefix/` so the final output set remains compact
without hiding reviewable run details.

The workflow now defaults IQ-TREE to `--iqtree-seed 1` and `--iqtree-threads 1`
so the checked output bundle is reproducible across reruns. Ultrafast bootstrap
support remains the supported branch-support backend here, which means
`--bootstrap-replicates` must be at least `1000`.

When `--resume` is enabled, the workflow reuses only stage outputs whose
manifest, input checksums, command, and detected engine version still match the
current run. The composite `prefix.manifest.json` now records one
`stage_fingerprints` ledger for raw-input validation, alignment, trimming,
model selection, inference, support, and final reporting so reviewers can see
which downstream stages were invalidated by changed inputs, changed bootstrap
settings, or changed engine binaries.

That same composite manifest is now the canonical bundle-export anchor for
portable review. `phylo bundle` copies the raw input when it is still present,
extracts the workflow config into one stable JSON file, keeps the reviewer
outputs together with the step-level engine artifacts, and writes one
reviewer-facing HTML report plus a bundle-local rerun ledger. The companion
`phylo validate-bundle` command fails when that handoff surface is incomplete.

Named MAFFT strategies are now first-class on both `adapter align` and
`adapter fasta-to-tree`. Use `--mode` or `--alignment-mode` with one of
`auto`, `linsi`, `ginsi`, `einsi`, or `fast`. The runtime expands each named
strategy into the explicit MAFFT arguments it runs, and the workflow manifest
captures both the resolved command and the detected MAFFT version for review.

The direct IQ-TREE adapter surface now preserves the native engine artifacts
that reviewers expect instead of collapsing everything into one opaque tree
file. `adapter model-select` retains `.iqtree`, `.log`, the native model
sidecar, and a generated `.model-candidates.tsv`; `adapter infer-ml` retains
`.treefile`, `.iqtree`, and `.log`; `adapter bootstrap` retains `.treefile`,
`.iqtree`, `.log`, `.ufboot`, and `.contree` when IQ-TREE emits them; and
`adapter consensus` retains the consensus `.contree` plus the matching
`.iqtree` and `.log`. The JSON report and persisted manifest for those runs
also expose the parsed `selected_model`, `selected_criterion`,
`candidate_model_count`, `best_model_aic`, `best_model_aicc`, `best_model_bic`,
`log_likelihood`, and support-value counts so downstream review does not have
to scrape IQ-TREE text again by hand.

Coding-sequence alignment is now first-class on `adapter align` through
`--codon-aware`. That workflow accepts raw coding DNA or RNA FASTA, excludes
frame-broken sequences and sequences with premature stop codons, aligns a
translated amino-acid guide with MAFFT, and back-translates guide gaps into
codon triplets so the final nucleotide alignment preserves reading-frame
boundaries. The workflow also writes audit artifacts for the translated guide
input, the aligned guide, and the excluded-sequence ledger.

Named trimAl strategies are also first-class on `adapter trim` and
`adapter fasta-to-tree`. Use `--mode` or `--trimming-mode` with one of
`gap-threshold`, `gappyout`, `strict`, `strictplus`, or `automated1`. The
workflow report and manifest record retained sites, removed sites, and
gap percentage before and after trimming alongside the explicit trimAl command
and detected version.

Raw sequence hygiene is now explicit. Use `alignment sequence-type` when you
need the raw FASTA type decision itself: compatible types, selected default,
confidence, and mixed or invalid blocking signals. Use `alignment validate-input`
to inspect duplicate identifiers, illegal sequence characters, empty records,
sequence-length outliers, and the same sequence-type report in one payload.
Use `alignment repair-input` when you want the runtime to normalize identifiers
or remove invalid records into a new FASTA. The same repair controls are
available on `adapter fasta-to-tree`; without them, the workflow now fails fast
instead of silently continuing on bad raw input. Mixed raw inputs must now
either be fixed or forced with an explicit `--sequence-type` choice before the
workflow continues.

The raw FASTA preflight surfaces now scan inputs linearly instead of
materializing whole record lists just to count duplicate identifiers, illegal
characters, empty records, and length outliers. That keeps
`alignment validate-input` practical on thousand-sequence inputs before any
external engine is invoked.

The checked real-dataset regression corpus for this end-to-end workflow now
lives under `packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree/`.
It currently pins reviewer-facing output bundles for:

- `gnathostome-ortholog-protein-benchmark`
- `gnathostome-ortholog-proteins`
- `gnathostome-ortholog-coding-sequences`
- `strnog-enog411bqtj-proteins`

Those governed workflow bundles are now compared by scientific equivalence
rather than by raw text alone. Tree outputs are checked by rooted clades plus
branch-length tolerance, support tables are checked by clade-aware numeric
tolerance, tabular ledgers are checked by schema plus stable row identity, HTML
reports are checked by headings and required linked artifacts, and manifest
payloads still carry exact output hashes when byte-identity matters.

For coding nucleotide phylogenetics, use `adapter align --codon-aware` first
and then run downstream inference steps such as `adapter model-select`,
`adapter infer-ml`, and `adapter bootstrap` on the codon alignment it writes.
The current `adapter fasta-to-tree` workflow remains a generic alignment and
trim pipeline, so it is not the codon-preserving entrypoint for coding DNA.

Internal alignment-quality scoring is available without any external engine.
Use `alignment quality` when you want one scored view that combines
per-sequence gap fractions, per-column gap fractions, invariant-site counts,
parsimony-informative-site counts, missing-data concentration, and a direct
suspicious-alignment verdict. The reviewer-facing alignment reports reuse the
same scoring contract so command-line JSON, HTML reports, and checked-in
reference fixtures stay aligned. The quality surface now reuses one loaded
alignment instead of re-reading the same matrix for each sub-diagnostic, and it
skips pairwise near-duplicate scans above the governed large-matrix threshold
with an explicit warning rather than silently doing quadratic work on
thousand-sequence alignments.

## Alignment Filter Profiles

The built-in alignment filtering profiles are `conservative`, `moderate`,
`aggressive`, `coding-safe`, and `phylogenomics-scale`.

Comparative workflows now include checked-in external reference validation and
a dedicated maturity audit for BM, OU, PGLS, residual diagnostics, and
leave-one-taxon-out sensitivity:

```bash
uv run bijux-phylogenetics parity --json
uv run bijux-phylogenetics parity --extended --json
uv run bijux-phylogenetics parity --reference-source ape-live --json
uv run bijux-phylogenetics parity --reference-source geiger-live --json
uv run bijux-phylogenetics parity --reference-source geiger-live --generated-report-out artifacts/geiger-parity-report.md --generated-report-json-out artifacts/geiger-parity-report.json --json
uv run bijux-phylogenetics comparative validate-reference --json
uv run bijux-phylogenetics comparative maturity tree.nwk traits.tsv --formula "height_cm ~ body_mass + habitat" --lambda-value 1.0 --json
uv run bijux-phylogenetics comparative pgls tree.nwk traits.tsv --formula "height_cm ~ log(body_mass) + habitat" --json
```

`parity` is the governed cross-surface reference suite. The default run stays
CI-sized, while `--extended` adds real primate comparative fits and the larger
posterior-tree bundle. The observation ledger now records expected failure
classification, overlap policy, and shared-versus-exclusive taxa so reviewers
can see whether a mismatch comes from topology, branch lengths, missing-taxa
handling, numerical tolerance, or model assumptions. The PGLS parity lane now
covers fixed-Brownian regression, treatment-coded categorical predictors,
treatment-coded interaction terms, and one governed estimated-lambda primate
regression against checked R `ape` plus `nlme` outputs for coefficients,
standard errors, p-values, likelihood, and AIC.

`parity --reference-source ape-live` is the live external parity harness. It
launches the checked-in R runner through `Rscript`, records the R version,
`ape` version, Bijux version, Bijux commit, function name, input fixture,
tolerance, pass or fail state, and mismatch reason for each governed case, and
writes reproducible failure or skip artifacts whenever the live lane disagrees
or `ape` is unavailable. The live observation table is structured rather than
string-based, so tree summaries, tip ledgers, normalized Newick outputs,
DNAbin state ledgers, DNA-state frequency tables, DNA-distance ledgers, and
translated amino-acid rows are compared as owned artifacts rather than scraped
console text. The governed live cases now cover `ape::read.tree`, `ape::write.tree`,
`ape::consensus`, `ape::prop.clades`, `ape::root`, `ape::unroot`, `ape::drop.tip`, `ape::keep.tip`, `ape::extract.clade`, `ape::getMRCA`, `ape::is.monophyletic`, `ape::cophenetic.phylo`, `ape::dist.topo`, `ape::vcv.phylo`, `ape::ace`, `ape::node.depth.edgelength`, `ape::branching.times`, `ape::is.ultrametric`, `ape::nj`, `ape::pic`, `ape::base.freq`, `ape::seg.sites`, `ape::dist.dna`, and `ape::trans` over shared tree, trait-table, and DNA
fixture ids. The lane now also covers governed `ape::rtree` and `ape::rcoal`
simulation-envelope cases over shared random-tree and coalescent simulation
fixtures. The tree, trait-table, DNA, and simulation inputs for that lane now come from
the governed shared fixture catalogs in
`tests/fixtures/metadata/shared_tree_fixture_catalog.json`,
`tests/fixtures/metadata/shared_trait_table_fixture_catalog.json`,
`tests/fixtures/metadata/shared_dna_alignment_fixture_catalog.json`, and
`tests/fixtures/metadata/shared_tree_simulation_fixture_catalog.json`, so Bijux
and `ape` resolve the same durable fixture identities instead of hand-picked
path lists. The `ape::read.tree` lane now checks structured clade rows rather
than only raw parse success, including rooted and unrooted trees, branch
lengths, internal node labels, support labels, quoted labels, one governed
multiple-tree input, and one governed malformed-Newick rejection case. Those
cases now flow through one owned native Newick parser and writer on top of
`PhyloTree` rather than an external tree reader. The
`ape::consensus` lane now covers majority-rule and strict consensus over
governed conflicting and posterior-style tree sets, writes one normalized
consensus Newick plus one clade-frequency TSV ledger per case, and fails
explicitly when the tree set does not share one exact taxon set. The
`ape::prop.clades` lane now covers reference-tree clade support mapping over
duplicate, reordered, posterior-style, and mismatched shared tree sets, and
the owned `tree-set support-map` surface writes one `reference-tree-support.tsv`
ledger keyed by descendant tip set rather than transient node index. The
runtime keeps the real `ape` edge case explicit too: unsupported
root-adjacent splits are left unscored instead of being mislabeled as zero
support.
The `ape::as.DNAbin` lane now covers clean, lowercase, gap-bearing, and
ambiguity-bearing DNA fixtures. On the owned Bijux side, DNA distance,
ape-style nucleotide composition, ape-style segregating-site review, and
aligned coding translation now all load through one DNAbin-compatible
nucleotide matrix that preserves taxon order and alignment length, normalizes
case, keeps gaps, ambiguity codes, and explicit missing states literal, writes
FASTA back without nucleotide-state loss, and rejects unsupported symbols
instead of silently degrading them to missing data. The same matrix is now
also reusable directly from Python through
`load_dna_bin_alignment(...)`,
`compute_alignment_base_frequency_report_from_dna_bin_alignment(...)`,
`compute_alignment_segregating_site_report_from_dna_bin_alignment(...)`, and
`compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(...)`. The
same owned matrix now also feeds
`inspect_coding_alignment_from_dna_bin_alignment(...)` and
`translate_coding_alignment_from_dna_bin_alignment(...)`, so review workflows
do not have to reparse one FASTA file separately for composition,
segregating-site, nucleotide-distance, aligned coding diagnostics, and
aligned translation inspection.

`parity --reference-source geiger-live` is the governed live `geiger`
execution harness for the owned continuous-mode fit surfaces. It launches the
checked-in `Rscript` runner through the same live parity contract as the `ape`
and `phytools` lanes, records the R version, `geiger` version, Bijux version,
Bijux commit, function name, fixture id, model name, invoked optimizer
settings, tolerance, pass or fail state, mismatch reason, and reproducible
artifact root for every governed case, and writes one summary TSV plus one
observation TSV instead of relying on console scraping. The initial registry is
intentionally narrow for this harness goal: it currently covers
`geiger::fitContinuous(model='BM')`,
`geiger::fitContinuous(model='white')`,
`geiger::fitContinuous(model='lambda')`,
`geiger::fitContinuous(model='kappa')`,
`geiger::fitContinuous(model='delta')`,
`geiger::fitContinuous(model='OU')`, and
`geiger::fitContinuous(model='EB')` over the four-taxon comparative smoke
fixture plus the governed shared `geiger` continuous-trait fixtures. Those
shared fixtures live in
`tests/fixtures/metadata/shared_geiger_continuous_fixture_catalog.json` and
cover twenty-four-taxon and one-hundred-twenty-eight-taxon ultrametric signal
surfaces, a rooted non-ultrametric control tree, OU and early-burst
known-truth traits, white-noise low-signal traits, missing-value pruning,
constant-trait blockers, an explicit outlier surface, one trend proxy for
later `fitContinuous` model expansion, and one governed per-taxon
standard-error review surface. Failed or skipped
cases always persist their case payload, structured summaries, parameter
ledgers, and mismatch reason under `artifacts/geiger-parity-failures/` so the
live parity lane stays reviewer-usable instead of collapsing into hidden
manual reruns. The same live report now also emits governed optimizer-disagreement
triage rows and an optional `--optimizer-triage-out` TSV. That triage keeps the
stored optimizer settings visible, records objective deltas, parameter-surface
agreement, boundary hits, trace row counts, local-optimum counts where the
owned or reference surface exposes a parameter trace, and classifies each case
as no algorithm mismatch, same-likelihood-different-parameters,
different-likelihood-same-parameters, boundary-solution review, or another
explicit review class instead of treating every failure as an owned bug or
waiving every mismatch away as optimizer noise.
The same report now also carries a governed parameterization registry plus an
optional `--parameterization-registry-out` TSV. That registry records the
reference and Bijux surface contracts, raw and canonical parameter names,
parameter conversions, bounds conversions, root-state and variance
parameterization notes, the current direct log-likelihood comparison policy, and
any expected divergence with explicit reproducible evidence. The current
accepted divergence is the continuous early-burst sign-and-bound convention:
raw `geiger` uses `a`, Bijux exposes `rate_change`, and the raw `geiger` bound
surface excludes exact zero even after sign conversion.
The same report now also carries a governed likelihood-constant policy table
plus an optional `--likelihood-policy-out` TSV. That table records whether a
case has one directly comparable raw `logLik` surface or only a ranked
row-level AIC/AICc comparison surface, checks that `AIC = 2k - 2logLik`
wherever a raw likelihood is supposed to be directly comparable, and keeps the
owned ranking guard explicit: `compare_fitcontinuous_model_ranking(...)` now
refuses to rank candidate modes if their likelihood-constant policies differ.
The same report now also carries a governed model-confidence table plus an
optional `--model-confidence-out` TSV. That table derives Akaike weights from
comparable delta AICc rows only, records delta AIC and delta AICc threshold
membership under one shared `2.0`-unit policy, keeps non-comparable models on
blank weight cells instead of pretending they received support, and preserves
the selected-model uncertainty class when the stored reference and owned
surfaces agree only up to confidence class rather than identical rounded
weight wording.
The same report now also carries a governed boundary-warning registry plus an
optional `--boundary-warning-out` TSV. That table records the affected
parameter, declared lower and upper bounds, explicit boundary hits,
near-boundary review, flat-likelihood-near-boundary review where trace data
exists, reviewer-facing warning kinds, and whether the reference and owned
surfaces still support a stable conclusion. Boundary-dominated cases such as
lambda at `0` or `1`, alpha near `0`, and EB, kappa, or delta fits that
collapse onto a limit surface now stay explicit instead of being narrated as
stable parameter recovery.
The same live run can now also generate one governed parity report with
`--generated-report-out` and `--generated-report-json-out`. That generated
report is not a hand-maintained checklist. It is built from the live geiger
parity observations plus governed artifact-backed summaries for the explicit
exclusion lanes, the `sim.char` envelope validator, the continuous and
discrete simulation-recovery panels, and the large-tree and real-dataset
macroevolution benchmarks. The benchmark rows are now loaded from packaged
governed benchmark bundles under `resources/benchmarks/` rather than rerunning
those slower benchmark surfaces inline during report generation. The generated
Markdown and JSON surfaces report the
goal-tranche inventory for goals `251` through `289`, live pass or fail or skip
counts, R and `geiger` versions, covered and excluded models, optimizer
mismatch categories, tolerance rules, model-boundary warning summaries, and
simulation-recovery counts without hiding boundary-dominated or review-only
surfaces behind blanket success language.
The owned Brownian and white-noise fits also now store that Gaussian
normalizing-constant policy directly in the fit result object, and the governed
tests pin both lanes against analytical known-answer likelihood calculations.
The BM lane now has three governed cases instead of one generic smoke check:
the four-taxon comparative example, a twenty-four-taxon sigma-recovery surface,
and a twenty-four-taxon missing-value pruning surface that must prove explicit
tip exclusion, missing-versus-nonnumeric classification, and the declared
no-standard-error policy before likelihood comparison is treated as valid.
The white lane now has three governed `fitContinuous(model='white')` cases:
a twenty-four-taxon strong-signal review surface that should fit worse than a
Brownian covariance model, a twenty-four-taxon weak-signal review surface where
independent variance is a plausible baseline, and a twenty-four-taxon
missing-value pruning review surface. Those cases govern the identity-covariance
owner fit, explicit no-phylogenetic-correlation warnings, closed-form
root-state and variance recovery, log-likelihood, AIC, AICc, and the declared
missing-value exclusion contract. The owned Bijux side exposes this through
`fit_continuous_evolutionary_mode(mode='white-noise')` so the live parity lane
compares a real repository-owned independent baseline instead of an ad hoc
reference-only shortcut.
Trend parity is explicitly excluded in this round. The shared fixture corpus
retains a directional trend-proxy surface, but the owned runtime does not claim
`fitContinuous(model='trend')` support because local `geiger` exposes distinct
`rate_trend` and `mean_trend` likelihood contracts and the shorthand `trend`
falls through to the rate-trend path by partial matching. Bijux now rejects
`trend`, `rate_trend`, and `mean_trend` requests explicitly instead of
pretending one interpretation is canonical.
Standard-error parity is also explicitly excluded in this round. Local
`geiger::fitContinuous` does expose the top-level `SE=` surface, so the
exclusion is explicit rather than silent: the owned Bijux continuous-mode fit
rejects `standard_error_trait`, the shared fixture corpus now retains one
governed per-taxon standard-error review surface for future parity work, and
the live `geiger` lane reports the durable policy string
`fitcontinuous-standard-error-explicitly-excluded-this-round`.
The lambda lane now has four governed `fitContinuous(model='lambda')` cases:
a twenty-four-taxon strong-signal review surface that lands on the Brownian
upper boundary, a twenty-four-taxon weak-signal review surface that collapses
to the zero-signal lower boundary, a twenty-four-taxon missing-value pruning
review surface, and a bounded-control review surface that keeps lambda inside
`[0.2, 0.6]` under an explicit reference `control` contract. Those cases
govern Pagel-lambda covariance transformation, explicit lambda bounds,
sigma-squared-backed rate, root-state recovery, log-likelihood, AIC, AICc,
reviewer-facing boundary warnings, and the durable mapping between live
`geiger` control settings and the owned Bijux search-control surface. The
owned Bijux side exposes this through
`fit_continuous_evolutionary_mode(mode='pagel-lambda')`, which follows the
`geiger::fitContinuous(model='lambda')` intercept-only likelihood contract
rather than reusing the repo's `phytools::phylosig`-style signal summary
surface. That distinction matters because the parity lane compares full
fit-likelihood, rate, and information-criterion behavior, not just one lambda
estimate.
The kappa lane now has three governed `fitContinuous(model='kappa')` cases:
a twenty-four-taxon strong-signal transformed-branch review surface, a
twenty-four-taxon weak-signal review surface that collapses to the lower
punctuational boundary, and a twenty-four-taxon missing-value pruning review
surface. Those cases govern direct branch-length power transformation,
explicit kappa bounds, sigma-squared-backed rate, root-state recovery,
log-likelihood, AIC, AICc, and reviewer-facing lower-boundary warnings. The
owned Bijux side exposes this through
`fit_continuous_evolutionary_mode(mode='pagel-kappa')` and
`rescale_tree_pagel_kappa(...)`, which follow the
`geiger::fitContinuous(model='kappa')` contract by transforming branch lengths
directly rather than by reusing the Pagel-lambda covariance scaling surface.
The delta lane now has three governed `fitContinuous(model='delta')` cases:
a twenty-four-taxon strong-signal temporal-concentration review surface, a
twenty-four-taxon weak-signal review surface that lands on the late-change
upper boundary, and a twenty-four-taxon missing-value pruning review surface.
Those cases govern transformed node-depth scaling, explicit delta bounds,
sigma-squared-backed rate, root-state recovery, log-likelihood, AIC, AICc,
and reviewer-facing upper-boundary warnings. The owned Bijux side exposes this
through `fit_continuous_evolutionary_mode(mode='pagel-delta')` and
`rescale_tree_pagel_delta(...)`, which follow the
`geiger::fitContinuous(model='delta')` contract by raising each node-depth
proportion to delta before recomputing branch lengths, rather than by applying
an edge-wise power transform like Pagel kappa.
The OU lane now has four governed cases with the same reviewer-facing rigor:
the twenty-four-taxon OU known-truth recovery surface, a missing-value pruning
review surface, a lower-boundary review surface on the rooted non-ultrametric
control tree, and a bounded-control review surface that keeps alpha inside
`[0.2, 1.0]` under explicit reference control settings. Those cases govern
alpha, sigma-squared-backed rate, optimum-or-root parameterization,
log-likelihood, AIC, AICc, declared parameter-bound policy, explicit
lower-boundary flags, and the durable control mapping between the live
reference lane and the owned Bijux search controls.
The Bijux side now records weak-identifiability and boundary warnings directly
from the owned continuous-mode fit report instead of hiding them inside
parity-only logic.
That owned report now also exposes a `boundary_assessment` for parameterized
continuous-mode fits, and `compare_fitcontinuous_model_ranking(...)` withholds
`stable_conclusion_supported` when the selected fit is dominated by a
boundary-hugging or flat boundary profile.
For parameterized continuous-mode fits, the owned side now also exposes
`ContinuousModeSearchControls` so bounded parity cases can declare coarse-grid
size, fine-grid size, and one explicit initial parameter value while the
governed diagnostics preserve the starting-value policy, first evaluated
likelihood, grid counts, and total evaluation count.
The EB lane now has two governed cases on the same public `rate_change`
surface: a twenty-four-taxon early-burst truth-recovery surface and a
Brownian-signal lower-boundary review surface. Bijux keeps a positive public
`rate_change` scale for the owned early-burst fit, while the live `geiger`
runner maps that surface to `fitContinuous(model='EB')` through the equivalent
negative `a` parameter and explicit bounded control. Those cases govern
sigma-squared-backed rate, root-state recovery, log-likelihood, AIC, AICc,
lower-boundary detection, and reviewer-facing weak-identifiability warnings
when the fitted rate change collapses back toward the Brownian boundary.
The same comparative owner package now also exposes a governed
`geiger::rescale` parity surface instead of treating branch transformation as
an internal side effect. Public wrappers now cover
`rescale_tree_pagel_lambda(...)`, `rescale_tree_pagel_kappa(...)`,
`rescale_tree_pagel_delta(...)`, `rescale_tree_early_burst(...)`, and
`rescale_tree_white_noise(...)`, and the shared transform engine is the same
one used by the owned continuous and discrete parity surfaces. The governed
reference contract compares transformed branch lengths and Brownian covariance
surfaces against local `geiger::rescale` on the same rooted example tree,
including the no-phylogeny white transform. The early-burst rescale contract
is explicit: the owned public wrapper keeps the positive `rate_change` scale,
while local `geiger::rescale(model='EB')` exposes the equivalent negative `a`
parameter, so parity compares those surfaces with the sign mapping stated
rather than hidden.
The live lane now also has a dedicated
`geiger::fitContinuous(model comparison)` surface over the governed
seven-model set `BM`, `white`, `lambda`, `kappa`, `delta`, `OU`, and `EB`.
Bijux exposes the owned side through `compare_fitcontinuous_model_ranking(...)`
while the live reference runner executes the same bounded model set per
fixture and records one full row per candidate model with AIC, AICc, deltas,
selected-model flags, and comparability notes. The governed comparison cases
cover Brownian truth, OU truth, early-burst truth, and white-noise truth on
the shared twenty-four-taxon fixtures. The parity contract intentionally
compares those model-comparison rows by model identity and numeric fit surface
rather than by derived whole-list ordering alone, because near-tied secondary
ordering can shift at machine precision even when the selected model and the
information-criterion surface still agree. The owned comparison report now also
records Akaike weights, lists the models that remain within the governed delta
threshold on both AIC and AICc, withholds weights from non-comparable models,
and emits explicit uncertainty language instead of reducing the result to one
bare selected-model name.
The same geiger-style fit surfaces now share one explicit tree and trait
alignment contract instead of rebuilding taxon overlap independently in each
loader. The owned side exposes that through
`align_tree_and_trait_table(...)`, which matches `geiger::treedata` on shared
taxon intersection, dropped tree-only tips, dropped trait-only rows, and final
tree-tip ordering, while keeping missing-value pruning as a second explicit
policy layered on top of the shared overlap report. Comparative readiness,
ancestral dataset loading, and discrete Mk parity summaries now all consume
that same alignment evidence so dropped-tip and dropped-row reporting stays
consistent across the owned and governed surfaces.
The same shared overlap policy now also exposes a first-class
`geiger::name.check` parity surface through
`check_tree_and_trait_taxon_names(...)` and the `traits name-check` command.
That contract keeps exact case-sensitive taxon matching, reports
`tree_not_data` and `data_not_tree` counts separately, rejects duplicate taxa
through the shared table loader, returns a clean `OK`-style outcome when the
tree and table agree, and can write a machine-readable mismatch table for
review. Dataset readiness preflight now consumes that same owned name-check
report so reviewer-facing blockers and warnings come from one durable mismatch
surface instead of ad hoc overlap logic.
The same comparative layer now also exposes a governed
`geiger::dtt` parity surface through
`summarize_disparity_through_time(...)` and the `comparative dtt` command.
That owned path computes geiger-style continuous clade disparity from the mean
squared Euclidean distance across all tip-trait pairs inside each internal
clade, carries those clade disparities into the same branching-time
aggregation rule that geiger DTT uses, supports one-column and multivariate
continuous trait matrices, writes curve, clade, exclusion, and optional
time-bin ledgers, and renders an SVG review figure. The current parity scope
is the observed DTT curve itself; simulated null envelopes and MDI-style
simulation summaries are still intentionally out of scope for this round.
That same owned clade kernel is also available directly through
`summarize_continuous_clade_disparity(...)` and the `comparative disparity`
command. The direct disparity surface reports one internal-clade disparity
row per ape-style internal node, keeps two-tip and larger clades on the same
formula, supports univariate and multivariate continuous matrices, writes a
dedicated summary ledger with the explicit average-squared-Euclidean method
formula, and is governed against local `geiger::disparity` references on the
simple four-tip clade partitions as well as the larger ultrametric DTT
fixtures. DTT now consumes that same owned clade-disparity surface rather
than reimplementing internal clade disparity separately.
The same governed discrete fixture layer now backs a live
`geiger::fitDiscrete(model='ER')` parity lane. The current ER registry covers
three reviewer-facing surfaces: a binary twenty-four-taxon known-truth panel,
a three-state missing-value pruning review panel, and a four-state
tree-versus-table tip-intersection review panel. Bijux exposes the owned side
through `fit_discrete_mk_model(...)`, and that owned comparative Mk fit now
uses the observed-root likelihood contract that local `geiger::fitDiscrete`
applies by default instead of the flatter equal-root prior used by the
ancestral reconstruction surface. The underlying governed catalog still lives
in `tests/fixtures/metadata/shared_geiger_discrete_fixture_catalog.json`,
keeping binary, three-state, four-state, and sparse six-state discrete
surfaces under durable fixture ids, explicit ER, SYM, and ARD
transition-matrix metadata, and missing-state, constant-state, and
tree-versus-table control panels so later discrete parity lanes do not fall
back to one-off state tables or hand-entered rate matrices.
The same governed discrete fixture layer now also backs a live
`geiger::fitDiscrete(model='ER', transform='lambda')` parity lane. The
current lambda registry covers three reviewer-facing surfaces: an equal-rates
strong-signal review, a curated near-zero lambda weak-signal review, and a
missing-value pruning review. Bijux exposes that owned side through
`fit_discrete_mk_model(model='equal-rates', transform='lambda')`, reuses the
same observed-root likelihood contract as the non-transformed ER lane, and
records the fitted lambda parameter separately from the directed transition
rates. The parity contract compares log-likelihood, information criteria,
pruning behavior, transform identity, and the fitted lambda value against
real local `geiger` output while still checking the directed rate rows on the
transformed branch surface for the stable strong-signal and missing-value
cases. The near-zero weak-signal review remains summary-level on the parity
surface, but the reason is now more specific than raw-rate instability alone:
live `geiger` can return materially different plateau lambda values across
repeated runs on that flat surface while still matching Bijux on
log-likelihood, AIC, AICc, transform identity, and pruning behavior. Bijux
resolves that same weak surface to the lower boundary, so the governed weak
lambda contract compares the stable objective surface rather than pretending
the plateau lambda parameter itself is deterministic.
The same owned surface now also supports
`fit_discrete_mk_model(..., transform='kappa')` and a governed live
`geiger::fitDiscrete(transform='kappa')` parity lane. That registry covers a
branch-length-sensitive ER review, a near-zero weak-signal ER review, and a
missing-value SYM review. For this discrete lane, Bijux now follows the
local `geiger` kappa contract rather than the continuous-model one: kappa is
searched within `[0, 1]`, the transformed tree is built by raising branch
lengths directly to the fitted kappa value, and near-boundary kappa fits
record explicit flattening or untransformed-branch warnings. The weak-signal
kappa review remains summary-level on the parity surface because the
likelihood surface still agrees even when the raw ER rate maximum becomes
unstable near the zero-contrast boundary.
The same owned surface now also supports
`fit_discrete_mk_model(..., transform='delta')` and a governed live
`geiger::fitDiscrete(transform='delta')` parity lane. That registry covers
an ER late-change boundary review, a SYM earliest-change review, and a
missing-value SYM earliest-change review. Bijux follows the local
`geiger` discrete-delta contract by searching delta within
`[exp(-5), 3]`, carrying the transform parameter in the discrete Mk
information-criterion surface, and recording whether the transformed fit
actually improves over the untransformed same-model baseline. The curated
binary time-sensitive delta surface remains a governed owner-side review
fixture rather than a row-level live parity claim, because both Bijux and
local `geiger` agree on the broad boundary behavior there while the interior
rate maximum is not stable enough to present as durable reviewer-facing
pairwise-rate parity.
The same owned surface now also supports
`fit_discrete_mk_model(..., transform='early-burst')` and a governed live
`geiger::fitDiscrete(transform='EB')` parity lane. This discrete lane follows
the raw local `geiger` `a` parameter directly rather than the positive
`rate_change` contract used by the owned continuous early-burst fit, because
the branch-rescaling contract is the same but the reviewer-facing parameter
surface is not. The registry covers one stable ER early-change review, one
weak-signal review, one late-change review, and one missing-value late-change
review. Only the stable early-change surface remains row-level live parity.
The weak-signal and late-change review surfaces stay summary-level because
local `geiger` and Bijux agree on the likelihood envelope, transform identity,
and pruning behavior there while the fitted `a` and raw ER rate maxima remain
flat enough to move between near-equivalent solutions.
The same live lane now also covers `geiger::fitDiscrete(model='SYM')` over
three governed multistate surfaces: a three-state known-truth panel, a
four-state known-truth panel, and a three-state missing-value pruning review
panel. Bijux exposes that owned side through
`fit_discrete_mk_model(model='symmetric')`, and the governed parity contract
compares the full directed rate table, log-likelihood, AIC, and AICc against
real local `geiger` outputs. The sparse six-state symmetric surface remains a
governed owner-side overparameterization review surface for now rather than a
row-level live parity claim, because its weakly identified multirate maximum
is useful for trust review but not yet stable enough to present as durable
pairwise-rate parity.
The same governed discrete lane now also covers
`geiger::fitDiscrete(model='ARD')` over three reviewer-facing surfaces: a
binary known-truth panel, a three-state missing-value pruning panel, and a
four-state weak-identification review panel. Bijux exposes that owned side
through `fit_discrete_mk_model(model='all-rates-different')`, and the stable
binary and missing-value ARD cases now compare the full directed rate table,
log-likelihood, AIC, and AICc against real local `geiger` outputs. The
four-state ARD surface remains in the live lane as a governed weakly
identified review case with a relaxed directional-rate tolerance, because
local `geiger` and Bijux reach near-equivalent maxima on that twelve-parameter
surface while still agreeing on the likelihood envelope and the need for
boundary and equal-rates-baseline caution. The sparse six-state ARD surface
still stays owner-side only, where it continues to govern overfit and
overparameterization warnings without pretending its full row-level maximum is
durable enough for live parity.
`geiger::fitDiscrete(model='meristic')` is explicitly excluded in this round.
Local `geiger` exposes a distinct meristic lane with an integer-state contract,
while Bijux currently offers generic ordered-state Mk support over user-named
categorical states. That ordered-state surface remains available for reviewer
work, but it is not claimed as meristic parity, and the runtime now rejects
`meristic` requests directly with that reason instead of letting them blur
into the ordered ER, SYM, or ARD paths.

`parity --reference-source phytools-live` is the governed live `phytools`
execution harness. It uses the same checked-in `Rscript` orchestration model
as the live `ape` lane, records the R version, `phytools` version, Bijux
version, Bijux commit, function name, input fixtures, tolerance, pass or fail
state, mismatch reason, and reproducible artifact root for every governed case,
and writes one summary TSV plus one observation TSV just like the other parity
surfaces. The live `phytools` runner itself now lives as a real package under
`parity/phytools/runner/`: registry selection and case serialization live with
the registry, `dispatch.py` owns Bijux payload routing, `comparison.py` owns
summary and row equivalence policy, `execution.py` owns orchestration, and
`__init__.py` is only the curated API gateway. The initial live `phytools`
registry is intentionally narrow for this
goal: it currently covers `phytools::phylosig(method='lambda')`,
`phytools::phylosig(method='K')`, `phytools::fitMk(model='ER')`,
`phytools::fitMk(model='SYM')`, `phytools::fitMk(model='ARD')`,
`phytools::make.simmap(model='ER')`, `phytools::make.simmap(model='SYM')`,
`phytools::make.simmap(model='ARD')`, `phytools::countSimmap`,
`phytools::densityMap`,
`phytools::describe.simmap`,
`phytools::sim.history`,
`phytools::fastBM`,
`phytools::sim.corrs`,
`phytools::pgls.SEy`,
`phytools::rerootingMethod`,
`phytools::fastAnc`, and `phytools::anc.ML` on governed strong-signal,
weak-signal, non-ultrametric, discrete-state, and missing-value comparative
fixtures. The live lambda lane includes one
non-ultrametric case that tracks the live `phytools` likelihood surface within
tolerance, and the governed checked-in parity surface now also covers strong
and weak twenty-four-taxon `phytools::phylosig` lambda and K references. The
owned signal surface also keeps its input policy explicit across all three
entrypoints: rooted non-ultrametric trees are accepted and reported, missing
overlapping trait values are pruned and listed, seeded permutation p-values
stay reproducible, and constant post-pruning trait vectors fail with one clear
comparative-method error instead of returning misleading scalar summaries. The
owned signal surface now also exposes fixed-lambda likelihood evaluation,
likelihood-ratio reporting against the zero-signal boundary, and optimizer
diagnostics instead of reducing Pagel's lambda to one opaque scalar. The owned
K-test surface now also keeps seeded permutation p-values plus explicit
null-distribution summaries. The owned discrete Mk surface now
also exposes one flat-root `fitMk`-style likelihood contract with ER rate
fitting, SYM pairwise-rate fitting, log-likelihood, AIC, AICc,
missing-value pruning audit, ER baseline comparison metrics, and one directed
rate-matrix ledger for binary and multistate traits. The live `fitMk` lane
now covers governed ER binary and multistate cases, governed SYM
multistate cases, and governed ARD binary plus weakly identified multistate
cases, including missing-value-pruned surfaces. The live `pgls.SEy` lane now
covers governed fixed-lambda Brownian covariance cases for one simple numeric
regression plus one categorical and one interaction-coded regression. That
boundary is explicit on purpose: installed `phytools 2.5.2` does not export a
general `phytools::pgls` surface, so the governed live lane stays on
`phytools::pgls.SEy` with `lambda = 1.0`, while the broader exact PGLS
contract for estimated lambda and full coefficient parity remains the
checked-in R `ape` plus `nlme` reference lane. The live
`make.simmap` lane now covers governed clean binary, clean multistate, and
missing-value-pruned binary ER cases; governed clean multistate and
missing-value-pruned multistate SYM cases; and governed binary plus
missing-value-pruned binary ARD cases at one governed seed and one governed
replicate count of 128 maps per case. It now also carries the owned fitted
discrete-Mk audit through the stochastic-mapping surface, including fitted
model identity, parameter count, log-likelihood, AIC, AICc, baseline-model
comparison, optimizer convergence, and weak-fit warnings. It compares
distributional envelopes only: excluded taxa, total-transition-count mean plus
interval, transition-count summary rows, and time-in-state summary rows. It
does not claim exact stochastic-history identity with `phytools`. Governed
multistate ARD cases stay on summary-envelope parity only when weakly
identified boundary rates make row-level transition summaries unstable across
optimizers.
The same owned stochastic-map surface now also exposes one `describe.simmap`-
style summary contract over saved map collections, including total-change
summary, transition-count rows, time-in-state rows, and per-branch
state-occupancy rows. The live `describe.simmap` lane now covers governed
clean binary, clean multistate, clean multistate SYM, and missing-value-pruned
binary cases at one governed seed and one governed replicate count of 128 maps
per case. It compares those summary rows, including branch occupancy, without
claiming exact stochastic-history identity.
The owned stochastic-map collection surface now also exposes one
`countSimmap`-style transition-count contract over saved map collections,
including one per-replicate total-transition ledger, one flat event table, one
aggregate transition matrix, and one per-branch directional transition summary.
The live `countSimmap` lane now covers governed clean binary, clean
multistate, clean multistate SYM, and missing-value-pruned binary cases at one
governed seed and one governed replicate count of 128 maps per case. It
compares total-transition envelopes plus directional transition-count rows,
including zero diagonal state pairs such as `0->0`, without claiming exact
stochastic-history identity.
The same owned stochastic-map collection surface now also exposes one
`densityMap`-style branch-probability contract over saved map collections,
including one branch-probability table, one branch-level probability envelope,
one slice-level probability table at governed resolution, and one report-ready
HTML or SVG artifact. The live `densityMap` lane is intentionally narrower
than the owned surface: it currently covers governed clean binary and
missing-value-pruned binary ER cases only. It compares per-branch posterior
probability summaries and branch-level uncertainty against live
`phytools::densityMap`, and does not claim pixel-perfect plotting parity.
The owned simulation surface now also exposes one fixed-tree
`simulate_discrete_histories(...)` contract over an explicit discrete rate
matrix. It supports binary and multistate states, fixed or probabilistic root
states, seeded replication, true branch-history segments, true transition
events, tip-state truth tables, node-state truth tables, branch-history truth
tables, and one parity-ready summary ledger over transition counts,
time-in-state totals, and tip-state frequencies. The live `sim.history` lane
now covers governed binary and multistate no-change plus high-rate fixtures on
fixed trees, plus one governed probabilistic binary root-prior fixture. It
compares distribution-summary envelopes against real
`phytools::sim.history`, including total-transition-count summaries,
transition-count rows, time-in-state rows, tip-state-frequency rows, and the
declared root-state policy, and does not claim exact simulated-history
identity across languages.
The same owned simulation surface now also exposes
`simulate_brownian_trait_collection(...)` plus one replicate trait-table writer
and one Brownian summary writer over tip distributions and tip covariances.
`simulate_brownian_traits(...)` and `simulate traits-brownian` now both accept
either `sigma` or explicit `sigma_squared`, preserve the resolved Brownian rate
parameter in the report, and keep seeded replicate collections on fixed trees
for downstream recovery work. The live `fastBM` lane now covers governed
low-variance, root-shift high-variance, and six-taxon Brownian cases. It
compares distribution summaries and tip-covariance rows against real
`phytools::fastBM` without claiming exact cross-language draws.
The same owned simulation layer now also exposes
`simulate_speciational_traits(...)` and
`simulate_speciational_trait_collection(...)` for the narrower
`geiger::sim.char(model='speciational')` contract, where every positive branch
contributes one Brownian step regardless of branch-length magnitude while
zero-length branches remain zero. The governed `geiger::sim.char` envelope now
covers three local reference cases on the fixed rooted non-ultrametric control
tree: one Brownian case, one speciational case, and one asymmetric binary
discrete rate-matrix case. Bijux exposes that review surface through
`validate_geiger_sim_char_reference_examples()` and
`simulate validate-sim-char-reference`, and the contract is explicit about its
scope: it compares tip-distribution, tip-covariance, and tip-state-frequency
summary envelopes against stored local `geiger::sim.char` outputs, keeps the
true generating parameters in the report, uses fixed seeds, and does not claim
exact cross-language random draws.
The same governed simulation stack now also exposes one packaged
`discrete_mode_recovery_panel` benchmark over explicit discrete Mk truth
matrices and stored local `geiger::fitDiscrete` references. It covers one
stable three-state ER recovery case, one stable three-state SYM recovery case,
one three-state ARD weak-identification review, and one five-state
overparameterized ARD failure review across the governed twelve-taxon and
twenty-four-taxon trees. Bijux exposes that benchmark through
`run_discrete_mode_recovery_panel_workflow()` and
`demo discrete-mode-recovery-panel`, while the packaged dataset and artifact
materialization surface also stays available through
`load_discrete_mode_recovery_panel_dataset()`,
`write_discrete_mode_recovery_panel_workflow_bundle(...)`, and
`run_discrete_mode_recovery_panel_demo(...)`. The benchmark writes paired
rate-recovery, rate-comparison, model-choice, and transform-parameter
ledgers, compares Bijux truth error against stored local `geiger` truth error
wherever the panel defines that truth, and keeps the scientific boundary
explicit: the ARD surfaces are review cases and are not counted as successful
recovery when the information-criterion surface prefers simpler models, when
fitted state pairs drop out entirely, or when transform-bearing cases are not
part of the governed panel.
The same owned simulation surface now also exposes
`simulate_correlated_brownian_trait_collection(...)` for two or more
continuous traits on one fixed tree from one explicit evolutionary covariance
matrix, plus writers for replicate tip ledgers and multivariate summary rows
over root states, evolutionary covariance, tip distributions, and tip
covariances. The CLI surface `simulate traits-brownian-correlated` accepts
either one covariance matrix directly or one correlation matrix plus per-trait
standard deviations, rejects invalid non-positive-definite covariance inputs
explicitly, and keeps the generating parameter matrix in the returned report.
The live `sim.corrs` lane now covers governed low-correlation,
negative-correlation root-shift, and three-trait six-taxon cases. It compares
distribution summaries, tip-covariance rows, and tip-correlation rows against
real `phytools::sim.corrs` without claiming exact cross-language draws.
The live
`rerootingMethod` lane now covers governed ER binary, governed ER multistate,
governed ER missing-value-pruned, governed SYM multistate, and governed SYM
missing-value-pruned node-probability cases under the same flat equal root
prior that `phytools::rerootingMethod` inherits from `fitMk`. Bijux reports
that boundary explicitly: ER and SYM runs with `--root-prior-mode equal` are
governed live `phytools::rerootingMethod` parity surfaces, while Fitch,
ordered-state, ARD, empirical-root-prior, and fixed-root-prior runs remain
owned Bijux review surfaces without a false live `phytools` parity claim.
The live
continuous ancestral lanes now compare stable node-signature rows, standard
errors, and 95% intervals against real `phytools` execution instead of only
checked-in expected JSON.
The `ape::nj` lane now covers one governed analytical three-taxon matrix plus
four-taxon ultrametric and non-ultrametric matrices. On the owned Bijux side,
neighbor joining no longer delegates through Biopython for that method: Bijux
now builds one deterministic NJ tree in-repo, validates zero-diagonal and
nonnegative matrix assumptions explicitly, produces branch lengths, and
resolves tied joins by stable taxon ordering so distance-tree recovery stays
reproducible. The same owned distance-tree core is now also reusable directly
from Python through `build_distance_tree_from_genetic_distance_matrix(...)`,
so one loaded `GeneticDistanceMatrix` no longer has to restart from a
path-based alignment or table wrapper before recovering one NJ, BIONJ,
single-linkage, complete-linkage, UPGMA, or WPGMA tree.
Within that clustering lane, UPGMA, WPGMA, single-linkage, and
complete-linkage now route through one shared agglomerative engine with
explicit method-specific update rules instead of carrying separate merge loops.
Fixed topologies can also be re-fitted against imported distance matrices and
scored by total fitted branch length through the owned minimum-evolution
surface, which ignores any pre-existing branch lengths on the input tree.
That fixed-topology lane now also includes owned classical
Fitch-Margoliash weighted least-squares fitting, which fits branch lengths by
distance-weighted residual minimization instead of reusing the unweighted OLS
solution when weighted fitting is requested.
The same owned fixed-topology lane now exposes ordinary least squares as a
first-class report surface with fitted branch lengths, a residual matrix, RSS,
matrix rank, condition number, and explicit negative-branch reporting instead
of silently clipping unconstrained solutions.
That same lane now also exposes nonnegative least-squares fitting as a
separate owned surface, which constrains every branch length to zero or above
and reports the active zero-constraint branches instead of reusing the
unconstrained OLS solution.
That fixed-topology distance lane now also exposes owned patristic residual
diagnostics, which compare one observed distance matrix against one supplied
tree's branch-length-implied tip distances and export ranked
`distance_residuals.tsv` rows for every unique taxon pair instead of only one
global RSS summary.
The alignment-side corrected-distance lane now also exposes explicit
pair-level saturation diagnostics through `alignment distance-saturation`.
JC69, K80, and TN93 corrected distances that become undefined or tend to
infinity now produce governed pair warnings and block distance-tree inference
before NJ, BIONJ, or linkage builders are called.
That same tree-building lane now also owns missing-distance handling through
one explicit policy engine. Alignment-derived and imported matrices can reject
incomplete pairs or impute them by mean distance, nearest valid distance, or
the tightest available triangle bound before distance-tree inference runs.
The distance diagnostics surface now also owns explicit ultrametricity testing
for raw matrices through the three-point condition, with governed violating
triple rows, the maximum observed deviation, and the applied tolerance instead
of relying on row-sum heuristics or on whether one later UPGMA tree happens to
be ultrametric after clustering.
That same diagnostics family now also owns quartet additivity review through
the four-point condition, with governed `four_point_violations.tsv` output that
records each violating quartet, the three competing pair-sum totals, the
best-supported split, and the violation magnitude instead of collapsing the
matrix to one opaque “valid” boolean.
Imported distance matrices now also support one leave-one-taxon-out taxon
influence diagnostic that rebuilds the chosen distance tree and recomputes
patristic residual RSS after each exclusion, then ranks taxa by RSS
improvement and rooted RF improvement against a supplied reference tree
instead of proxying influence by raw missingness alone.
Imported distance matrices now also support one owned leave-one-taxon-out
jackknife that prunes the baseline inferred tree, rebuilds the reduced tree for
each removed taxon, and exports rooted RF distance, residual change, affected
clades, and rebuilt Newick trees for every jackknife row.
Imported distance matrices now also support one owned distance-method
comparison workflow that builds NJ, BIONJ, UPGMA, and WPGMA trees on the same
matrix, scores those exact topologies by patristic residual RSS, balanced
minimum evolution, and OLS refits, and exports a rooted RF matrix plus
assumption-warning ledger instead of ranking methods by only one criterion.
Alignment-side distance bootstrap now also records the exact resampled site
indices and rebuilt replicate tree for each seeded replicate, so the owned
bootstrap workflow can prove it resamples alignment columns with replacement
before recomputing distances and rebuilding trees rather than bootstrapping one
already-computed distance matrix.
Imported distance matrices can also start from one owned NJ or BIONJ tree and
be hill-climbed by rooted NNI under the owned balanced minimum-evolution
objective, with governed search traces that record each accepted objective
decrease until the final local optimum.
The owned tree-manipulation core now also lives directly on `PhyloTree`.
Outgroup rooting, unrooting, keep-tip pruning, drop-tip pruning, clade
extraction, MRCA lookup, and monophyly review all run through one native
topology surface with explicit rootedness and branch-length policy instead of
splitting those transforms across separate tree backends.
That topology lane now also lives as a real `bijux_phylogenetics.core.topology`
package instead of one overgrown flat module. Shared report dataclasses live
in `core.topology.models`, branch-preserving change summaries and collapse
workflow live in `core.topology.transformation`, subtree extraction plus MRCA
and monophyly review live in `core.topology.subtree`, deterministic child-order
transforms live in `core.topology.ordering`, and outgroup rooting, midpoint
rerooting, unrooting, and rooting-ledger writing live in
`core.topology.rooting`. The package root is now the curated gateway so
downstream tree, pruning, parity, and dataset surfaces import one stable
topology API without depending on a foundational kitchen-sink file.
The taxon-audit core now follows that same package shape under
`bijux_phylogenetics.core.taxonomy`. Shared taxon dataclasses live in
`core.taxonomy.models`, label-collision heuristics live in
`core.taxonomy.identity`, downstream-safe label rewriting lives in
`core.taxonomy.normalization`, synonym-table parsing and accepted-name export
live in `core.taxonomy.synonyms`, namespace and rank heuristics live in
`core.taxonomy.classification`, duplicate-identity and mapping-conflict review
live in `core.taxonomy.conflicts`, and the final reviewer-facing audit
composition lives in `core.taxonomy.audit`. The package root is now a curated
gateway with explicit exports instead of another foundational spill file.
The owned simulation baseline now also exposes direct single-tree surfaces
through `simulate_random_tree(...)` and `simulate_coalescent_tree(...)`, so a
caller that needs one governed baseline tree with its summary record no longer
has to go through one batch-only wrapper.
That simulation lane now also lives as a real Python package instead of one
overgrown `simulation/__init__.py` file. The package root is a curated lazy
public API surface, shared dataclasses live in `simulation.models`, summary
math lives in `simulation.statistics`, trait and Brownian propagation lives in
`simulation.propagation`, discrete state validation lives in
`simulation.discrete_policy`, and shared stochastic jump counting lives in
`simulation.stochastic`, so dataset, parity, and comparative consumers no
longer depend on package-root implementation spillover.
The `ape::pic` lane now covers balanced rooted ultrametric, pectinate rooted
non-ultrametric, and six-taxon clean comparative fixtures through one governed
trait-table catalog. On the owned Bijux side, `comparative contrasts
--contrasts-out <table.tsv>` writes one `independent-contrasts.tsv` ledger with
stable ape-style `node_id` values, left-versus-right descendant partitions,
standardized contrasts, and expected variances. The JSON report now also
preserves one explicit input audit with rootedness and ultrametric reporting,
minimum and maximum root-to-tip depths, the owned missing-value pruning policy,
and the pruned taxon list. Missing trait values therefore remain an explicit
owned pruning surface, while negative branch lengths are rejected explicitly as
an invalid comparative-analysis boundary instead of being pushed through to
live `ape::pic`. The same owned Brownian comparative core is now also reusable
directly from Python through
`summarize_brownian_covariance_from_tree(...)` on one native `PhyloTree` and
`compute_phylogenetic_independent_contrasts_from_dataset(...)` on one loaded
`ComparativeDataset`, so covariance review and PIC no longer need to restart
from path-based loaders once the in-memory comparative inputs already exist.
The governed `phytools::phyl.resid` lane now covers six-taxon Brownian,
estimated-lambda, and missing-value-pruned allometry fixtures. On the owned
Bijux side, `comparative phylogenetic-residuals` preserves that same
response-predictor review in native Python while keeping one explicit
excluded-taxa ledger and one coefficient ledger. The governed
`phytools::phylANOVA` lane now covers six-taxon unequal-group and
missing-value-pruned group-effect fixtures. On the owned Bijux side,
`comparative phylogenetic-anova` carries the same seeded null-distribution
surface, one pairwise Holm-adjusted comparison ledger, and one low-sample
warning contract without relying on opaque plotting-side behavior.
The `ape::ace` continuous lane now covers balanced rooted ultrametric,
pectinate rooted non-ultrametric, six-taxon Brownian, and pruned missing-value
fixtures through that same governed trait-table catalog. On the owned Bijux
side, `ancestral continuous` now carries one explicit Brownian fit diagnostic
surface with ultrametric state, root-to-tip depth bounds, covariance rank and
conditioning, solver regularization status, and one GLS likelihood summary.
Its Brownian estimator surface is now explicit: the default `ace-pic`
estimator preserves the governed `ape::ace(type='continuous', method='pic',
CI=TRUE)` parity lane, `--estimator anc-ml` exposes the governed live
`phytools::anc.ML` lane with Brownian log-likelihood, fitted sigma-squared,
closed-form optimizer diagnostics, stable node-signature tables, standard
errors, 95% intervals, and explicit missing-value pruning, while
`--estimator fast-anc` exposes the governed live `phytools::fastAnc` lane
with the same node-signature review surface. The same owned continuous
ancestral runtime is now also reusable directly from Python through
`reconstruct_continuous_ancestral_states_from_dataset(...)` once one
`AncestralContinuousDataset` has already been loaded.
The `ape::ace` discrete lane now covers governed ER, SYM, and ARD fixtures,
including balanced, pectinate, six-taxon, and pruned missing-value cases
through that same shared trait-table catalog. On the owned Bijux side,
`ancestral discrete` now emits fitted transition-rate tables, log-likelihood,
parameter count, AIC, weak-fit warnings, and ER baseline-comparison data
alongside node probabilities, and it supports owned `equal`, `empirical`, and
`fixed` root-prior policies. The live parity lane is scoped honestly to
`ape::ace(type='discrete', model='ER'|'SYM'|'ARD')`, so root-prior controls
remain an explicit Bijux-owned review surface rather than a false live `ape`
parity claim. The same summary surface now also reports whether a given
likelihood run is comparable to live `phytools::rerootingMethod`: ER and SYM
with the equal root prior are governed rerooting-parity surfaces, while ARD,
Fitch, ordered-state, empirical-root-prior, and fixed-root-prior runs are
flagged explicitly as non-comparable. The same owned discrete ancestral runtime is now also reusable
directly from Python through
`reconstruct_discrete_ancestral_states_from_dataset(...)` once one
`AncestralDiscreteDataset` has already been loaded.
The `ape::dist.dna` lane now covers raw nucleotide distance, JC69, K80, F81,
and TN93 distance over governed clean, gapped pairwise-deletion, gapped
complete-deletion, ambiguity-bearing, identical-sequence, high-divergence,
missing-data, and unequal-length-invalid fixtures. Bijux accepts the
ape-compatible `raw`, `jc69`, `k80`, `f81`, and `tn93` model aliases on the
owned distance surface, keeps `p-distance`, `jukes-cantor`,
`kimura-2-parameter`, `felsenstein-81`, and `tamura-nei-93` as the canonical
internal labels, and reports saturated corrected-distance pairs explicitly as
either undefined or infinite instead of hiding them. F81 and TN93 estimate
alignment-wide resolved nucleotide frequencies, write one `--parameters-out`
ledger for reviewer-facing model coefficients, and warn explicitly when the
resolved composition breaks TN93 assumptions instead of silently falling back
to a simpler model. `alignment distance-matrix` can also write one
`--components-out` TSV ledger with pairwise mismatch, transition,
transversion, AG-transition, CT-transition, ambiguity, and saturation fields
alongside the distance matrix. Unequal-length alignments still fail
explicitly instead of deferring that failure until later matrix handling.
The `ape::base.freq` lane now covers lowercase, ambiguity-bearing,
missing-data, and all-gap-or-missing alignments. On the owned Bijux side,
`alignment composition --base-frequency-out <table.tsv>` writes one combined
alignment-plus-sequence TSV ledger with `scope`, `identifier`, `state`,
`count`, and `frequency` columns, returns the same literal-state frequencies
in JSON output, and reports composition outlier sequences beside those base
frequency rows. Ambiguity codes, gaps, and explicit missing states are counted
as literal states to match `ape::base.freq`, and all-gap or missing inputs
warn explicitly instead of fabricating canonical A/C/G/T content.
The `ape::seg.sites` lane now covers lowercase, invariant,
one-variable-site, gap-bearing, ambiguity-bearing, missing-data, and
all-gap-or-missing DNA alignments. On the owned Bijux side,
`alignment segregating-sites --site-table-out <table.tsv>` writes one
reviewer-facing `segregating-sites.tsv` ledger with site positions plus
literal and ape-normalized state summaries. Leading and trailing gaps are
normalized to `N` to match live `ape::seg.sites`, explicit missing states do
not create segregating sites by themselves, and incompatible ambiguity states
or internal gaps remain visible as real segregating-site evidence.
The `ape::trans` lane now covers valid-reading-frame, ambiguous-codon,
internal-stop, terminal-stop, frame-truncation, and vertebrate-mitochondrial
genetic-code fixtures. On the owned Bijux side, `alignment translate
--codon-validation-out <table.tsv> --excluded-sequences-out <table.tsv>`
writes one amino-acid FASTA plus a codon-level validation ledger. The aligned
translation surface now matches live `ape::trans` by truncating trailing
partial codons with an explicit warning instead of hard-failing, while the
stricter `prepare_coding_sequences_for_alignment` surface still owns the
pre-alignment exclusion policy for frame errors, ambiguous codons, and
premature stop codons in serious codon-aware workflows.
The
`ape::root` lane now uses the same shared tree catalog for single-tip
outgroups, monophyletic multi-tip outgroups, already-rooted trees, missing
outgroups, and non-monophyletic outgroups, with rooted clades and branch
lengths compared against live `ape::root` and ambiguous rooting rejected
explicitly on the Bijux side. That surface now reroots trees natively on
`PhyloTree` rather than delegating outgroup rooting through Biopython. The
`ape::unroot` lane now compares rooted-tree unrooting, post-outgroup-rooting
unrooting, already-unrooted inputs, and malformed tree failures against live
`ape::unroot`, and it exposes the root-edge redistribution policy explicitly:
Bijux now merges the removed root-edge length into the retained sibling branch
to match `ape::unroot` rather than moving that length into the expanded clade.
The `ape::drop.tip` lane now compares single-tip exclusion, multi-tip
exclusion, rooted and unrooted root-state changes, and unknown-tip handling
against live `ape::drop.tip`, while keeping one explicit product-safety rule:
Bijux rejects pruning requests that would leave fewer than two retained taxa
instead of emitting one-tip trees into downstream scientific workflows.
The `ape::keep.tip` lane now compares valid keep-set pruning against live
`ape::keep.tip`, including selected-tip order differences and rootedness
changes after pruning. Bijux keeps two workflow-facing extensions outside that
live parity subset: tree and trait matching paths still report absent
requested taxa instead of failing immediately, and they stop clearly when
fewer than two retained taxa would remain.
The `ape::extract.clade` lane now compares rooted subtree extraction against
live `ape::extract.clade` for root and internal-node cases on the governed
internal-label fixture, keeps internal labels where the extracted source node
had them, and treats tip-node or out-of-bounds requests as explicit extraction
errors. Alongside that live parity lane, Bijux now exposes one owned
descendant-taxa extractor so callers can resolve the same subtree by stable
taxon identity instead of only by ape-style node number.
The `ape::getMRCA` lane now compares stable ape-style internal node ids and
matched descendant tip sets for two-tip, many-tip, root-clade, duplicate-tip,
rooted-polytomy, and already-rooted-outgroup cases, while keeping one explicit
workflow-side rule on the Bijux surface: missing requested taxa fail clearly
instead of bubbling through as an opaque reference-side parser condition.
The `ape::is.monophyletic` lane now compares rooted and unrooted monophyly
calls against live `ape` with explicit reroot policy, full-tip-set behavior,
singleton and mixed-missing handling, rooted-polytomy behavior, post-rooting
behavior, and all-missing reroot failures. Bijux also records the matched
MRCA node and any extra descendant taxa that make a direct clade
non-monophyletic, so the parity lane remains reviewer-usable instead of only
returning a boolean.
The `ape::cophenetic.phylo` lane now compares rooted and unrooted pairwise
tip-distance matrices against live `ape`, keeps the taxon order explicit, and
matches one governed long-form distance ledger instead of only checking a
printed matrix. On the owned Bijux side, missing branch lengths now fail
explicitly for tip-distance calculations unless the caller opts into an
explicit unit-length fallback policy.
The `ape::dist.topo` lane now compares identical rooted trees, rooted
child-order rotations, one-conflict rooted pairs, rooted tree-versus-polytomy
pairs, one governed unrooted split conflict, and one governed 128-tip rooted
pair against live `ape`. It matches one explicit RF-style split ledger rather
than only a scalar distance and keeps rooted-versus-unrooted policy explicit
per case. On the owned Bijux side, `adapter compare --split-table-out` now
writes the same split ledger directly for review and downstream automation.
Those distance and support workflows now share one native clade-set core for
canonical rooted clades, canonical unrooted splits, RF metrics, partial-taxon
scope handling, and polytomy reporting, so parity, topology comparison,
tree-set support mapping, and posterior clade summaries all use the same
split identity contract.
The `ape::vcv.phylo` lane now compares Brownian shared-ancestry covariance on
rooted ultrametric, rooted non-ultrametric, unrooted branch-length, and
singular zero-branch trees against live `ape`, and it persists the compared
covariance tables automatically whenever parity fails. On the owned Bijux
side, `summarize_brownian_covariance(...)` now exposes the same matrix with
explicit tip order, rejects missing or negative branch lengths instead of
silently regularizing them away, and reports singular-versus-near-singular
state directly from the raw covariance matrix.
The `ape::node.depth.edgelength` lane now compares rooted ultrametric,
rooted non-ultrametric, zero-branch-length, and post-outgroup-rooting trees
against live `ape`, using stable ape-style node ids for both tips and
internal nodes. On the owned Bijux side, `compute_tree_node_depths(...)`
records root-to-node branch-length depth for every tip and internal node,
writes one deterministic node-depth ledger, and rejects incomplete branch
lengths instead of guessing edge-count or zero-length fallbacks.
The `ape::branching.times` lane now compares rooted ultrametric trees,
internal-node labels, one medium ultrametric tree, and one zero-length
internal-branch ultrametric tree against live `ape`, again preserving
ape-style internal node ids. On the owned Bijux side,
`compute_tree_branching_times(...)` writes one deterministic internal-node
branching-time ledger, reports the root age explicitly, and intentionally
rejects non-ultrametric trees instead of returning the misleading negative or
inconsistent node ages that `ape::branching.times` will still emit on invalid
inputs.
The `ape::is.ultrametric` lane now compares exact ultrametric, governed
near-ultrametric, tight-tolerance near-ultrametric, and clearly
non-ultrametric trees against live `ape`. On the owned Bijux side,
`assess_tree_ultrametricity(...)` reports the ape-style criterion name,
criterion value, tolerance, maximum tip-depth deviation, offending taxa, and
one deterministic `ultrametric-diagnostics.tsv` ledger. That same surface is
now reused before rooted Brownian, OU, and diversification workflows claim
time-tree compatibility, so those methods no longer rely on separate hidden
ultrametric checks with mismatched tolerances.
The
`ape::write.tree` lane now roundtrips Bijux-written Newick through live `ape`
for rooted, unrooted, internal-label, support-label, quoted-label, and
multiple-tree cases, while the Bijux writer rejects unnamed tips, empty tree
sets, and non-finite branch lengths instead of writing malformed Newick
silently. Tree IO parity in that lane is structural rather than text-based:
equivalent child reorderings pass, while rootedness, tip-set, clade or split,
branch-length, and internal-label drift produce specific mismatch reasons.
The live `ape` runner itself is now packaged by owned boundaries under
`parity/ape/runner/`: `dispatch.py` maps governed cases to Bijux payload
builders, `reference_payloads.py` loads live runner outputs, `comparison.py`
and `mismatch_policy.py` own equivalence rules, `failure_artifacts.py` owns
reproducibility bundles, `execution.py` owns orchestration, and
`__init__.py` stays as the public API gateway instead of remaining the real
module.

The smaller ancestral review lanes now use the same governed shared trait-table
fixture catalog in `tests/fixtures/metadata/shared_trait_table_fixture_catalog.json`.
That corpus covers continuous, binary discrete, multistate discrete,
missing-value, mismatch, duplicate-row, constant-trait, categorical-predictor,
and misordered-row cases, and the owned discrete-reference checks plus the
live `ape::ace` spot checks resolve those trait tables by durable fixture id
instead of keeping separate ad hoc path lists for Bijux and `ape`.

`report release-truth` is the governed pre-release summary surface. It consumes
actual pytest JUnit XML reports for the full test lane and the real-engine test
lane, reruns the owned workflow-validation, release-gate, parity, and
stress-suite checks, and writes one HTML report plus one machine manifest that
states total tests passed, failed, and skipped; real-engine tests passed,
failed, and skipped; supported and experimental workflow surfaces; available
flagship datasets; reference-parity case coverage; stress-suite workload
coverage; and aggregated known limitations.

Tree-distance hard cases now have their own governed validation surface:

```bash
uv run bijux-phylogenetics topology distance-reference --json
```

That suite locks rooted RF, unrooted RF, normalized RF, and branch-score
distance against checked DendroPy 5.0.8 references, including polytomies,
star-tree collapse, branch-length-only disagreement, rooting-only disagreement,
and shared-taxa pruning cases. It also keeps the owned
`--taxon-overlap-policy require-identical` rejection behavior explicit for both
RF and branch-score comparison.

Branch-support parsing now has its own governed validation surface as well:

```bash
uv run bijux-phylogenetics topology support-reference --json
```

That suite locks IQ-TREE UFBoot labels, composite SH-aLRT/UFBoot labels,
FastTree local support, and posterior clade frequencies against checked
fixtures, proves that support comparison still maps values by clade instead of
node order, and flags bootstrap-versus-posterior topology mismatch explicitly
when one support surface contains clades the other does not.

The Python tree-comparison surface now also lives as a real
`bijux_phylogenetics.compare.topology` package instead of one overgrown module.
Shared report dataclasses and enums live in `compare.topology.models`,
rooted-versus-unrooted RF workflow and shared-taxa policy live in
`compare.topology.comparison`, branch-length and branch-score review live in
`compare.topology.branch_lengths`, clade-overlap and shared-taxa pruning live
in `compare.topology.overlap`, support pairing and conflict classification live
in `compare.topology.support`, and reviewer-facing TSV writers live in
`compare.topology.tables`. The package root is now the curated public gateway,
so downstream code imports one stable comparison surface without reaching
through a single kitchen-sink file.

## Read this next

- package docs: [Runtime package docs](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
- source directory: [Runtime source directory](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)
- changelog: [Runtime package changelog](https://github.com/bijux/bijux-phylogenetics/blob/main/packages/bijux-phylogenetics/CHANGELOG.md)
- security policy: [Security policy](https://github.com/bijux/bijux-phylogenetics/blob/main/SECURITY.md)
