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
macroevolution analysis, evidence bundle creation, and HTML report generation.

## Install

`bijux-phylogenetics` supports Python 3.11 and newer.

```bash
python3.11 -m pip install bijux-phylogenetics
bijux-phylogenetics --help
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
- reconstruct continuous ancestral states under Brownian or OU-style trait models
- reconstruct discrete ancestral states under Fitch parsimony or likelihood-style ER, SYM, and ARD models with explicit ambiguity and low-confidence reporting, root-prior controls, and fitted transition-rate ledgers
- validate discrete ancestral likelihood surfaces against governed `ape::ace` references plus root-prior, ordered, irreversible, and ambiguity policy checks
- compare continuous ancestral reconstructions across two supported models, summarize ancestral sensitivity across model, tree, pruning, or coding choices, and package publication-ready ancestral figures
- validate discrete geographic state coding, detect incomplete ordered vocabularies, estimate ancestral node states under ordered or unordered assumptions, compare equal-rates, symmetric, and all-rates-different models, export node and transition tables, highlight model-sensitive ancestral regions, simulate approximate stochastic maps, and render discrete-state HTML reports
- audit alignment inference readiness, validate model-selection outputs against engine artifacts, verify inferred-tree taxa against the alignment, inspect metadata-group clustering, classify inference failures, and validate bootstrap tree sets before interpreting engine outputs
- estimate lineage-through-time curves, simple Yule or birth-death diversification rates, sampling-aware corrections, clade outlier summaries, and trait-linked diversification tables for rooted ultrametric trees
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
- build Neighbor-Joining or UPGMA trees from computed distance matrices, bootstrap site-resampled trees, summarize clade support, and write reproducibility bundles
- validate imported long-form distance matrices, detect nonmetric violations, and build trees from imported distances
- audit NJ and UPGMA method assumptions, including explicit UPGMA ultrametric-clock violations for computed or imported distance matrices
- load posterior tree sets, compute consensus trees, and export clade-frequency or pairwise tree-distance summaries
- cluster identical rooted topologies, detect unstable taxa or clades, and compare two posterior tree sets
- simulate birth-death or coalescent trees, Brownian, OU, or early-burst continuous traits, discrete traits, and DNA or protein alignments
- ship governed internal recovery panels that check whether Brownian, OU, and early-burst trait models recover known simulation truth with explicit parameter-error and warning ledgers
- benchmark validation, tree comparison, and alignment diagnostics across increasing problem sizes
- root trees on explicit outgroups or reroot them by midpoint
- audit rooting, ordering, clade extraction, and pruning transforms with before/after summaries and retained-versus-removed taxon reasoning
- validate tree roundtrips across Newick, Nexus, and phyloXML formats with topology-preservation checks, support-label audits, and semantic-loss warnings
- audit ambiguous taxon identities, synonym candidates, namespace mixing, workflow taxon loss, and cross-run taxon stability before downstream comparison or linkage
- produce HTML reports and file-level evidence manifests

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
bijux-phylogenetics alignment distance-suitability alignment.fasta --model jukes-cantor --json
bijux-phylogenetics alignment distance-assumptions alignment.fasta --model p-distance --json
bijux-phylogenetics alignment build-tree alignment.fasta --method upgma --out upgma-tree.nwk
bijux-phylogenetics alignment compare-distance-to-tree alignment.fasta inferred-tree.nwk --method neighbor-joining --json
bijux-phylogenetics alignment bootstrap-tree alignment.fasta --method neighbor-joining --replicates 200 --support-out artifacts/distance-support.tsv --tree-set-out artifacts/distance-bootstrap.trees --json
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
bijux-phylogenetics report taxonomy --tree tree.nwk --synonym-table taxonomy.tsv --metadata metadata.tsv --traits traits.tsv --alignment alignment.fasta --reported-taxa reviewer-table.tsv --out artifacts/taxonomy-report.html --json
bijux-phylogenetics taxonomy synonyms tree.nwk --synonym-table synonyms.tsv --json
bijux-phylogenetics taxonomy resolve-synonyms tree.nwk --synonym-table synonyms.tsv --out normalized-tree.nwk --mapping-out synonym-map.tsv --json
bijux-phylogenetics taxonomy namespaces tree.nwk --json
bijux-phylogenetics taxonomy loss tree.nwk --metadata metadata.csv --traits traits.csv --alignment alignment.fasta --filtered-alignment filtered.fasta --inference-tree inferred.nwk --reported-taxa reported.csv --json
bijux-phylogenetics taxonomy stability --run tree=tree.nwk --run alignment=alignment.fasta --run filtered=filtered.fasta --json
bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative brownian tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative compare-models tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json
bijux-phylogenetics comparative pgls tree.nwk traits.tsv --formula "height_cm ~ body_mass * habitat" --json
bijux-phylogenetics comparative covariance-audit tree.nwk traits.tsv --analysis pgls --formula "height_cm ~ body_mass + habitat" --summary-out artifacts/covariance-audit-summary.tsv --candidates-out artifacts/covariance-audit-candidates.tsv --excluded-taxa-out artifacts/covariance-audit-excluded.tsv --json
bijux-phylogenetics comparative multiple-testing tree.nwk traits.tsv --responses height_cm range_km --predictors body_mass log_range --json
bijux-phylogenetics comparative report tree.nwk traits.tsv --formula "height_cm ~ body_mass + habitat" --out artifacts/comparative-report.html --json
bijux-phylogenetics comparative compare-trees tree-a.nwk tree-b.nwk traits.tsv --response height_cm --predictors body_mass log_range --json
bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --json
bijux-phylogenetics ancestral discrete tree.nwk traits.tsv --trait habitat --model equal-rates --root-prior-mode empirical --summary-out artifacts/ancestral-discrete-summary.tsv --probabilities-out artifacts/ancestral-discrete-probabilities.tsv --transitions-out artifacts/ancestral-discrete-transitions.tsv --json
bijux-phylogenetics ancestral sensitivity tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --json
bijux-phylogenetics ancestral report tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --out artifacts/ancestral-report.html
bijux-phylogenetics ancestral package tree.nwk traits.tsv --trait habitat --kind discrete --model symmetric --state-ordering ordered --ordered-states low,medium,high --out-dir artifacts/ancestral-package --json
bijux-phylogenetics discrete-evolution model tree.nwk geography.tsv --trait region --model symmetric --state-ordering ordered --ordered-states north,south,island --node-table-out artifacts/node-states.tsv --transitions-out artifacts/transitions.tsv --json
bijux-phylogenetics discrete-evolution stochastic-map tree.nwk geography.tsv --trait region --model symmetric --replicates 200 --collection-out artifacts/geography-maps.json --summary-out artifacts/geography-stochastic-summary.tsv --json
bijux-phylogenetics discrete-evolution summarize-maps artifacts/geography-maps.json --summary-out artifacts/geography-stochastic-summary.tsv --json
bijux-phylogenetics discrete-evolution report tree.nwk geography.tsv --trait region --compare-model symmetric --out artifacts/geography-report.html
bijux-phylogenetics diversification estimate tree.nwk --metadata sampling.tsv --model birth-death --json
bijux-phylogenetics diversification report tree.nwk --metadata sampling.tsv --traits traits.tsv --trait habitat --out artifacts/diversification-report.html
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
bijux-phylogenetics tree-set report posterior.trees --out artifacts/tree-uncertainty-report.html --max-tree-count 5000 --max-report-table-rows 50 --memory-warning-threshold-bytes 134217728
bijux-phylogenetics demo gnathostome-ortholog-protein-benchmark --out artifacts/gnathostome-ortholog-protein-benchmark --json
bijux-phylogenetics demo rabies-method-sensitivity-panel --out artifacts/rabies-method-sensitivity-panel --json
bijux-phylogenetics demo rabies-cross-host-geography-panel --out artifacts/rabies-cross-host-geography-panel --config src/bijux_phylogenetics/resources/datasets/pathogens/rabies_cross_host_geography_panel/workflow-config.json --json
bijux-phylogenetics demo catarrhine-data-quality-stress-panel --out artifacts/catarrhine-data-quality-stress-panel --json
bijux-phylogenetics demo known-answer-reference-panel --out artifacts/known-answer-reference-panel --json
bijux-phylogenetics simulate tree-birth-death --tree-count 5 --tip-count 16 --out simulated.trees
bijux-phylogenetics simulate traits-early-burst tree.nwk --root-state 1.0 --sigma 0.5 --rate-change 4.0 --out simulated-early-burst.tsv --json
bijux-phylogenetics demo continuous-mode-recovery-panel --out artifacts/continuous-mode-recovery-panel --json
bijux-phylogenetics simulate alignment-dna tree.nwk --sequence-length 500 --out simulated-alignment.fasta
bijux-phylogenetics benchmark tree-comparison --replicates 3 --json
bijux-phylogenetics diagnose assumptions tree.nwk --metadata metadata.tsv --json
bijux-phylogenetics alignment translate coding.fasta --out translated.fasta --codon-validation-out artifacts/codon-validation.tsv --excluded-sequences-out artifacts/translation-exclusions.tsv
bijux-phylogenetics report dataset tree.nwk metadata.tsv traits.tsv --alignment alignment.fasta --tip-dates tip-dates.tsv --calibrations calibrations.tsv --out artifacts/dataset-report.html --json
bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk
bijux-phylogenetics phylo preflight --workflow fasta-to-tree --json
bijux-phylogenetics phylo run workflow-config.yaml --json
```

`demo rabies-cross-host-geography-panel` is the repository's flagship public
biological workflow surface. In addition to the dataset and workflow
subdirectories, the demo now writes one reviewer-facing package layer at the
output root:

- `dataset/source-accessions.tsv` for machine-readable accession provenance
- `rabies-cross-host-geography-overview.html` for one direct public handoff
- `rabies-cross-host-geography-package.manifest.json` for the biological
  question, short answer, config provenance, output checksums, and key metrics

Its packaged `workflow-config.json` now also carries explicit resource-budget
controls alongside the scientific settings: `iqtree_threads`,
`timeout_seconds`, `max_bootstrap_tree_count`, `max_report_table_rows`, and
`memory_warning_threshold_bytes`. The runtime records observed workflow and
bootstrap-review runtime or memory where available and emits structured budget
failures or warnings instead of silently overrunning those limits.

Its JSON metrics also surface the same `biological_question` and
`short_answer` directly so reviewers do not need to open the nested HTML report
to understand the intended scientific claim.

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
Its reviewer-facing HTML report is now intentionally compact: it surfaces one
summary card set plus explicit links to the governed TSV and JSON ledgers
instead of embedding those tables directly, and the bundle now includes
`workflow/report-artifacts/rabies-method-sensitivity-report.manifest.json`
with linked-artifact checksums and byte counts. The JSON metrics for the demo
also report `report_linked_artifact_count`, `report_html_size_bytes`,
`report_linked_artifact_bytes`, and `report_total_output_bytes`.

`tree-set report` now follows the same scaling contract. The HTML report keeps
top-level uncertainty summaries in the page body, writes large reviewer tables
to one sibling `*.artifacts/` directory as TSV or JSON, links those artifacts
explicitly, records report mode as either `full-review` or `scaled-summary`,
and reports `html_size_bytes`, `linked_artifact_bytes`, and
`total_output_bytes`. Tree sets with `1,000+` trees stay reviewable by
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

`comparative signal` keeps its input policy explicit in JSON output: rooted
trees are accepted whether or not they are ultrametric, overlapping missing
trait values are pruned and reported, permutation rows are reproducible from
the supplied seed, and constant post-pruning trait vectors fail with a typed
comparative-method error instead of returning misleading scalar summaries.

The BEAST adapter surface now makes its evidence state explicit. `adapter
beast-prepare` only prepares XML, `adapter beast-log`, `beast-trees`, and
`beast-consensus` parse existing posterior outputs, and reviewer-facing
diagnostics such as `adapter bayesian-methods` state when they are only
summarizing a prepared XML versus parsing existing log/tree outputs. When a
matching `analysis.manifest.json` from `adapter beast-run` is present, those
diagnostics identify the posterior log and tree set as outputs from a recorded
prior BEAST inference rather than implying the report executed BEAST itself.

The Bayesian runtime controls are intentionally strict. `adapter beast-run`
and `adapter mrbayes-run` leave an explicit `.incomplete.json` marker not only
for timeouts and nonzero exits but also when the engine exits yet the emitted
posterior files fail validation. `--resume` reuses only one verified completed
manifest, `--incomplete-run-policy clean` is the governed way to discard that
partial state, and a missing executable stops before any incomplete-run marker
is written because no engine run started.

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

The external-engine trust surface now has two distinct verification lanes.
Fast `engine_contract` tests keep fake-executable and parser behavior stable in
routine verification, while `tests/real_local` carries the governed
`engine_real` lane for installed MAFFT, trimAl, IQ-TREE, FastTree, and
MrBayes executables plus the checked-in real BEAST XML/log/tree corpus. Those
real-local tests write one external validation matrix JSON artifact per lane
that records reviewer-facing engine names, validation modes, executable paths,
version text, commands, exit codes, runtime, output paths, and output hashes.

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
scientific result for that workflow surface.

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
fixture ids. The tree, trait-table, and DNA inputs for that lane now come from
the governed shared fixture catalogs in
`tests/fixtures/metadata/shared_tree_fixture_catalog.json`,
`tests/fixtures/metadata/shared_trait_table_fixture_catalog.json`, and
`tests/fixtures/metadata/shared_dna_alignment_fixture_catalog.json`, so Bijux
and `ape` resolve the same durable fixture identities instead of hand-picked
path lists. The `ape::read.tree` lane now checks structured clade rows rather
than only raw parse success, including rooted and unrooted trees, branch
lengths, internal node labels, support labels, quoted labels, one governed
multiple-tree input, and one governed malformed-Newick rejection case. The
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
instead of silently degrading them to missing data.
The `ape::nj` lane now covers one governed analytical three-taxon matrix plus
four-taxon ultrametric and non-ultrametric matrices. On the owned Bijux side,
neighbor joining no longer delegates through Biopython for that method: Bijux
now builds one deterministic NJ tree in-repo, validates zero-diagonal and
nonnegative matrix assumptions explicitly, produces branch lengths, and
resolves tied joins by stable taxon ordering so distance-tree recovery stays
reproducible.
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
live `ape::pic`.
The `ape::ace` continuous lane now covers balanced rooted ultrametric,
pectinate rooted non-ultrametric, six-taxon Brownian, and pruned missing-value
fixtures through that same governed trait-table catalog. On the owned Bijux
side, `ancestral continuous` now carries one explicit Brownian fit diagnostic
surface with ultrametric state, root-to-tip depth bounds, covariance rank and
conditioning, solver regularization status, and one GLS likelihood summary,
while the live parity lane remains scoped honestly to
`ape::ace(type='continuous', method='pic', CI=TRUE)` because that is the
governed shared closed-form Brownian reference surface.
The `ape::ace` discrete lane now covers governed ER, SYM, and ARD fixtures,
including balanced, pectinate, six-taxon, and pruned missing-value cases
through that same shared trait-table catalog. On the owned Bijux side,
`ancestral discrete` now emits fitted transition-rate tables, log-likelihood,
parameter count, AIC, weak-fit warnings, and ER baseline-comparison data
alongside node probabilities, and it supports owned `equal`, `empirical`, and
`fixed` root-prior policies. The live parity lane is scoped honestly to
`ape::ace(type='discrete', model='ER'|'SYM'|'ARD')`, so root-prior controls
remain an explicit Bijux-owned review surface rather than a false live `ape`
parity claim.
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
explicitly on the Bijux side. The
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

## Read this next

- package docs: [Runtime package docs](https://bijux.io/bijux-phylogenetics/public/phylogenetics/)
- source directory: [Runtime source directory](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)
- changelog: [Runtime package changelog](https://github.com/bijux/bijux-phylogenetics/blob/main/packages/bijux-phylogenetics/CHANGELOG.md)
- security policy: [Security policy](https://github.com/bijux/bijux-phylogenetics/blob/main/SECURITY.md)
