# bijux-phylogenetics

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/bijux-phylogenetics/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/workflows/repo%20/%20verify/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://github.com/bijux/bijux-phylogenetics/workflows/release-pypi/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://github.com/bijux/bijux-phylogenetics/workflows/release-ghcr/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://github.com/bijux/bijux-phylogenetics/workflows/release-github/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)
[![Release](https://img.shields.io/github/v/release/bijux/bijux-phylogenetics?display_name=tag&label=release)](https://github.com/bijux/bijux-phylogenetics/releases)
[![GHCR packages](https://img.shields.io/badge/ghcr-2%20packages-181717?logo=github)](https://github.com/bijux?tab=packages&repo_name=bijux-phylogenetics)
[![Published packages](https://img.shields.io/badge/published%20packages-2-2563EB)](https://github.com/bijux/bijux-phylogenetics/tree/main/packages)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)
[![phylogenetics](https://img.shields.io/pypi/v/phylogenetics?label=phylogenetics&logo=pypi)](https://pypi.org/project/phylogenetics/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)
[![phylogenetics](https://img.shields.io/badge/phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetics)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![phylogenetics docs](https://img.shields.io/badge/docs-phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

`bijux-phylogenetics` is a reproducible phylogenetics workbench for tree
inspection, validation, comparison, metadata linkage, evidence capture, and
publishable reporting.

## Repository Layout

The repository keeps these durable top-level surfaces:

- `packages/` for published runtime, alias, and maintainer packages
- `docs/` for documentation source
- `examples/` for tracked workflow examples
- `datasets/` for durable repository-owned reference inputs
- `reports/` for tracked proof outputs and reviewer-facing artifacts
- `tests/` for cross-package test assets as the repository grows

The repository follows the same Bijux Python workspace pattern used by
`bijux-canon`, `bijux-proteomics`, and `bijux-pollenomics`: shared `.bijux`
assets, standardized `makes/` and `configs/`, a root `uv` workspace, a primary
runtime package, a compatibility alias package, and a repository-owned
maintainer package.




This repository publishes `2` public packages. The current runtime focuses on
reproducible tree hygiene, taxon normalization, tree and trait pruning,
metadata linkage, shared-clade comparison, alignment-quality diagnostics,
alignment readiness classification, raw-versus-aligned FASTA detection,
sliding-window region diagnostics, ambiguity-aware uncertainty accounting,
alignment filtering profiles, cleaned-alignment comparison, transparent
alignment quality scoring, alignment forensic reporting, one-shot dataset
readiness audits, dataset crosswalk tables, completeness matrices, exclusion
tables, ordering-drift detection, pruning step retention summaries,
readiness-level reporting, alignment trimming, coding-sequence translation,
identity-matrix export, DNA distance-matrix analysis, distance-tree
construction, explicit rooting transforms, comparative trait readiness,
phylogenetic independent contrasts, phylogenetic signal estimation,
standalone Brownian and OU trait-model fitting, comparative-model
comparison, formula-driven phylogenetic generalized least-squares,
comparative multiple-testing correction, comparative audit and influence
reporting, continuous and discrete
ancestral-state reconstruction, ancestral uncertainty reporting, ancestral sensitivity
summaries, supplement-style ancestral reporting, ancestral tree
rendering, ordered versus unordered discrete-state modeling, discrete-state geographic transition modeling and reporting,
stochastic-map approximation and uncertainty summarization for discrete-state
transitions, model-sensitive ancestral-region comparison,
governed external-engine orchestration for alignment, trimming, model
selection, tree inference, bootstrap retention, inference-readiness auditing,
model-selection validation, inferred-tree taxon validation, metadata-group
clustering audits, and initial MrBayes posterior analysis, tree-set
consensus and posterior uncertainty analysis, tree and alignment simulation,
scientific benchmarking, deterministic SVG tree rendering, publication figure
packaging, evidence manifests, and HTML report generation rather than full
likelihood or Bayesian inference engines.
Distance workflows now also validate built-in reference examples, support
Kimura 2-parameter and amino-acid p-distance models, handle ambiguity codes
through explicit policies, report saturation and low-information diagnostics,
bootstrap site-resampled distance trees, summarize clade support, and write
reproducibility bundles.
Diversification and macroevolution workflows now also estimate lineage-through-time
curves, simple Yule or birth-death rates, sampling-aware corrections, clade
outlier summaries, and trait-linked diversification tables for rooted
ultrametric trees.

Recent tree diagnostics now also classify internal-node child counts, missing
internal versus terminal branch lengths, singleton internal nodes, branch-length
outlier nodes, support-like versus name-like internal labels, metadata-declared
branch-length units, and explicit time-tree versus substitution-tree
compatibility assumptions.

Taxon auditing now also classifies mixed naming levels across species, genus,
sample, accession, and population-style labels, exports accepted-name mappings
through synonym tables, flags duplicate biological identity candidates, and
produces reviewer-readable taxon conflict audits.

Dataset auditing now also writes cross-surface mismatch tables, transparent
risk components, minimal-fix plans, and reviewer checklists that summarize what
still blocks downstream analyses.

Tree rendering now supports rectangular cladograms, rectangular phylograms with
scale bars, circular trees, support-value labels, categorical and continuous tip
traits, collapsed named clades, metadata strips, trait heatmaps, and
publication-style figure bundles.

## Start Here

- read the docs home: [Documentation home](https://bijux.io/bijux-phylogenetics/)
- inspect the runtime package source: [packages/bijux-phylogenetics](packages/bijux-phylogenetics)
- inspect the alias package source: [packages/phylogenetics](packages/phylogenetics)
- inspect maintainer tooling: [packages/bijux-phylogenetics-dev](packages/bijux-phylogenetics-dev)

## What This Repository Produces

Today, the checked-in repository produces these durable outcomes:

- a Python runtime package for phylogenetic validation, inspection, comparison, and reporting
- a compatibility alias distribution for the shorter `phylogenetics` command
- a repository-owned maintainer package for docs, release, and quality automation
- a MkDocs documentation site that builds into `artifacts/root/docs/site/`

## Common Workflows

- `make install` syncs the editable environment from the tracked `uv.lock`
- `make check` runs the main repository verification pass: lock check, lint, tests, docs, and distribution verification
- `make docs-serve` serves the docs locally at `http://127.0.0.1:8000/`
- `make package-verify` runs wheel, sdist, and smoke-install package proof targets
- `make sync-badges` renders the shared badge catalog into managed README surfaces

## Example Commands

- `uv run bijux-phylogenetics alignment classify sequences.fasta --json`
- `uv run bijux-phylogenetics alignment profiles --json`
- `uv run bijux-phylogenetics alignment windows alignment.fasta --window-size 50 --step-size 10 --json`
- `uv run bijux-phylogenetics alignment readiness alignment.fasta --json`
- `uv run bijux-phylogenetics alignment low-information alignment.fasta --json`
- `uv run bijux-phylogenetics alignment duplicate-policy alignment.fasta --identity-threshold 0.99 --json`
- `uv run bijux-phylogenetics alignment ambiguous-columns alignment.fasta --threshold 0.5 --json`
- `uv run bijux-phylogenetics alignment sequence-ranking alignment.fasta --json`
- `uv run bijux-phylogenetics alignment length-outliers sequences.fasta --json`
- `uv run bijux-phylogenetics alignment forensic alignment.fasta --json`
- `uv run bijux-phylogenetics alignment filter alignment.fasta --profile coding-safe --out cleaned.fasta --json`
- `uv run bijux-phylogenetics alignment compare original.fasta cleaned.fasta --json`
- `uv run bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4`
- `uv run bijux-phylogenetics alignment coding coding-alignment.fasta --json`
- `uv run bijux-phylogenetics alignment translate coding-alignment.fasta --out translated.fasta`
- `uv run bijux-phylogenetics alignment identity-matrix alignment.fasta --out identity.tsv`
- `uv run bijux-phylogenetics alignment distance-matrix alignment.fasta --model jukes-cantor --gap-handling complete-deletion --out distances.tsv`
- `uv run bijux-phylogenetics alignment distance-matrix proteins.fasta --model amino-acid-p-distance --ambiguity-policy partial-match --out protein-distances.tsv`
- `uv run bijux-phylogenetics alignment distance-quality alignment.fasta --model kimura-2-parameter --json`
- `uv run bijux-phylogenetics alignment build-tree alignment.fasta --method neighbor-joining --out nj-tree.nwk`
- `uv run bijux-phylogenetics alignment compare-distance-trees alignment.fasta --json`
- `uv run bijux-phylogenetics alignment bootstrap-tree alignment.fasta --method neighbor-joining --replicates 200 --support-out artifacts/distance-support.tsv --tree-set-out artifacts/distance-bootstrap.trees --json`
- `uv run bijux-phylogenetics alignment distance-bundle alignment.fasta --method neighbor-joining --replicates 200 --out-dir artifacts/distance-bundle --json`
- `uv run bijux-phylogenetics distance validate exported-distances.tsv --json`
- `uv run bijux-phylogenetics distance build-tree exported-distances.tsv --method upgma --out imported-upgma.nwk`
- `uv run bijux-phylogenetics distance report exported-distances.tsv --out artifacts/distance-report.html`
- `uv run bijux-phylogenetics distance reference --json`
- `uv run bijux-phylogenetics report alignment --alignment alignment.fasta --out artifacts/alignment-report.html --json`
- `uv run bijux-phylogenetics report dataset tree.nwk metadata.tsv traits.tsv --alignment alignment.fasta --tip-dates tip-dates.tsv --calibrations calibrations.tsv --out artifacts/dataset-report.html --json`
- `uv run bijux-phylogenetics report dataset tree.nwk metadata.tsv traits.tsv --alignment alignment.fasta --out artifacts/dataset-review.html --json`
- `uv run bijux-phylogenetics report taxonomy --tree tree.nwk --synonym-table taxonomy.tsv --metadata metadata.tsv --traits traits.tsv --alignment alignment.fasta --reported-taxa reviewer-table.tsv --out artifacts/taxonomy-report.html --json`
- `uv run bijux-phylogenetics taxonomy rank-consistency tree.nwk --json`
- `uv run bijux-phylogenetics taxonomy accepted-names tree.nwk --synonym-table taxonomy.tsv --out accepted-names.tsv --json`
- `uv run bijux-phylogenetics taxonomy audit tree.nwk --synonym-table taxonomy.tsv --json`
- `uv run bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative contrasts tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative brownian tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative ou tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative compare-models tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json`
- `uv run bijux-phylogenetics comparative pgls tree.nwk traits.tsv --formula "height_cm ~ body_mass * habitat" --json`
- `uv run bijux-phylogenetics comparative multiple-testing tree.nwk traits.tsv --responses height_cm range_km --predictors body_mass log_range --json`
- `uv run bijux-phylogenetics comparative report tree.nwk traits.tsv --formula "height_cm ~ body_mass + habitat" --out artifacts/comparative-report.html --json`
- `uv run bijux-phylogenetics comparative influence tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json`
- `uv run bijux-phylogenetics comparative compare-trees tree-a.nwk tree-b.nwk traits.tsv --response height_cm --predictors body_mass log_range --json`
- `uv run bijux-phylogenetics comparative compare-pruning tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --drop-taxa OutlierTaxon --json`
- `uv run bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --json`
- `uv run bijux-phylogenetics ancestral discrete tree.nwk traits.tsv --trait habitat --model symmetric --state-ordering ordered --ordered-states low,medium,high --json`
- `uv run bijux-phylogenetics ancestral compare tree.nwk traits.tsv --trait height_cm --left-model brownian --right-model ou --json`
- `uv run bijux-phylogenetics ancestral sensitivity tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --json`
- `uv run bijux-phylogenetics ancestral report tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --out artifacts/ancestral-report.html`
- `uv run bijux-phylogenetics ancestral package tree.nwk traits.tsv --trait habitat --kind discrete --model symmetric --state-ordering ordered --ordered-states low,medium,high --out-dir artifacts/ancestral-package --json`
- `uv run bijux-phylogenetics discrete-evolution validate-coding tree.nwk geography.tsv --trait region --allowed-states north,south,island --json`
- `uv run bijux-phylogenetics discrete-evolution model tree.nwk geography.tsv --trait region --model symmetric --state-ordering ordered --ordered-states north,south,island --node-table-out artifacts/node-states.tsv --transitions-out artifacts/transitions.tsv --json`
- `uv run bijux-phylogenetics discrete-evolution stochastic-map tree.nwk geography.tsv --trait region --model symmetric --replicates 200 --collection-out artifacts/geography-maps.json --summary-out artifacts/geography-stochastic-summary.tsv --json`
- `uv run bijux-phylogenetics discrete-evolution summarize-maps artifacts/geography-maps.json --summary-out artifacts/geography-stochastic-summary.tsv --json`
- `uv run bijux-phylogenetics discrete-evolution render tree.nwk geography.tsv --trait region --out artifacts/geography.svg --json`
- `uv run bijux-phylogenetics discrete-evolution report tree.nwk geography.tsv --trait region --compare-model all-rates-different --out artifacts/geography-report.html --json`
- `uv run bijux-phylogenetics diversification ltt tree.nwk --out artifacts/ltt.tsv --json`
- `uv run bijux-phylogenetics diversification estimate tree.nwk --metadata sampling.tsv --model birth-death --json`
- `uv run bijux-phylogenetics diversification clades tree.nwk --out artifacts/clades.tsv --json`
- `uv run bijux-phylogenetics diversification report tree.nwk --metadata sampling.tsv --traits traits.tsv --trait habitat --out artifacts/diversification-report.html --json`
- `uv run bijux-phylogenetics adapter align unaligned.fasta --out aligned.fasta --json`
- `uv run bijux-phylogenetics adapter model-select alignment.fasta --out-dir artifacts/model-select --prefix mammals --json`
- `uv run bijux-phylogenetics adapter infer-ml alignment.fasta --out-dir artifacts/ml --model GTR+G --prefix mammals --json`
- `uv run bijux-phylogenetics adapter bootstrap alignment.fasta --out-dir artifacts/bootstrap --model GTR+G --replicates 1000 --prefix mammals --json`
- `uv run bijux-phylogenetics adapter consensus artifacts/bootstrap/mammals.ufboot --out-dir artifacts/consensus --prefix mammals --json`
- `uv run bijux-phylogenetics adapter infer-fast alignment.fasta --out artifacts/fasttree.nwk --json`
- `uv run bijux-phylogenetics adapter compare --fast-tree artifacts/fasttree.nwk --ml-tree artifacts/ml/mammals.treefile --out artifacts/engine-comparison.html --json`
- `uv run bijux-phylogenetics adapter mrbayes-prepare alignment.fasta --out artifacts/mrbayes/analysis.nex --ngen 20000 --samplefreq 100 --json`
- `uv run bijux-phylogenetics adapter mrbayes-run artifacts/mrbayes/analysis.nex --resume --json`
- `uv run bijux-phylogenetics adapter mrbayes-summarize artifacts/mrbayes/analysis.run1.t --burnin-fraction 0.25 --json`
- `uv run bijux-phylogenetics adapter mrbayes-traces artifacts/mrbayes/analysis.run1.p --json`
- `uv run bijux-phylogenetics adapter mrbayes-ess artifacts/mrbayes/analysis.run1.p --json`
- `uv run bijux-phylogenetics adapter mrbayes-convergence artifacts/mrbayes/analysis.run1.p --ess-threshold 200 --json`
- `uv run bijux-phylogenetics adapter mrbayes-report artifacts/mrbayes/analysis.run1.t --traces artifacts/mrbayes/analysis.run1.p --out artifacts/mrbayes/posterior-report.html --json`
- `uv run bijux-phylogenetics adapter beast-prepare alignment.fasta --out artifacts/beast/analysis.xml --tree tree.nwk --calibrations calibrations.tsv --tip-dates tip-dates.tsv --clock-model relaxed-lognormal --tree-prior birth-death --json`
- `uv run bijux-phylogenetics adapter beast-calibrations tree.nwk calibrations.tsv --json`
- `uv run bijux-phylogenetics adapter beast-tip-dates tree.nwk tip-dates.tsv --alignment alignment.fasta --json`
- `uv run bijux-phylogenetics adapter beast-log artifacts/beast/run.log --json`
- `uv run bijux-phylogenetics adapter beast-convergence artifacts/beast/run.log --ess-threshold 200 --json`
- `uv run bijux-phylogenetics adapter beast-calibration-report tree.nwk calibrations.tsv --tip-dates tip-dates.tsv --alignment alignment.fasta --out artifacts/beast/calibration-audit.html --json`
- `uv run bijux-phylogenetics adapter bayesian-evidence --out-dir artifacts/bayesian-bundle --inputs alignment.fasta calibrations.tsv tip-dates.tsv --configs artifacts/beast/analysis.xml --trees tree.nwk --logs artifacts/beast/run.log --diagnostics diagnostics.json --reports artifacts/beast/calibration-audit.html --json`
- `uv run bijux-phylogenetics adapter report artifacts/mrbayes/analysis.manifest.json --out artifacts/mrbayes/inference-report.html --json`
- `uv run bijux-phylogenetics tree-set inspect posterior.trees --json`
- `uv run bijux-phylogenetics tree-set consensus posterior.trees --out consensus.nwk`
- `uv run bijux-phylogenetics tree-set compare posterior-a.trees posterior-b.trees --json`
- `uv run bijux-phylogenetics tree-set report posterior.trees --out artifacts/tree-uncertainty-report.html`
- `uv run bijux-phylogenetics simulate tree-birth-death --tree-count 10 --tip-count 32 --out simulated.trees`
- `uv run bijux-phylogenetics simulate traits-brownian tree.nwk --sigma 0.5 --out simulated-traits.tsv`
- `uv run bijux-phylogenetics simulate alignment-dna tree.nwk --sequence-length 500 --out simulated-alignment.fasta`
- `uv run bijux-phylogenetics benchmark tree-validation --replicates 3 --json`
- `uv run bijux-phylogenetics diagnose assumptions tree.nwk --metadata metadata.tsv --json`
- `uv run bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk`
- `uv run bijux-phylogenetics topology reroot-midpoint tree.nwk --out midpoint-rooted.nwk`
- `uv run bijux-phylogenetics render tree.nwk --layout phylogram --support-labels --metadata metadata.tsv --label-column species --metadata-strip-columns location --traits traits.tsv --categorical-column habitat --continuous-column height_cm --heatmap-columns height_cm,status --package-dir artifacts/tree-figure --out artifacts/tree.svg`
- `uv run bijux-phylogenetics inspect tree-with-support.nwk --json`

## Working Rules


- treat `.bijux/shared/`, `makes/`, and `configs/` as managed standard surfaces
- keep generated outputs under `artifacts/` unless the task explicitly governs another tracked destination
- update README and docs claims only when the current repository state actually supports them
