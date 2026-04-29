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
alignment trimming, coding-sequence translation, identity-matrix export,
DNA distance-matrix analysis, distance-tree construction, explicit rooting
transforms, comparative trait readiness, phylogenetic independent contrasts,
phylogenetic signal estimation, phylogenetic generalized least-squares,
tree-set consensus and posterior uncertainty analysis, tree and alignment
simulation, scientific benchmarking, deterministic SVG tree rendering,
publication figure packaging, evidence manifests, and HTML report generation
rather than likelihood or Bayesian tree inference.

Recent tree diagnostics now also classify internal-node child counts, missing
internal versus terminal branch lengths, singleton internal nodes, branch-length
outlier nodes, support-like versus name-like internal labels, metadata-declared
branch-length units, and explicit time-tree versus substitution-tree
compatibility assumptions.

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

- `uv run bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4`
- `uv run bijux-phylogenetics alignment coding coding-alignment.fasta --json`
- `uv run bijux-phylogenetics alignment translate coding-alignment.fasta --out translated.fasta`
- `uv run bijux-phylogenetics alignment identity-matrix alignment.fasta --out identity.tsv`
- `uv run bijux-phylogenetics alignment distance-matrix alignment.fasta --model jukes-cantor --gap-handling complete-deletion --out distances.tsv`
- `uv run bijux-phylogenetics alignment build-tree alignment.fasta --method neighbor-joining --out nj-tree.nwk`
- `uv run bijux-phylogenetics alignment compare-distance-trees alignment.fasta --json`
- `uv run bijux-phylogenetics distance validate exported-distances.tsv --json`
- `uv run bijux-phylogenetics distance build-tree exported-distances.tsv --method upgma --out imported-upgma.nwk`
- `uv run bijux-phylogenetics distance report exported-distances.tsv --out artifacts/distance-report.html`
- `uv run bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative contrasts tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json`
- `uv run bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json`
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
