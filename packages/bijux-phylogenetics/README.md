# bijux-phylogenetics

<!-- bijux-phylogenetics-badges:generated:start -->
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-3776AB?logo=python&logoColor=white)](https://pypi.org/project/bijux-phylogenetics/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-0F766E)](https://github.com/bijux/bijux-phylogenetics/blob/main/LICENSE)
[![Verify](https://github.com/bijux/bijux-phylogenetics/workflows/repo%20/%20verify/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/verify.yml?query=branch%3Amain)
[![Release PyPI](https://github.com/bijux/bijux-phylogenetics/workflows/release-pypi/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-pypi.yml)
[![Release GHCR](https://github.com/bijux/bijux-phylogenetics/workflows/release-ghcr/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-ghcr.yml)
[![Release GitHub](https://github.com/bijux/bijux-phylogenetics/workflows/release-github/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/release-github.yml)
[![Docs](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/bijux/bijux-phylogenetics/actions/workflows/deploy-docs.yml)

[![bijux-phylogenetics](https://img.shields.io/pypi/v/bijux-phylogenetics?label=bijux--phylogenetics&logo=pypi)](https://pypi.org/project/bijux-phylogenetics/)
[![phylogenetics](https://img.shields.io/pypi/v/phylogenetics?label=phylogenetics&logo=pypi)](https://pypi.org/project/phylogenetics/)

[![bijux-phylogenetics](https://img.shields.io/badge/bijux--phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fbijux-phylogenetics)
[![phylogenetics](https://img.shields.io/badge/phylogenetics-ghcr-181717?logo=github)](https://github.com/bijux/bijux-phylogenetics/pkgs/container/bijux-phylogenetics%2Fphylogenetics)

[![bijux-phylogenetics docs](https://img.shields.io/badge/docs-bijux--phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
[![phylogenetics docs](https://img.shields.io/badge/docs-phylogenetics-2563EB?logo=materialformkdocs&logoColor=white)](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
<!-- bijux-phylogenetics-badges:generated:end -->

Runtime package for the bijux-phylogenetics repository.

This package provides the Python API and CLI for tree validation, inspection,
comparison, metadata linkage, evidence bundle creation, and HTML report
generation.

## Install

`bijux-phylogenetics` supports Python 3.11 and newer.

```bash
python3.11 -m pip install bijux-phylogenetics
bijux-phylogenetics --help
```

## Current Scope

- parse Newick trees and FASTA alignments
- inspect tree shape and branch-length health
- inspect internal child counts, singleton nodes, missing internal versus terminal branch lengths, branch-length outlier nodes, support normalization, and tree-assumption compatibility
- normalize unsafe taxon labels and audit normalization collisions
- prune trees from explicit taxa, exclusions, traits, or metadata tables
- classify internal node labels as support-like or name-like and detect suspicious or mixed support scales
- compare shared clades, clade changes, and shared-split branch lengths between trees
- validate trait and metadata linkage against tree tips
- export joined metadata rows and missing trait-value diagnostics
- inspect alignment alphabets, composition, GC content, duplicates, composition outliers, coding stop codons, and frameshift-like sequence lengths
- trim all-gap or all-missing columns and remove high-missingness sequences
- translate coding nucleotide alignments to amino-acid alignments and export pairwise identity matrices
- compute p-distance or Jukes-Cantor DNA distance matrices with pairwise-deletion or complete-deletion gap handling
- build Neighbor-Joining or UPGMA trees from computed DNA distance matrices and compare their topologies
- validate imported long-form distance matrices, detect nonmetric violations, and build trees from imported distances
- load posterior tree sets, compute consensus trees, and export clade-frequency or pairwise tree-distance summaries
- cluster identical rooted topologies, detect unstable taxa or clades, and compare two posterior tree sets
- simulate birth-death or coalescent trees, Brownian or OU continuous traits, discrete traits, and DNA or protein alignments
- benchmark validation, tree comparison, and alignment diagnostics across increasing problem sizes
- root trees on explicit outgroups or reroot them by midpoint
- produce HTML reports and file-level evidence manifests

## Example CLI Runs

```bash
bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4
bijux-phylogenetics alignment distance-matrix alignment.fasta --model p-distance --out distances.tsv
bijux-phylogenetics alignment build-tree alignment.fasta --method upgma --out upgma-tree.nwk
bijux-phylogenetics distance validate distances.tsv --json
bijux-phylogenetics tree-set inspect posterior.trees --json
bijux-phylogenetics tree-set consensus posterior.trees --out consensus.nwk
bijux-phylogenetics tree-set report posterior.trees --out artifacts/tree-uncertainty-report.html
bijux-phylogenetics simulate tree-birth-death --tree-count 5 --tip-count 16 --out simulated.trees
bijux-phylogenetics simulate alignment-dna tree.nwk --sequence-length 500 --out simulated-alignment.fasta
bijux-phylogenetics benchmark tree-comparison --replicates 3 --json
bijux-phylogenetics diagnose assumptions tree.nwk --metadata metadata.tsv --json
bijux-phylogenetics alignment translate coding.fasta --out translated.fasta
bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk
```

## Read this next

- package docs: [Runtime package docs](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- source directory: [Runtime source directory](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)
- changelog: [Runtime package changelog](https://github.com/bijux/bijux-phylogenetics/blob/main/packages/bijux-phylogenetics/CHANGELOG.md)
- security policy: [Security policy](https://github.com/bijux/bijux-phylogenetics/blob/main/SECURITY.md)
