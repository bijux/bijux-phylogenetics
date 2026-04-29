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

## Current Scope

- parse Newick trees and FASTA alignments
- inspect tree shape and branch-length health
- inspect internal child counts, singleton nodes, missing internal versus terminal branch lengths, branch-length outlier nodes, support normalization, rootedness confidence, and tree-assumption compatibility
- classify tree validity, biological safety, unsafe external labels, node-label conflicts, and downstream forensic readiness for topology, time-tree, comparative, visualization, and publication use
- normalize unsafe taxon labels and audit normalization collisions
- prune trees from explicit taxa, exclusions, traits, or metadata tables
- classify internal node labels as support-like or name-like and detect suspicious or mixed support scales
- compare shared clades, clade changes, and shared-split branch lengths between trees
- validate trait and metadata linkage against tree tips
- check comparative readiness for rooted trees and numeric traits
- compute phylogenetic independent contrasts, Blomberg's K, Pagel's lambda, and permutation-based signal tests
- fit phylogenetic generalized least-squares models with one or more numeric predictors
- reconstruct continuous ancestral states under Brownian or OU-style trait models
- reconstruct discrete ancestral states under Fitch parsimony with explicit ambiguity reporting
- compare continuous ancestral reconstructions across two supported models and render annotated ancestral trees
- validate discrete geographic state coding, detect state imbalance, estimate ancestral node states, compare equal-rates and all-rates-different models, export node and transition tables, and render discrete-state HTML reports
- estimate lineage-through-time curves, simple Yule or birth-death diversification rates, sampling-aware corrections, clade outlier summaries, and trait-linked diversification tables for rooted ultrametric trees
- run governed MAFFT-, trimAl-, IQ-TREE-, and FastTree-style external workflows with captured commands, versions, logs, and warning summaries
- prepare and run deterministic MrBayes analyses, summarize posterior trees after burn-in filtering, parse parameter traces, and compute per-parameter ESS values
- compare fast approximate and maximum-likelihood trees through the same deterministic tree-comparison report surface
- resume inference only when saved manifests, inputs, and outputs still match, and render standalone HTML inference workflow reports from those manifests
- export joined metadata rows and missing trait-value diagnostics
- inspect alignment alphabets, composition, GC content, duplicates, raw-sequence length outliers, sliding-window quality, suspicious alignment regions, coding stop codons, and frameshift-like sequence lengths
- classify FASTA inputs as aligned, raw-sequence, or equal-length-but-shape-ambiguous and report method-specific alignment readiness
- detect mixed coding versus noncoding behavior inside the same nucleotide dataset
- define named alignment-filtering profiles, generate cleaned alignments, compare original versus cleaned versions, and warn when filtering removes signal or biases taxon groups
- score alignment quality with transparent components and emit one-shot alignment forensic reports
- audit tree, metadata, traits, alignment, tip dates, and calibrations together through one-shot dataset readiness decisions
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
- audit rooting, ordering, clade extraction, and pruning transforms with before/after summaries and retained-versus-removed taxon reasoning
- validate tree roundtrips across Newick, Nexus, and phyloXML formats with topology-preservation checks, support-label audits, and semantic-loss warnings
- audit ambiguous taxon identities, whitespace or underscore collisions, and suspicious near-duplicate labels before downstream comparison or linkage
- produce HTML reports and file-level evidence manifests

## Example CLI Runs

```bash
bijux-phylogenetics alignment classify sequences.fasta --json
bijux-phylogenetics alignment profiles --json
bijux-phylogenetics alignment windows alignment.fasta --window-size 50 --step-size 10 --json
bijux-phylogenetics alignment readiness alignment.fasta --json
bijux-phylogenetics alignment length-outliers sequences.fasta --json
bijux-phylogenetics alignment forensic alignment.fasta --json
bijux-phylogenetics alignment filter alignment.fasta --profile moderate --out cleaned.fasta --json
bijux-phylogenetics alignment compare alignment.fasta cleaned.fasta --json
bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4
bijux-phylogenetics alignment distance-matrix alignment.fasta --model p-distance --out distances.tsv
bijux-phylogenetics alignment build-tree alignment.fasta --method upgma --out upgma-tree.nwk
bijux-phylogenetics distance validate distances.tsv --json
bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json
bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json
bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --json
bijux-phylogenetics ancestral discrete tree.nwk traits.tsv --trait habitat --json
bijux-phylogenetics ancestral report tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --out artifacts/ancestral-report.html
bijux-phylogenetics discrete-evolution model tree.nwk geography.tsv --trait region --node-table-out artifacts/node-states.tsv --transitions-out artifacts/transitions.tsv --json
bijux-phylogenetics discrete-evolution report tree.nwk geography.tsv --trait region --compare-model all-rates-different --out artifacts/geography-report.html
bijux-phylogenetics diversification estimate tree.nwk --metadata sampling.tsv --model birth-death --json
bijux-phylogenetics diversification report tree.nwk --metadata sampling.tsv --traits traits.tsv --trait habitat --out artifacts/diversification-report.html
bijux-phylogenetics adapter align unaligned.fasta --out aligned.fasta --json
bijux-phylogenetics adapter model-select alignment.fasta --out-dir artifacts/model-select --prefix mammals --json
bijux-phylogenetics adapter infer-ml alignment.fasta --out-dir artifacts/ml --model GTR+G --prefix mammals --json
bijux-phylogenetics adapter bootstrap alignment.fasta --out-dir artifacts/bootstrap --model GTR+G --replicates 1000 --prefix mammals --json
bijux-phylogenetics adapter consensus artifacts/bootstrap/mammals.ufboot --out-dir artifacts/consensus --prefix mammals --json
bijux-phylogenetics adapter infer-fast alignment.fasta --out artifacts/fasttree.nwk --json
bijux-phylogenetics adapter compare --fast-tree artifacts/fasttree.nwk --ml-tree artifacts/ml/mammals.treefile --out artifacts/engine-comparison.html --json
bijux-phylogenetics adapter mrbayes-prepare alignment.fasta --out artifacts/mrbayes/analysis.nex --ngen 20000 --samplefreq 100 --json
bijux-phylogenetics adapter mrbayes-run artifacts/mrbayes/analysis.nex --resume --json
bijux-phylogenetics adapter mrbayes-summarize artifacts/mrbayes/analysis.run1.t --burnin-fraction 0.25 --json
bijux-phylogenetics adapter mrbayes-traces artifacts/mrbayes/analysis.run1.p --json
bijux-phylogenetics adapter mrbayes-ess artifacts/mrbayes/analysis.run1.p --json
bijux-phylogenetics adapter report artifacts/mrbayes/analysis.manifest.json --out artifacts/mrbayes/inference-report.html --json
bijux-phylogenetics tree-set inspect posterior.trees --json
bijux-phylogenetics tree-set consensus posterior.trees --out consensus.nwk
bijux-phylogenetics tree-set report posterior.trees --out artifacts/tree-uncertainty-report.html
bijux-phylogenetics simulate tree-birth-death --tree-count 5 --tip-count 16 --out simulated.trees
bijux-phylogenetics simulate alignment-dna tree.nwk --sequence-length 500 --out simulated-alignment.fasta
bijux-phylogenetics benchmark tree-comparison --replicates 3 --json
bijux-phylogenetics diagnose assumptions tree.nwk --metadata metadata.tsv --json
bijux-phylogenetics alignment translate coding.fasta --out translated.fasta
bijux-phylogenetics report dataset tree.nwk metadata.tsv traits.tsv --alignment alignment.fasta --tip-dates tip-dates.tsv --calibrations calibrations.tsv --out artifacts/dataset-report.html --json
bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk
```

## Alignment Filter Profiles

The built-in alignment filtering profiles are `conservative`, `moderate`,
`aggressive`, `coding-safe`, and `phylogenomics-scale`.

## Read this next

- package docs: [Runtime package docs](https://bijux.io/bijux-phylogenetics/01-bijux-phylogenetics/)
- source directory: [Runtime source directory](https://github.com/bijux/bijux-phylogenetics/tree/main/packages/bijux-phylogenetics)
- changelog: [Runtime package changelog](https://github.com/bijux/bijux-phylogenetics/blob/main/packages/bijux-phylogenetics/CHANGELOG.md)
- security policy: [Security policy](https://github.com/bijux/bijux-phylogenetics/blob/main/SECURITY.md)
