---
title: Repository Overview
audience: public
type: handbook
status: active
owner: bijux-phylogenetics
last_reviewed: 2026-04-29
---

# Repository Overview

`bijux-phylogenetics` provides a governed Python surface for tree validation,
inspection, comparison, metadata linkage, alignment trimming, coding-sequence
translation, explicit rooting transforms, comparative trait analysis,
ancestral-state reconstruction, external engine orchestration, deterministic
tree rendering, publication figure packaging, evidence bundles, and HTML
report generation.

The repository intentionally does not reimplement inference engines. Its
current product surface is the reproducible orchestration and evidence layer
around trees, alignments, and trait tables.

## Current CLI Surface

- `bijux-phylogenetics validate tree.nwk`
- `bijux-phylogenetics inspect tree.nwk`
- `bijux-phylogenetics compare left.nwk right.nwk`
- `bijux-phylogenetics annotate tree.nwk --metadata traits.tsv`
- `bijux-phylogenetics render tree.nwk --layout phylogram --support-labels --metadata samples.tsv --label-column species --metadata-strip-columns location --traits traits.tsv --categorical-column habitat --continuous-column height_cm --heatmap-columns height_cm,status --package-dir artifacts/tree-figure --out artifacts/tree.svg`
- `bijux-phylogenetics bundle run/ --out evidence-pack/`
- `bijux-phylogenetics report --tree tree.nwk --alignment alignment.fasta --traits traits.tsv --metadata samples.tsv --out phylo-report.html`
- `bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4`
- `bijux-phylogenetics alignment coding coding.fasta --json`
- `bijux-phylogenetics alignment translate coding.fasta --out translated.fasta`
- `bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative contrasts tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json`
- `bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --json`
- `bijux-phylogenetics ancestral discrete tree.nwk traits.tsv --trait habitat --json`
- `bijux-phylogenetics ancestral compare tree.nwk traits.tsv --trait height_cm --left-model brownian --right-model ou --json`
- `bijux-phylogenetics ancestral report tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --out artifacts/ancestral-report.html`
- `bijux-phylogenetics adapter align unaligned.fasta --out aligned.fasta --json`
- `bijux-phylogenetics adapter model-select alignment.fasta --out-dir artifacts/model-select --prefix mammals --json`
- `bijux-phylogenetics adapter infer-ml alignment.fasta --out-dir artifacts/ml --model GTR+G --prefix mammals --json`
- `bijux-phylogenetics adapter bootstrap alignment.fasta --out-dir artifacts/bootstrap --model GTR+G --replicates 1000 --prefix mammals --json`
- `bijux-phylogenetics adapter consensus artifacts/bootstrap/mammals.ufboot --out-dir artifacts/consensus --prefix mammals --json`
- `bijux-phylogenetics adapter infer-fast alignment.fasta --out artifacts/fasttree.nwk --json`
- `bijux-phylogenetics adapter compare --fast-tree artifacts/fasttree.nwk --ml-tree artifacts/ml/mammals.treefile --out artifacts/engine-comparison.html --json`
- `bijux-phylogenetics alignment identity-matrix alignment.fasta --out identity.tsv`
- `bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk`
- `bijux-phylogenetics topology reroot-midpoint tree.nwk --out midpoint-rooted.nwk`
- `bijux-phylogenetics inspect tree-with-support.nwk --json`

## Tree Diagnostics Highlights

- internal-node child counts are reported explicitly for every internal node
- missing branch lengths are separated into internal-branch and terminal-branch diagnostics
- singleton internal nodes are detected instead of being silently treated as ordinary branching structure
- long and short nonzero branch outliers are reported as concrete affected nodes
- numeric internal labels are classified as likely support values, textual internal labels as likely clade names
- support-like labels are checked for suspicious ranges and mixed probability-versus-percentage scales

## Tree Figure Highlights

- rectangular cladograms, branch-length phylograms, and circular tree layouts are supported from one render surface
- phylogram output includes a branch-length scale bar instead of leaving absolute lengths visually ambiguous
- internal support values can be rendered directly on branch junctions for review-ready trees
- tip annotations can combine categorical markers, continuous bars, metadata strips, and aligned heatmap columns
- named clades can be collapsed into figure summaries without mutating the input tree
- figure bundles can be emitted as a reusable package containing `figure.svg`, a JSON manifest, a caption draft, and aligned tip-annotation tables

## Comparative Analysis Highlights

- rooted trees and numeric traits can be checked explicitly for comparative readiness before modeling
- phylogenetic independent contrasts are available as one deterministic internal-node table per trait
- trait signal can be summarized with Blomberg's K, Pagel's lambda, and a permutation-based significance surface
- phylogenetic generalized least-squares accepts one or more numeric predictors and rejects categorical predictors or branch-length-incomplete trees explicitly

## Ancestral-State Highlights

- continuous ancestral-state reconstruction supports Brownian and OU-style trait models over a rooted pruned analysis tree
- discrete ancestral-state reconstruction supports Fitch parsimony with explicit ambiguous state sets and node-level probability summaries
- uncertainty is surfaced directly through continuous confidence intervals and discrete state-probability tables instead of hidden heuristics
- ancestral trees can be rendered with internal-node labels, exported as deterministic tables, compared across supported continuous models, and bundled into standalone HTML reports

## External Engine Highlights

- MAFFT-, trimAl-, IQ-TREE-, and FastTree-style workflows are exposed through one governed adapter surface instead of ad hoc shell snippets
- every engine workflow records the resolved executable, exact command vector, version output, stdout, stderr, and extracted warning lines in a deterministic manifest
- model selection emits a stable selected-model artifact, ML inference emits validated Newick trees, and bootstrap workflows retain both support trees and bootstrap tree sets
- fast approximate and ML trees can be compared through the same topology, clade, support, and branch-length report surface already used for native tree comparison
