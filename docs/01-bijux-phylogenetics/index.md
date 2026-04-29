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
translation, explicit rooting transforms, deterministic tree rendering,
publication figure packaging, evidence bundles, and HTML report generation.

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
