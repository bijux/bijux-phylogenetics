---
title: Repository Overview
audience: public
type: handbook
status: active
owner: bijux-phylogenetics
last_reviewed: 2026-04-28
---

# Repository Overview

`bijux-phylogenetics` provides a governed Python surface for tree validation,
inspection, comparison, metadata linkage, alignment trimming, coding-sequence
translation, explicit rooting transforms, evidence bundles, and HTML report
generation.

The repository intentionally does not reimplement inference engines. Its
current product surface is the reproducible orchestration and evidence layer
around trees, alignments, and trait tables.

## Current CLI Surface

- `bijux-phylogenetics validate tree.nwk`
- `bijux-phylogenetics inspect tree.nwk`
- `bijux-phylogenetics compare left.nwk right.nwk`
- `bijux-phylogenetics annotate tree.nwk --metadata traits.tsv`
- `bijux-phylogenetics render tree.nwk --metadata traits.tsv --out report.html`
- `bijux-phylogenetics bundle run/ --out evidence-pack/`
- `bijux-phylogenetics report --tree tree.nwk --alignment alignment.fasta --traits traits.tsv --metadata samples.tsv --out phylo-report.html`
- `bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4`
- `bijux-phylogenetics alignment coding coding.fasta --json`
- `bijux-phylogenetics alignment translate coding.fasta --out translated.fasta`
- `bijux-phylogenetics alignment identity-matrix alignment.fasta --out identity.tsv`
- `bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk`
- `bijux-phylogenetics topology reroot-midpoint tree.nwk --out midpoint-rooted.nwk`
