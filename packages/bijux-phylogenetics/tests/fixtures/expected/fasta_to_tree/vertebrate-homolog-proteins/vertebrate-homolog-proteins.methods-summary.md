# Tree Inference Methods Summary

The workflow manifest `vertebrate-homolog-proteins.manifest.json` records one Bijux `fasta-to-tree` run from raw FASTA input `vertebrate-homolog-proteins.fasta` through alignment, trimming, model selection, maximum-likelihood inference, and bootstrap-supported tree finalization. The selected substitution model was `Q.insect+I`, the sequence type used for inference was `protein`, and the final delivered tree was `vertebrate-homolog-proteins.tree`.

## Input And Alignment Preparation

- raw input path: `/Users/bijan/bijux/bijux-phylogenetics/packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/inputs/vertebrate-homolog-proteins.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input path: `/Users/bijan/bijux/bijux-phylogenetics/packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/inputs/vertebrate-homolog-proteins.fasta`
- validated sequence count: `9`
- total raw residue count: `1560`
- inferred raw sequence alphabet: `protein`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/vertebrate-homolog-proteins.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `185` of `185`
- removed alignment sites: `0`
- trimmed alignment path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/vertebrate-homolog-proteins.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `283`
- governing information criterion: `BIC`
- selected substitution model: `Q.insect+I`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/engine-artifacts/vertebrate-homolog-proteins/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `Q.insect+I`
- maximum-likelihood log-likelihood: `-1704.768`
- unannotated maximum-likelihood tree artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/engine-artifacts/vertebrate-homolog-proteins/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/engine-artifacts/vertebrate-homolog-proteins/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `Q.insect+I`
- bootstrap replicates: `1000`
- supported internal nodes: `6` of `6`
- minimum/median/maximum support: `45.0` / `79.5` / `93.0`
- weakly supported clade count: `2`
- bootstrap-supported tree artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/engine-artifacts/vertebrate-homolog-proteins/bootstrap-support/bootstrap-support.contree`
- bootstrap tree-set artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/engine-artifacts/vertebrate-homolog-proteins/bootstrap-support/bootstrap-support.ufboot`
- support manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/engine-artifacts/vertebrate-homolog-proteins/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/vertebrate-homolog-proteins.tree`
- reviewer-facing model table: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/vertebrate-homolog-proteins.model.tsv`
- reviewer-facing support table: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/vertebrate-homolog-proteins.support.tsv`
- workflow manifest path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/vertebrate-recheck/vertebrate-homolog-proteins.manifest.json`

## Workflow Warnings

- WARNING: Normalizing state frequencies so that sum of them equals to 1
- one or more internal clades remain weakly supported
- input contains sequence length outliers
