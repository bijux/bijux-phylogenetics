# Tree Inference Methods Summary

The workflow manifest `strnog-enog411bqtj-proteins.manifest.json` records one Bijux `fasta-to-tree` run from raw FASTA input `strnog-enog411bqtj-proteins.fasta` through alignment, trimming, model selection, maximum-likelihood inference, and bootstrap-supported tree finalization. The selected substitution model was `JTTDCMut+G4`, the sequence type used for inference was `protein`, and the final delivered tree was `strnog-enog411bqtj-proteins.tree`.

## Input And Alignment Preparation

- raw input path: `/Users/bijan/bijux/bijux-phylogenetics/packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/inputs/strnog-enog411bqtj-proteins.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input path: `/Users/bijan/bijux/bijux-phylogenetics/packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/inputs/strnog-enog411bqtj-proteins.fasta`
- validated sequence count: `31`
- total raw residue count: `4314`
- inferred raw sequence alphabet: `protein`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/strnog-enog411bqtj-proteins.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `145` of `324`
- removed alignment sites: `179`
- trimmed alignment path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/strnog-enog411bqtj-proteins.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `283`
- governing information criterion: `BIC`
- selected substitution model: `JTTDCMut+G4`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/engine-artifacts/strnog-enog411bqtj-proteins/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `JTTDCMut+G4`
- maximum-likelihood log-likelihood: `-1589.6619`
- unannotated maximum-likelihood tree artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/engine-artifacts/strnog-enog411bqtj-proteins/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/engine-artifacts/strnog-enog411bqtj-proteins/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `JTTDCMut+G4`
- bootstrap replicates: `1000`
- supported internal nodes: `28` of `28`
- minimum/median/maximum support: `0.0` / `80.5` / `100.0`
- weakly supported clade count: `12`
- bootstrap-supported tree artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/engine-artifacts/strnog-enog411bqtj-proteins/bootstrap-support/bootstrap-support.contree`
- bootstrap tree-set artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/engine-artifacts/strnog-enog411bqtj-proteins/bootstrap-support/bootstrap-support.ufboot`
- support manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/engine-artifacts/strnog-enog411bqtj-proteins/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/strnog-enog411bqtj-proteins.tree`
- reviewer-facing model table: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/strnog-enog411bqtj-proteins.model.tsv`
- reviewer-facing support table: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/strnog-enog411bqtj-proteins.support.tsv`
- workflow manifest path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/fasta-to-tree-golden-refresh/strnog-recheck/strnog-enog411bqtj-proteins.manifest.json`

## Workflow Warnings

- WARNING: Normalizing state frequencies so that sum of them equals to 1
- one or more internal clades remain weakly supported
- input contains sequence length outliers
