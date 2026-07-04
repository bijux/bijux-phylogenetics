# Tree Inference Methods Summary

This summary describes one governed Bijux `fasta-to-tree` run recorded in `gnathostome-ortholog-proteins.manifest.json`. The workflow starts from raw FASTA input `inputs/gnathostome-ortholog-proteins.fasta`, runs alignment, trimming, model selection, maximum-likelihood inference, and bootstrap support, then delivers the final tree `gnathostome-ortholog-proteins.tree`. The selected substitution model was `Q.insect+I`, and the reviewed sequence type was `protein`.

## Input And Alignment Preparation

- raw input artifact: `inputs/gnathostome-ortholog-proteins.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input artifact: `inputs/gnathostome-ortholog-proteins.fasta`
- validated sequence count: `9`
- total raw residue count: `1560`
- inferred raw sequence alphabet: `protein`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output artifact: `gnathostome-ortholog-proteins.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `185` of `185`
- removed alignment sites: `0`
- trimmed alignment artifact: `gnathostome-ortholog-proteins.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `283`
- governing information criterion: `BIC`
- selected substitution model: `Q.insect+I`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest artifact: `engine-artifacts/gnathostome-ortholog-proteins/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `Q.insect+I`
- maximum-likelihood log-likelihood: `-1704.768`
- unannotated maximum-likelihood tree artifact: `engine-artifacts/gnathostome-ortholog-proteins/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest artifact: `engine-artifacts/gnathostome-ortholog-proteins/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `Q.insect+I`
- bootstrap replicates: `1000`
- supported internal nodes: `6` of `7`
- minimum/median/maximum support: `45.0` / `79.5` / `93.0`
- weakly supported clade count: `2`
- bootstrap-supported tree artifact: `engine-artifacts/gnathostome-ortholog-proteins/bootstrap-support/bootstrap-support.contree`
- bootstrap tree-set artifact: `engine-artifacts/gnathostome-ortholog-proteins/bootstrap-support/bootstrap-support.ufboot`
- support manifest artifact: `engine-artifacts/gnathostome-ortholog-proteins/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree artifact: `gnathostome-ortholog-proteins.tree`
- reviewer-facing model table: `gnathostome-ortholog-proteins.model.tsv`
- reviewer-facing support table: `gnathostome-ortholog-proteins.support.tsv`
- workflow manifest artifact: `gnathostome-ortholog-proteins.manifest.json`

## Workflow Warnings

- WARNING: Normalizing state frequencies so that sum of them equals to 1
- one or more internal nodes did not expose numeric support labels
- one or more internal clades remain weakly supported
- input contains sequence length outliers
