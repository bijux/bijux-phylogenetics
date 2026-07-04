# Tree Inference Methods Summary

This summary describes one governed Bijux `fasta-to-tree` run recorded in `strnog-enog411bqtj-proteins.manifest.json`. The workflow starts from raw FASTA input `inputs/strnog-enog411bqtj-proteins.fasta`, runs alignment, trimming, model selection, maximum-likelihood inference, and bootstrap support, then delivers the final tree `strnog-enog411bqtj-proteins.tree`. The selected substitution model was `JTTDCMut+G4`, and the reviewed sequence type was `protein`.

## Input And Alignment Preparation

- raw input artifact: `inputs/strnog-enog411bqtj-proteins.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input artifact: `inputs/strnog-enog411bqtj-proteins.fasta`
- validated sequence count: `31`
- total raw residue count: `4314`
- inferred raw sequence alphabet: `protein`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output artifact: `strnog-enog411bqtj-proteins.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `145` of `324`
- removed alignment sites: `179`
- trimmed alignment artifact: `strnog-enog411bqtj-proteins.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `283`
- governing information criterion: `BIC`
- selected substitution model: `JTTDCMut+G4`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest artifact: `engine-artifacts/strnog-enog411bqtj-proteins/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `JTTDCMut+G4`
- maximum-likelihood log-likelihood: `-1589.6619`
- unannotated maximum-likelihood tree artifact: `engine-artifacts/strnog-enog411bqtj-proteins/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest artifact: `engine-artifacts/strnog-enog411bqtj-proteins/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `JTTDCMut+G4`
- bootstrap replicates: `1000`
- supported internal nodes: `28` of `28`
- minimum/median/maximum support: `0.0` / `80.5` / `100.0`
- weakly supported clade count: `12`
- bootstrap-supported tree artifact: `engine-artifacts/strnog-enog411bqtj-proteins/bootstrap-support/bootstrap-support.contree`
- bootstrap tree-set artifact: `engine-artifacts/strnog-enog411bqtj-proteins/bootstrap-support/bootstrap-support.ufboot`
- support manifest artifact: `engine-artifacts/strnog-enog411bqtj-proteins/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree artifact: `strnog-enog411bqtj-proteins.tree`
- reviewer-facing model table: `strnog-enog411bqtj-proteins.model.tsv`
- reviewer-facing support table: `strnog-enog411bqtj-proteins.support.tsv`
- workflow manifest artifact: `strnog-enog411bqtj-proteins.manifest.json`

## Workflow Warnings

- WARNING: Normalizing state frequencies so that sum of them equals to 1
- one or more internal clades remain weakly supported
- input contains sequence length outliers
