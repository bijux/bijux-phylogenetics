# Tree Inference Methods Summary

This summary describes one governed Bijux `fasta-to-tree` run recorded in `influenza-a-ha-reference-panel.manifest.json`. The workflow starts from raw FASTA input `sequences.fasta`, runs alignment, trimming, model selection, maximum-likelihood inference, and bootstrap support, then delivers the final tree `influenza-a-ha-reference-panel.tree`. The selected substitution model was `TPM3u+F+R2`, and the reviewed sequence type was `dna`.

## Input And Alignment Preparation

- raw input artifact: `sequences.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input artifact: `sequences.fasta`
- validated sequence count: `6`
- total raw residue count: `10561`
- inferred raw sequence alphabet: `dna`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output artifact: `influenza-a-ha-reference-panel.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `1820` of `1820`
- removed alignment sites: `0`
- trimmed alignment artifact: `influenza-a-ha-reference-panel.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `71`
- governing information criterion: `BIC`
- selected substitution model: `TPM3u+F+R2`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest artifact: `engine-artifacts/influenza-a-ha-reference-panel/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `TPM3u+F+R2`
- maximum-likelihood log-likelihood: `-9379.188`
- unannotated maximum-likelihood tree artifact: `engine-artifacts/influenza-a-ha-reference-panel/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest artifact: `engine-artifacts/influenza-a-ha-reference-panel/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `TPM3u+F+R2`
- bootstrap replicates: `1000`
- supported internal nodes: `3` of `4`
- minimum/median/maximum support: `98.0` / `99.0` / `99.0`
- weakly supported clade count: `0`
- bootstrap-supported tree artifact: `engine-artifacts/influenza-a-ha-reference-panel/bootstrap-support/bootstrap-support.treefile`
- bootstrap tree-set artifact: `engine-artifacts/influenza-a-ha-reference-panel/bootstrap-support/bootstrap-support.ufboot`
- support manifest artifact: `engine-artifacts/influenza-a-ha-reference-panel/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree artifact: `influenza-a-ha-reference-panel.tree`
- reviewer-facing model table: `influenza-a-ha-reference-panel.model.tsv`
- reviewer-facing support table: `influenza-a-ha-reference-panel.support.tsv`
- workflow manifest artifact: `influenza-a-ha-reference-panel.manifest.json`

## Workflow Warnings

- one or more internal nodes did not expose numeric support labels
- automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone
