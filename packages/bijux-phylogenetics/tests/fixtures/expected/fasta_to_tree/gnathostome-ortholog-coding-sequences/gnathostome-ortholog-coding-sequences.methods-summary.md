# Tree Inference Methods Summary

This summary describes one governed Bijux `fasta-to-tree` run recorded in `gnathostome-ortholog-coding-sequences.manifest.json`. The workflow starts from raw FASTA input `inputs/gnathostome-ortholog-coding-sequences.fasta`, runs alignment, trimming, model selection, maximum-likelihood inference, and bootstrap support, then delivers the final tree `gnathostome-ortholog-coding-sequences.tree`. The selected substitution model was `TNe+G4`, and the reviewed sequence type was `dna`.

## Input And Alignment Preparation

- raw input artifact: `inputs/gnathostome-ortholog-coding-sequences.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input artifact: `inputs/gnathostome-ortholog-coding-sequences.fasta`
- validated sequence count: `9`
- total raw residue count: `4707`
- inferred raw sequence alphabet: `dna`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output artifact: `gnathostome-ortholog-coding-sequences.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `558` of `558`
- removed alignment sites: `0`
- trimmed alignment artifact: `gnathostome-ortholog-coding-sequences.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `50`
- governing information criterion: `BIC`
- selected substitution model: `TNe+G4`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest artifact: `engine-artifacts/gnathostome-ortholog-coding-sequences/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `TNe+G4`
- maximum-likelihood log-likelihood: `-3440.0991`
- unannotated maximum-likelihood tree artifact: `engine-artifacts/gnathostome-ortholog-coding-sequences/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest artifact: `engine-artifacts/gnathostome-ortholog-coding-sequences/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `TNe+G4`
- bootstrap replicates: `1000`
- supported internal nodes: `6` of `7`
- minimum/median/maximum support: `59.0` / `85.0` / `97.0`
- weakly supported clade count: `2`
- bootstrap-supported tree artifact: `engine-artifacts/gnathostome-ortholog-coding-sequences/bootstrap-support/bootstrap-support.treefile`
- bootstrap tree-set artifact: `engine-artifacts/gnathostome-ortholog-coding-sequences/bootstrap-support/bootstrap-support.ufboot`
- support manifest artifact: `engine-artifacts/gnathostome-ortholog-coding-sequences/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree artifact: `gnathostome-ortholog-coding-sequences.tree`
- reviewer-facing model table: `gnathostome-ortholog-coding-sequences.model.tsv`
- reviewer-facing support table: `gnathostome-ortholog-coding-sequences.support.tsv`
- workflow manifest artifact: `gnathostome-ortholog-coding-sequences.manifest.json`

## Workflow Warnings

- one or more internal nodes did not expose numeric support labels
- one or more internal clades remain weakly supported
- input contains sequence length outliers
- automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone
