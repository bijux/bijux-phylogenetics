# Tree Inference Methods Summary

This summary describes one governed Bijux `fasta-to-tree` run recorded in `pleistocene-bear-cytb-fragments.manifest.json`. The workflow starts from raw FASTA input `sequences.fasta`, runs alignment, trimming, model selection, maximum-likelihood inference, and bootstrap support, then delivers the final tree `pleistocene-bear-cytb-fragments.tree`. The selected substitution model was `HKY+F`, and the reviewed sequence type was `dna`.

## Input And Alignment Preparation

- raw input artifact: `sequences.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input artifact: `sequences.fasta`
- validated sequence count: `5`
- total raw residue count: `4540`
- inferred raw sequence alphabet: `dna`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output artifact: `pleistocene-bear-cytb-fragments.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `1140` of `1140`
- removed alignment sites: `0`
- trimmed alignment artifact: `pleistocene-bear-cytb-fragments.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `71`
- governing information criterion: `BIC`
- selected substitution model: `HKY+F`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest artifact: `engine-artifacts/pleistocene-bear-cytb-fragments/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `HKY+F`
- maximum-likelihood log-likelihood: `-1872.2463`
- unannotated maximum-likelihood tree artifact: `engine-artifacts/pleistocene-bear-cytb-fragments/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest artifact: `engine-artifacts/pleistocene-bear-cytb-fragments/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `HKY+F`
- bootstrap replicates: `1000`
- supported internal nodes: `2` of `3`
- minimum/median/maximum support: `80.0` / `90.0` / `100.0`
- weakly supported clade count: `0`
- bootstrap-supported tree artifact: `engine-artifacts/pleistocene-bear-cytb-fragments/bootstrap-support/bootstrap-support.treefile`
- bootstrap tree-set artifact: `engine-artifacts/pleistocene-bear-cytb-fragments/bootstrap-support/bootstrap-support.ufboot`
- support manifest artifact: `engine-artifacts/pleistocene-bear-cytb-fragments/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree artifact: `pleistocene-bear-cytb-fragments.tree`
- reviewer-facing model table: `pleistocene-bear-cytb-fragments.model.tsv`
- reviewer-facing support table: `pleistocene-bear-cytb-fragments.support.tsv`
- workflow manifest artifact: `pleistocene-bear-cytb-fragments.manifest.json`

## Workflow Warnings

- WARNING: 1 sequences contain more than 50% gaps/ambiguity
- Warning! Some parameters hit the boundaries
- one or more internal nodes did not expose numeric support labels
- input contains sequence length outliers
- automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone
