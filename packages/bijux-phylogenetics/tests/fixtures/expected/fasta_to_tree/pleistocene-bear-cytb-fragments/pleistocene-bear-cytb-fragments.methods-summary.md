# Tree Inference Methods Summary

The workflow manifest `pleistocene-bear-cytb-fragments.manifest.json` records one Bijux `fasta-to-tree` run from raw FASTA input `sequences.fasta` through alignment, trimming, model selection, maximum-likelihood inference, and bootstrap-supported tree finalization. The selected substitution model was `HKY+F`, the sequence type used for inference was `dna`, and the final delivered tree was `pleistocene-bear-cytb-fragments.tree`.

## Input And Alignment Preparation

- raw input path: `/Users/bijan/bijux/bijux-phylogenetics/packages/bijux-phylogenetics/src/bijux_phylogenetics/resources/datasets/ancient_dna/pleistocene_bear_cytb_fragments/sequences.fasta`
- input FASTA was used directly without identifier or record repair
- prepared input path: `/Users/bijan/bijux/bijux-phylogenetics/packages/bijux-phylogenetics/src/bijux_phylogenetics/resources/datasets/ancient_dna/pleistocene_bear_cytb_fragments/sequences.fasta`
- validated sequence count: `5`
- total raw residue count: `4540`
- inferred raw sequence alphabet: `dna`
- alignment engine: `mafft`
- alignment mode: `auto`
- aligned output path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/pleistocene-bear-cytb-fragments.aln`
- trimming engine: `trimal`
- trimming mode: `gap-threshold`
- trimming gap threshold: `0.1`
- retained alignment length: `1140` of `1140`
- removed alignment sites: `0`
- trimmed alignment path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/pleistocene-bear-cytb-fragments.trimmed.aln`

## Model Selection

- model-selection engine: `iqtree`
- candidate substitution models reviewed: `71`
- governing information criterion: `BIC`
- selected substitution model: `HKY+F`
- iqtree random seed: `1`
- iqtree threads: `1`
- model-selection manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/engine-artifacts/pleistocene-bear-cytb-fragments/model-selection/model-selection.manifest.json`

## Maximum-Likelihood Inference

- inference engine: `iqtree`
- inference model: `HKY+F`
- maximum-likelihood log-likelihood: `-1872.2463`
- unannotated maximum-likelihood tree artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/engine-artifacts/pleistocene-bear-cytb-fragments/maximum-likelihood/maximum-likelihood.treefile`
- inference manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/engine-artifacts/pleistocene-bear-cytb-fragments/maximum-likelihood/maximum-likelihood.manifest.json`

## Branch Support And Final Tree

- support engine: `iqtree`
- support workflow: ultrafast bootstrap support on the same trimmed alignment under `HKY+F`
- bootstrap replicates: `1000`
- supported internal nodes: `2` of `3`
- minimum/median/maximum support: `80.0` / `90.0` / `100.0`
- weakly supported clade count: `0`
- bootstrap-supported tree artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/engine-artifacts/pleistocene-bear-cytb-fragments/bootstrap-support/bootstrap-support.treefile`
- bootstrap tree-set artifact: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/engine-artifacts/pleistocene-bear-cytb-fragments/bootstrap-support/bootstrap-support.ufboot`
- support manifest: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/engine-artifacts/pleistocene-bear-cytb-fragments/bootstrap-support/bootstrap-support.manifest.json`

## Tree Processing And Traceability

- the final delivered tree is copied from the bootstrap-supported inference artifact so branch support remains attached to the reviewer-facing tree
- the workflow records a separate unannotated maximum-likelihood tree under engine-artifacts for audit and comparison
- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest
- final tree path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/pleistocene-bear-cytb-fragments.tree`
- reviewer-facing model table: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/pleistocene-bear-cytb-fragments.model.tsv`
- reviewer-facing support table: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/pleistocene-bear-cytb-fragments.support.tsv`
- workflow manifest path: `/Users/bijan/bijux/bijux-phylogenetics/artifacts/tmp4s0lj80m/pleistocene-bear-cytb-fragments/pleistocene-bear-cytb-fragments.manifest.json`

## Workflow Warnings

- WARNING: 1 sequences contain more than 50% gaps/ambiguity
- Warning! Some parameters hit the boundaries
- one or more internal nodes did not expose numeric support labels
- input contains sequence length outliers
- automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone
