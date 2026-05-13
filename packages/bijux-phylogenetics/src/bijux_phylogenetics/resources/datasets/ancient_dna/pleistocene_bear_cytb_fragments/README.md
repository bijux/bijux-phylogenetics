# Pleistocene Bear CYTB Fragment Panel

This packaged dataset contains a small real ancient-DNA-style mitochondrial
panel built for missing-data-aware sequence-to-tree workflow review. It pairs
modern bear cytochrome b references with real cave-bear cytochrome b
sequences reduced to short fragment-style inputs with explicit internal
missing-data blocks.

## Source

Source accessions:

- `OQ318974.1` brown bear complete mitochondrial genome
- `NC_003428.1` polar bear complete mitochondrial genome
- `OQ318956.1` American black bear complete mitochondrial genome
- `KX641337.1` cave bear ancient mitochondrial sequence
- `KX641335.1` cave bear ancient mitochondrial sequence

All packaged sequences derive from the cytochrome b coding region of those
published accessions.

## Fragment Construction

Modern reference sequences keep the full extracted cytochrome b region.

Ancient-style degraded sequences are derived from the real cave-bear
cytochrome b sequences with explicit sparse-fragment rules:

- `cave_bear_ud1838_fragment`: `KX641337.1` cytochrome b positions `121..760`
  plus one internal missing-data block of `24` `N` characters
- `cave_bear_wk01_fragment`: `KX641335.1` cytochrome b positions `281..760`
  plus one internal missing-data block of `18` `N` characters

This keeps the packaged panel honest about two missingness modes at once:
short degraded sequence length and explicit unresolved internal sites.

## Contents

- `sequences.fasta`: raw unaligned bear cytochrome b panel
- `expected/`: governed reference outputs from the owned FASTA-to-tree runtime

## Governed Workflow

The packaged workflow reruns the following owned sequence-to-tree chain:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE model selection
- IQ-TREE maximum-likelihood inference
- IQ-TREE bootstrap support estimation
- explicit missingness cleanup review over the aligned output

The governed workflow uses:

- sequence type: `dna`
- site missingness threshold: `0.15`
- sequence missingness threshold: `0.15`
- IQ-TREE seed: `1`
- IQ-TREE threads: `1`
- bootstrap replicates: `1000`
