# Gnathostome Ortholog Protein Benchmark

This packaged dataset contains one small real amino-acid benchmark for the
owned FASTA-to-tree workflow. It starts from unaligned ortholog protein
sequences and supports MAFFT alignment, trimAl trimming, IQ-TREE protein model
selection, maximum-likelihood inference, and bootstrap support review through a
public runtime surface.

## Source

The packaged FASTA derives from the repository's governed trimAl reference
corpus entry `example.009.AA`.

Transformation:

- removed alignment gap characters from the aligned source FASTA
- removed placeholder missing-data marks
- kept the resulting raw amino-acid sequences as the benchmark input

The benchmark contains 9 protein records and keeps the original stable
identifiers from that governed corpus.

## Contents

- `sequences.fasta`: raw unaligned protein FASTA input
- `expected/`: governed benchmark outputs from the owned protein workflow

## Governed Workflow

The packaged workflow reruns the following owned sequence-to-tree chain:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE protein model selection
- IQ-TREE maximum-likelihood inference
- IQ-TREE bootstrap support estimation

The governed workflow uses:

- sequence type: `protein`
- IQ-TREE sequence keyword: `AA`
- IQ-TREE seed: `1`
- IQ-TREE threads: `1`
- bootstrap replicates: `1000`

This is a protein benchmark, not a translated coding-DNA workflow. No codon
translation, codon-position partitioning, or nucleotide-specific model
assumptions are part of this surface.
