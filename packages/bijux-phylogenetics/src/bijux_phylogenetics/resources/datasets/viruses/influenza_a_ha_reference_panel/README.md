# Influenza A Hemagglutinin Reference Panel

This packaged dataset contains a small real viral FASTA panel built for the
owned FASTA-to-tree workflow. It starts from unaligned Influenza A
hemagglutinin segment-4 sequences and supports alignment, trimming,
maximum-likelihood inference, and bootstrap support review through the public
runtime surface.

## Source

Source accessions:

- `NC_002017.1` Influenza A virus `A/Puerto Rico/8/1934(H1N1)` segment 4
- `CY033655.1` Influenza A virus `A/Texas/36/1991(H1N1)` segment 4
- `CY046787.1` Influenza A virus `A/Wisconsin/629-D02276/2009(H1N1)` segment 4
- `NC_007366.1` Influenza A virus `A/New York/392/2004(H3N2)` segment 4
- `NC_007374.1` Influenza A virus `A/Korea/426/1968(H2N2)` segment 4
- `AY653200.1` Influenza A virus `A/chicken/Jilin/9/2004(H5N1)` segment 4

All sequences were fetched from the NCBI nucleotide database and normalized to
stable taxon identifiers in the packaged FASTA file.

## Contents

- `sequences.fasta`: raw unaligned viral nucleotide sequences
- `expected/`: governed reference outputs from the owned FASTA-to-tree runtime

## Governed Workflow

The packaged workflow reruns the following owned sequence-to-tree chain:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE model selection
- IQ-TREE maximum-likelihood inference
- IQ-TREE bootstrap support estimation

The governed workflow uses:

- sequence type: `dna`
- IQ-TREE seed: `1`
- IQ-TREE threads: `1`
- bootstrap replicates: `1000`
