# Catarrhine Mitogenome Five-Locus Panel

This packaged dataset contains a small real multi-locus panel built for
concatenation, partitioning, occupancy review, and partitioned inference. It
uses five mitochondrial coding-gene alignments extracted from complete
catarrhine mitochondrial genomes and ships with one explicit taxon ledger plus
one locus directory.

## Source

Source accessions:

- `NC_012920.1` `Homo sapiens`
- `NC_001643.1` `Pan troglodytes`
- `NC_011120.1` `Gorilla gorilla gorilla`
- `NC_001646.1` `Pongo pygmaeus`
- `NC_002082.1` `Hylobates lar`
- `NC_005943.1` `Macaca mulatta`

All packaged loci were extracted from complete mitochondrial genome RefSeq
records. The panel keeps one old-world-monkey taxon (`Macaca mulatta`) plus
five hominoid taxa so the partitioned inference surface has a stable internal
topology and one clear external branch.

## Locus Design

The dataset ships five aligned mitochondrial coding loci:

- `mt-cox1`
- `mt-cox2`
- `mt-cox3`
- `mt-cytb`
- `mt-nd2`

Every locus is already aligned and keyed by the same six `taxon` identifiers.
That makes the packaged workflow focus on multi-locus assembly and partitioned
inference rather than raw read cleanup or unaligned-sequence preprocessing.

## Contents

- `taxa.csv`: taxon and source-accession ledger
- `loci/*.fasta`: aligned per-locus FASTA inputs
- `expected/`: governed multi-locus workflow outputs regenerated from the
  owned runtime surface

## Governed Workflow

The packaged workflow reruns the following surfaces:

- concatenation of the five aligned loci into one supermatrix
- partition-file generation over named loci
- locus-occupancy summary, per-taxon coverage, per-locus coverage, and matrix
  ledgers
- IQ-TREE partitioned model selection on the concatenated supermatrix
- IQ-TREE partitioned bootstrap-support inference on the same supermatrix and
  partition file

This packaged surface is intentionally honest about scale. It is a small real
multi-gene panel built to exercise the repository’s partitioned phylogenomics
workflow contract, not a claim of broad genome-scale taxon sampling.
