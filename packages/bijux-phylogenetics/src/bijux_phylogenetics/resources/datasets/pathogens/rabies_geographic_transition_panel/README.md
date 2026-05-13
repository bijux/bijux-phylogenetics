# Rabies Geographic Transition Panel

This packaged dataset contains a small real pathogen panel built for
geographic transition review. It uses rabies virus nucleoprotein coding
sequences from five macroregions and ships with one rooted tree plus one
region metadata table.

## Source

Source accessions:

- `MG458305` rabies virus isolate `RV108`
- `MG458304` rabies virus isolate `RV50`
- `PV641713` rabies virus isolate `PKDOW/PDB011`
- `PX845689` rabies virus isolate `4092/26/2025/Novosibirsk`
- `OQ693985` rabies virus isolate `RV_2277140521J_Mazovie_Raccoon_dog_Pol_2021`
- `PX845683` rabies virus isolate `4092/20/2025/Novosibirsk`
- `PX845681` rabies virus isolate `4092/18/2025/Novosibirsk`
- `PX845678` rabies virus isolate `4092/15/2025/Novosibirsk`
- `PX845676` rabies virus isolate `4092/13/2025/Novosibirsk`

All packaged sequences derive from the rabies virus nucleoprotein coding
region. Country labels come from the GenBank source qualifiers shipped with
those records.

## Geographic Coding

The metadata table carries both source-country provenance and one grouped
workflow state:

- `country`: raw country or locality label from the accession source feature
- `region_group`: grouped macroregion used by the packaged geography demo

The packaged workflow uses `region_group` rather than the raw country field
because this compact panel contains multiple singleton localities. Grouped
regions keep the biogeography review surface interpretable without pretending
that one small demonstration panel can resolve every locality-level movement
cleanly.

The grouped regions in this panel are:

- `south_america`
- `north_america`
- `europe`
- `south_asia`
- `north_asia`

## Tree Provenance

The packaged rooted tree was derived from the shipped nucleoprotein FASTA
panel with the owned sequence-to-tree workflow:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE model selection and maximum-likelihood inference
- IQ-TREE bootstrap support estimation
- outgroup rooting on `bat_chile_rv108` (`MG458305`)

The rooted tree is packaged directly so the public geography demo does not
depend on external executables.

## Contents

- `sequences.fasta`: raw rabies nucleoprotein coding sequences
- `tree.nwk`: rooted rabies tree used by the geography workflow
- `regions.csv`: geographic metadata keyed by `taxon`
- `expected/`: governed geography-transition outputs regenerated from the
  owned runtime surface

## Governed Workflow

The packaged workflow reruns the following review surfaces:

- geographic state modeling over `region_group`
- discrete ancestral model `ard`
- internal ancestral region probability ledgers
- geographic transition-rate and branch-event ledgers
- branchwise geographic migration-event extraction
