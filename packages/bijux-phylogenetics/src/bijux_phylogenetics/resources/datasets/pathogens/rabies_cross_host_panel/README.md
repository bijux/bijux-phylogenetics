# Rabies Cross-Host Panel

This packaged dataset contains a small real pathogen panel built for
host-switching workflow review. It uses rabies virus nucleoprotein coding
sequences from bat, canid, and livestock hosts and ships with one rooted tree
plus one host metadata table.

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
region. Host labels come from the GenBank source qualifiers shipped with those
records.

## Host Coding

The metadata table carries both exact host species and one grouped host state:

- `host_species`: raw host label from the accession source feature
- `host_group`: grouped workflow state used by the packaged host-switching demo

The packaged workflow uses `host_group` rather than `host_species` because
this compact panel contains multiple singleton host species. Grouped host
states keep the host-switching review surface interpretable without pretending
that one small demonstration panel can resolve every species-level transition
cleanly.

## Tree Provenance

The packaged rooted tree was derived from the shipped nucleoprotein FASTA
panel with the owned sequence-to-tree workflow:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE model selection and maximum-likelihood inference
- IQ-TREE bootstrap support estimation
- outgroup rooting on `bat_chile_rv108` (`MG458305`)

The rooted tree is packaged directly so the public host-switching demo does
not depend on external executables.

## Contents

- `sequences.fasta`: raw rabies nucleoprotein coding sequences
- `tree.nwk`: rooted rabies tree used by the host-switching workflow
- `hosts.csv`: host metadata keyed by `taxon`
- `expected/`: governed host-switching review outputs regenerated from the
  owned runtime surface

## Governed Workflow

The packaged workflow reruns the following review surface:

- host-switching analysis over `host_group`
- discrete ancestral model `ard`
- summary, node, branch, count, fit, unsupported-claim, and exclusion ledgers
