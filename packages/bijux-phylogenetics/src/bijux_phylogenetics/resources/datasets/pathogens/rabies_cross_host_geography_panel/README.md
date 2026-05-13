# Rabies Cross-Host Geography Panel

This packaged dataset contains one real rabies virus nucleoprotein panel built
for full sequence-to-result workflow review. It starts from raw coding
sequences and one combined metadata table, then reruns alignment, trimming,
maximum-likelihood tree inference, bootstrap support estimation, host-state
reconstruction, geographic transition analysis, migration-event extraction,
and a final integrated report.

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
region. Host and locality labels come from the accession source qualifiers
shipped with those records.

## Metadata Coding

The combined metadata table carries both raw provenance columns and grouped
workflow states:

- `host_species`: raw host label from the accession source feature
- `host_group`: grouped host state used for host-switching analysis
- `country`: raw country or locality label from the accession source feature
- `region_group`: grouped macroregion used for geographic transition analysis

The packaged workflow uses grouped host and region states because this compact
panel contains several singleton species and localities. The grouped coding
keeps the scientific claims interpretable without pretending that one small
demonstration panel can resolve every species-level host jump or every
locality-level geographic move cleanly.

The bundled centroid table provides one explicit map placement for each grouped
region so the final geographic review can render a stable self-contained map.

## Governed Workflow

The packaged workflow reruns the following owned surfaces from raw-ish inputs:

- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- IQ-TREE model selection
- IQ-TREE maximum-likelihood inference
- IQ-TREE bootstrap support estimation
- outgroup rooting on `bat_chile_rv108` (`MG458305`)
- host-switching analysis over `host_group`
- geographic transition analysis over `region_group`
- migration-event extraction
- integrated HTML report with tree, tables, and map outputs

## Contents

- `sequences.fasta`: raw rabies nucleoprotein coding sequences
- `metadata.csv`: combined host and geography metadata keyed by `taxon`
- `region-centroids.csv`: explicit region centroids for geographic map review
- `expected/`: governed stable workflow outputs regenerated from the owned
  runtime surface
