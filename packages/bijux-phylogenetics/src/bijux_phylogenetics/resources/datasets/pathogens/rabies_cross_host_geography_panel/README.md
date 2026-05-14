# Rabies Cross-Host Geography Panel

This packaged dataset contains one real rabies virus nucleoprotein panel built
for full sequence-to-result workflow review. It starts from raw coding
sequences and one combined metadata table, then reruns alignment, trimming,
alignment-quality review, maximum-likelihood tree inference, bootstrap support
estimation, bootstrap topology summary, clade extraction, host-state
reconstruction, geographic transition analysis, migration-event extraction,
one comparative model, and a final integrated report.

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

- FASTA validation with explicit sequence-type detection
- MAFFT multiple-sequence alignment
- trimAl alignment trimming
- alignment quality scoring and sequence ranking
- IQ-TREE model selection
- IQ-TREE maximum-likelihood inference
- IQ-TREE bootstrap support estimation
- bootstrap tree-set topology review and consensus summary
- outgroup rooting on `bat_chile_rv108` (`MG458305`)
- rooted-tree clade extraction with host and region metadata summaries
- host-switching analysis over `host_group`
- geographic transition analysis over `region_group`
- comparative regression on derived regional longitude with explicit
  branch-length repair when the rooted tree carries nonpositive branches
- migration-event extraction
- integrated HTML report with tree, tables, bootstrap review, comparative
  outputs, and map outputs

The packaged `workflow-config.json` file is the governed reproduction entry
point for this workflow. It resolves the input files relative to the dataset
directory and records the comparative formula, clade columns, bootstrap review
thresholds, and engine-facing workflow defaults used by the public runtime
surface.

The governed workflow bundle now also includes:

- `workflow-config-audit.tsv`: config checks over input presence, taxon
  crosswalks, centroid coverage, and comparative-response availability
- `workflow-config.resolved.json`: resolved workflow parameters plus input
  checksums
- `bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv`: one stable
  ML-versus-consensus topology summary
- `bootstrap-review/rooted-tree-vs-bootstrap-consensus.comparison.tsv`: split
  comparison ledger across topology, support, and branch-length surfaces
- `scientific-findings.tsv`: reviewer-facing biological findings with linked
  source artifacts and cautions

## Contents

- `sequences.fasta`: raw rabies nucleoprotein coding sequences
- `metadata.csv`: combined host and geography metadata keyed by `taxon`
- `region-centroids.csv`: explicit region centroids for geographic map review
- `workflow-config.json`: one config file that reproduces the full workflow
- `expected/`: governed stable workflow outputs regenerated from the owned
  runtime surface

When the public demo is materialized, it also writes
`dataset/source-accessions.tsv` so the packaged accession set remains available
as one machine-readable provenance ledger alongside the copied dataset files.
