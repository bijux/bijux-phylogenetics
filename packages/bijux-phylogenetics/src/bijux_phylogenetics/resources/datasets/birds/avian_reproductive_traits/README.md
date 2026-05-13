# Avian Reproductive Trait Dataset

This packaged dataset contains a rooted bird phylogeny pruned to the observed
trait table plus a cleaned comparative trait matrix for reproductive and
developmental review.

## Contents

- `tree.nwk`: rooted 94-taxon bird tree used by the packaged workflow
- `traits.csv`: comparative trait table keyed by `species`
- `expected/`: governed reference workflow outputs regenerated from the owned
  runtime surface

## Trait Columns

- `species`: stable taxon identifier used across the tree and trait table
- `common_name`: human-readable bird name
- `mating_system`: observed mating-system category
- `testes_mass`: continuous reproductive trait
- `body_mass`: continuous body-size trait
- `multiple_paternity_percentage`: continuous percentage trait
- `development_mode`: categorical developmental state

## Governed Workflow

The packaged workflow bundle reruns the following review surfaces:

- PGLS with `testes_mass ~ body_mass`
- Brownian trait evolution for `testes_mass`
- OU trait evolution for `testes_mass`
- phylogenetic signal for `testes_mass`
- continuous ancestral reconstruction for `testes_mass`
- discrete ancestral reconstruction for `mating_system`
- clade-specific trait summaries for `mating_system`
