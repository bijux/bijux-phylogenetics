# Central European Seashore Flora Dataset

This packaged dataset contains a rooted plant phylogeny pruned to a published
saltwater-and-seashore subset of the Central European flora plus a cleaned
comparative trait matrix for non-animal workflow review.

## Source

Source dataset:

- Vandelook, Filip; Janssens, Steven B.; Matthies, Diethart (2018). Data from:
  Ecological niche and phylogeny explain distribution of seed mass in the
  Central European flora. Dryad. `10.5061/dryad.0st06f0`

Packaged subset rule:

- retain taxa whose published `habitat` class is `saltwater_and_seashore`

## Contents

- `tree.nwk`: rooted 42-taxon plant tree used by the packaged workflow
- `traits.csv`: comparative trait table keyed by `species`
- `expected/`: governed reference workflow outputs regenerated from the owned
  runtime surface

## Trait Columns

- `species`: stable taxon identifier used across the tree and trait table
- `seed_mass`: dispersal-unit dry mass in milligrams
- `plant_height`: maximum plant height in centimeters
- `lifeform`: descriptive lifeform state from the published coding
- `longevity`: descriptive longevity state from the published coding
- `light`: Ellenberg indicator value for light
- `temperature`: Ellenberg indicator value for temperature
- `continentality`: Ellenberg indicator value for continentality
- `humidity`: Ellenberg indicator value for humidity
- `reaction`: Ellenberg indicator value for reaction
- `nitrogen`: Ellenberg indicator value for nitrogen
- `salt`: Ellenberg indicator value for salt
- `habitat`: descriptive habitat state from the published coding

## Governed Workflow

The packaged workflow bundle reruns the following review surfaces:

- PGLS with `seed_mass ~ plant_height`
- Brownian trait evolution for `seed_mass`
- OU trait evolution for `seed_mass`
- phylogenetic signal for `seed_mass`
- continuous ancestral reconstruction for `seed_mass`
- discrete ancestral reconstruction for `lifeform`
- clade-specific trait summaries for `lifeform`
