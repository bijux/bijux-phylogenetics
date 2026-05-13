## Primate Comparative Dataset

This packaged dataset is the first real mammal comparative dataset shipped with
`bijux-phylogenetics`.

- `tree.nwk` is a governed trimmed primate tree with 75 tip taxa.
- `traits.csv` is a governed primate comparative table keyed by `species`.
- `expected/` contains checked reference outputs from the owned Bijux
  comparative workflow surface.

The packaged inputs are derived from the repository-owned primate study inputs
under `evidence-book/studies/primate-longevity-signal/datasets/`. They are
promoted here so users can run comparative workflows without depending on the
evidence-book layout.

The reference workflow uses:

- continuous trait: `longevity`
- PGLS predictor: `social_group_size`
- discrete trait: `mating_system`

Continuous traits present in the table include `body_mass`, `gestation`,
`home_range`, `longevity`, and `social_group_size`.

Categorical traits present in the table include `family`, `sex_dimorphism`, and
`mating_system`.
