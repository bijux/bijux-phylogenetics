# Evidence Book

This directory holds checked-in comparative trust studies and curated evidence
bundles for `bijux-phylogenetics`.

These files are intentionally different from transient local outputs under
`artifacts/`:

- `evidence-book/` contains reviewed, versioned evidence summaries that are
  meant to stay in the repository
- `artifacts/` contains scratch execution products, caches, and local build
  outputs that should not define repository trust on their own

Current tracked studies:

- [`primate-longevity-signal`](./studies/primate-longevity-signal/): block-by-block
  comparative parity against established R phylogenetics packages
  Study manifest: [`study.json`](./studies/primate-longevity-signal/study.json)
  Current bundle: [`evidence-001`](./studies/primate-longevity-signal/evidence-001/README.md)
- [`primate-pgls-and-signal`](./studies/primate-pgls-and-signal/): governed parity
  study for the Lund regression and phylogenetic-signal lecture sections, with
  explicit coverage boundaries for the still-open EB and ancestral fragments
  Study manifest: [`study.json`](./studies/primate-pgls-and-signal/study.json)
  Current bundle: [`evidence-001`](./studies/primate-pgls-and-signal/evidence-001/README.md)
- [`taxon-trust`](./studies/taxon-trust/): reviewer-facing evidence for taxon identity,
  normalization, crosswalk, exclusion, loss, and stability behavior
  Study manifest: [`study.json`](./studies/taxon-trust/study.json)
  Current bundle: [`evidence-001`](./studies/taxon-trust/evidence-001/README.md)
