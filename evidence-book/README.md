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

- [`phylogenetic-101`](./phylogenetic-101/): block-by-block comparative
  validation against established R phylogenetics packages
  Registry: [`examples/registry.json`](./phylogenetic-101/examples/registry.json)
  Current example: [`primate-longevity-signal-workflow`](./phylogenetic-101/examples/primate-longevity-signal-workflow/accepted-tool-parity/README.md)
- [`taxon-trust`](./studies/taxon-trust/): reviewer-facing evidence for taxon identity,
  normalization, crosswalk, exclusion, loss, and stability behavior
  Study manifest: [`study.json`](./studies/taxon-trust/study.json)
  Current bundle: [`evidence-001`](./studies/taxon-trust/evidence-001/README.md)
