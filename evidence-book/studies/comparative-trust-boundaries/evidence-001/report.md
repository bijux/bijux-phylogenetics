# Comparative input rejection bundle

Governed expected-failure cases for rootedness, complete branch lengths, and numeric-response requirements in PGLS workflows.

- study: `comparative-trust-boundaries`
- evidence: `evidence-001`
- comparison mode: `bijux_native_reinterpretation`
- expected verdict: `matched`

## Local Artifacts

- `reference.R`: r-reference-program
- `analysis.py`: python-analysis-program
- `checks.json`: machine-check-contract
- `report.md`: human-report
- `provenance.json`: provenance-record

## Claims

- `comparative-input-rejection-governed`

## Source Basis

- `evidence-book/studies/comparative-trust-boundaries/provenance/runtime-sources.json`
- `packages/bijux-phylogenetics/tests/fixtures/trees/example_tree_unrooted.nwk`
- `packages/bijux-phylogenetics/tests/fixtures/trees/example_tree_no_lengths.nwk`
- `packages/bijux-phylogenetics/tests/fixtures/metadata/example_traits_comparative.tsv`

## Governed Primary Outputs

- `evidence-book/studies/comparative-trust-boundaries/evidence-001/expected-failure-cases.json`

## Results Directory

- `evidence-book/studies/comparative-trust-boundaries/evidence-001/results/README.md`
- `evidence-book/studies/comparative-trust-boundaries/evidence-001/results/manifest.json`

## Limits

- This bundle audits representative rejection cases, not every possible malformed comparative input.
