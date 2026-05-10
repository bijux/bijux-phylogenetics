# Taxon workflow review bundle

Fixture-backed taxon trust evidence spanning identity, crosswalk, exclusion, loss, and stability behavior.

- study: `taxon-trust`
- evidence: `evidence-001`
- comparison mode: `bijux_native_reinterpretation`
- expected verdict: `matched`

## Claims

- `taxon-spelling-variant-audit`
- `taxonomic-synonym-candidate-detection`
- `controlled-synonym-resolution`
- `ambiguous-synonym-rejection`
- `taxon-namespace-classification`
- `mixed-namespace-detection`
- `taxon-crosswalk-table`
- `taxon-exclusion-reasoning`
- `workflow-taxon-loss-report`
- `taxon-stability-report`

## Source Basis

- `packages/bijux-phylogenetics/tests/fixtures`

## Primary Outputs

- `evidence-book/studies/taxon-trust/evidence-001/taxonomy_report_machine_manifest.json`

## Limits

- Covers one governed taxon-workflow bundle rather than the full repository surface.
