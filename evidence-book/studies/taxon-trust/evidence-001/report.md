# Taxon workflow review bundle

Fixture-backed taxon trust evidence spanning identity, crosswalk, exclusion, loss, and stability behavior.

- study: `taxon-trust`
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

## Governed Primary Outputs

- `evidence-book/studies/taxon-trust/evidence-001/taxonomy_report_machine_manifest.json`

## Results Directory

- `evidence-book/studies/taxon-trust/evidence-001/results/README.md`
- `evidence-book/studies/taxon-trust/evidence-001/results/manifest.json`

## Limits

- Covers one governed taxon-workflow bundle rather than the full repository surface.
