# OU identifiability warning bundle

Governed audit of built-in reference cases that should trigger OU identifiability warnings.

- study: `comparative-trust-boundaries`
- evidence: `evidence-003`
- comparison mode: `bijux_native_reinterpretation`
- expected verdict: `matched`

## Local Artifacts

- `reference.R`: r-reference-program
- `analysis.py`: python-analysis-program
- `checks.json`: machine-check-contract
- `report.md`: human-report
- `provenance.json`: provenance-record

## Claims

- `ou-identifiability-boundaries-detected`

## Source Basis

- `evidence-book/studies/comparative-trust-boundaries/provenance/runtime-sources.json`
- `packages/bijux-phylogenetics/tests/fixtures`
- `packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/models.py`

## Governed Primary Outputs

- `evidence-book/studies/comparative-trust-boundaries/evidence-003/ou-identifiability-audit.json`

## Results Directory

- `evidence-book/studies/comparative-trust-boundaries/evidence-003/results/README.md`
- `evidence-book/studies/comparative-trust-boundaries/evidence-003/results/manifest.json`

## Limits

- The audit covers known reference warning families and does not claim that every possible OU pathology is enumerated here.
