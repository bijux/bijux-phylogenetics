# Primate missing-data accounting parity bundle

Governed evidence for how missing values are handled before the processed primate table is used downstream.

- study: `primate-longevity-signal`
- evidence: `evidence-004`
- comparison mode: `direct_parity`
- expected verdict: `matched`

## Local Artifacts

- `reference.R`: r-reference-program
- `analysis.py`: python-analysis-program
- `checks.json`: machine-check-contract
- `report.md`: human-report
- `provenance.json`: provenance-record

## Claims

- `pcm1-missing-data-accounting-parity`

## Source Basis

- `evidence-book/studies/primate-longevity-signal/evidence-001/reference_primate.csv`
- `evidence-book/studies/primate-longevity-signal/evidence-004/missing-data-accounting-parity.json`

## Governed Primary Outputs

- `evidence-book/studies/primate-longevity-signal/evidence-004/missing-data-accounting-parity.json`

## Results Directory

- `evidence-book/studies/primate-longevity-signal/evidence-004/results/README.md`
- `evidence-book/studies/primate-longevity-signal/evidence-004/results/manifest.json`

## Limits

- The raw workbook remains external, so the bundle records the governed missing-value contract at the processed artifact boundary.
