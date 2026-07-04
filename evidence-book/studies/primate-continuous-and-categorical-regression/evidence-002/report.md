# Primate baseline regression structure bundle

Governed study slice for the PCM3 baseline GLS fit with one continuous predictor and one categorical predictor.

- study: `primate-continuous-and-categorical-regression`
- evidence: `evidence-002`
- comparison mode: `direct_parity`
- expected verdict: `not_comparable`

## Local Artifacts

- `reference.R`: r-reference-program
- `analysis.py`: python-analysis-program
- `checks.json`: machine-check-contract
- `report.md`: human-report
- `provenance.json`: provenance-record

## Claims

- `pcm3-baseline-gls-with-categorical-predictor`

## Source Basis

- `evidence-book/studies/primate-continuous-and-categorical-regression/provenance/lund-course-sources.json`
- `evidence-book/studies/primate-continuous-and-categorical-regression/datasets/registry.json`
- `evidence-book/studies/primate-continuous-and-categorical-regression/reference/README.md`
- `packages/bijux-phylogenetics/tests/test_pgls_interaction_coefficients.py`
- `packages/bijux-phylogenetics/tests/test_regression_model_selection.py`

## Governed Primary Outputs


## Results Directory

- `evidence-book/studies/primate-continuous-and-categorical-regression/evidence-002/results/README.md`
- `evidence-book/studies/primate-continuous-and-categorical-regression/evidence-002/results/manifest.json`

## Limits

- The repository exposes nearby regression surfaces, but this bundle does not yet claim lecture-specific coefficient parity.
