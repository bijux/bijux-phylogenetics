# Primate baseline regression structure bundle

Governed study slice for the PCM3 baseline GLS fit with one continuous predictor and one categorical predictor.

- study: `Primate Continuous And Categorical Regression`
- study id: `primate-continuous-and-categorical-regression`
- evidence id: `evidence-002`
- comparison mode: `direct_parity`
- current verdict: `not_comparable`

## Why This Bundle Exists

This bundle is checked in and reviewable, but lecture-specific PCM3 parity for this slice is not yet claimed by the repository.

## Source Fragments

- `baseline-gls-fit`
- `categorical-encoding-contract`

## Current Repository Posture

- `comparative_regression_surface_available_but_lecture_parity_not_yet_governed`

## Source Basis

- `evidence-book/studies/primate-continuous-and-categorical-regression/provenance/lund-course-sources.json`
- `evidence-book/studies/primate-continuous-and-categorical-regression/datasets/registry.json`
- `evidence-book/studies/primate-continuous-and-categorical-regression/reference/README.md`
- `packages/bijux-phylogenetics/tests/test_pgls_interaction_coefficients.py`
- `packages/bijux-phylogenetics/tests/test_regression_model_selection.py`

## Current Limits

- The repository exposes nearby regression surfaces, but this bundle does not yet claim lecture-specific coefficient parity.

## Next Step

- Check in coefficient and model-fit ledgers for the baseline categorical regression surface.
