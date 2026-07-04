# Primate scaling and intercept semantics bundle

Governed study slice for PCM3 predictor scaling and the intercept-semantics changes that follow standardization.

- study: `Primate Continuous And Categorical Regression`
- study id: `primate-continuous-and-categorical-regression`
- evidence id: `evidence-006`
- comparison mode: `direct_parity`
- current verdict: `not_comparable`

## Why This Bundle Exists

This bundle is checked in and reviewable, but lecture-specific PCM3 parity for this slice is not yet claimed by the repository.

## Source Fragments

- `predictor-scaling`
- `standardized-intercept-semantics`

## Current Repository Posture

- `comparative_regression_surface_available_but_lecture_parity_not_yet_governed`

## Source Basis

- `evidence-book/studies/primate-continuous-and-categorical-regression/provenance/lund-course-sources.json`
- `evidence-book/studies/primate-continuous-and-categorical-regression/datasets/registry.json`
- `evidence-book/studies/primate-continuous-and-categorical-regression/reference/README.md`
- `packages/bijux-phylogenetics/tests/test_pgls_interaction_coefficients.py`
- `packages/bijux-phylogenetics/tests/test_regression_model_selection.py`

## Current Limits

- This bundle keeps the standardization story explicit without yet claiming checked-in parity for the lecture’s scaled outputs.

## Next Step

- Add governed before-and-after ledgers for scaled predictors, standardized coefficients, and intercept reinterpretation.
