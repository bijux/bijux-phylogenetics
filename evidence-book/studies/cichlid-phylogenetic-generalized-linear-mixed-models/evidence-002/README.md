# Cichlid posterior phylogenetic signal bundle

Governed study slice for the PCM5 posterior-derived phylogenetic signal calculation on the Gaussian mixed model.

- study: `Cichlid Phylogenetic Generalized Linear Mixed Models`
- study id: `cichlid-phylogenetic-generalized-linear-mixed-models`
- evidence id: `evidence-002`
- comparison mode: `direct_parity`
- current verdict: `not_comparable`

## Why This Bundle Exists

This bundle is checked in and reviewable, but the repository does not yet claim governed parity for the PCM5 generalized posterior workflow.

## Source Fragments

- `posterior-phylogenetic-signal`

## Current Repository Posture

- `generalized_phylogenetic_mixed_model_parity_open_with_partial_binary_regression_support`

## Source Basis

- `evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/provenance/lund-course-sources.json`
- `evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/datasets/registry.json`
- `evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/reference/README.md`
- `packages/bijux-phylogenetics/tests/test_phylogenetic_logistic.py`

## Current Limits

- This bundle keeps the posterior signal calculation explicit without claiming repository-owned parity for the variance-component workflow.

## Next Step

- Add governed variance-component and posterior-signal ledgers with the lecture denominator rules preserved.
