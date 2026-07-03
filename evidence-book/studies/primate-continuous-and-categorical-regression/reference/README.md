# PCM3 Runtime Reference Mapping

The repository already exposes several runtime surfaces that line up with the
PCM3 lecture family, even though a lecture-specific governed parity bundle has
not yet been checked in.

Relevant owned surfaces include:

- `compare_comparative_regression_models(...)` for model ranking and shared
  exclusion accounting
- `summarize_pgls_categorical_contrasts(...)` for categorical encoding and
  reference-level interpretation
- `summarize_pgls_interaction_coefficients(...)` for explicit interaction-term
  ledgers
- `summarize_phylogenetic_residuals(...)` for reviewer-facing diagnostic
  summaries
- `summarize_ou_covariance_pgls(...)` and `summarize_brownian_covariance_pgls(...)`
  for covariance-structure review around comparative regression fits

Existing tests that demonstrate nearby runtime support include:

- `packages/bijux-phylogenetics/tests/test_pgls_interaction_coefficients.py`
- `packages/bijux-phylogenetics/tests/test_reference_parity.py`
- `packages/bijux-phylogenetics/tests/test_regression_model_selection.py`
- `packages/bijux-phylogenetics/tests/test_comparative_cli.py`

The missing piece is not a total lack of runtime support. The missing piece is
the governed lecture-parity packaging that turns those owned surfaces into one
PCM3 evidence program with checked-in study outputs.
