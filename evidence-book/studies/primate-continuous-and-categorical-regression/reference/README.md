# PCM3 Runtime Reference Mapping

The repository already exposes several runtime surfaces that line up with the
PCM3 lecture family, and the evidence book now checks in governed parity bundle
IDs for each PCM3 evidence slice. What remains open is the lecture-specific
parity verdict for those slices.

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
the lecture-specific coefficient, diagnostic, and hypothesis-test parity that
would turn the deployed governed parity bundles into a closed PCM3 parity
program.

Checked-in study-local lecture source copies:

- `lund_pcm3_contcat_source.R`
- `../datasets/primate_raw.xlsx`
- `../datasets/primate.RData`
- `../datasets/primatetree.nex`
- `../datasets/reference_primate.csv`
- `../datasets/reference_trimmed_primatetree.nwk`
