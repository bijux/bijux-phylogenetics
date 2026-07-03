# Primate Continuous And Categorical Regression

This study captures PCM3, the Lund lecture on multivariable comparative
regression with continuous predictors, categorical predictors, interactions,
and phylogenetically corrected regression.

It is deployed here as a complete study dossier: source provenance, dataset
registry, evidence decomposition, and runtime-surface mapping are all checked
in. What is not yet deployed is a full governed lecture-parity bundle for each
evidence slice.

## Study Metadata

- study id: `primate-continuous-and-categorical-regression`
- lecture family: `PCM3_contcat`
- deployment state: `study_dossier_deployed`
- current runtime posture:
  `comparative_regression_surface_available_but_lecture_parity_not_yet_governed`
- evidence rule: keep each lecture claim separate enough to show whether the
  repository reproduces the data contract, the fitted coefficients, the
  diagnostics, or only the surrounding workflow

## Study Story

PCM3 shifts from single-predictor examples into multivariable model design.
The lecture asks whether primate longevity can be modeled with:

- one continuous predictor
- one categorical predictor
- explicit interaction terms
- standardized predictors that change the intercept interpretation
- a phylogenetically corrected covariance structure through `corMartins`

That makes PCM3 the bridge between the earlier primate signal studies and the
later mixed-model lectures.

## Evidence Map

| Evidence ID | Focus | Current repository state |
| --- | --- | --- |
| `primate-continuous-and-categorical-regression/evidence-001` | Source intake, primate workspace reload, and shared dataset alignment | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-002` | Baseline GLS with one continuous and one categorical predictor | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-003` | Log-response repair, residual diagnostics, and marginal versus sequential tests | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-004` | Intercept-free parameterization for group-specific estimates | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-005` | Continuous-by-categorical interaction coefficients and interpretation | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-006` | Predictor scaling and intercept semantics after standardization | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-007` | Phylogenetic regression with categorical predictors and `corMartins` covariance | Dossier defined |
| `primate-continuous-and-categorical-regression/evidence-008` | Post-correction diagnostics and hypothesis testing | Dossier defined |

## Ideal Comparative Checks

- Confirm that the shared primate inputs match the lecture workspace and the
  earlier PCM1 and PCM2 governed derivatives.
- Reproduce categorical encoding, reference-level handling, and interaction
  expansion exactly enough to compare coefficients term by term.
- Track the effect of the log transform and the scaled predictor on both
  diagnostics and coefficient interpretation.
- Verify that the phylogenetic correction changes the regression story for the
  right reasons, not because taxa were silently excluded or columns were
  encoded differently.

## Elite Comparative Checks

- Keep coefficient ledgers, categorical-contrast ledgers, and interaction
  ledgers alongside the main parity report.
- Compare baseline GLS, Pagel-lambda PGLS, and Martins-covariance PGLS as one
  continuous model family instead of isolated runs.
- Preserve the exact lecture distinction between workflow-level plots and
  claim-worthy numeric surfaces.
- Trace each evidence slice to the corresponding runtime surface so later
  parity bundles can be implemented without rediscovering the study structure.

## Checked-In Study Inputs

- [datasets](./datasets/)
- [provenance](./provenance/)
- [runtime reference mapping](./reference/README.md)
