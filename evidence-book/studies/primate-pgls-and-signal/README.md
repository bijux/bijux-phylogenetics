# Primate PGLS And Signal

This study turns the regression, transformed-tree, evolutionary-mode fit,
likelihood-ratio, and ancestral sections of the Lund primate comparative
lecture into governed evidence slices backed by checked-in R reference outputs
and repository reproductions.

It is the deeper primate methods study: comparative model fitting, PGLS,
diagnostics, and explicit open boundaries all live here.

## Study Metadata

- study id: `primate-pgls-and-signal`
- lecture family: `PCM2_modes_pgls`
- deployment state: `governed_parity_deployed`
- evidence rule: keep every major lecture claim in its own bundle, and treat
  any uncovered fragment as an explicit boundary rather than silent scope

## Study Story

PCM2 asks whether the repository can stay aligned with the lecture once the
workflow becomes statistically richer:

- transformed trees change the covariance structure
- baseline GLS becomes phylogenetically corrected regression
- evolutionary mode fitting and likelihood ratios become central
- signal and residual review must still remain interpretable
- ancestral results need to stay tied to the same model family

## Evidence Map

| Evidence ID | Focus |
| --- | --- |
| `primate-pgls-and-signal/evidence-001` | Reload semantics |
| `primate-pgls-and-signal/evidence-002` | Baseline GLS parity |
| `primate-pgls-and-signal/evidence-003` | Pagel-lambda regression parity |
| `primate-pgls-and-signal/evidence-004` | Phylogenetic signal parity |
| `primate-pgls-and-signal/evidence-005` | Residual diagnostics parity |
| `primate-pgls-and-signal/evidence-006` | Transformed-tree parity |
| `primate-pgls-and-signal/evidence-007` | Evolutionary mode-fit parity |
| `primate-pgls-and-signal/evidence-008` | Likelihood-ratio parity |
| `primate-pgls-and-signal/evidence-009` | Ancestral-mode parity |
| `primate-pgls-and-signal/evidence-010` | Intercept-sweep coverage boundary |

## Ideal Comparative Checks

- Compare baseline and phylogenetically corrected regressions on shared taxa
  and shared response definitions.
- Keep transformed-tree, model-fit, and likelihood-ratio evidence separate
  enough that disagreements remain localized.
- Preserve diagnostics as evidence, not as an afterthought to the coefficient
  table.

## Elite Comparative Checks

- Make the one remaining uncovered lecture fragment visible as its own bundle
  instead of burying it in caveats.
- Keep ancestral reconstruction tied to the same reviewed model family rather
  than treating it as a detached plotting step.
- Reuse the governed primate derivatives from PCM1 so cross-study agreement is
  based on explicit shared inputs.

## Study Surfaces

- [datasets](./datasets/)
- [provenance](./provenance/)
- [reference adapters](./reference/)
