# Cichlid Phylogenetic Generalized Linear Mixed Models

This study captures PCM5, the Lund lecture on phylogenetic generalized linear
mixed models for count data, binomial response data, repeated measures, and
posterior-based phylogenetic signal calculations.

The repository now carries a complete study dossier for this lecture family.
That means the lecture sources, datasets, evidence decomposition, and runtime
boundaries are all explicit. It does not mean the repository already owns a
governed `MCMCglmm`-equivalent PGLMM surface.

## Study Metadata

- study id: `cichlid-phylogenetic-generalized-linear-mixed-models`
- lecture family: `PCM5_PGLMM`
- deployment state: `study_dossier_deployed`
- current runtime posture:
  `generalized_phylogenetic_mixed_model_parity_open_with_partial_binary_regression_support`
- evidence rule: isolate phylogenetic signal, ancestral-state extraction,
  Poisson count models, multinomial infection models, probability-scale
  back-transforms, and repeated-measure random effects into separate review
  slices

## Study Story

PCM5 extends the cichlid mixed-model story in three directions:

- phylogenetic signal becomes a posterior-derived quantity rather than a
  likelihood-ratio or lambda-only estimate
- ancestral-state reconstruction is expressed as posterior node effects
- generalized models are fitted for Poisson counts and binomial-style
  infection outcomes, including a repeated-measure species effect

That combination makes PCM5 the most demanding comparative lecture family in
the evidence book. It combines phylogeny, non-Gaussian likelihoods, posterior
uncertainty, and repeated measures in one workflow.

## Evidence Map

| Evidence ID | Focus | Current repository state |
| --- | --- | --- |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-001` | Fish workspace reload and repeated-measure data contract | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-002` | Posterior-derived phylogenetic signal for the Gaussian mixed model | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-003` | Ancestral-state posterior extraction and export rules | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-004` | Tree-linked ancestral plotting contract | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-005` | Poisson phylogenetic count model for parasite abundance | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-006` | Multinomial infection model and logit-scale summaries | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-007` | Probability-scale back-transformation for intercepts and slopes | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-008` | Repeated-measure random effects for `species + treetip` | Dossier defined |
| `cichlid-phylogenetic-generalized-linear-mixed-models/evidence-009` | Variance partition between phylogeny, species, and residual components | Dossier defined |

## Ideal Comparative Checks

- Keep the Gaussian, Poisson, and multinomial workflows separate enough that
  link-function changes remain visible.
- Preserve the lecture's distinction between posterior summaries on the logit
  scale and derived summaries on the probability scale.
- Track repeated-measure handling as its own evidence slice instead of letting
  it disappear inside one omnibus PGLMM run.
- Record how phylogenetic signal is computed from posterior variance
  components, including the different denominator used for binomial-style
  models.

## Elite Comparative Checks

- Add governed ledgers for posterior node estimates, probability-scale slope
  transforms, and variance partitions.
- Compare observation-level and species-level models without silently changing
  the phylogenetic denominator or the random-effect story.
- Distinguish clearly between repository-owned binary comparative regression
  approximations and the external-engine generalized mixed-model lecture path.
- Keep the study readable to a reviewer who wants to understand the scientific
  workflow without already knowing `MCMCglmm`.

## Checked-In Study Inputs

- [datasets](./datasets/)
- [provenance](./provenance/)
- [runtime reference mapping](./reference/README.md)
