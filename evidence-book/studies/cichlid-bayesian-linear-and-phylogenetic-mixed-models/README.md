# Cichlid Bayesian Linear And Phylogenetic Mixed Models

This study captures PCM4, the Lund lecture that moves from non-phylogenetic
Gaussian regression into Bayesian linear mixed models and then into a
phylogenetic random-effect formulation on cichlid parasite data.

The repository now carries a complete study dossier for this lecture family
with governed parity evidence bundles for each lecture slice: dataset
registry, source provenance, evidence decomposition, explicit runtime
boundaries, and bundle-local review artifacts are all checked in. What it does
not claim is full governed parity with the `MCMCglmm` posterior workflow.

## Study Metadata

- study id: `cichlid-bayesian-linear-and-phylogenetic-mixed-models`
- lecture family: `PCM4_PLMM`
- deployment state: `governed_parity_bundle_program_deployed_with_external_posterior_boundary`
- governed evidence ids: `evidence-001` through `evidence-007`
- current runtime posture:
  `bayesian_mixed_model_parity_open_with_pgls_comparator_surface_only`
- evidence rule: separate data-tree alignment, Gaussian baseline fits,
  posterior diagnostics, phylogenetic variance components, and comparator
  models into distinct evidence slices

## Study Story

PCM4 introduces a new biological system and a new modeling logic:

- raw cichlid observations must be aligned to the tree
- repeated observations are summarized to one species-level table
- a standard Gaussian regression provides the baseline
- `MCMCglmm` adds posterior sampling and random-effect uncertainty
- a phylogenetic random effect is introduced through the inverse covariance
  matrix of the tree
- the resulting Bayesian fit is compared back to a Brownian PGLS baseline

That lecture is central because it is the first place where the evidence book
must distinguish between a workflow the repository can already reproduce and a
workflow whose main scientific claim still lives in an external engine.

## Evidence Map

| Evidence ID | Focus | Current repository state |
| --- | --- | --- |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-001` | Cichlid data-tree alignment and species-level summarization | Governed parity bundle deployed; posterior verdict still open |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-002` | Non-phylogenetic Gaussian baseline fit | Governed parity bundle deployed; posterior verdict still open |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-003` | Bayesian Gaussian regression chain diagnostics and posterior summaries | Governed parity bundle deployed; posterior verdict still open |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-004` | Tree inverse construction and phylogenetic random-effect setup | Governed parity bundle deployed; posterior verdict still open |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-005` | Bayesian phylogenetic mixed-model fixed-effect posteriors | Governed parity bundle deployed; posterior verdict still open |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-006` | Phylogenetic variance-component posteriors | Governed parity bundle deployed; posterior verdict still open |
| `cichlid-bayesian-linear-and-phylogenetic-mixed-models/evidence-007` | Comparison back to Brownian PGLS | Governed parity bundle deployed; posterior verdict still open |

## Ideal Comparative Checks

- Verify that cichlid species names match the tree before any summarization or
  posterior run is accepted.
- Keep the non-phylogenetic Gaussian baseline separate from the Bayesian
  posterior workflow so disagreements remain interpretable.
- Record priors, chain settings, convergence checks, and credible intervals as
  first-class study outputs rather than prose summaries.
- Compare the phylogenetic mixed model against the Brownian PGLS baseline on
  shared response and predictor definitions.

## Elite Comparative Checks

- Preserve posterior summaries and variance-component ledgers in a format that
  can be re-audited without rerunning the lecture from scratch.
- Track the construction of the inverse phylogenetic covariance as its own
  evidence slice rather than letting it disappear inside one MCMC run.
- Distinguish repository-owned comparator surfaces from external-engine
  posterior surfaces so the study never over-claims native support.
- Add later parity bundles only when the repository can explain posterior
  diagnostics and uncertainty as clearly as the lecture does.

## Checked-In Study Inputs

- [datasets](./datasets/)
- [provenance](./provenance/)
- [runtime reference mapping](./reference/README.md)
- governed parity bundles: `evidence-001` to `evidence-007`
