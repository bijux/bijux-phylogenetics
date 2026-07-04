# PCM4 Runtime Reference Mapping

Governed parity bundle IDs are now checked in for each PCM4 evidence slice,
but no governed repository adapter currently reproduces the core `MCMCglmm`
posterior workflow used in PCM4.

What the repository can already support nearby:

- tree loading and tip-alignment validation
- standard comparative regression and diagnostic summaries
- Brownian covariance PGLS surfaces that match the lecture's comparator model

What remains outside the current repository-owned surface:

- posterior sampling for Gaussian mixed models
- credible-interval summaries for mixed-model fixed and random effects
- chain convergence diagnostics equivalent to the lecture's `MCMCglmm` review
- phylogenetic random-effect posterior ledgers

This dossier exists so the evidence-book tells the truth cleanly: PCM4 now has
source-backed study records and governed parity bundle IDs, but native or
governed parity for its Bayesian mixed-model core is not yet present.

Checked-in study-local lecture source copies:

- `lund_pcm4_mixed_model_source.R`
- `../datasets/cichlid.csv`
- `../datasets/cichlidtree.nex`
- `../datasets/fishdata.RData`
