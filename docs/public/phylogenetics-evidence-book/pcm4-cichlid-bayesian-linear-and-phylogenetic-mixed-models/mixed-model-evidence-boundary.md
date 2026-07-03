---
title: PCM4 Mixed-Model Evidence Boundary
audience: public
type: study-guide
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# PCM4 Mixed-Model Evidence Boundary

PCM4 deserves a direct boundary page because it is easy to over-read. The
lecture uses a Bayesian mixed-model workflow whose main scientific outputs are
posterior summaries, convergence checks, and phylogenetic variance estimates.

## What Is Source-Backed And Deployed

- the cichlid dataset and tree provenance
- the species-level summarization step used by the lecture
- the decomposition of the lecture into smaller evidence slices
- the Brownian PGLS comparator named at the end of the lecture

## What Is Not Yet A Repository-Owned Parity Claim

- Gaussian `MCMCglmm` posterior sampling
- chain-convergence diagnostics equivalent to the lecture
- governed posterior ledgers for phylogenetic random effects
- governed credible-interval summaries for the mixed-model parameters

## Why The Boundary Is Useful

Saying this directly helps readers separate two very different questions:

- Has the repository documented the lecture family well enough to review it?
- Has the repository already reproduced the lecture's Bayesian posterior
  outputs under governed evidence?

For PCM4 today, the answer is yes to the first question and no to the second.
