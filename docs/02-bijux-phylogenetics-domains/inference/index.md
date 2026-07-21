---
title: Inference
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Inference

Inference connects observed characters to parameters, ancestral assignments,
and topology. The runtime distinguishes probabilistic objectives from
parsimony cost models and preserves the search or sampling record needed to
interpret either one.

| Family | Objective | Essential review evidence |
| --- | --- | --- |
| maximum likelihood | maximize data likelihood under a named model | starts, moves, convergence, fitted parameters, support method |
| Bayesian inference | characterize a posterior under likelihood and priors | proposals, acceptance, burn-in, ESS, chains, topology distribution |
| parsimony | minimize declared character-change cost | state/cost model, search trace, equal-best trees, consensus, support |

Read [Likelihood And Inference](likelihood-and-inference.md) for finite-state
likelihood, optimization, model comparison, and posterior sampling. Read
[Parsimony](parsimony.md) for Fitch, Wagner, Sankoff, Dollo, Camin–Sokal,
topology search, reconstruction, and support.

## Shared Limits

A reported optimum is conditional on the model and explored search space. A
posterior is conditional on priors, proposals, and retained samples. A most
parsimonious tree is conditional on coding, costs, and search coverage. None of
these labels alone proves a historically certain topology.
