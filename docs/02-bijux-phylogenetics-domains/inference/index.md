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

## Select The Inferential Contract

| If the question is… | Required declaration | Do not omit |
| --- | --- | --- |
| Which parameterized model has the strongest support? | likelihood model, candidate set, fitting policy | comparable likelihood convention and parameter counts |
| Which topology maximizes likelihood? | substitution model, partitions, starting trees, moves | search coverage, convergence, branch-length constraints |
| What posterior uncertainty remains? | likelihood, priors, proposals, chains, retention policy | burn-in, acceptance, ESS, topology diversity |
| Which tree minimizes character-change cost? | coding, cost model, ancestral constraints | equal-best trees, topology moves, search limits |
| How stable is the selected tree? | resampling or posterior support method | replicate count, support meaning, failed replicates |

The same tree can be an input, a search state, a point estimate, or a sample
from uncertainty. Preserve that role in result and artifact labels.

## Review The Search, Not Only The Winner

A final topology or score is insufficient when different starts, moves, or
chains explored materially different regions. Retain accepted and rejected
moves where supported, convergence and termination reasons, equal or near-equal
solutions, support replicates, and warnings. These records distinguish a
well-supported optimum from the best result found under a narrow search.

## Shared Limits

A reported optimum is conditional on the model and explored search space. A
posterior is conditional on priors, proposals, and retained samples. A most
parsimonious tree is conditional on coding, costs, and search coverage. None of
these labels alone proves a historically certain topology.

Inference should refuse interpretation when the model is undefined for the
input, optimization terminates without a usable solution, chains lack adequate
mixing, or support generation is incomplete. A partial artifact may remain
valuable for diagnosis, but its status must remain partial.
