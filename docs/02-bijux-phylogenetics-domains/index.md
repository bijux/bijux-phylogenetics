---
title: Scientific Domains
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Scientific Domains

The runtime is organized around scientific responsibilities rather than a
single end-to-end pipeline. Tree and alignment contracts support likelihood,
inference, comparative, ancestral, and parsimony methods. Simulation and
datasets provide controlled inputs. Validation, benchmarks, and reporting make
results inspectable.

```mermaid
flowchart TB
    inputs["Taxa · trees · alignments<br/>traits · partitions"]
    foundations["Topology · branch lengths<br/>state spaces · covariance"]
    likelihood["Likelihood, distance<br/>and Bayesian inference"]
    comparative["Comparative and ancestral methods"]
    histories["Ecology and phylogeography"]
    parsimony["Parsimony and topology search"]
    context["Simulation and datasets"]
    review["Validation · benchmarks<br/>reports · evidence"]

    inputs --> foundations
    foundations --> likelihood
    foundations --> comparative
    comparative --> histories
    foundations --> parsimony
    context --> likelihood
    context --> comparative
    context --> parsimony
    likelihood --> review
    comparative --> review
    histories --> review
    parsimony --> review
```

## Domain Map

| Domain | Core scientific objects | Representative outputs |
| --- | --- | --- |
| [Trees, alignments, and taxa](data-and-structure/trees-alignments-and-taxa.md) | topology, rootedness, branch lengths, support, sequences, partitions, identifiers | validated trees, reconciled alignments, topology and support ledgers |
| [Likelihood and inference](inference/likelihood-and-inference.md) | state spaces, rate matrices, likelihoods, priors, search and sampling | fitted models, inferred trees, search traces, posterior samples |
| [Distance-based inference](inference/distance-based-inference.md) | pairwise distances, matrix assumptions, joining and least-squares objectives | inferred trees, residuals, influence ledgers, bootstrap support |
| [Bayesian inference](inference/bayesian-inference.md) | likelihoods, priors, proposals, chains, clocks and topology distributions | run manifests, checkpoints, posterior samples, diagnostics and predictive summaries |
| [Comparative and ancestral methods](evolutionary-analysis/comparative-and-ancestral-methods.md) | trait covariance, signal, regression, evolutionary modes, node states | coefficients, model comparisons, uncertainty and ancestral ledgers |
| [Ecology and phylogeography](evolutionary-analysis/ecology-and-phylogeography.md) | host, niche, region and coordinate histories on trees | transition, movement, outlier, exclusion and map ledgers |
| [Parsimony](inference/parsimony.md) | character costs, ancestral sets, topology moves, resampling | scores, reconstructions, consensus and support summaries |
| [Simulation and datasets](scientific-assurance/simulation-and-datasets.md) | generative models, known truths, reference and stress corpora | simulated trees/traits/alignments, dataset manifests |
| [Validation, benchmarks, and reporting](scientific-assurance/validation-benchmarks-and-reporting.md) | invariants, tolerances, metrics, provenance, publication contracts | diagnostics, benchmark reports, manifests, tables and figures |
| [Diagnostics and sensitivity](scientific-assurance/diagnostics-and-sensitivity.md) | tree assumptions, temporal signal, alternative trees, models and policies | findings, stability rows, refusal reasons and decision records |

## Choose By Scientific Question

| Question | Begin in | Continue into |
| --- | --- | --- |
| Are the taxa, topology, lengths, states, and partitions coherent? | [data and structure](data-and-structure/index.md) | the method that consumes the validated objects |
| Which topology or model parameters best explain observed characters? | [inference](inference/index.md) | support, sensitivity, and assurance review |
| How do traits vary with ancestry, predictors, or internal nodes? | [evolutionary analysis](evolutionary-analysis/index.md) | model comparison and tree/model sensitivity |
| Is an implementation recovering known truth or matching a reference? | [scientific assurance](scientific-assurance/index.md) | benchmark or Evidence Book artifacts |

Scientific questions commonly cross domains, but ownership remains local. A
taxon mismatch is resolved before model fitting; a convergence failure is not
recast as a reporting concern; and a visually persuasive figure cannot replace
an uncertainty or sensitivity record.

## Domain Handoffs

```mermaid
flowchart LR
    intake["Data and structure<br/>validated objects"]
    inference["Inference<br/>trees and parameters"]
    evolution["Evolutionary analysis<br/>traits and histories"]
    assurance["Scientific assurance<br/>checks and evidence"]

    intake --> inference
    intake --> evolution
    inference --> evolution
    inference --> assurance
    evolution --> assurance
    assurance -. failures and limits .-> intake
    assurance -. failures and limits .-> inference
    assurance -. failures and limits .-> evolution
```

Each handoff carries identities and assumptions, not only values. Tree scale,
rooting, taxon order, state encoding, model parameterization, search or sampler
configuration, and software provenance remain attached to downstream results.

## Cross-Domain Handoff Record

| Field | Why the next domain needs it |
| --- | --- |
| object identity | prevents a tree, matrix, alignment, or dataset revision from being substituted silently |
| taxon population and order | preserves joins, pruning decisions, matrix row meaning, and denominators |
| units and scale | distinguishes time, substitutions, distance, rates, transformed traits, and probabilities |
| model and state conventions | fixes alphabets, ambiguity, partitions, covariance, transition rules, and parameterization |
| execution state | distinguishes complete, partial, refused, failed, and externally unavailable work |
| diagnostics and exclusions | prevents downstream presentation from retaining only successful values |
| provenance | identifies implementation or engine owner, version, configuration, and source inputs |

A downstream method may add context, but it must not erase upstream identity or
reinterpret an unresolved condition as an ordinary value.

## Recognize When A Handoff Changes The Question

Some transformations produce a new scientific object or estimand rather than
merely preparing the old one. The downstream record must name that change so a
result is not attributed to the original question unchanged.

| Handoff decision | Original question | New question that must be declared |
| --- | --- | --- |
| prune to shared taxa | relationship on the declared taxa | relationship on the retained intersection |
| trim sites or exclude sequences | inference from the admitted character matrix | inference from a selected matrix under a removal policy |
| root, reroot, or rescale a tree | topology or process under the supplied tree identity | topology, direction, or rate under a transformed tree |
| standardize, log, aggregate, or impute traits | effect on the observed response and population | effect on a transformed or reconstructed quantity |
| select a model after screening candidates | fit of one predeclared model | selection-conditional result over a candidate set |
| summarize across trees, chains, or replicates | result for one execution state | distribution or consensus over a declared population |
| normalize outputs across implementations | engine- or method-specific quantity | correspondence under an explicit mapping rule |

The transformation can be scientifically valid. The defect is allowing the
new population, scale, conditioning event, or estimand to inherit the old
label without a parent-child record and an interpretation update.

## Attribute Failure To The Deciding Domain

Parsing and identity defects belong to data and structure. Undefined
likelihoods, search exhaustion, or poor chain behavior belong to inference.
Model sensitivity and unidentified covariance or transition assumptions belong
to evolutionary analysis. Missing independent truth, reference incompatibility,
or incomplete denominators belong to scientific assurance. Reports carry these
states; they do not become their owner.

## Shared Scientific Invariants

- Taxon identity must remain stable across every joined object.
- Tree rootedness and branch-length meaning must match the chosen method.
- Character or sequence state space must be explicit.
- Model parameterization and transformation conventions must be retained with
  estimates.
- Search, optimization, and sampling diagnostics qualify final values.
- Structured results are canonical; figures and prose are projections.
- Validation scope is named rather than inherited from a neighboring method.

## Native And External Ownership

The package owns substantial tree, likelihood, inference, comparative,
ancestral, parsimony, simulation, validation, and reporting implementations.
It also integrates MAFFT, trimAl, IQ-TREE2, FastTree, MrBayes, and BEAST where
those programs remain the computational owner. Both paths converge on typed
results and artifacts, but their installation, provenance, and failure
boundaries remain different.

## Connecting Domains To Evidence

Scientific-domain coverage describes what the runtime implements. The
[Evidence Book](../03-bijux-phylogenetics-evidence-book/index.md) records
bounded study comparisons. A domain can be mature in code and tests without a
course-study bundle; an evidence bundle can depend on an external engine
without claiming native ownership.

For cross-implementation questions, use the
[verification matrix](../04-bijux-phylogenetics-parities/correspondence/verification-matrix.md)
to distinguish live parity, checked fixtures, and governed external execution.
