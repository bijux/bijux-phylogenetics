---
title: R Ecosystem Correspondence
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# R Ecosystem Correspondence

Three native parity harnesses execute checked-in R runners against live APE,
GEIGER, and phytools installations. Additional study fixtures retain narrower
R-derived comparisons, including `ape` plus `nlme` PGLS outputs and course
workflows that use packages such as `MCMCglmm` without claiming a live
drop-in-equivalence surface.

| Reference | Strongest governed relationship | Main scientific territory |
| --- | --- | --- |
| [APE](ape.md) | live case registry plus checked fixtures | trees, topology, distances, DNA utilities, covariance, contrasts, ancestral states |
| [GEIGER](geiger.md) | live model registry plus generated review report | continuous and discrete evolutionary-model fitting |
| [phytools](phytools.md) | live case registry plus stochastic-summary policies | signal, Mk models, simulation, stochastic maps, ancestral estimates |
| [`nlme`](study-references.md#ape-plus-nlme) | checked reference outputs with APE covariance | GLS/PGLS coefficients and model comparison |
| [`treeio`, `tidytree`, `ggtree`](study-references.md#treeio-tidytree-and-ggtree) | governed structural evidence or plot-only provenance | study-specific tree data and visualization surfaces |
| [`MCMCglmm`](study-references.md#mcmcglmm) | external study boundary | Bayesian phylogenetic mixed models in PCM4 and PCM5 |

## Choose The Reference From The Operation

| Operation needing correspondence | Reference owner | Bijux surface represented in the registry |
| --- | --- | --- |
| tree parsing, clades, rooting, pruning, MRCA, tree distance | APE | topology, I/O, pruning, and comparison contracts |
| DNA encoding, base frequencies, segregating sites, sequence distance | APE | alignment and native genetic-distance records |
| covariance, contrasts, foundational ancestral reconstruction | APE | comparative and ancestral structured results |
| BM, OU, early-burst, lambda, kappa, delta, or white-noise fits | GEIGER | continuous evolutionary-mode fitting and ranking |
| ER, SYM, or ARD discrete fits and tree transforms | GEIGER | finite-state Mk fitting with explicit state/rate identity |
| signal, specialized ancestral estimates, residuals, or phylogenetic ANOVA | phytools | signal, reconstruction, PGLS/residual, and group-summary records |
| stochastic maps and continuous/discrete simulation | phytools | event/history and collection summaries under distribution rules |
| study PGLS using APE covariance and `nlme::gls` | checked study fixture | claim-scoped coefficients, likelihoods, model comparisons, diagnostics |
| cichlid phylogenetic mixed-model posterior | `MCMCglmm` source boundary | no substitute surface; PCM4/PCM5 remain `not_comparable` |

The mapping is operation-specific. An APE tree case does not validate a
GEIGER optimizer, and a phytools stochastic summary does not validate an
`MCMCglmm` posterior. When several packages appear in one source workflow,
retain which package owns each calculated quantity.

## Correspondence Key

Every observation joins `reference function + Bijux function + case + fixture
+ normalization + rule`. Function names alone are insufficient: defaults,
state order, root treatment, likelihood convention, and missing-data behavior
can change the quantity while leaving the printed method name unchanged.

## Registry Lifecycle

```mermaid
flowchart LR
    registry["Versioned cases<br/>and shards"]
    selection["Reference source<br/>filters and denominator"]
    reference["Independent R<br/>execution"]
    bijux["Owned Bijux<br/>execution"]
    observations["Exact · structural<br/>tolerance · distribution"]
    report["Summary · observations<br/>failure artifacts"]

    registry --> selection
    selection --> reference
    selection --> bijux
    reference --> observations
    bijux --> observations --> report
```

A registry entry is the unit of eligibility. A shard is a durable grouping by
owned Python surface and reference operation. A run selects cases from that
registry; it does not discover successful examples after execution.

## Execution Boundary

Live parity is opt-in through `RUN_LIVE_PARITY=1` and requires `Rscript` plus
the named package. The command emits structured summary and observation rows;
failed and skipped cases can be preserved under a failure-artifact directory.
The default parity command uses checked fixtures and does not imply that R ran
on the current machine.

## Read A Harness Result

| Field | Interpretation |
| --- | --- |
| reference source | proves whether the case used live APE, GEIGER, phytools, or a checked fixture |
| registry and shard identity | bounds the operation, model, and fixtures that were eligible for selection |
| selected count | denominator requested by the command after filters |
| passed count | cases whose registered exact, structural, tolerance, or distribution rules all held |
| failed count | executed cases with a Bijux/reference failure, missing output, or comparison mismatch |
| skipped count | selected cases not executed because the live reference environment was unavailable |
| mismatch reason and failure artifacts | preserve why a case did not pass and what is needed to reproduce it |

Report all four counts—selected, passed, failed, and skipped—with the selected
shards. Do not describe only the successful subset as the registry result.

## Execution Identity

Retain the Bijux version, registry case and shard IDs, reference-source label,
R executable identity, R version, package version, fixture hashes, filters,
environment gate, command arguments, start and finish state, and output hashes.
For stochastic cases, also retain the replicate population, seed policy,
summary definition, and invalid-replicate count.

The package version is observed from the live reference process. A current
installation does not retroactively identify an older checked fixture, and a
fixture cannot be relabelled as live because the corresponding package is now
installed.

## Treat The R Session As Part Of The Reference

An R function is not fully identified by package and function name. Live
correspondence must retain the session conditions that can change dispatch,
defaults, optimization, ordering, or stochastic output.

| Session component | Why it matters |
| --- | --- |
| R and package versions | define exported behavior, object structures, defaults and numerical dependencies |
| namespace-qualified call | prevents search-path masking and ambiguous method ownership |
| loaded dependency versions | can change matrix, optimizer, tree, serialization and compiled behavior |
| global options, contrasts and locale | affect factor coding, printing, parsing, ordering and numerical presentation |
| random-number kind and seed | define stochastic case generation and reference draws |
| platform and numeric libraries | can affect optimization paths, tolerances and floating-point results |
| warnings, errors and session metadata | preserve diagnostic behavior rather than only returned values |

Capture these fields from the reference process itself. Shell environment
claims or a later `sessionInfo()` from another process are insufficient. When
a session field changes materially, retain a new execution observation even
if every compared scalar still passes.

## Treat Reference Drift As A New Observation

An R, package, dependency, or platform change can alter defaults, object
representation, optimizer behavior, warnings, and numerical output. Registry
case identity remains stable, but its execution observation is versioned.

| Detected change | Required response |
| --- | --- |
| reference package or R version | execute the selected cases again and retain the prior observation |
| default, warning, or return structure | inspect the native R output before changing normalization or parsers |
| optimizer or stochastic behavior | keep seeds, starts, replicate population, convergence, and failure counts comparable |
| fixture or registry case | create a new case/input identity; do not overwrite the old denominator |
| tolerance or distributional rule | record a new rule revision and rationale before adjudicating new output |
| indirect dependency or platform | retain environment identity and qualify differences that cannot be isolated |

A prior passing observation remains evidence for its recorded environment. It
does not automatically fail or pass the new environment, and refreshing a
fixture without live execution does not establish current correspondence.

## Programmatic Entry Points

Each live registry exposes the same review shape from its public module:

| Module | Discover | Execute | Persist |
| --- | --- | --- | --- |
| `bijux_phylogenetics.parity.ape` | `list_ape_parity_cases`, `list_ape_parity_shards` | `run_ape_parity_cases` | summary and observation table writers |
| `bijux_phylogenetics.parity.geiger` | `list_geiger_parity_cases`, `list_geiger_parity_shards` | `run_geiger_parity_cases` | summary, observation, triage, policy, parameterization, and generated-report writers |
| `bijux_phylogenetics.parity.phytools` | `list_phytools_parity_cases`, `list_phytools_parity_shards` | `run_phytools_parity_cases` | summary and observation table writers |

The CLI composes these surfaces. Direct Python use does not weaken the same
selection, environment, status, and artifact requirements.

## Why Several Reference Packages Matter

No single R package owns every convention used in phylogenetics. APE is the
structural baseline, GEIGER is a model-fitting reference, and phytools provides
specialized comparative and stochastic-history behavior. Agreement with one
does not silently grant agreement with another.

The three harnesses share `passed`, `failed`, and `skipped` case status, but
their observation rules remain tool- and shard-specific. That common status
does not erase differences among exact tree structure, tolerance-based model
parameters, and distributional stochastic summaries.

[Study-specific R references](study-references.md) explains the narrower
relationships that enter through Evidence Book workflows rather than a live
package registry.
