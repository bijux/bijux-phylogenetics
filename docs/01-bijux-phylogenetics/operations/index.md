---
title: Runtime Operations
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Runtime Operations

Operating scientific software means preserving the reasoning around a run,
not merely obtaining a zero exit status. These guides connect environment,
inputs, method choice, execution, diagnostics, artifacts, and review.

```mermaid
flowchart LR
    install["Install and identify<br/>versions"] --> prepare["Validate inputs<br/>and assumptions"]
    prepare --> execute["Execute native method<br/>or adapter"]
    execute --> inspect["Inspect diagnostics<br/>and sensitivity"]
    inspect --> preserve["Preserve manifests<br/>and review artifacts"]
```

## Start By Outcome

| Outcome | Guide |
| --- | --- |
| create a reliable environment | [Installation and setup](installation-and-setup.md) |
| move from input data to a reviewable result | [Common workflows](common-workflows.md) |
| run owned maximum-likelihood inference | [Native maximum-likelihood workflows](native-maximum-likelihood-workflows.md) |
| run a supported Bayesian DNA path | [Native Bayesian workflows](native-bayesian-workflows.md) |
| interpret benchmark records | [Native benchmark review](native-benchmark-review.md) |
| understand what automation cannot decide | [Operational boundaries](operational-boundaries.md) |

## Establish The Run Identity First

Before selecting a model, record the input set and the runtime that will own
the result:

```bash
bijux-phylogenetics --version
bijux-phylogenetics env inspect --json
bijux-phylogenetics inspect dataset/tree.nwk --json
bijux-phylogenetics alignment validate-input dataset/alignment.fasta --json
```

The last two commands are independent intake checks; use only the objects
present in the analysis. For trait work, validate the tree–table join and
retain the reconciliation output before fitting a model. For external engines,
capture executable discovery and version output before execution.

## Follow The Operational State Machine

```mermaid
stateDiagram-v2
    [*] --> Intake
    Intake --> Refused: invalid identity or assumptions
    Intake --> Ready: validated inputs and configuration
    Ready --> Running
    Running --> Partial: some declared components unavailable
    Running --> Failed: computation or infrastructure failure
    Running --> Complete: required outputs and diagnostics present
    Partial --> Review
    Failed --> Review
    Complete --> Review
    Refused --> Review
```

`Review` does not mean success. It means that the terminal state, diagnostics,
and output denominator are preserved well enough to understand what happened.

## Select One Accepted Execution

Retries, alternate models, and sensitivity runs form branches. They must not
write into a shared result identity and later appear to be one uninterrupted
successful execution.

```mermaid
flowchart LR
    parent["Validated inputs<br/>and declared question"]
    first["Execution A<br/>terminal record"]
    retry["Execution B<br/>changed seed or recovery"]
    alternate["Execution C<br/>changed model or data"]
    decision["Acceptance decision<br/>selected scope and reason"]

    parent --> first --> decision
    parent --> retry --> decision
    parent --> alternate --> decision
```

The acceptance decision identifies every candidate execution, the scientific
differences among them, the selected record, the rule used to select it, and
the conclusions excluded by failed or rejected candidates. A retry with
identical identity may demonstrate recoverability; a changed seed, model,
engine, tolerance, input, or incomplete-run policy is a new branch.

Preserve non-selected branches when they affect failure rates, sensitivity,
search coverage, or the credibility of the selection rule. Deleting an
unfavorable attempt changes the denominator and is not operational cleanup.

## The Operational Contract

A defensible run records:

- installed runtime and, where applicable, external-engine versions;
- immutable or checksummed input identities;
- taxon reconciliation and validation outcomes;
- model, partition, search, sampling, and seed configuration;
- warnings, convergence state, support, uncertainty, or sensitivity diagnostics;
- a declared output inventory and the command or call that produced it.

The exact fields vary by method. The invariant is that a reviewer should not
have to guess which assumptions created the result.

## Native And External Execution

Native workflows fail within the Python runtime and can return rich diagnostic
state directly. Adapter workflows add executable discovery, subprocess exit,
captured stdout/stderr, output parsing, and engine-version provenance. Plan for
those additional failure and reproducibility boundaries.

## Reproducibility Levels

| Level | What it establishes |
| --- | --- |
| command recorded | invocation can be inspected |
| environment recorded | software identity is known |
| inputs and configuration fixed | the run can be reconstructed |
| deterministic or seeded replay passes | the declared execution reproduces within its contract |
| reference or evidence checks pass | bounded comparison claims gain additional support |

Do not collapse these levels into “reproducible.” State which one applies.

## Before Interpretation

Review validation warnings, optimizer or sampler diagnostics, support and
uncertainty, missing-data behavior, taxon coverage, and sensitivity to material
choices. A result object existing is an execution fact; it is not an automatic
scientific verdict.

Use a refusal or incomplete status as information. A taxon mismatch, missing
engine, unconverged optimizer, low effective sample size, incomplete support
replicates, or absent required artifact should stop the corresponding claim.
Changing presentation or omitting the failed check does not make the run
complete.

## Minimum Handoff Record

Preserve enough state for a second person to answer these questions without
terminal history:

- Which exact inputs and taxa were analyzed?
- Which runtime and external-engine versions executed?
- Which model, partitions, transformations, priors, costs, and seeds applied?
- Which diagnostics passed, warned, failed, or remained unavailable?
- Which files are canonical machine output and which are presentations?
- Which comparison or Evidence Book claim, if any, supports interpretation?

## Generated References

The exhaustive [CLI surface](../interfaces/cli-surface.md) and generated parts
of [common workflows](common-workflows.md) are lookup references. Use this
operations section to decide *why* and *how* to run; use generated references
to confirm exact arguments and available methods.
