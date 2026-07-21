---
title: Evidence Honesty And Limits
audience: public
type: explanation
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Evidence Honesty And Limits

Evidence earns trust when readers can see both what it supports and where it
stops. The project therefore treats unexplained mismatches, non-comparable
claims, open follow-through, and external blockers as useful information rather
than wording problems.

## Claim Anatomy

```mermaid
flowchart LR
    claim["Bounded claim"] --> source["Source and rationale"]
    source --> inputs["Identified inputs"]
    inputs --> execution["Recorded method and versions"]
    execution --> checks["Explicit checks and tolerances"]
    checks --> verdict["Comparison verdict"]
```

If a link is absent, say so. A narrative that leaps from source material to a
passing verdict without inputs and checks is not a complete evidence chain.

## Governed Verdict Meanings

| Evidence verdict | Meaning |
| --- | --- |
| `matched` | every declared comparable observation satisfies its exact or structural rule |
| `matched_with_tolerance` | declared numerical differences remain within explicit tolerances |
| `mismatch_explained` | a mismatch is retained and attributed to a governed assumption or implementation difference |
| `mismatch_unexplained` | a required comparable observation disagrees without a resolved cause |
| `not_comparable` | execution, assumptions, ownership, or outputs do not support the proposed comparison |

Resolution workflow is a separate axis. A mismatch or non-comparable bundle
may have open follow-through, a named blocker, or an accepted bounded state.
Do not replace the evidence verdict with that workflow label. Verdicts apply to
registered claims and bundles, not to an entire package or study by contagion.

## Observation Status And Bundle Verdict

A bundle can contain several observation statuses and one governed summary
verdict. The summary must account for mismatches, plot-only fragments,
unsupported cases, and missing comparisons; it cannot report only the passing
subset. Read the observation ledger before interpreting the bundle headline.

## Common Overclaims

- **Parity** without naming compared outputs, conventions, and tolerance.
- **Reproduced** when only the Python-side preparation ran.
- **Native** for a workflow whose computation is owned by an external engine.
- **Validated** because nearby methods have evidence.
- **Current** when the bundle predates material runtime or input changes.
- **Scientifically correct** because an execution is deterministic.

## External Boundaries

When an engine, posterior workflow, source dataset, or lecture-specific
interpretation remains external, record that boundary in the claim and bundle.
An adapter can be fully implemented and tested while the end-to-end study
verdict remains `not_comparable`. Those are compatible facts.

## Freshness

Evidence is tied to code, configuration, inputs, and sources. A change to any
material dependency can make a prior verdict stale without making it
historically false. Freshness reports should identify that state; documentation
should not keep presenting an old verdict as current.

## Negative Results

Failed comparisons, fragile examples, and sensitivity are part of the
scientific record. Preserve the configuration and observed difference, explain
whether the cause is understood, and define the next closure criterion. Do not
delete a difficult case from an aggregate without a declared exclusion rule.

## How To Write A Public Claim

Prefer: “For the named PCM study bundle, the implementation agrees with the
reference output within the declared tolerance; the posterior extension remains
external and non-comparable.” Avoid: “The comparative workflow has full parity.”

The first statement lets a reader verify strength and limitation. The second
borrows confidence from unspecified work.

## Reviewer Standard

Trace the claim ID to its bundle, source provenance, inputs, runtime versions,
checks, and verdict. Confirm it is fresh. Then read exclusions and open
boundaries before the headline result. If the chain cannot be followed, treat
the statement as documentation, not evidence.
