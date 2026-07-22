---
title: Claim Verdicts And Resolution
audience: public
type: evidence-guide
status: canonical
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-22
---

# Claim Verdicts And Resolution

A verdict reports the outcome of a registered comparison. Resolution reports
what happens next. Keeping them separate prevents open work from weakening a
real mismatch and prevents a detailed plan from masquerading as completed
evidence.

## Locate The Status Level

```mermaid
flowchart TB
    observation["Observation status<br/>one scalar, row or structure"] --> claim["Claim verdict<br/>one bounded statement"]
    claim --> bundle["Bundle summary<br/>registered claim collection"]
    bundle --> index["Repository indexes<br/>coverage and freshness"]
    claim --> resolution["Resolution record<br/>owner, blocker and closure"]
```

Never lift a favorable status upward without checking its denominator. A
matched observation does not make its claim matched; a matched claim does not
make every observation in the bundle matched; a repository readiness result
does not make every scientific method validated.

## Verdict Vocabulary

| Verdict | Required meaning | Honest downstream language |
| --- | --- | --- |
| `matched` | all claim-required comparable observations satisfy the registered exact or structural rule | the bounded claim matches under the recorded contract |
| `matched_with_tolerance` | all claim-required comparable observations satisfy the predeclared numerical budget | the bounded claim agrees within the named tolerance |
| `mismatch_explained` | a comparable difference exists and its registered cause and consequence are established | the claim differs for the recorded reason; state its impact |
| `mismatch_unexplained` | a comparable difference violates the rule without an accepted explanation | the claim does not currently correspond |
| `not_comparable` | the required estimand, identity, execution, or output alignment is unavailable | no valid comparison is currently owned |

`not_comparable` is not a polite match and `mismatch_unexplained` is not merely
“needs review.” Neither may be converted into `matched_with_tolerance` by
widening the tolerance after observing the result.

## Derive A Claim Verdict From Its Required Population

Start from the claim dependency set, not from whichever rows are easiest to
compare:

```mermaid
flowchart TD
    required["Enumerate required observations"] --> aligned{"Every required observation<br/>has aligned identities and outputs?"}
    aligned -- no --> incomparable["not_comparable"]
    aligned -- yes --> rules{"Apply every predeclared rule"}
    rules -- all exact --> exact["matched"]
    rules -- all within budgets --> tolerance["matched_with_tolerance"]
    rules -- violation with accepted cause --> explained["mismatch_explained"]
    rules -- violation without accepted cause --> unexplained["mismatch_unexplained"]
```

The denominator includes every observation the claim declares, including
missing, excluded, failed, and non-comparable members. A verdict may summarize
that population only when the aggregation rule is defined in advance and the
machine record retains each member. “Most rows matched” is a descriptive count,
not one of the governed claim verdicts.

When exact and tolerance rules coexist, report the strongest vocabulary that
is true for the complete claim population. Do not label the claim `matched`
because its exact subset passed while required tolerance-based observations
were merely within budget.

## Apply Verdicts To The Present Studies

The same repository can legitimately retain several evidence states because
each state belongs to a different claim and observation population.

| Present record | Correct interpretation | Incorrect compression |
| --- | --- | --- |
| PCM1 focused preparation bundle marked `matched` | the named preparation invariant corresponds for its registered inputs | every PCM1 output is exact or every figure is equivalent |
| PCM1 aggregate bundle marked `matched_with_tolerance` | all claim-required comparable observations meet their predeclared budgets | numerical identity or proof on unrelated primate data |
| PCM2 focused matches beside 12 broader `mismatch_unexplained` scalar rows | bounded bundle claims and unresolved broader observations coexist | PCM2 “passes overall” because most bundles are favorable |
| PCM2 intercept sweep marked `not_comparable` | no governed aligned observation currently supports that comparison | the two implementations were compared and disagreed |
| PCM3–PCM5 contract-only bundles marked `not_comparable` | the intended claim, inputs, source, and missing execution boundary are governed | wrappers, nearby methods, or complete manifests constitute numerical evidence |

These examples are not exceptions to the vocabulary. They demonstrate why a
verdict must travel with its claim identifier, population, method,
observations, and limitations.

## Keep Other Statuses In Their Own Records

Execution states such as blocked, failed, partial, and complete belong to run
records. Freshness states belong to dependency review. Coverage gaps describe
missing proof. Plot-only, workflow-only, artifact-only, and seeded-input-only
classify tracked surfaces that are not analytical correspondence. Preserve
these dimensions beside the verdict instead of forcing them into it.

## Resolve By Earliest Owning Cause

| Problem | Owning layer | Closure evidence |
| --- | --- | --- |
| source or dataset unavailable | provenance and dataset registry | governed locator, checksum, license and admitted population |
| reference program never executed | reference execution | versioned invocation, terminal record and primary outputs |
| Bijux program never executed | Bijux execution | runtime identity, invocation, diagnostics and primary outputs |
| estimands or identities differ | comparison design | aligned parameterization and stable observation keys |
| registered tolerance fails | observation | corrected implementation or accepted causal explanation plus consequence |
| claim lacks required observations | claim dependency map | complete eligible observation denominator |
| material dependency changed | freshness record | regenerated bundle and current indexes |

Resolve the earliest broken edge. Editing a report cannot repair a missing
execution; rerunning a program cannot repair an ill-defined estimand; marking
an item closed cannot repair a failed tolerance.

## Preserve Contradictions

When bundle, claim, observation, report, and generated index disagree, retain
all conflicting records and stop the affected claim. Identify the owning
source, rebuild projections, and review whether the contradiction changed a
published interpretation. Do not choose the most favorable status or delete
the row that prevents a clean summary.

## Record Resolution Without Rewriting History

A resolution record names the affected claim and observations, current
verdict, blocker or causal explanation, responsible owner, required execution
or evidence, acceptance rule, and closure artifacts. When work closes, create
new governed outputs and freshness state. Historical mismatches and
non-comparable states remain part of the audit trail.

A corrected implementation may produce a matching new comparison, but it does
not mean the earlier mismatch never occurred. Report both when the history
matters to released artifacts or scientific conclusions.

Closure requires a new observation chain, not only a resolution-state edit.
The closure record should link the corrected input, execution, result,
comparison, claim verdict, and freshness projection. If an accepted causal
explanation leaves the numerical difference in place, use
`mismatch_explained`; explanation makes the discrepancy interpretable, not
equal.

## Write The Public Statement

Use this pattern:

> Claim `<claim_id>` is `<verdict>` for `<population and method>` under
> `<comparison rule>`. The record excludes `<limits>`. Freshness is `<state>`
> at `<repository revision>`. Resolution state is `<state and owner>`.

Omit no material qualifier. If a statement cannot name its claim identifier,
population, method, comparison rule, and current verdict, it is broader than
the Evidence Book record supports.
