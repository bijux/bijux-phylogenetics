---
title: Validation Lanes
audience: mixed
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Validation Lanes

Use the narrowest truthful lane name when describing confidence. Lanes are
orthogonal proof families; passing one never silently grants another.

| Lane | Core question | Typical proof | Does not prove |
| --- | --- | --- | --- |
| structural | Is the input or artifact well formed? | parsers, schemas, manifests | scientific appropriateness |
| contract | Does the public interface obey its documented behavior? | unit and API tests | numerical reference agreement |
| numerical | Does computation satisfy invariants or tolerances? | known answers, properties, reference values | workflow completeness |
| integration | Do components compose and preserve state? | workflow and artifact tests | external scientific parity |
| benchmark | How does behavior compare on a named corpus? | metrics with environment and thresholds | universal superiority |
| evidence | Do governed records support named study claims? | provenance, checks, verdicts | claims outside the study |
| release | Are packages, docs, schemas, and evidence aligned? | build and publication gates | future or untested behavior |

## Lane Relationships

```mermaid
flowchart LR
    input["Structural"] --> contract["Contract"]
    contract --> integration["Integration"]
    numerical["Numerical"] --> integration
    integration --> benchmark["Benchmark"]
    benchmark --> evidence["Evidence"]
    integration --> release["Release"]
    evidence --> release
```

This is a common dependency shape, not a claim that all capabilities must reach
every node. Release checks may legitimately publish a well-tested feature with
no study evidence if documentation says exactly that.

## Select A Lane From The Claim

| Proposed statement | Required lane | Missing proof if only the previous lane passes |
| --- | --- | --- |
| “the file is interpretable” | structural | public behavior and scientific values remain unchecked |
| “the documented command behaves as specified” | contract | numerical correctness remains unchecked |
| “the method returns the declared value within tolerance” | numerical | composed workflow and artifact behavior remain unchecked |
| “the workflow preserves inputs, diagnostics, and outputs” | integration | reference agreement and scale behavior remain unchecked |
| “behavior meets a target on this corpus” | benchmark | transfer beyond the corpus remains unchecked |
| “this study claim has the recorded verdict” | evidence | claims outside the registered scope remain unchecked |
| “this revision is publishable under repository policy” | release | future revisions and undeclared scientific scope remain unchecked |

Start from the statement being made, not the check that is easiest to run.

## Structural Lane

Checks formats, required fields, taxon uniqueness, manifest shape, and schema
compatibility. Structural validity establishes that later stages can interpret
the input, not that its biological content is correct.

## Contract Lane

Checks types, public imports, documented errors, CLI exit behavior, and output
shape. Contract tests defend user expectations while leaving numerical and
scientific questions to their own lanes.

## Numerical Lane

Checks invariants, known-answer fixtures, tolerances, and property behavior.
Tolerance must reflect scale and algorithmic expectations; unexplained loose
tolerances are not honest robustness.

## Integration Lane

Checks end-to-end composition, adapter capture, result normalization, and
artifact completeness. Integration success can coexist with a diagnostic that
makes a scientific result unusable; preserve both outcomes.

## Benchmark Lane

Checks defined metrics on named corpora under recorded environments. Always
report failures and exclusions. See [native benchmark review](../operations/native-benchmark-review.md).

For every aggregate, retain selected, executed, successful, failed, skipped,
and missing counts. A benchmark that reports only completed cases cannot
support a corpus-wide claim.

## Evidence Lane

Checks claim IDs against sources, inputs, executions, and verdict rules in the
Evidence Book. Open follow-through and named blockers are valid resolution
states and must not be rewritten as matched because adjacent claims pass.

Evidence Book verdicts use their governed comparison vocabulary:
`matched`, `matched_with_tolerance`, `mismatch_explained`,
`mismatch_unexplained`, and `not_comparable`. Workflow resolution such as open
follow-through is recorded separately; do not substitute one namespace for the
other.

## Release Lane

Checks that distributions build, public metadata agrees, docs render, schemas
validate, and declared evidence is fresh enough for publication. It protects
the declared release surface; it cannot broaden the underlying proofs.

## Claim Template

State: **subject + lane + scope + acceptance rule + verdict + revision**.
For example: “The native inference result passed topology agreement tolerance
T on corpus C under configuration K at revision R.” This is more useful than
“native inference is validated.”
