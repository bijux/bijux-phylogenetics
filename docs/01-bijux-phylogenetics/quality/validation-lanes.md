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

## Evidence Lane

Checks claim IDs against sources, inputs, executions, and verdict rules in the
Evidence Book. Open and blocked verdicts are valid records and must not be
rewritten as passing because adjacent claims pass.

## Release Lane

Checks that distributions build, public metadata agrees, docs render, schemas
validate, and declared evidence is fresh enough for publication. It protects
the release narrative; it cannot broaden the underlying proofs.

## Claim Template

State: **subject + lane + scope + acceptance rule + verdict + revision**.
For example: “The native inference result passed topology agreement tolerance
T on corpus C under configuration K at revision R.” This is more useful than
“native inference is validated.”
