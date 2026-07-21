---
title: Cross-Surface Contract Maintenance
audience: maintainer
type: handbook
status: canonical
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-21
---

# Cross-Surface Contract Maintenance

One runtime capability can appear as a Python import, CLI command, API schema,
artifact field, package description, workflow example, benchmark subject, and
evidence dependency. These are projections of one contract and must change
together when semantics change.

```mermaid
flowchart TD
    contract["Owned scientific contract"]
    python["Python API"]
    cli["CLI"]
    schema["OpenAPI and artifacts"]
    packages["Distribution metadata"]
    guides["Usage and method guides"]
    proof["Tests, benchmarks,<br/>evidence"]

    contract --> python
    contract --> cli
    contract --> schema
    contract --> packages
    contract --> guides
    contract --> proof
```

## Canonical Sources And Projections

| Concern | Canonical source | Projections to inspect |
| --- | --- | --- |
| Python behavior | curated runtime implementation and tests | API guide, examples, serialization |
| command behavior | command registration and tests | generated CLI reference, operations guides |
| service contract | pinned OpenAPI schema | API guide, clients, schema hash |
| artifact meaning | result and manifest contracts | reports, consumer guide, evidence bundles |
| package identity | package metadata | root/package READMEs, badges, release records |
| evidence verdict | bundle claims and checks | indexes, study guides, release claims |

Never edit a projection to conceal disagreement with its canonical source.

## Change Impact Review

For each material change, record:

1. the owned contract and compatibility class;
2. affected Python, CLI, schema, and artifact surfaces;
3. affected package descriptions and examples;
4. benchmarks whose configuration or interpretation changed;
5. evidence bundles made stale by code, input, or configuration changes;
6. verification commands and observed outcomes.

## Contract Change Classes

| Change | Required review |
| --- | --- |
| additive public surface | discovery route, reference generation, package exposure, and bounded quality statement |
| changed default or interpretation | examples, method guidance, result semantics, compatibility, and evidence freshness |
| removed or renamed surface | deprecation path, alias behavior, migration guidance, generated references, and built-artifact smoke |
| schema or artifact field change | versioning decision, producer/consumer compatibility, migration, and evidence consumers |
| numerical implementation change | known-answer and invariant checks, benchmark comparison, affected bundle freshness |
| ownership change | native/external attribution across code, docs, metadata, diagnostics, and reports |

A textual correction can be contract-relevant when it changes what users are
promised. Classify by semantic effect, not by file type or diff size.

## Completeness Audit

For the changed subject, search public imports, command registration, OpenAPI
and artifact schemas, package metadata, README files, handbook guides,
benchmarks, and evidence records. Record each occurrence as one of:

- canonical owner;
- generated projection;
- maintained explanatory projection;
- historical or obsolete reference requiring removal;
- unaffected use with a documented reason.

This inventory makes an omission reviewable. A list of files edited does not
prove that every projection was found.

## Close The Projection Dependency Record

For each public statement, retain a dependency record from the canonical owner
to every material projection and proof. Each edge names the consumer, meaning
carried across the edge, source revision or digest, freshness rule, and review
state.

```mermaid
flowchart LR
    owner["Canonical contract<br/>and revision"]
    interface["Installed interface"]
    artifact["Artifact and schema"]
    guide["Reader guidance"]
    proof["Test · benchmark<br/>evidence"]
    release["Release decision"]

    owner --> interface --> artifact
    owner --> guide
    owner --> proof
    artifact --> proof
    guide --> release
    proof --> release
```

Classify every expected edge as current, unaffected with a reason, stale,
blocked, or not executed. Missing edges are not implicitly unaffected. A
projection is current only when its material inputs resolve to the reviewed
owner; matching wording or nearby timestamps do not establish freshness.

Close the record from both directions: every material consumer must name its
owner, and every changed owner must account for its consumers. This catches
orphaned schemas, examples, summaries, and evidence dependencies that a
source-only search cannot discover.

## Drift Patterns

- A new native method exists but only an internal module path exposes it.
- The canonical package shows an example the alias cannot execute.
- CLI help changes while the generated command reference stays unchanged.
- A result field changes meaning without an artifact schema change.
- A study guide summarizes a verdict no longer present in the indexes.
- Package metadata calls an adapter-backed workflow native.
- A report renders successfully while omitting warnings from its structured
  result.

## Confidence Vocabulary

Use `implemented`, `tested`, `benchmarked`, `evidence-backed`, and `released`
only with the scope defined by the corresponding contract. Avoid umbrella
labels such as “validated” when the actual proof is one unit test or one study
bundle.

## Status Ownership

Status values are part of their record contracts, not interchangeable
sentiment labels.

| Record | Status namespace | Forbidden inference |
| --- | --- | --- |
| live parity case | `passed`, `failed`, `skipped` | a pass is not an Evidence Book match |
| Evidence Book claim | `matched`, `matched_with_tolerance`, `mismatch_explained`, `mismatch_unexplained`, `not_comparable` | `not_comparable` is not a failed runtime execution |
| resolution record | open or blocked follow-through | open work is not an “open” verdict |
| release readiness | ready or blocked | ready does not imply scientific equivalence |

When a status crosses a surface, carry the source record and translate its
meaning explicitly. Do not copy the word into a different namespace.

## Release Language Decision

Derive publication language from the narrowest current record:

1. identify the exact capability or scientific claim;
2. confirm that its computational owner and installed invocation are current;
3. trace reported values to structured output and retained inputs;
4. read the claim-level verdict and any separate resolution work;
5. state only the supported scope, denominator, tolerance, and limitation;
6. block publication if any required projection contradicts that record.

A declared `not_comparable` verdict may be published as a bounded limitation
when readiness policy permits it. It must not be rewritten as a match, a
failure, or unfinished execution.

## Publication Gate

A change is ready for publication when canonical sources and projections
agree, compatibility effects are explicit, generated outputs have been
regenerated rather than hand-edited, evidence freshness has been evaluated,
and the repository's package and documentation checks pass from a clean
environment. Open or blocked resolution work may coexist with a release only
when its governing verdict and limitation remain explicit. An unacknowledged
mismatch between runtime behavior and its published contracts cannot.

Block publication when the canonical owner is unclear, compatibility is
undeclared, a generated projection is stale, a public example cannot execute
through the installed surface, warnings or missing cases vanish across a
boundary, or an evidence claim no longer resolves to a current governed
record.
