---
title: Trust Surfaces
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-21
---

# Trust Surfaces

Trust is a traceability property. A public statement is supportable only when
its computational owner, installed interface, structured output, validation
record, evidence verdict, and release language agree. No surface may silently
borrow authority from another.

```mermaid
flowchart LR
    owner["Computational owner"] --> interface["Installed interface"]
    interface --> result["Structured result"]
    result --> validation["Validation record"]
    validation --> claim["Governed claim"]
    claim --> language["Permitted release language"]

    limits["Failures · exclusions<br/>missing outputs"] --> validation
    limits --> claim
    limits --> language
```

## Questions Owned By Each Surface

| Surface | Question it answers | Record that carries authority |
| --- | --- | --- |
| computation | Who performed the numerical work? | native implementation or named external engine |
| interface | How can an installed user invoke it? | curated Python export, CLI registration, or service schema |
| result | What was produced, omitted, or refused? | typed result, serialized payload, manifest, and warnings |
| validation | What behavior was exercised under which denominator? | test, benchmark, or live parity record |
| evidence | What claim follows from governed observations? | claim-scoped Evidence Book bundle and verdict |
| publication | What may a release say now? | readiness decision reconciled with current limitations |

A passing test does not create a parity verdict. A live parity case does not
create an Evidence Book verdict. A ready release does not turn an explicit
`not_comparable` boundary into equivalence.

## Detect Cross-Surface Drift

| Drift | Detection question | Required correction |
| --- | --- | --- |
| ownership drift | Is an adapter or external posterior described as native? | restore the actual computational owner everywhere |
| interface drift | Do examples use undocumented imports, commands, or fields? | move examples to supported surfaces or change the contract |
| artifact drift | Can a report value be traced to structured output? | repair the producer/consumer contract before presentation |
| evidence drift | Does prose exceed the bundle verdict, scope, or freshness? | narrow the statement and link the governing record |
| status drift | Has a runtime state been restated as an evidence or release status? | restore the vocabulary owned by each record |
| package drift | Do canonical and alias distributions imply different capability? | align metadata, README, install smoke, and release records |
| release drift | Does publication proceed with failed or missing governed output? | stop the gate and preserve the blocking status |

## Trust Decision Packet

For a material scientific or public-interface change, retain enough context to
reconstruct the decision without relying on CI history or maintainer memory:

- repository revision and resolved package versions;
- computational owner and supported invocation surface;
- input, configuration, seed, model, and external-tool identity;
- selected-case denominator and produced/missing artifact inventory;
- exact validation commands and observed outcomes;
- affected claims, bundle revisions, verdicts, and resolution records;
- permitted user-facing statement and named limitations;
- release-readiness result and decision owner.

The packet is incomplete when a skipped optional lane disappears into an
aggregate pass, a report cannot be reconstructed from structured output, or a
claim lacks a current governed record.

## Maintained Guides

- [Cross-surface contract maintenance](public-trust-surface-maintenance.md)
  maps runtime, evidence, package, documentation, and release changes across
  their canonical owners and public projections.
- [Public interface and artifact maintenance](public-interface-and-artifact-maintenance.md)
  covers Python and CLI selection, structured results, manifests, schemas,
  and report consumers.
- [Evidence claim and study-guide maintenance](trust-lane-and-study-guide-maintenance.md)
  governs claim promotion, status ownership, evidence freshness, and the
  language exposed by PCM study guides.

## Release Review Questions

- Can every capability claim identify native or external ownership?
- Can every example execute through a supported installed surface?
- Can every reported value be traced to typed or machine-readable output?
- Can every evidence statement identify a current claim and verdict?
- Are failures, skips, `not_comparable` cases, and incomplete artifacts visible?
- Do both distribution names, site navigation, and release metadata agree?
- Could a reviewer reconstruct the decision from retained records alone?
