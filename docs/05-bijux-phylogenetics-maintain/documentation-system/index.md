---
title: Documentation System
audience: maintainer
type: handbook
status: active
owner: bijux-phylogenetics-dev
last_reviewed: 2026-07-21
---

# Documentation System

Every maintained statement has one canonical owner and may have several
projections. Runtime behavior is owned by source and executable contracts;
scientific verdicts are owned by Evidence Book records; package identity is
owned by distribution metadata; generated references are owned by their
generator inputs. Handwritten pages explain those authorities without
replacing them.

```mermaid
flowchart TD
    change["Behavior, schema, evidence,<br/>or package change"]
    owner{"Canonical owner"}
    generated["Generated references<br/>and indexes"]
    handbook["Reader-facing explanation<br/>and examples"]
    contracts["Semantic and build<br/>validation"]

    change --> owner
    owner --> generated
    owner --> handbook
    generated --> contracts
    handbook --> contracts
```

## Resolve Ownership Before Editing

| Changed fact | Canonical owner | Maintained projections |
| --- | --- | --- |
| Python signature or typed result | public source surface and contract tests | API guide, examples, generated references |
| CLI command, option, or exit behavior | command registry and CLI tests | CLI reference, workflow guide, README example |
| artifact field or manifest relationship | schema/producer and consumer tests | artifact contract and consumption guide |
| numerical ownership or scientific boundary | implementation, numerical tests, parity/evidence records | domain, workflow, quality, and parity pages |
| evidence verdict or freshness | bundle manifests and generated indexes | Evidence Book study guides and landing ledger |
| distribution name, requirement, or URL | package metadata | root and package READMEs, installation guide |

Edit the owner first. Regenerate owned derivatives with repository commands.
Then revise handwritten explanations whose reader-visible meaning changed.
Never hand-edit a generated page to mask disagreement with its input.

## Separate Reader And Maintainer Language

Reader-facing pages explain the scientific decision: which method applies,
what input contract it accepts, what result it returns, how failure appears,
and what evidence bounds the claim. They should not discuss the writing plan,
page quality, desired impression, or internal delivery sequence.

Maintainer pages explain ownership and propagation: which source establishes a
fact, which projections repeat it, what generator owns derived content, and
which checks block publication. A maintenance instruction belongs here only
when it changes how repository facts are kept coherent.

| Content | Reader-facing location | Maintainer location |
| --- | --- | --- |
| choosing ML, Bayesian, parsimony, simulation, or an adapter | domain and workflow guide | affected-surface inventory and validation route |
| interpreting a result, warning, refusal, or limitation | interface, operations, or quality guide | producer/consumer ownership and compatibility obligations |
| understanding a study verdict | Evidence Book study guide | freshness, regeneration, and publication gate |
| changing an API, command, schema, or artifact | migration guidance when consumers are affected | canonical source, projection set, checks, and release decision |

## Propagate A Contract Change

```mermaid
flowchart LR
    classify["Classify changed meaning"]
    owner["Change canonical owner<br/>and contract tests"]
    generate["Regenerate owned<br/>projections"]
    explain["Revise affected<br/>reader decisions"]
    evidence["Evaluate evidence<br/>and benchmark freshness"]
    verify["Build and validate<br/>the publication set"]

    classify --> owner --> generate --> explain --> evidence --> verify
```

1. State the changed behavior or meaning without referring to filenames.
2. Locate the implementation, schema, registry, package metadata, or evidence
   record that owns it.
3. Identify every generated and handwritten projection by consumer decision.
4. Update the owner and its contract tests before changing explanatory prose.
5. Regenerate repository-owned derivatives and inspect their diff separately.
6. Re-evaluate benchmarks and Evidence Book freshness when scientific behavior,
   inputs, defaults, or normalization changed.
7. Run the narrow contract tests, documentation truth contract, and strict site
   build; then inspect the complete publication diff for contradictory claims.

## Documentation Review Packet

A reviewable documentation change records:

- the changed contract and its canonical owner;
- the affected reader decisions and maintained projections;
- generated outputs, with the generator command and source identity;
- examples executed through the installed or explicitly named development
  surface;
- evidence and benchmark records evaluated for freshness;
- checks run, observed results, and any publication blocker;
- deliberate exclusions with the reason they are semantically unaffected.

The packet is an accountability record, not a claim of completeness by volume.
A short correction can require broad review when it changes public meaning; a
large generated diff can require no prose change when semantics are unchanged.

## Maintained Guides

- [Documentation ownership](documentation-ownership.md) explains which
  source owns a fact, how canonicality is resolved, and which duplication is
  intentional.
- [Native surface documentation](native-surface-documentation.md) explains how
  an owned scientific implementation is connected to API, workflow,
  diagnostics, benchmark, and evidence boundaries.

## Completion Criteria

A documentation change is complete when the canonical owner and every affected
projection agree, generated assets are current, examples use installed public
surfaces, links resolve in a strict site build, scientific strength is bounded
by executable evidence, and no maintainer-only delivery commentary appears in
public scientific guidance.
