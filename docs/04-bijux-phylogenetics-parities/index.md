---
title: Parities And Correspondence
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-22
---

# Parities And Correspondence

Bijux Phylogenetics is compared with established phylogenetic software at
explicit, reviewable boundaries. The comparison may be a live call to an R
package, a checked reference fixture, a benchmark against a named
implementation, or an adapter that executes an external engine. Those are
different relationships and the documentation keeps them separate.

```mermaid
flowchart LR
    question["Named scientific<br/>operation"]
    bijux["Bijux result"]
    reference["Live tool or<br/>checked fixture"]
    normalize["Aligned inputs,<br/>units and conventions"]
    compare["Declared fields<br/>and tolerances"]
    outcome["Passed · failed<br/>or skipped case"]

    question --> bijux
    question --> reference
    bijux --> normalize
    reference --> normalize
    normalize --> compare --> outcome
```

## Reference Landscape

| Relationship | Tools | What the repository owns | What that record does not establish |
| --- | --- | --- | --- |
| live R parity | APE, GEIGER, phytools | governed cases, R runners, normalization, tolerances, mismatch artifacts | package-wide interchangeability or unselected operations |
| checked R-derived correspondence | `nlme`, selected course workflows and fixtures | frozen inputs/outputs, comparison rules, provenance and limits | a current live execution of the reference package |
| external alignment and tree execution | MAFFT, trimAl, IQ-TREE2, FastTree | adapters, configuration, execution capture, parsing and normalized artifacts | that Bijux owns or independently reproduces the external algorithm |
| external posterior execution and review | MrBayes, BEAST | preparation, execution adapters, trace/tree parsing, diagnostics and reports | posterior agreement without a separate aligned comparison |

An adapter is not automatically a parity proof. It establishes that the
external program can participate in a governed Bijux workflow. Numerical or
structural equivalence requires a separate comparison with aligned assumptions
and explicit observation rules.

## Choose The Relationship From The Question

| Question | Required relationship | Minimum retained record |
| --- | --- | --- |
| Does one registered operation agree with APE, GEIGER, or phytools in this environment? | live R parity case | tool versions, case identity, both outputs, normalization, observations, tolerance, status |
| Does Bijux reproduce a frozen reference result without R? | checked-fixture correspondence | fixture provenance, frozen output identity, Bijux output, comparison rule |
| Can an external executable participate in a governed workflow? | external-engine adapter execution | executable version, invocation, streams, native files, parser result, workflow status |
| Does one native result agree with an external engine? | explicit cross-engine comparison | aligned models, inputs, meanings, observations, tolerances, both denominators |
| Does repository evidence support a bounded study claim? | Evidence Book bundle | claim ID, source, inputs, execution, checks, observations, verdict, freshness |

One run can contribute to several records, but it does not merge their status
names. Select the relationship before choosing a command or interpreting a
green result.

## Correspondence Is Not Substitutability

Two surfaces can agree on selected observations while remaining unsuitable as
drop-in replacements. A substitution claim needs a broader contract than a
numerical parity claim.

| Dimension | Correspondence asks | Substitution additionally requires |
| --- | --- | --- |
| input | were the selected scientific objects aligned? | accepted formats, validation, defaults, missing-data policy, and refusal behavior align for intended users |
| method | were the compared estimands and parameterizations aligned? | the replacement exposes every required model, constraint, and configuration |
| result | did registered observations satisfy their rules? | result fields, ordering, units, uncertainty, omissions, and status semantics remain usable |
| operation | did the selected case execute? | failure, retry, resource, concurrency, and artifact behavior satisfy the consuming workflow |
| assurance | is the named case supported? | coverage and independent evidence match the full intended substitution scope |

The repository publishes operation-scoped correspondence. Unless a specific
consumer contract says otherwise, a passing case does not authorize replacing
APE, GEIGER, phytools, MAFFT, IQ-TREE2, BEAST, or another owner across an
analysis pipeline.

## Three Records, Three Meanings

| Record | Status owner | What it establishes |
| --- | --- | --- |
| parity harness case | live or checked-fixture runner: `passed`, `failed`, `skipped` | whether one registered case executed and satisfied its comparisons |
| external-engine bundle | adapter and workflow result state | whether a named external computation was prepared, executed, parsed, or refused |
| Evidence Book bundle | governed verdict vocabulary | whether identified scientific claims matched, mismatched, or remained non-comparable |

These records may be linked but are not interchangeable. A skipped R case has
no Evidence Book verdict. A successful adapter run is not a parity case. An
open follow-through item describes remaining resolution work; it does not
replace a bundle's recorded verdict.

## Reading Paths

- [R ecosystem](r-ecosystem/index.md) describes direct live parity and checked
  R-derived references.
- [External engines](external-engines/index.md) describes correspondence with
  executables that remain computational owners.
- [Interpreting correspondence](correspondence/index.md) explains case
  selection, normalization, outcomes, artifacts, and claim limits.
- [Verification matrix](correspondence/verification-matrix.md) distinguishes
  live parity, checked fixtures, external execution, and study provenance.
- [Artifacts and reproduction](correspondence/artifacts-and-reproduction.md)
  shows how a comparison remains inspectable after execution.

## The Public Contract

Every parity claim names the external tool and version, Bijux surface, fixture,
parameterization, compared fields, tolerance, execution mode, and case status.
Unavailable R packages and unsupported reference cases remain visible as
skips or out-of-scope boundaries rather than being counted as passes.

A complete claim can be read as:

> For **operation**, Bijux **version and surface** was compared with
> **reference version and surface** on **fixture or corpus**, after
> **normalization**, using **observations and rules**. Of **selected cases**,
> **passed**, **failed**, and **skipped** counts were retained at **artifact
> location**. The conclusion is limited to **named scope**.

Omitting the denominator, reference identity, or comparison rule turns a
reviewable correspondence statement into an unverifiable similarity claim.
