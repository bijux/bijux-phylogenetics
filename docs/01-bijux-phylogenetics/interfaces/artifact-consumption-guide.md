---
title: Artifact Consumption Guide
audience: public
type: how-to
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Artifact Consumption Guide

Use this guide when the repository has already produced files and the next
question is how to consume them correctly.

## Start By Artifact Kind

- I need one reviewer-facing report or figure
  Start with the HTML or SVG artifact and its paired manifest or JSON summary.
- I need one machine-readable workflow result
  Start with JSON and the matching manifest or stable TSV ledgers.
- I need one downstream table for scripting or analysis
  Start with the governed TSV families.
- I need to know whether one file name or schema is stable enough to depend on
  Start with [artifact contracts](artifact-contracts.md).

## Practical Reading Rules

- do not treat one HTML report alone as the whole contract when a manifest or
  JSON companion exists
- do not scrape terminal text when the workflow already emits JSON or TSV
- prefer manifest-backed bundles when one report links multiple outputs
- prefer the documented stable families instead of inferring stability from one
  incidental example file

## Common Stable Families

The public artifact families most worth depending on are:

- workflow and report manifests such as `prefix.manifest.json`
- model-selection tables such as `prefix.model.tsv`
- support ledgers such as `prefix.support.tsv`
- clade, branch, and event ledgers
- comparative trait and comparative summary ledgers
- benchmark and parity JSON or TSV summaries where the surface documents them

## When To Use Bundle Outputs

Prefer bundle outputs when:

- one result includes multiple linked files
- one review needs both human-readable and machine-readable outputs
- one publication or audit trail needs stable provenance inside the same output
  surface

That is the normal case for serious report and evidence-oriented workflows in
this repository.

## When To Escalate To Quality Or Evidence Docs

Artifact shape is not the same thing as claim trust.

If the next question becomes:

- can I trust this benchmark difference
- does this release look publishable
- does this study close a broader claim

then stop reading interface docs alone and move to the quality or evidence
surfaces.
