---
title: Artifact Contracts
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-03
---

# Artifact Contracts

The runtime emits review artifacts, not just terminal text. Those artifacts are
part of the public surface when the documentation names them as governed output
families.

## Main Public Artifact Families

Common governed families include:

- SVG tree figures
- HTML reports
- JSON summaries
- bundle directories with manifests and supporting files
- reviewer-facing TSV ledgers

## Stable Families Worth Depending On

The most important schema-governed public families include:

- FASTA-to-tree model tables such as `prefix.model.tsv`
- FASTA-to-tree support tables such as `prefix.support.tsv`
- clade review tables such as `clade-table.tsv`
- branch review tables such as `host-switch-branches.tsv`
- derived comparative trait tables such as `comparative-traits.tsv`
- comparative model-output tables such as `comparative-summary.tsv`
- event ledgers such as `biogeography/event-table.tsv`
- workflow and report manifests such as `prefix.manifest.json` and
  `comparative-report.manifest.json`

When the repository documents one family at that level, downstream code should
anchor on the stable headers or top-level keys instead of scraping ad hoc
examples.

## Contract Reading Rule

- use HTML and SVG outputs for human review
- use JSON and TSV outputs for downstream automation
- use manifests when one workflow emits multiple linked files
- use benchmark and parity summaries through their documented output surfaces
  rather than inferring structure from one test fixture

## What This Page Does Not Claim

- it does not claim that every output file ever written by the repository is a
  stable public contract
- it does not claim that every artifact is equally important for trust or
  publication review
- it does not turn one output schema into proof that the underlying science is
  already closed

Those outputs matter because they let downstream readers inspect the result
instead of trusting a one-line success message.

## Matching Practical Guide

For day-to-day output handling, start with the
[artifact consumption guide](artifact-consumption-guide.md).
