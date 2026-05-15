---
title: API Surface
audience: reader
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-11
---

# API Surface

The stable public API of this repository is deliberately small. Most users
should treat the CLI, Python import surface, and checked-in evidence files as
the main public contract. The HTTP-facing surface is frozen as a checked-in
OpenAPI bundle, not as a live server promise.

## The Main Public Surfaces

- the `bijux-phylogenetics` and `phylogenetic` console scripts
- stable Python imports under `bijux_phylogenetics`
- checked-in evidence studies under `evidence-book/studies/`
- the frozen public API contract under `apis/bijux-phylogenetics/v1/`

## Python Workflow Entry Points

For notebook, script, and pipeline use, the stable workflow-oriented Python
surface lives under `bijux_phylogenetics.api`.

The named entry points are:

- `run_fasta_validation_workflow`
- `run_alignment_workflow`
- `run_tree_inference_workflow`
- `run_sequence_to_tree_workflow`
- `run_tree_comparison_workflow`
- `run_comparative_model_workflow`
- `run_ancestral_reconstruction_workflow`
- `render_report_workflow`
- `run_configured_phylo_workflow`

These functions do not invent a second result contract. They return the same
typed runtime report objects already used by the CLI surfaces:

- `FastaInputValidationReport`
- `EngineWorkflowReport`
- `FastaToTreeWorkflowReport`
- `TreeComparisonReport`
- `PGLSResult`
- `DiscreteAncestralReport`
- `ReportBuildResult`
- `WorkflowConfigRunReport`

For most users, that is enough. The repository prefers a reviewable frozen
contract over pretending that every internal helper is part of the promise.

## Frozen HTTP Contract

The checked-in OpenAPI bundle describes a future-compatible read surface for:

- repository health and release identity
- runtime package ownership summaries
- evidence-study catalog access
- evidence-study detail summaries for the governed PCM1 and PCM2 studies

The freeze bundle lives here:

- `apis/bijux-phylogenetics/v1/schema.yaml`
- `apis/bijux-phylogenetics/v1/pinned_openapi.json`
- `apis/bijux-phylogenetics/v1/schema.hash`

## Practical Rule

If you need an integration contract, start from the CLI, the named Python
surface, or the frozen OpenAPI bundle. If an internal module or helper is not
named in those surfaces, assume it can move.
