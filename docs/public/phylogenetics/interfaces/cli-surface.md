---
title: CLI Surface
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-11
---

# CLI Surface

The CLI is the primary operational surface for most users.

## Major Command Families

- `validate`, `inspect`, `compare`, `annotate`, `render`
- `alignment ...`
- `comparative ...`
- `ancestral ...`
- `adapter ...`
- `bundle` and `report`

The public rule is simple: commands should produce explicit, reviewable outputs
and should not hide important assumptions behind silent defaults.

The alignment family includes matrix-audit commands such as `alignment occupancy`
for concatenated multi-locus FASTA plus partition inputs. That workflow reports
per-taxon coverage, per-locus coverage, low-coverage flags, TSV tables, and an
optionally filtered retained matrix with remapped partitions.

The adapter family also includes `adapter fasta-to-tree`, which is the
supported end-to-end inference entrypoint for raw FASTA inputs. It emits a
reviewable aligned matrix, trimmed matrix, selected-model table, supported
tree, support summary table, run log, and manifest in one command instead of
forcing users to stitch separate adapter steps together by hand.

For raw input hygiene before alignment, the alignment family now includes
`alignment sequence-type`, `alignment validate-input`, and
`alignment repair-input`. Those commands expose the same raw sequence-type,
duplicate-ID, illegal-character, empty-record, length-outlier, and
identifier-normalization contract that `adapter fasta-to-tree` uses internally,
including the rule that mixed raw inputs must be forced with an explicit
`--sequence-type` before the workflow can continue.
