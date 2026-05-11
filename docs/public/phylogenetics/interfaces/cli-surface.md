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
