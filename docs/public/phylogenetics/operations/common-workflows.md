---
title: Common Workflows
audience: public
type: how-to
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-11
---

# Common Workflows

Typical public workflows include:

- validate and inspect a tree before downstream use
- trim and inspect an alignment before reporting
- run one command from raw FASTA to a supported inference bundle
- audit a concatenated multi-locus matrix before inference
- run a comparative model and capture JSON plus report artifacts
- generate a review-ready figure package for a tree

The public workflow contract is that important outputs should be inspectable
after the command finishes.

For concatenated phylogenomics inputs, run `alignment occupancy` against an
aligned FASTA plus partition file before tree inference. The command can emit
per-taxon and per-locus TSV tables, a taxon-by-locus occupancy matrix, and a
retained alignment plus remapped partition file after applying explicit
coverage thresholds.

## Raw FASTA To Tree

Use `adapter fasta-to-tree` when you need one governed command from unaligned
FASTA to a reviewable inference bundle.

```bash
bijux-phylogenetics alignment validate-input raw-sequences.fasta --json
bijux-phylogenetics alignment repair-input raw-sequences.fasta \
  --out artifacts/raw-sequences.repaired.fasta \
  --normalize-identifiers \
  --remove-invalid-records \
  --json
bijux-phylogenetics adapter fasta-to-tree raw-sequences.fasta \
  --out-dir artifacts/fasta-to-tree \
  --prefix mammals \
  --normalize-identifiers \
  --remove-invalid-records \
  --bootstrap-replicates 1000 \
  --json
```

The workflow accepts DNA and protein FASTA inputs. It writes these durable
user-facing outputs:

- `mammals.aln`
- `mammals.trimmed.aln`
- `mammals.tree`
- `mammals.log`
- `mammals.model.tsv`
- `mammals.support.tsv`
- `mammals.manifest.json`

It also retains step-specific engine artifacts under
`artifacts/fasta-to-tree/engine-artifacts/mammals/` for auditability.

Use `alignment validate-input` when you need an explicit report of duplicate
identifiers, illegal sequence characters, empty records, and raw-sequence
length outliers before any engine is invoked.

Use `alignment repair-input` or the matching `adapter fasta-to-tree`
`--normalize-identifiers` and `--remove-invalid-records` controls when you want
the runtime to prepare a repaired FASTA explicitly rather than proceeding on
silent assumptions. Without those repair flags, `adapter fasta-to-tree` now
fails fast on duplicate identifiers, empty sequences, or illegal characters.

The checked repository examples that currently exercise this workflow are:

- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_sequences_raw.fasta`
- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_alignment.fasta`
- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_alignment_protein.fasta`
- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_sequences_invalid_input.fasta`

Those checks cover raw DNA input, already aligned DNA input, and protein input,
and they verify that the workflow emits the aligned matrix, trimmed matrix,
tree, log, model table, and support table with stable names.
