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

When the matrix is already aligned and you need to validate the partition file
itself, run `alignment partition-summary` first. That command reports assigned
and unassigned sites, mixed declared datatypes, and one row per locus, and it
can write the review table directly as TSV.

## Partitioned Multi-Locus Inference

Use the partition summary command before sending a concatenated matrix into the
adapter inference surface, then pass the same partition file into the adapter
step that needs it.

```bash
bijux-phylogenetics alignment partition-summary multilocus.aln.fasta \
  multilocus.partitions \
  --out artifacts/multilocus.partition-summary.tsv \
  --json
bijux-phylogenetics adapter model-select multilocus.aln.fasta \
  --partitions multilocus.partitions \
  --out-dir artifacts/multilocus-model \
  --prefix multilocus \
  --json
bijux-phylogenetics adapter infer-ml multilocus.aln.fasta \
  --partitions multilocus.partitions \
  --out-dir artifacts/multilocus-ml \
  --model GTR+G \
  --prefix multilocus \
  --json
```

If every partition is DNA or every partition is protein, the adapter passes a
normalized partition scheme to IQ-TREE on the original aligned matrix. If the
partition file mixes DNA and protein loci, the adapter writes one extracted
alignment per partition and a generated NEXUS scheme before invoking IQ-TREE.
For that mixed-datatype path, do not force one fixed single model across every
partition. Use a model-selection keyword such as `MF`, `MFP`, `TEST`, or
`TESTMERGE` so the engine can choose partition-appropriate models honestly.

## Coding DNA Alignment

Use `adapter align --codon-aware` when you need a nucleotide alignment that
preserves codon triplets for downstream phylogenetic inference.

```bash
bijux-phylogenetics adapter align coding-cds.fasta \
  --out artifacts/coding.aligned.fasta \
  --mode linsi \
  --codon-aware \
  --json
bijux-phylogenetics adapter model-select artifacts/coding.aligned.fasta \
  --out-dir artifacts/coding-model \
  --prefix coding \
  --sequence-type dna \
  --json
bijux-phylogenetics adapter infer-ml artifacts/coding.aligned.fasta \
  --out-dir artifacts/coding-ml \
  --model GTR+G \
  --prefix coding \
  --sequence-type dna \
  --json
```

The codon-aware alignment workflow excludes frame-broken sequences and
sequences with premature stop codons before MAFFT runs. It then aligns a
translated amino-acid guide and projects guide gaps back as nucleotide triplets,
so the final alignment length stays divisible by three and codon boundaries are
retained. The workflow also writes the translated guide input, the aligned
guide, and a TSV ledger of excluded sequences for review.

## Raw FASTA To Tree

Use `adapter fasta-to-tree` when you need one governed command from unaligned
FASTA to a reviewable inference bundle.

```bash
bijux-phylogenetics alignment sequence-type raw-sequences.fasta --json
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

Use `alignment sequence-type` when you need the raw FASTA type decision before
any engine is invoked. It reports compatible types, the selected default, the
confidence level, and mixed or invalid blocking signals.

Use `alignment validate-input` when you need one broader report of duplicate
identifiers, illegal sequence characters, empty records, raw-sequence length
outliers, and sequence-type detection before any engine is invoked.

Use `alignment repair-input` or the matching `adapter fasta-to-tree`
`--normalize-identifiers` and `--remove-invalid-records` controls when you want
the runtime to prepare a repaired FASTA explicitly rather than proceeding on
silent assumptions. Without those repair flags, `adapter fasta-to-tree` now
fails fast on duplicate identifiers, empty sequences, or illegal characters.
Mixed raw inputs also fail fast unless you declare a compatible
`--sequence-type` explicitly and remove incompatible records.

The checked repository examples that currently exercise this workflow are:

- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_sequences_raw.fasta`
- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_alignment.fasta`
- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_alignment_protein.fasta`
- `packages/bijux-phylogenetics/tests/fixtures/alignments/example_sequences_invalid_input.fasta`

Those checks cover raw DNA input, already aligned DNA input, and protein input,
and they verify that the workflow emits the aligned matrix, trimmed matrix,
tree, log, model table, and support table with stable names.

For coding DNA, prefer the codon-aware alignment workflow above and then run
the downstream adapter steps explicitly. The current `adapter fasta-to-tree`
path is still a generic trimAl-based pipeline rather than the codon-preserving
entrypoint.
