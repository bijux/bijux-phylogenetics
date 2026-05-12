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
- assemble aligned loci into a concatenated supermatrix
- audit a concatenated multi-locus matrix before inference
- run a comparative model and capture JSON plus report artifacts
- generate a review-ready figure package for a tree

The public workflow contract is that important outputs should be inspectable
after the command finishes.

When your starting point is one aligned FASTA per locus, run
`alignment concatenate` first. That workflow writes the concatenated alignment,
the remapped partition file, and the taxon-by-locus occupancy matrix in one
step while preserving taxon identifiers and inserting `?` blocks for absent
taxa.

```bash
bijux-phylogenetics alignment concatenate loci/alpha-dna.fasta \
  loci/beta-protein.fasta \
  loci/gamma-dna.fasta \
  --data-type DNA \
  --data-type PROTEIN \
  --data-type DNA \
  --out artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions-out artifacts/mixed-locus-supermatrix.partitions.txt \
  --matrix-out artifacts/mixed-locus-supermatrix.matrix.tsv \
  --json
```

Use repeated `--data-type` flags when short or ambiguity-rich loci would make a
DNA locus and a protein locus look identical by characters alone. The
concatenation workflow records those explicit datatypes in the partition file
so downstream partitioned inference uses the honest locus contract.

For concatenated phylogenomics inputs, run `alignment occupancy` against an
aligned FASTA plus partition file before tree inference. The command can emit
per-taxon and per-locus TSV tables, a taxon-by-locus occupancy matrix, and a
retained alignment plus remapped partition file after applying explicit
coverage thresholds.

Use `--minimum-locus-occupancy` when partial locus fragments should count as
absent for taxon/locus thresholding. This keeps the retained matrix honest on
datasets where one or two recovered characters should not be treated as
meaningful locus presence. The TSV outputs include `site_coverage_fraction`
columns so the binary coverage calls can be reviewed against overall retained
signal.

```bash
bijux-phylogenetics alignment occupancy \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  artifacts/mixed-locus-supermatrix.partitions.txt \
  --taxon-coverage-threshold 0.6 \
  --locus-coverage-threshold 0.6 \
  --minimum-locus-occupancy 0.75 \
  --taxa-out artifacts/occupancy/taxa.tsv \
  --loci-out artifacts/occupancy/loci.tsv \
  --matrix-out artifacts/occupancy/matrix.tsv \
  --filtered-alignment-out artifacts/occupancy/filtered.fasta \
  --filtered-partitions-out artifacts/occupancy/filtered-partitions.txt \
  --json
```

When the matrix is already aligned and you need to validate the partition file
itself, run `alignment partition-summary` first. That command reports assigned
and unassigned sites, mixed declared datatypes, and one row per locus, and it
can write the review table directly as TSV.

## Partitioned Multi-Locus Inference

Use the concatenation workflow first, then run the partition summary command
before sending the resulting matrix into the adapter inference surface. Pass
the same partition file into the adapter step that needs it.

```bash
bijux-phylogenetics alignment partition-summary \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  artifacts/mixed-locus-supermatrix.partitions.txt \
  --out artifacts/multilocus.partition-summary.tsv \
  --json
bijux-phylogenetics adapter model-select \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions artifacts/mixed-locus-supermatrix.partitions.txt \
  --out-dir artifacts/multilocus-model \
  --prefix multilocus \
  --json
bijux-phylogenetics adapter infer-ml \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions artifacts/mixed-locus-supermatrix.partitions.txt \
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

Those direct adapter runs now preserve IQ-TREE's own `.iqtree` and `.log`
artifacts, plus the model-selection sidecar, a generated
`.model-candidates.tsv`, and `.treefile`, `.ufboot`, or `.contree` where the
invoked step produces them. The emitted JSON and manifests also include the
parsed `selected_model`, `selected_criterion`, `candidate_model_count`,
`best_model_aic`, `best_model_aicc`, `best_model_bic`, `log_likelihood`, and
support-value counts so reviewers can verify the ML result against structured
fields instead of manually scraping engine text files.

For ultrafast bootstrap review specifically, `adapter bootstrap` now writes
three reviewer-facing TSV artifacts alongside the native IQ-TREE files:

- `PREFIX.support.tsv` for every supported internal branch
- `PREFIX.low-support.tsv` for the subset of branches below the governed weak-support threshold
- `PREFIX.support-histogram.tsv` for the support distribution buckets used in reports

Those artifacts map support values back onto explicit descendant-taxon clades,
flag low-support branches directly, and preserve the same `lt50`, `50to69`,
`70to89`, and `ge90` buckets that appear in the manifest and HTML workflow
report.

Use `adapter sh-alrt` when you need SH-aLRT support alongside ultrafast
bootstrap support on the same supported tree.

```bash
bijux-phylogenetics adapter sh-alrt \
  aligned-matrix.fasta \
  --out-dir artifacts/sh-alrt-support \
  --model GTR+G \
  --alrt-replicates 1000 \
  --bootstrap-replicates 1000 \
  --prefix mammals \
  --json
```

That workflow runs IQ-TREE with both `-alrt` and `-bb`, retains the native
`.treefile`, `.ufboot`, `.iqtree`, and `.log` outputs, writes
`mammals.support.tsv` with both support measures on each branch, and writes
`mammals.conflicting-support.tsv` for branches where SH-aLRT and UFBoot imply
different confidence postures under the governed thresholds. The structured
manifest and JSON also expose parsed SH-aLRT minima, maxima, annotated-branch
counts, and conflicting-signal counts directly.

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
  --iqtree-seed 1 \
  --iqtree-threads 1 \
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

Within that engine-artifact directory, the IQ-TREE stages preserve the native
`.iqtree`, `.log`, model-selection sidecar, `.model-candidates.tsv`,
`.treefile`, `.ufboot`, and `.contree` files where each stage produces them.
The bootstrap-support stage also exports `.support.tsv`, `.low-support.tsv`,
and `.support-histogram.tsv` so the supported tree can be reviewed without
re-parsing Newick labels manually. The corresponding manifest entries expose
the parsed selected model, selected criterion, AIC/AICc/BIC winners,
candidate-model counts, log-likelihood, support summary, and weak-backbone
summary directly.

The dedicated SH-aLRT support workflow follows the same pattern: it retains the
native IQ-TREE files, exports `.support.tsv` and
`.conflicting-support.tsv`, and records the parsed combined SH-aLRT/UFBoot
branch summary in the manifest.

The IQ-TREE part of the workflow now defaults to deterministic execution with
`--iqtree-seed 1` and `--iqtree-threads 1`. Ultrafast bootstrap support is the
governed support workflow here, so `--bootstrap-replicates` must be at least
`1000`.

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

The checked real-dataset workflow corpus now lives under:

- `packages/bijux-phylogenetics/tests/fixtures/fasta_to_tree/real/`
- `packages/bijux-phylogenetics/tests/fixtures/expected/fasta_to_tree/`

Those checks pin reviewer-facing output bundles for:

- `gnathostome-ortholog-proteins`
- `gnathostome-ortholog-coding-sequences`
- `strnog-enog411bqtj-proteins`

They verify that the workflow emits the aligned matrix, trimmed matrix, tree,
log, model table, and support table with stable names on real DNA and protein
inputs.

For coding DNA, prefer the codon-aware alignment workflow above and then run
the downstream adapter steps explicitly. The current `adapter fasta-to-tree`
path is still a generic trimAl-based pipeline rather than the codon-preserving
entrypoint.
