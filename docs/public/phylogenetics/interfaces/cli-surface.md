---
title: CLI Surface
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-12
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

The alignment family includes matrix-assembly and matrix-audit commands for
concatenated multi-locus inputs. `alignment concatenate` assembles one
supermatrix from aligned per-locus FASTA inputs, preserves taxon identities,
inserts `?` blocks for absent taxa, writes the remapped partition file, and can
emit the taxon-by-locus occupancy matrix in the same run. `alignment occupancy`
then audits an existing concatenated FASTA plus partition file, reports
per-taxon coverage, per-locus coverage, low-coverage flags, explicit
`site_coverage_fraction` summaries, TSV tables, and an optionally filtered
retained matrix with remapped partitions. Use
`--minimum-locus-occupancy` when partial fragments should count as absent for
thresholding instead of being treated as covered from a single observed site.

For partitioned inference preparation, the alignment family also includes
`alignment partition-summary`. It validates one partition file against an
aligned matrix, reports assigned versus unassigned sites, detects mixed declared
datatypes, and can write a stable TSV summary for review before any engine is
invoked.

The adapter family also includes `adapter fasta-to-tree`, which is the
supported end-to-end inference entrypoint for raw FASTA inputs. It emits a
reviewable aligned matrix, trimmed matrix, selected-model table, supported
tree, support summary table, run log, and manifest in one command instead of
forcing users to stitch separate adapter steps together by hand. The workflow
also exposes `--iqtree-seed` and `--iqtree-threads` so the checked inference
bundle can be reproduced exactly when the same engine versions are available.

For coding nucleotide inputs, `adapter align --codon-aware` is the supported
alignment entrypoint. It excludes frame-broken sequences and sequences with
premature stop codons, aligns a translated amino-acid guide, and back-translates
guide gaps as nucleotide triplets so the resulting alignment stays codon-safe
for downstream inference steps.

For aligned multi-locus matrices, `adapter model-select`, `adapter infer-ml`,
and `adapter bootstrap` now accept `--partitions`. On single-datatype matrices
they pass a normalized partition scheme directly to IQ-TREE. On mixed
DNA/protein matrices they materialize one partition alignment per locus and a
generated NEXUS scheme. Fixed single-model requests are not accepted for mixed
DNA/protein runs; use a model-selection keyword such as `MF`, `MFP`, `TEST`, or
`TESTMERGE` instead.

The direct IQ-TREE adapter commands also preserve the native engine artifacts
that correspond to each run. `adapter model-select` keeps `.iqtree`, `.log`,
the native model sidecar, and a generated `.model-candidates.tsv`; `adapter infer-ml`
keeps `.treefile`, `.iqtree`, and `.log`; `adapter bootstrap` keeps `.treefile`,
`.iqtree`, `.log`, `.ufboot`, `.contree`, `.support.tsv`, `.low-support.tsv`,
and `.support-histogram.tsv` when those artifacts exist; `adapter sh-alrt`
keeps `.treefile`, `.iqtree`, `.log`, `.ufboot`, `.support.tsv`, and
`.conflicting-support.tsv`; and `adapter consensus` keeps the consensus
`.contree` with the matching `.iqtree` and `.log`. Their JSON summaries expose
parsed `selected_model`,
`selected_criterion`, `candidate_model_count`, `best_model_aic`,
`best_model_aicc`, `best_model_bic`, `log_likelihood`, support-value counts,
support minima and maxima, weak-support counts, weak-backbone counts, and the
governed support histogram so review surfaces can rely on structured engine
outputs instead of re-parsing free text. The SH-aLRT command also exposes
annotated-branch counts, SH-aLRT minima and maxima, and conflicting-signal
counts for the combined SH-aLRT/UFBoot review surface.

`adapter infer-fast` is the governed FastTree surface for aligned matrices when
speed matters more than fully optimized ML search. It keeps the inferred tree
plus `.support.tsv`, `.low-support.tsv`, and `.support-histogram.tsv` sidecars,
and its JSON summary exposes `approximate_method`, `support_label_kind`,
`support_scale`, `annotated_node_count`, local-support minima and maxima, and
weakly supported clade counts. The public interpretation rule is explicit:
FastTree support labels are SH-like local-support proportions on a `0..1`
scale, and the workflow should be treated as approximate evidence rather than
as a silent substitute for the IQ-TREE ML workflows.

For raw input hygiene before alignment, the alignment family now includes
`alignment sequence-type`, `alignment validate-input`, and
`alignment repair-input`. Those commands expose the same raw sequence-type,
duplicate-ID, illegal-character, empty-record, length-outlier, and
identifier-normalization contract that `adapter fasta-to-tree` uses internally,
including the rule that mixed raw inputs must be forced with an explicit
`--sequence-type` before the workflow can continue.
