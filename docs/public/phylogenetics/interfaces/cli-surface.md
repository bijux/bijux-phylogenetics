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
- `topology ...`
- `adapter ...`
- `bundle` and `report`

The public rule is simple: commands should produce explicit, reviewable outputs
and should not hide important assumptions behind silent defaults.

The topology family provides direct tree-transformation commands for already
inferred trees. `topology root-outgroup` accepts one outgroup taxon or one
expected outgroup clade, writes the rooted tree, and can emit a one-row TSV
report with `--report-out`. That report records requested taxa, matched taxa,
absent taxa, ingroup taxa, whether the matched outgroup is monophyletic in the
input tree, the matched outgroup MRCA, any extra MRCA taxa that break
monophyly, the taxa isolated on the rooted outgroup side, and any rooting
warnings. Its JSON metrics also expose matched, absent, ingroup, rooted
outgroup, rooted ingroup, MRCA spillover, and warning counts so pipelines can
detect non-monophyletic or incomplete outgroup requests without scraping text.

`topology reroot-midpoint` is the exploratory rooted-tree surface when no
explicit outgroup is available. With `--report-out`, it writes a one-row TSV
that records the anchor tip pair used for the selected midpoint path, the
tip-to-tip path length, the midpoint distance from the anchor tip, whether the
midpoint landed on an original node or within an original branch, the taxa on
the anchor side of the new root, the taxa on the opposite side, and whether
the input tree was suitable for straightforward midpoint interpretation. Its
JSON metrics expose the same placement fields plus `midpoint_suitable` and
warning counts so exploratory midpoint-rooted trees can be filtered or flagged
without re-parsing the written TSV.

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

`adapter mrbayes-prepare` is the governed Bayesian input-generation surface for
one aligned FASTA file. It writes a MrBayes NEXUS analysis file with the data
matrix, model block, and MCMC settings, and it also accepts `--partitions` for
same-datatype partition files. When partitions are present, the generated
NEXUS includes named charsets plus one active partition declaration inside the
MrBayes block, and the JSON summary exposes `partitioned`, `partition_count`,
and `partition_warning_count` so review surfaces can separate flat versus
partitioned Bayesian preparation without scraping the written NEXUS text.

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

`adapter infer-large` is the governed aligned-matrix FastTree surface for
larger inputs when you need streamed preflight validation, direct resource
reporting, and resumable output checks. It keeps the inferred tree plus
`.support.tsv`, `.low-support.tsv`, `.support-histogram.tsv`, `.resources.tsv`,
`.log`, and `.manifest.json` sidecars. Its JSON summary exposes sequence
count, alignment length, total site cells, resolved sequence type, resumed
status, timeout budget, and the maximum observed peak-memory measurement across
the recorded workflow stages.

`adapter compare-engines` is the governed side-by-side inference mode for one
aligned matrix. It runs IQ-TREE model selection, IQ-TREE ultrafast bootstrap
support inference, and FastTree approximate inference on the same input, then
emits the two inferred trees, an HTML comparison report, a flat comparison
table, a shared-clade ledger, a conflicting-clade ledger, and a manifest in
one command. Its JSON summary exposes the selected model, shared-taxon count,
Robinson-Foulds distance, shared-clade count, conflicting-clade count, and the
count of support disagreements detected after fraction normalization.

The support-normalization rule is public and narrow by design: FastTree
SH-like local support and IQ-TREE UFBoot support are both rendered as fractions
for side-by-side review, but that normalization does not claim the two support
families are biologically or statistically interchangeable.

`adapter reproducibility` is the governed rerun-check surface for one aligned
matrix when you need to detect deterministic versus unstable supported IQ-TREE
outputs. It runs model selection once to choose a fixed model, reruns the same
bootstrap-supported inference settings multiple times, and emits
`.runs.tsv`, `.comparisons.tsv`, `.support-deltas.tsv`, and `.manifest.json`
artifacts in one command. Its JSON summary exposes the selected model,
`overall_status`, repeat count, unstable-comparison count, and
equivalent-comparison count so automation can separate exact repeatability from
acceptable equivalence and from genuine instability.

For raw input hygiene before alignment, the alignment family now includes
`alignment sequence-type`, `alignment validate-input`, and
`alignment repair-input`. Those commands expose the same raw sequence-type,
duplicate-ID, illegal-character, empty-record, length-outlier, and
identifier-normalization contract that `adapter fasta-to-tree` uses internally,
including the rule that mixed raw inputs must be forced with an explicit
`--sequence-type` before the workflow can continue.
