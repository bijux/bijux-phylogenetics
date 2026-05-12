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
- `tree-set ...`
- `topology ...`
- `adapter ...`
- `bundle` and `report`

The public rule is simple: commands should produce explicit, reviewable outputs
and should not hide important assumptions behind silent defaults.

`comparative pgls` is the governed regression surface for continuous trait
association under phylogenetic covariance. Its JSON metrics now report
`coefficient_count`, `confidence_interval_count`,
`residual_degrees_of_freedom`, and `coefficient_inference_distribution` so
review tooling can distinguish a minimally identified model from one with
meaningful residual support.

The same command also owns the public formula-expansion contract. Use
`--formula` when the scientific hypothesis is easier to express as one formula
than as separate `--response` plus `--predictors` flags. The formula surface
supports:
- continuous predictors
- categorical predictors
- interaction expansion through `*` and `:`
- intercept-free formulas through `0 + ...` or `... - 1`

Each coefficient row under `data.model.coefficients` preserves the durable
uncertainty contract directly:
- `estimate`
- `standard_error`
- `test_statistic`
- `p_value`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`
- `degrees_of_freedom`
- `inference_distribution`

The coefficient-level inference is explicit. Bijux uses Student-t coefficient
tests and 95% confidence intervals with the fitted residual degrees of freedom,
not a silent large-sample normal approximation. That matters most on smaller
comparative datasets, where a visually large coefficient can still carry wide
intervals and modest nominal support once phylogenetic covariance and limited
taxon counts are taken seriously.

When `--model-matrix-out` is supplied, `comparative pgls` writes the encoded
design matrix as CSV or TSV. Its JSON metrics then also report
`intercept_included`, `model_matrix_row_count`, and
`model_matrix_column_count`, while `data.inputs.model_matrix` preserves the
encoded column names and one taxon-level row per analyzed observation. This
surface exists so reviewers can inspect the actual fitted predictors instead of
inferring the design matrix indirectly from coefficient names.

The `compare` family includes direct topology-distance review for two existing
trees. `compare LEFT RIGHT` now exposes `--rf-mode rooted|unrooted` and
`--taxon-overlap-policy prune-to-shared|require-identical`. The default RF mode
is rooted, which means root placement contributes to the reported
Robinson-Foulds distance. Switch to `--rf-mode unrooted` when you want the
distance to ignore root placement and compare only the recovered splits.

The overlap policy is explicit for partial taxon overlap. By default
`prune-to-shared` computes RF distance only on the taxa present in both trees
and reports the shared, left-only, and right-only taxa in JSON. Use
`--taxon-overlap-policy require-identical` when any taxon-set mismatch should
stop the comparison instead of being pruned away for review.

`compare prune` is the governed shared-taxon pruning surface for two trees. It
reduces both trees to their exact shared taxon set, keeps the richer per-tree
pruning audits from the core pruning layer, and includes one
`post_pruning_comparison` report so the retained trees are compared
immediately instead of leaving that as a manual follow-up step.

When `--out` is supplied, `compare prune` expects an output directory rather
than a single file. The command then writes a stable bundle:
- `left-shared.nwk` and `right-shared.nwk` for the retained trees
- `shared-taxa-pruning.tsv` with one summary row per input tree
- `shared-taxa-removed.tsv` with one row per removed taxon and reason
- `shared-taxa-comparison.tsv` with the retained-tree comparison ledger

The JSON metrics make the main review boundaries explicit: shared taxon count,
removed taxon counts per side, whether the retained trees still match in
topology, and the retained-tree Robinson-Foulds distance.

The same overlap policy also governs `compare branch-lengths`. That command now
keeps the existing shared-clade branch-length table and adds a governed
branch-score summary under `data.branch_score`. The summary reports whether the
trees shared the same taxon set, how many branch-score splits were shared or
unique to one side, and the final Felsenstein branch-score distance when every
matched split has a numeric branch length.

The branch-score contract is explicit. Missing splits contribute as
zero-length branches on the opposite side, so topology disagreement increases
the distance directly. Missing branch lengths on a split that is present do not
silently become zero: the split is counted in `missing_length_split_count` and
the final branch-score distance is reported as unavailable until the missing
length is resolved.

`compare support` is the governed support-aware tree-comparison surface for two
trees over their shared taxon set. It keeps one shared-clade row for clades
present in both trees, normalizes support values onto `0..1` fractions for
comparison, and flags support disagreements when the normalized support delta is
at least `0.15`. It also emits one conflicting-clade row for clades present in
only one tree, then classifies that topology conflict by the strongest observed
support on the present side.

The conflict classification is explicit. A conflicting clade with normalized
support at or above `0.9` is reported as `high_support_conflict`. Conflicting
clades below `0.7` are reported as `low_support_disagreement`. Conflicting
clades between `0.7` and `0.9` are preserved separately as
`moderate_support_disagreement` instead of being overstated as either strong or
weak. If no support label was available on the present side, the row is marked
`support_unavailable`.

When `--out` is supplied, `compare support` writes a flat TSV ledger that keeps
both shared-clade support rows and conflicting-clade severity rows in one file.
That ledger is intended for reviewer-facing conflict triage, not just raw
support extraction.

`compare influence` is the governed leave-one-taxon-out comparison surface for
two trees. It starts from the shared taxon set, excludes one shared taxon at a
time from both trees, reruns the topology comparison and the support-aware
conflict comparison, and ranks taxa by how much those disagreement surfaces
change.

The JSON payload keeps the baseline topology and baseline support reports under
`data.baseline_topology` and `data.baseline_support`, then one ranked row per
excluded taxon under `data.rows`. Each row preserves the retained taxon set,
rooted and unrooted Robinson-Foulds deltas, changes in support disagreement
count, changes in conflicting-clade count, changes in high-support-conflict
count, and the final influence score.

When `--out` is supplied, `compare influence` writes one flat TSV ledger with
one row per excluded taxon. That ledger is intended for reviewer-facing taxon
triage: a high rank means that excluding the taxon materially changes the
disagreement surface, not that the taxon is automatically erroneous or should
always be removed.

`compare clades` is the governed overlap surface for two or more trees. It
takes two required tree paths plus additional trees through repeated `--tree`
flags, computes rooted clade overlap on the shared taxon set, and reports
which clades are present in every tree versus conflicting across trees. Its
JSON payload includes one tree summary per input plus one clade row per
observed clade, and each clade row keeps the per-tree presence flag and parsed
support value when one was available from the tree labels.

When `--out` is supplied, `compare clades` writes a flat TSV ledger with one
row per clade-per-tree observation instead of a wide dynamic matrix. That
format keeps the table stable as tree count changes while still preserving
which tree carried which clade and support value.

`tree-set bootstrap-summary` is the governed summary surface for bootstrap
replicate tree files. It reads one tree-set file, computes clade frequencies,
builds a consensus tree at the chosen threshold, summarizes topology diversity,
and writes a dedicated unstable-branch ledger for consensus branches that fall
below the robust bootstrap threshold or conflict with alternative clades across
the replicates.

The command writes a stable artifact bundle under `--out-dir`:
- `.summary.tsv` with one row of tree-count, topology-diversity, threshold, and consensus summary fields
- `.consensus.nwk` with the consensus topology labeled by bootstrap support percentages
- `.clade-frequencies.tsv` with one row per informative clade frequency
- `.unstable-branches.tsv` with one row per non-robust consensus branch
- `.unstable-clades.tsv` with the broader conflicting-clade ledger across the full replicate set
- `.distance-matrix.tsv` and `.topology-clusters.tsv` for direct topology-variation review

The unstable-branch contract is explicit. A consensus branch is omitted from the
unstable-branch ledger only when its replicate frequency reaches the robust
threshold and no conflicting alternative clade is present. That keeps a
bootstrap summary from overstating a majority-rule consensus as if every branch
were equally stable.

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

`topology clades` is the governed clade-extraction surface for one tree. It
writes one row per node-derived clade, including tips, internal clades, and
the root, and preserves member taxa, parsed support labels, incoming branch
length, root depth, descendant tip depths, and `node_age` when branch lengths
are complete and the tree is ultrametric. When `--metadata` plus repeated
`--metadata-column` flags are supplied, the command also flattens taxon-keyed
metadata into stable per-clade review fields such as matched taxa, missing
taxa, per-taxon values, and distinct observed values for each requested
column.

`topology shape` is the governed tree-shape surface for one tree. It writes one
summary row with Sackin imbalance, Colless imbalance where the tree is strictly
binary, cherry count, topological tree height, branch-length tree height where
all root-to-tip distances are available, and the stable shape summary
`balanced`, `skewed`, or `ladderized`. Its JSON payload also preserves whether
the tree is star-like, comb-like, or unusually imbalanced so review tooling can
filter strongly ladderized or star-topology cases directly.

`topology branch-lengths` is the governed branch-distribution surface for one
tree. It writes one row per non-root branch with the branch length, root depth,
descendant tip count, and explicit flags for missing, zero, negative, long, or
short branches. Its aggregate JSON metrics report branch-count totals plus
minimum, maximum, mean, and median branch length so odd scale shifts are visible
without scraping the full ledger by hand.

`tree-set clades` applies the same clade-table contract to every tree in one
tree-set file. It preserves the one-based source tree index for each row and
requires every tree in the set to carry the same taxon set before clades are
tabulated. That keeps the resulting table reviewable across posterior or
bootstrap samples instead of silently merging incompatible tree contents.

`tree-set shape` applies the same metrics to every tree in one tree-set file
and writes one summary row per sampled tree. Its aggregate JSON metrics count
how many trees are balanced, ladderized, star-like, or comb-like and summarize
mean cherry count, mean Sackin imbalance, and mean tree height. That keeps
shape variation reviewable across posterior or bootstrap samples instead of
collapsing shape into one representative tree too early.

`tree-set branch-lengths` applies the same branch ledger to every tree in one
tree-set file and preserves the one-based source tree index for each branch row.
Its aggregate JSON metrics then summarize set-wide branch-count totals,
zero-length and negative-length counts, long-outlier counts, and overall branch
length minima, maxima, means, and medians. That keeps one pathological sampled
tree from disappearing into a tree-set average or consensus summary.

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

`adapter beast-prepare` is the governed BEAST2 input-generation surface for
one aligned FASTA file plus optional dating metadata. It writes a real BEAST2
XML document rather than a placeholder summary file: the XML includes the
alignment block, a provided-tree or UPGMA starting tree, an HKY or JTT site
model chosen from the inferred sequence alphabet, a strict or uncorrelated
lognormal clock, a Yule or birth-death tree prior, explicit MCMC loggers, and
one MRCA prior per validated calibration target. Its JSON metrics expose
`taxon_count`, `calibration_count`, `tip_date_count`, `warning_count`,
`starting_tree_source`, and `beast_data_type`.

When `--calibrations` or `--tip-dates` are supplied, `--tree` is required.
That rule keeps the dated template anchored to the same guide topology and
named clades that were validated during preparation. Calibration translation is
also explicit in the JSON payload: bounded calibrations are preserved as hard
uniform bounds, while lower-bound-only calibrations are translated into offset
parametric priors with reviewable warnings instead of being copied as if BEAST2
accepted a one-sided hard uniform interval directly.

The warning surface also carries one dated-tree limitation explicitly: if
`--tip-dates` is combined with the standard `birth-death` prior, the XML still
validates, but the JSON warnings mark that combination as exploratory because
BEAST reports that the standard birth-death prior is not serial-sampling
aware.

`adapter beast-log` is the governed parser surface for one BEAST posterior log.
It accepts BEAST's native `Sample` header as well as lowercase `state`
variants, and it can apply `--burnin-fraction` before reporting summaries. Its
JSON metrics expose `row_count`, `column_count`, `burnin_fraction`,
`kept_row_count`, `posterior_parameter_count`, `likelihood_parameter_count`,
`prior_parameter_count`, `clock_parameter_count`, and `tree_parameter_count`.
When `--summary-out` is provided, the command also writes a TSV table with one
row per retained parameter containing effective sample size, mean, median,
sample standard deviation, 95% HPD interval, min/max, first-half mean,
second-half mean, standardized mean shift, and the retained state window.

`adapter beast-parameters` is the governed posterior-parameter diagnostics
surface for one BEAST posterior log. It applies the same burn-in handling as
`adapter beast-log` but reports only the burn-in-aware parameter summaries
instead of the raw parsed rows. Its JSON metrics expose `burnin_fraction`,
`kept_row_count`, `parameter_count`, and `posterior_parameter_count`. When
`--summary-out` is provided, the written TSV contains the same per-parameter
posterior diagnostics table used by the log parser surface.

`adapter beast-convergence` applies the same burn-in handling to ESS and drift
warnings directly. Its JSON metrics expose `warning_count`, `converged`,
`burnin_fraction`, and the post-burn-in `sample_count`, while the payload
lists each warning by parameter and warning code so convergence review does not
depend on re-reading the BEAST log manually.

`adapter beast-burnin-sensitivity` is the governed cross-fraction review
surface for one posterior tree set plus an optional BEAST log. By default it
tests burn-in fractions `0.05`, `0.1`, `0.25`, and `0.5`, though
`--burnin-fractions` can override that set explicitly. Its JSON metrics expose
`slice_count`, `parameter_shift_count`, `unstable_parameter_count`,
`clade_shift_count`, and `unstable_clade_count`. Optional `--slice-out`,
`--parameter-out`, and `--clade-out` ledgers write the per-fraction summary,
the cross-fraction parameter comparison table, and the cross-fraction clade
comparison table respectively.

The instability contract is explicit: a parameter is flagged when the tested
95% HPD intervals do not share a common overlap, and a clade is flagged when
its posterior probability crosses the majority-rule threshold across the
tested burn-in fractions.

`adapter beast-trees` is the governed parser surface for one BEAST posterior
tree file. It accepts native `.trees` NEXUS files, records the sampled
`STATE_*` generations, applies `--burnin-fraction`, extracts clade-frequency
rows from the retained trees, and can emit `--tree-set-out` as normalized
Newick for generic tree-set tooling. Its JSON metrics expose
`total_tree_count`, `kept_tree_count`, `rooted_tree_count`,
`burnin_fraction`, `clade_count`, and `sampled_state_count`. The payload keeps
the retained sampled states, normalized Newick trees, sorted tip set, and
clade-frequency table so posterior-tree review can stay on structured data
instead of manual NEXUS scraping.

`adapter beast-subsample` is the governed retained-subset surface for one BEAST
posterior tree file after optional burn-in handling. It requires
`--method evenly-spaced` plus `--thinning-interval`, or `--method random` plus
`--sample-count`. Random retained subsets become reproducible when `--seed` is
set explicitly. Optional `--tree-set-out` writes the retained normalized Newick
set, and optional `--sample-table-out` writes the retained metadata ledger.

Its JSON metrics expose `total_tree_count`, `burnin_tree_count`,
`pre_subsampling_tree_count`, `retained_tree_count`, `selection_method`, and
`retained_state_count`. The TSV ledger preserves the retained source index,
post-burn-in index, tree name, sampled state, rooted flag, and the governing
selection parameters used to retain that tree.

`adapter beast-consensus` is the governed summary surface for one BEAST
posterior tree file after burn-in handling. It builds a majority-rule consensus
tree from the retained posterior trees, writes that consensus as canonical
Newick through `--out`, can copy the retained normalized posterior tree set
through `--tree-set-out`, and can write a retained clade-frequency ledger
through `--clade-table-out`. Its JSON metrics expose `total_tree_count`,
`kept_tree_count`, `annotated_node_count`, `clade_frequency_count`, and
`burnin_fraction`.

The consensus tree labels use posterior clade probabilities on the `0..1`
scale. The clade-frequency ledger intentionally contains all informative
retained clades, including alternative groupings that do not survive as
majority clades in the final consensus topology, so downstream review can
distinguish a strongly resolved consensus from a superficially simple one.

`adapter beast-diversity` is the governed topology-dispersion surface for one
BEAST posterior tree file after burn-in handling. It can copy the retained
posterior tree set through `--tree-set-out`, write the full pairwise RF ledger
through `--distance-out`, write rooted topology clusters through
`--topology-out`, and write non-unanimous clade instability evidence through
`--unstable-clade-out`. Its JSON metrics expose `total_tree_count`,
`kept_tree_count`, `rooted_topology_count`, `dominant_topology_frequency`,
`pair_count`, `unstable_clade_count`, and `burnin_fraction`.

The command is intentionally broader than a consensus summary. The topology
cluster ledger records how many distinct rooted topologies remain after
burn-in and how often the dominant one appears. The unstable-clade ledger
preserves conflicting alternatives and support classifications for clades that
fail unanimity, so posterior uncertainty review does not collapse immediately
to one representative tree.

`adapter mrbayes-run` is the governed execution surface for one prepared
MrBayes NEXUS file. Its workflow manifest and JSON output now keep the native
posterior tree file (`.run1.t`), parameter trace table (`.run1.p`), MCMC
diagnostics table (`.mcmc`), and consensus tree (`.con.tre`) together so the
downstream review surface can stay on durable engine outputs rather than on
copied snippets.

The matching parser commands expose those artifacts directly:

- `adapter mrbayes-traces` reports trace `row_count` and `column_count` from `.run1.p`
- `adapter mrbayes-parameters` reports `burnin_fraction`, `kept_row_count`, and `parameter_count` from `.run1.p` after applying the requested burn-in cut
- `adapter mrbayes-burnin-sensitivity` reports `slice_count`, `parameter_shift_count`, `unstable_parameter_count`, `clade_shift_count`, and `unstable_clade_count` from `.run1.t` plus optional `.run1.p`
- `adapter mrbayes-trees` reports `tree_count`, `rooted_tree_count`, and `sampled_generation_count` from `.run1.t`
- `adapter mrbayes-subsample` reports `total_tree_count`, `burnin_tree_count`, `pre_subsampling_tree_count`, `retained_tree_count`, `selection_method`, and `retained_generation_count` from `.run1.t`
- `adapter mrbayes-mcmc` reports `row_count`, `column_count`, and `comment_count` from `.mcmc`
- `adapter mrbayes-consensus` reports `tip_count`, `annotated_node_count`, and `maximum_posterior_probability` from `.con.tre`

The `mrbayes-parameters` table and JSON payload expose posterior mean, median,
sample standard deviation, 95% HPD interval, effective sample size, and the
retained generation window for every retained parameter, so trace review does
not have to be reconstructed from the raw `.run1.p` table manually.

`adapter mrbayes-burnin-sensitivity` uses the same default burn-in fractions
`0.05`, `0.1`, `0.25`, and `0.5` unless overridden. Its optional
`--slice-out`, `--parameter-out`, and `--clade-out` outputs mirror the BEAST
burn-in workflow: one per-fraction summary ledger, one parameter-shift ledger,
and one clade-probability-shift ledger. Parameter instability is defined by
non-overlapping 95% HPD intervals, while clade instability is defined by
posterior probabilities that cross the majority-rule threshold across the
tested burn-in fractions.

`adapter mrbayes-subsample` is the governed retained-subset surface for one
MrBayes posterior tree file after optional burn-in handling. Like the BEAST
variant, it supports either evenly spaced thinning or seeded random
subsampling, can write the retained normalized tree set through
`--tree-set-out`, and can write the retained metadata ledger through
`--sample-table-out`.

Its JSON payload keeps the retained source indices and the retained tree
records directly. The TSV ledger preserves the retained source index,
post-burn-in index, tree name, sampled generation, rooted flag, and the
selection parameters so reviewers can reproduce or audit the retained subset
without reopening the full posterior tree file.

The consensus parser exists because MrBayes writes posterior-probability and
branch-length summaries as inline bracket annotations inside the NEXUS tree
text. The governed CLI strips those annotations only after parsing their
probability fields into structured metrics, which keeps the public review
contract honest about what the engine produced.

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
