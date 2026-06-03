---
title: CLI Surface
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-16
---

# CLI Surface

The CLI is the primary operational surface for most users.

For notebook and pipeline users, the workflow-level Python counterparts now
live under `bijux_phylogenetics.api`. They return typed workflow result
objects that wrap the same CLI-grade runtime reports used internally and add
stable JSON plus TSV export helpers.

## Major Command Families

- `validate`, `inspect`, `compare`, `annotate`, `render`
- `demo ...`
- `alignment ...`
- `comparative ...`
- `ancestral ...`
- `biogeography ...`
- `host-association ...`
- `ecological-niche ...`
- `phylogeography ...`
- `tree-set ...`
- `topology ...`
- `adapter ...`
- `bundle` and `report`

The public rule is simple: commands should produce explicit, reviewable outputs
and should not hide important assumptions behind silent defaults.

`demo rabies-cross-host-geography-panel` is the flagship public biological
workflow. In addition to the dataset export, workflow rerun, overview HTML,
and package manifest, it now writes one
`workflow/conclusion-stability/` directory with:

- `conclusion-stability-summary.tsv`
- `key-clade-stability.tsv`
- `support-value-stability.tsv`
- `ancestral-state-stability.tsv`
- `comparative-coefficient-stability.tsv`
- `conclusion-stability-report.html`

Its JSON metrics also report:

- `conclusion_stable_count`
- `conclusion_weak_count`
- `conclusion_unstable_count`

## Runtime Method Tiers

Serious workflow and report commands now publish one explicit method-tier
contract in JSON output so users can distinguish validated inference from
approximate, advisory, or parser-only surfaces.

The governed tier values are:

- `supported`
- `experimental`
- `advisory`
- `parser-only`

When present, the JSON metrics expose:

- `method_tier`
- `method_inference_mode`
- `method_validation_basis`
- `method_approximation`

Tier meaning is strict:

- `supported` requires reference parity or real-engine validation
- `experimental` emits a clear warning and names its approximation when one is used
- `advisory` is review output and should not be read as new inference
- `parser-only` means the command parsed external-engine artifacts and does not claim Bijux ran the inference itself

For automation and downstream notebooks, the canonical reviewer-facing TSV and
JSON outputs now have stable tested schemas. That governed contract covers the
major `.model.tsv`, `.support.tsv`, clade-table, branch-table,
comparative-traits, comparative-summary, event-table, and manifest artifacts
written by the public workflow surfaces.

`parity` is the governed reference-parity validation surface. By default it
checks the repository's core numerical methods against checked-in outputs from
established external tools on small fixtures, with one optional extended suite
for governed primate comparative fits and larger posterior-tree validation.
Its JSON metrics report:

- `all_passed`
- `case_count`
- `method_count`
- `failed_case_count`
- `reference_source`
- `extended`

The command can write:

- `reference-parity-summary.tsv`
- `reference-parity-observations.tsv`

The summary ledger contains one row per method with suite, case counts, pass
counts, fail counts, and contributing reference tools. The observation ledger
contains one row per checked case with:

- `input_fixtures`
- `reference_tool`
- `reference_version`
- `reference_source`
- `tolerance`
- `tolerance_reason`
- `expected_failure_mode`
- `taxon_overlap_policy`
- `shared_taxa`
- `left_only_taxa`
- `right_only_taxa`
- `passed`
- `mismatch_kind`
- `expected_output`
- `observed_output`

The core suite covers RF distance, branch-score distance, PGLS, Pagel's
lambda, Brownian and OU trait models, PIC, Blomberg's K, posterior clade
frequencies, and consensus tree generation. `mismatch_kind` is intentionally a
scientific review surface, not a generic failure flag. It distinguishes
topology, branch length, missing-taxa policy, numerical tolerance, and
model-assumption mismatches.

`parity --reference-source ape-live` switches to the live `ape` execution
harness. That lane runs the checked-in R parity runner through `Rscript`,
compares structured JSON, TSV, and normalized Newick outputs against the Bijux
runtime, and records for each governed case:

- `r_version`
- `ape_version`
- `bijux_version`
- `bijux_commit`
- `function_name`
- `input_fixture`
- `tolerance`
- `passed`
- `mismatch_reason`
- `reproducible_artifact_root`

Its JSON metrics report:

- `all_passed`
- `case_count`
- `function_count`
- `failed_case_count`
- `skipped_case_count`
- `reference_source`

Skipped live cases are not hidden. When `ape` or `Rscript` is unavailable, the
observation ledger records the skip reason and the harness writes one small
reproducible artifact bundle for that case.

The governed live `ape` cases now span shared tree, DNA, and simulation
fixtures. Today that lane covers `ape::read.tree`, `ape::write.tree`,
`ape::consensus`, `ape::prop.clades`, `ape::root`, `ape::unroot`, `ape::drop.tip`, `ape::keep.tip`, `ape::extract.clade`, `ape::getMRCA`, `ape::is.monophyletic`, `ape::cophenetic.phylo`, `ape::dist.topo`, `ape::vcv.phylo`, `ape::node.depth.edgelength`, `ape::branching.times`, `ape::is.ultrametric`, `ape::nj`, `ape::rcoal`, `ape::rtree`, `ape::base.freq`, `ape::seg.sites`, `ape::dist.dna`, and `ape::trans`, with durable inputs
resolved from `shared_tree_fixture_catalog.json`,
`shared_dna_alignment_fixture_catalog.json`, and
`shared_tree_simulation_fixture_catalog.json`. The `ape::read.tree` portion now
compares structured clade rows and covers branch lengths, internal labels,
support labels, quoted labels, one governed multiple-tree Newick input, and
one governed malformed-Newick rejection case. Those cases now exercise one
owned native Newick parser and writer on top of `PhyloTree`, including
location-aware parse failures, rather than routing tree reads through an
external parser. The `ape::root` portion now
uses the same shared tree catalog for one-tip outgroups, monophyletic
multi-tip outgroups, already-rooted trees, missing outgroups, and
non-monophyletic outgroups, and it compares rooted clades plus branch lengths
against live `ape::root` instead of only checking that a rooted flag changed.
On the owned Bijux side, outgroup rooting now runs through the same native
`PhyloTree` manipulation core as unrooting, pruning, clade extraction, MRCA
lookup, and monophyly review rather than delegating that reroot step through
Biopython.
The `ape::unroot` portion now covers rooted binary trees, post-outgroup-rooting
trees, already-unrooted inputs, and malformed input failures, and it makes the
root-edge policy explicit by matching `ape::unroot` branch-length
redistribution instead of silently moving the removed root-edge length into the
expanded clade.
The `ape::drop.tip` portion now covers rooted and unrooted exclusion cases,
unknown excluded tip names, and root-state changes after pruning, and it keeps
one explicit Bijux safety boundary: requests that would leave fewer than two
retained taxa fail clearly instead of producing one-tip workflow outputs.
The `ape::consensus` portion now covers majority-rule and strict consensus
over governed conflicting and posterior-style tree sets, compares one
normalized consensus topology plus one clade-frequency ledger per case, and
fails explicitly when the input tree set does not share one exact taxon set.
The `ape::prop.clades` portion now covers reference-tree clade support mapping
over duplicate, reordered, posterior-style, and mismatched shared tree sets.
It compares one `reference-tree-support.tsv` ledger keyed by descendant tip
set instead of transient node number, and the owned `tree-set support-map`
surface keeps one real `ape` edge case explicit: unsupported root-adjacent
splits are left unscored instead of being mislabeled as zero support.
On the owned Bijux side, `tree-set inspect`, `tree-set consensus`,
`tree-set support-map`, `tree-set compare`, and the posterior tree-set review
commands now read one native `PhyloTree` per Newick record instead of routing
tree-set loading through an external tree object model. Plain `.nwk` and
plain `.trees` inputs both work when the file content is one Newick record per
tree. Strict consensus and support commands keep exact-taxon-set validation as
an explicit hard stop, while tolerant review commands surface malformed-record
skips in their processing metrics instead of failing silently.
The `ape::as.DNAbin` portion now covers clean, lowercase, gap-bearing, and
ambiguity-bearing DNA fixtures. On the owned Bijux side there is no separate
CLI command for that matrix, but the same DNAbin-compatible nucleotide surface
now sits underneath DNA distance, ape-style nucleotide composition,
ape-style segregating-site review, and aligned coding translation. It
preserves taxon order and alignment length, normalizes case, keeps gaps,
ambiguity codes, and explicit missing states literal, writes FASTA back
without nucleotide-state loss, and rejects unsupported symbols explicitly.
The `ape::nj` portion now covers one governed analytical three-taxon matrix
plus four-taxon ultrametric and non-ultrametric matrices. On the owned Bijux
side, `alignment build-tree --method neighbor-joining` and `distance-matrix
build-tree --method neighbor-joining` now use one in-repo deterministic NJ
builder that validates zero-diagonal and nonnegative matrix assumptions and
breaks tied joins by stable taxon ordering instead of delegating the NJ method
through Biopython.

`parity --reference-source phytools-live` switches to the governed live
`phytools` execution harness. It runs one checked-in R parity runner through
`Rscript`, writes structured JSON and TSV outputs, and records for each
governed case:

- `r_version`
- `phytools_version`
- `bijux_version`
- `bijux_commit`
- `function_name`
- `input_fixtures`
- `tolerance`
- `passed`
- `mismatch_reason`
- `reproducible_artifact_root`

Its JSON metrics report:

- `all_passed`
- `case_count`
- `function_count`
- `failed_case_count`
- `skipped_case_count`
- `reference_source`

The initial live `phytools` registry is intentionally narrow for this goal. It
currently covers `phytools::phylosig(method='lambda')`,
`phytools::phylosig(method='K')`, `phytools::fitMk(model='ER')`,
`phytools::fitMk(model='SYM')`, `phytools::fitMk(model='ARD')`,
`phytools::make.simmap(model='ER')`, `phytools::make.simmap(model='SYM')`,
`phytools::make.simmap(model='ARD')`, `phytools::countSimmap`,
`phytools::densityMap`,
`phytools::describe.simmap`,
`phytools::sim.history`,
`phytools::fastBM`,
`phytools::sim.corrs`,
`phytools::pgls.SEy`,
`phytools::rerootingMethod`,
`phytools::fastAnc`, and `phytools::anc.ML` on governed twenty-four-taxon
comparative fixtures from the shared
`shared_phytools_comparative_fixture_catalog.json` corpus. The live lambda
lane now includes one non-ultrametric strong-signal fixture plus one
ultrametric weak-signal fixture so the harness proves both a near-boundary
high-signal fit and a near-zero-signal fit against real `phytools`
likelihood output instead of only one easy interior case. The live K lane now
includes strong-signal and weak-signal seeded permutation cases and compares
the observed K scalar, permutation p-value, and null-distribution summary
under one governed replicate count. The live `fitMk` lane now includes clean
binary, clean multistate, binary missing-value, and multistate missing-value
ER cases plus clean multistate and missing-value-pruned multistate SYM cases.
It now also includes clean and missing-value-pruned binary ARD cases at full
rate-row parity plus clean and missing-value-pruned multistate ARD cases at
summary parity when the optimizer reports weakly identified boundary rates. It
compares flat-root log-likelihood, AIC, AICc, excluded taxa, one explicit
ER-versus-SYM-versus-ARD model identity summary, and the directed rate matrix
when the governed case is identifiable against real `phytools` output.
The same live registry now also includes one governed `phytools::pgls.SEy`
lane for fixed-lambda Brownian covariance PGLS over one simple numeric
regression plus one categorical and one interaction-coded regression. That
claim stays deliberately narrow: installed `phytools 2.5.2` does not export a
general `phytools::pgls` surface, so the live lane proves `pgls.SEy` with
`lambda = 1.0`, while the broader exact PGLS contract for estimated lambda
and full coefficient parity remains the checked-in `ape` plus `nlme`
reference suite.
The live `make.simmap` lane now includes clean binary, clean multistate, and
missing-value-pruned binary ER cases; clean multistate and
missing-value-pruned multistate SYM cases; and binary plus
missing-value-pruned binary ARD cases. It also keeps governed multistate ARD
cases on summary-envelope parity only when weakly identified boundary rates
make row-level comparison untrustworthy across optimizers. The owned
`discrete-evolution stochastic-map` surface now also reports fitted-model
identity, parameter count, log-likelihood, AIC, AICc, baseline-model
comparison, optimizer convergence, and weak-fit warnings alongside the seeded
simulation output. The parity lane compares distributional envelopes only:
excluded taxa, total-transition-count mean plus interval,
transition-count summary rows, and time-in-state summary rows. It does not
claim exact stochastic-history identity with real `phytools`. It also writes
one flat branch-segment TSV, one per-state time-summary TSV, and one per-branch
state-occupancy TSV on the owned Bijux side through
`discrete-evolution stochastic-map`.
The owned CLI now also exposes one `discrete-evolution count-maps` surface over
saved stochastic-map collections. It writes one per-replicate count matrix,
one aggregate transition matrix, one per-branch directional transition table,
and one flat event ledger. The live `countSimmap` lane now covers clean
binary, clean multistate, clean multistate SYM, and missing-value-pruned
binary cases. It compares total-transition envelopes plus directional
transition-count rows, including zero diagonal state pairs, without claiming
exact stochastic-history identity.
The same owned CLI now also exposes one `discrete-evolution density-maps`
surface over saved stochastic-map collections. It writes one branch-probability
table, one branch-level density envelope, one slice-level probability table at
the requested resolution, and one report-ready HTML or SVG artifact. The live
`densityMap` lane is intentionally narrower than the owned CLI surface: it
currently covers binary ER collections only, including one
missing-value-pruned case. It compares per-branch posterior probability
summaries and branch-level uncertainty against real `phytools::densityMap`
without claiming pixel-perfect plotting parity.
The live `describe.simmap` lane now covers clean binary, clean multistate,
clean multistate SYM, and missing-value-pruned binary cases. It compares the
owned summary surface directly, including total changes, transition rows,
time-in-state rows, and per-branch state-occupancy rows.
The owned CLI now also exposes one `simulate history-discrete` surface for
fixed-tree discrete-history simulation from an explicit rate matrix. It writes
one tip-state truth table, one node-state truth table, one branch-history
truth table, one transition-event ledger, one branch-segment ledger, and one
parity summary table with transition-count, time-in-state, and tip-state
frequency rows. The live `sim.history` lane now covers governed binary and
multistate no-change plus high-rate fixed-tree cases and compares those
distribution-summary envelopes against real `phytools::sim.history` without
claiming exact history identity across languages.
The owned `simulate traits-brownian` surface now accepts either `--sigma` or
`--sigma-squared` and reports the resolved Brownian rate parameter in JSON
output. The same simulation family now also owns one Brownian replicate-review
surface over tip distributions and tip covariances, which underlies the live
`phytools::fastBM` lane for governed low-variance, root-shift high-variance,
and six-taxon fixed-tree cases. That lane compares summary envelopes and
tip-covariance rows against real `phytools::fastBM` without claiming
cross-language draw identity.
The same simulation family now also owns one
`simulate traits-brownian-correlated` surface for two or more continuous
traits on one fixed tree from one explicit evolutionary covariance matrix. It
accepts either repeated `--covariance-row` values directly or repeated
`--correlation-row` values plus one `--trait-standard-deviation` per trait,
writes one long-form replicate tip-trait ledger plus one optional summary
ledger, and reports trait count, replicate count, and the generating
covariance contract in JSON output. The live `phytools::sim.corrs` lane now
covers governed low-correlation, negative-correlation root-shift, and
three-trait six-taxon cases. It compares summary envelopes, tip-covariance
rows, and tip-correlation rows against real `phytools::sim.corrs` without
claiming exact cross-language draw identity.
The live
`rerootingMethod` lane now includes governed ER binary, ER multistate, ER
missing-value-pruned, SYM multistate, and SYM missing-value-pruned cases. It
compares one flat node-probability ledger keyed by stable node signature and
state label against real `phytools` output. That claim is intentionally
narrow: `phytools::rerootingMethod` is only governed here for ER or SYM under
the equal root prior inherited from `fitMk`, while ARD, Fitch, ordered-state,
empirical-root-prior, and fixed-root-prior ancestral runs remain owned Bijux
review surfaces without a live `phytools` parity claim. The live
`fastAnc` lane now includes ultrametric strong-signal, ultrametric
weak-signal, non-ultrametric strong-signal, and missing-value pruning cases,
and compares stable node-signature rows plus standard errors against real
`phytools` output. The live `anc.ML` lane now covers the same four fixture
shapes and compares stable node-signature rows, standard errors, 95%
intervals, Brownian log-likelihood, and fitted sigma-squared against real
`phytools` output.
For this round, `bionj` is explicitly excluded. The distance-tree CLI surfaces
therefore accept `--method bionj` only so the owned runtime can return one
structured out-of-scope error naming `ape::bionj`, rather than failing with
one generic parser-choice message.
The `ape::dist.dna` portion now covers raw nucleotide distance, JC69, K80,
F81, and TN93 distance over governed clean, gapped pairwise-deletion, gapped
complete-deletion, ambiguity-bearing, identical-sequence, high-divergence,
missing-data, and unequal-length-invalid fixtures. On the owned Bijux side,
`alignment distance-matrix --model raw`, `--model jc69`, `--model k80`,
`--model f81`, and `--model tn93` accept the ape-compatible aliases while
keeping `p-distance`, `jukes-cantor`, `kimura-2-parameter`,
`felsenstein-81`, and `tamura-nei-93` as the canonical internal labels.
Saturated JC69, K80, F81, and TN93 pairs are reported explicitly as either
undefined or infinite, `--components-out` writes one pairwise component
ledger for review, `--parameters-out` writes one model-parameter ledger for
reviewer-facing base-frequency and coefficient inspection, and unequal-length
alignments fail explicitly before any matrix is written. TN93 also warns
explicitly when the resolved alignment composition omits a nucleotide instead
of silently degrading to JC69 or K80.
The `ape::base.freq` portion now covers lowercase, ambiguity-bearing,
missing-data, and all-gap-or-missing alignments. On the owned Bijux side,
`alignment composition --base-frequency-out <table.tsv>` writes one combined
alignment-plus-sequence TSV ledger with `scope`, `identifier`, `state`,
`count`, and `frequency` columns, returns the same literal-state frequencies
in JSON output, and reports composition outlier sequences beside those base
frequency rows. Ambiguity codes, gaps, and explicit missing states are counted
as literal states to match `ape::base.freq`, and all-gap or missing inputs
warn explicitly instead of fabricating canonical A/C/G/T content.
The `ape::seg.sites` portion now covers lowercase, invariant,
one-variable-site, gap-bearing, ambiguity-bearing, missing-data, and
all-gap-or-missing alignments. On the owned Bijux side,
`alignment segregating-sites --site-table-out <table.tsv>` writes one
`segregating-sites.tsv` ledger with site positions plus literal and
ape-normalized state summaries. Leading and trailing gaps are normalized to
`N` to match live `ape::seg.sites`, explicit missing states do not create
segregating sites by themselves, and incompatible ambiguity states or internal
gaps remain governed live parity cases instead of being flattened into a total
count only.
The `ape::trans` portion now covers valid-reading-frame, ambiguous-codon,
internal-stop, terminal-stop, frame-truncation, and vertebrate-mitochondrial
genetic-code fixtures. On the owned Bijux side, `alignment translate
--codon-validation-out <table.tsv> --excluded-sequences-out <table.tsv>`
writes one amino-acid FASTA plus codon-level validation rows. The aligned
translation surface truncates trailing partial codons with the same explicit
warning as live `ape::trans`, while the stricter codon-preparation surface
still owns pre-alignment sequence exclusion for frame errors, ambiguous
codons, and premature stop codons.
The `ape::keep.tip` portion now covers valid rooted and unrooted keep-set
cases, selected-tip order differences, and rootedness changes after pruning.
Bijux keeps the workflow-facing absent-requested-taxon report and minimum-two
retained-taxa stop as explicit product extensions rather than pretending those
paths are live `ape::keep.tip` parity.
The `ape::extract.clade` portion now covers rooted root-clade and internal-node
subtree extraction plus explicit tip-node and out-of-bounds failures. Bijux
keeps one adjacent owned surface outside the live `ape` call shape too:
callers can extract the same subtree by exact descendant-taxa identity instead
of only by ape-style node number.
The `ape::getMRCA` portion now covers stable internal-node identity for
two-tip, many-tip, full-tip-set, duplicate-tip, rooted-polytomy, and
already-rooted-outgroup cases. Bijux keeps one adjacent workflow-side rule
outside the live `ape` call shape too: missing requested taxa fail clearly
instead of surfacing as a low-level parser-side error.
The `ape::is.monophyletic` portion now covers rooted and unrooted monophyly
calls with explicit reroot policy, full-tip-set behavior, singleton and
mixed-missing requests, rooted-polytomy behavior, post-rooting behavior, and
all-missing reroot failures. Bijux uses that same lane to expose matched MRCA
node identity and extra descendant taxa when a direct clade is not cleanly
monophyletic.
The `ape::cophenetic.phylo` portion now covers rooted and unrooted branch-length
trees, compares one governed long-form tip-distance ledger rather than only a
printed matrix, and keeps the taxon order explicit in the summary payload. On
the owned Bijux side, tip-distance calculations now reject missing branch
lengths unless the caller opts into one explicit unit-length fallback policy.
The `ape::dist.topo` portion now covers identical rooted trees, rooted
child-order rotations, one-conflict rooted pairs, rooted tree-versus-polytomy
pairs, one governed unrooted split conflict, and one governed 128-tip rooted
pair. It compares one explicit RF-style split ledger rather than only a
scalar distance, keeps rooted-versus-unrooted policy explicit per case, and
aligns directly with the owned `adapter compare --split-table-out` review
surface. Those same split rows now come from one native clade-set core shared
with support comparison, tree-set support mapping, posterior clade summaries,
and live `ape::dist.topo` parity rather than parallel helper
implementations.
The `ape::vcv.phylo` portion now covers rooted ultrametric, rooted
non-ultrametric, unrooted branch-length, and singular zero-branch trees. It
compares one governed long-form Brownian shared-ancestry covariance ledger,
persists the compared covariance tables automatically when parity fails, and
keeps the taxon order explicit in the summary payload. On the owned Bijux
side, `summarize_brownian_covariance(...)` now rejects missing or negative
branch lengths explicitly and reports singular-versus-near-singular state from
the raw covariance matrix instead of silently regularizing it away.
The `ape::node.depth.edgelength` portion now covers rooted ultrametric,
rooted non-ultrametric, zero-branch-length, and post-outgroup-rooting trees.
It compares one governed node-depth table keyed by stable ape-style node ids,
and the owned Bijux surface `compute_tree_node_depths(...)` rejects incomplete
branch lengths instead of substituting edge counts or implied zeros.
The `ape::branching.times` portion now covers rooted ultrametric trees with
and without internal labels, one medium ultrametric tree, and one zero-length
internal-branch ultrametric tree. It compares one governed internal-node
branching-time table keyed by stable ape-style node ids, and the owned Bijux
surface `compute_tree_branching_times(...)` rejects non-ultrametric trees
instead of forwarding the invalid negative or inconsistent node ages that
`ape::branching.times` can still produce on those inputs.
The `ape::gammaStat` portion now covers rooted ultrametric trees with and
without internal labels, one medium ultrametric tree, and one zero-internal-
branch ultrametric tree. It compares one governed one-row
`gamma-statistic.tsv` ledger, and the owned Bijux surface
`compute_diversification_gamma_statistic(...)` keeps two workflow-side
boundaries explicit instead of inheriting `ape`'s looser behavior: the tree
must stay fully bifurcating, and incomplete sampling remains a warning surface
rather than an implicit correction.
The `ape::is.ultrametric` portion now covers exact ultrametric,
near-ultrametric, tight-tolerance near-ultrametric, and clearly
non-ultrametric trees. It compares one governed tip-depth diagnostic table,
and the owned Bijux surface `assess_tree_ultrametricity(...)` reports the
criterion name, criterion value, tolerance, maximum tip-depth deviation,
offending taxa, and a deterministic `ultrametric-diagnostics.tsv` ledger.
That same ape-style surface is now reused before rooted Brownian, OU, and
diversification workflows claim time-tree compatibility.
The `ape::write.tree` portion
roundtrips Bijux-written Newick through live `ape` for rooted, unrooted,
internal-label, support-label, quoted-label, and multiple-tree cases. The DNA
cases include lowercase input, ambiguity, missing data, identical sequences,
high-divergence distances, and valid, ambiguous-codon, internal-stop,
terminal-stop, frame-truncation, or alternate-genetic-code coding translation
rows. Unequal-length DNA fixtures still stay on the diagnostic side of the
contract for distance workflows, but frame-error coding fixtures now stay in
the governed `ape::trans` parity-pass registry because the owned aligned
translation surface truncates trailing partial codons with the same explicit
warning that live `ape::trans` emits.

Bijux does not silently serialize malformed trees in that lane. Unnamed tips,
empty tree sets, and non-finite branch lengths fail on the Bijux side before
the live `ape` comparison is attempted.

Tree IO equality in that lane is structural rather than string-based. The
governed comparison accepts reordered-but-equivalent children and emits
specific mismatch reasons when rootedness, tip labels, clades or splits,
branch lengths, or internal labels differ.

The governed PGLS lane is not limited to one intercept-plus-slope example. The
core suite now includes one fixed-Brownian numeric regression, one
treatment-coded categorical regression, and one treatment-coded interaction
regression checked against R `ape` plus `nlme` outputs for coefficients,
standard errors, p-values, likelihood, AIC, and encoded model-matrix rows. The
extended suite adds one governed estimated-lambda primate regression against
the same external tool chain.

`benchmark stress-suite` is the governed large-dataset resource review surface.
It executes five owned workload families on one selected tier:

- large alignment inference
- multi-locus supermatrix assembly
- posterior or bootstrap tree-set consensus
- comparative independent contrasts
- tree-annotation table generation

Its JSON metrics report:

- `observation_count`
- `tier`

Each observation records:

- `input_size_bytes`
- `sequence_count` when applicable
- `alignment_length` when applicable
- `tree_count` when applicable
- `taxon_count` when applicable
- `locus_count` when applicable
- `runtime_seconds`
- `peak_memory_bytes`
- `memory_observation_kind`
- `output_row_count`

Use `--tier small` for the routine stress lane and `--tier heavy` for the
optional `1,000+` sequence and `1,000+` tree pressure check. The surface is
explicit about scope: it measures the repository's owned workflows, not a
synthetic micro-benchmark disconnected from user-facing outputs.

`report release-truth` is the governed pre-release summary surface. It
consumes actual pytest JUnit XML reports for the full test lane and the
real-engine lane, reruns the owned workflow-validation, release-gate, parity,
and stress-suite checks, and writes one HTML report plus one machine manifest.
Its JSON metrics report:

- `total_tests`
- `total_tests_passed`
- `total_tests_failed`
- `total_tests_skipped`
- `real_engine_tests`
- `real_engine_tests_passed`
- `real_engine_tests_failed`
- `real_engine_tests_skipped`
- `supported_workflow_count`
- `experimental_workflow_count`
- `flagship_dataset_count`
- `reference_parity_case_count`
- `stress_workload_count`

`report tree-package` is the governed full tree review surface. It takes one
tree and materializes a richer review directory than the older `report tree`
diagnostic. Its JSON metrics report:

- `tip_count`
- `supported_branch_count`
- `rendered_support_count`
- `long_outlier_count`

The command writes:

- `tree-report.html`
- `tree-image.svg`
- `support-table.tsv`
- `clade-table.tsv`
- `branch-stats.tsv`
- `tree-report.manifest.json`

The HTML report embeds the owned SVG tree image directly and includes reviewer
summary, support, clade, and branch-stat sections. The TSV ledgers remain the
durable flat review contract for downstream inspection and automation. Use
`report tree` when only the lightweight structural and forensic HTML audit is
needed; use `report tree-package` when the image and tabular review outputs are
required together. Its JSON and HTML surfaces now mark the package as
`advisory` rather than inference.

`demo primate-comparative` is the governed packaged mammal dataset surface. It
materializes the shipped primate comparative dataset into one output directory
and reruns the owned comparative workflow bundle over those packaged inputs.
Its JSON metrics report:

- `artifact_count`
- `dataset_taxon_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/tree.nwk`
- `dataset/traits.csv`
- `dataset/expected/*.tsv`
- `workflow/workflow-summary.tsv`
- `workflow/pgls-lambda-profile.tsv`
- `workflow/brownian-summary.tsv`
- `workflow/ou-summary.tsv`
- `workflow/signal-summary.tsv`
- `workflow/signal-permutations.tsv`
- `workflow/continuous-ancestral-summary.tsv`
- `workflow/continuous-ancestral-uncertainty.tsv`
- `workflow/discrete-ancestral-summary.tsv`
- `workflow/discrete-ancestral-probabilities.tsv`
- `overview.md`

The packaged dataset is keyed by `species`, carries both continuous and
categorical traits, and uses the following governed comparative workflow
choices:

- PGLS response `longevity`
- PGLS predictor `social_group_size`
- continuous ancestral trait `longevity`
- discrete ancestral trait `mating_system`

This command is intentionally a public data-and-workflow entrypoint rather than
another evidence-book wrapper. It gives users a real mammal comparative dataset
without requiring them to know the repository’s internal study layout.

`demo avian-reproductive-traits` is the governed packaged bird dataset surface.
It materializes the shipped avian reproductive dataset into one output
directory and reruns the owned comparative workflow bundle over those packaged
inputs. Its JSON metrics report:

- `artifact_count`
- `dataset_taxon_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/tree.nwk`
- `dataset/traits.csv`
- `dataset/expected/*.tsv`
- `workflow/workflow-summary.tsv`
- `workflow/pgls-lambda-profile.tsv`
- `workflow/brownian-summary.tsv`
- `workflow/ou-summary.tsv`
- `workflow/signal-summary.tsv`
- `workflow/signal-permutations.tsv`
- `workflow/continuous-ancestral-summary.tsv`
- `workflow/continuous-ancestral-uncertainty.tsv`
- `workflow/discrete-ancestral-summary.tsv`
- `workflow/discrete-ancestral-probabilities.tsv`
- `workflow/clade-trait-summary.tsv`
- `workflow/clade-trait-clades.tsv`
- `overview.md`

The packaged dataset is keyed by `species`, carries both continuous and
categorical reproductive traits, and uses the following governed comparative
workflow choices:

- PGLS response `testes_mass`
- PGLS predictor `body_mass`
- continuous ancestral trait `testes_mass`
- discrete ancestral trait `mating_system`
- clade summary trait `mating_system`

This command is intentionally a public bird comparative entrypoint rather than
an internal teaching-data wrapper. It gives users a real bird dataset that can
exercise trait-evolution and clade-pattern workflows immediately.

`demo central-european-seashore-flora` is the governed packaged plant dataset
surface. It materializes the shipped Central European flora subset into one
output directory and reruns the owned comparative workflow bundle over those
packaged inputs. Its JSON metrics report:

- `artifact_count`
- `dataset_taxon_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/tree.nwk`
- `dataset/traits.csv`
- `dataset/expected/*.tsv`
- `workflow/workflow-summary.tsv`
- `workflow/pgls-lambda-profile.tsv`
- `workflow/brownian-summary.tsv`
- `workflow/ou-summary.tsv`
- `workflow/signal-summary.tsv`
- `workflow/signal-permutations.tsv`
- `workflow/continuous-ancestral-summary.tsv`
- `workflow/continuous-ancestral-uncertainty.tsv`
- `workflow/discrete-ancestral-summary.tsv`
- `workflow/discrete-ancestral-probabilities.tsv`
- `workflow/clade-trait-summary.tsv`
- `workflow/clade-trait-clades.tsv`
- `overview.md`

The packaged dataset is keyed by `species`, carries both continuous and
categorical plant traits, and uses the following governed comparative workflow
choices:

- PGLS response `seed_mass`
- PGLS predictor `plant_height`
- continuous ancestral trait `seed_mass`
- discrete ancestral trait `lifeform`
- clade summary trait `lifeform`

This command is intentionally a public non-animal comparative entrypoint rather
than a generic flora dump. It exposes one documented published subset with a
fully rerunnable workflow contract.

`demo influenza-a-ha-reference-panel` is the governed packaged viral
sequence-to-tree surface. It materializes the shipped Influenza A
hemagglutinin FASTA panel into one output directory and reruns the owned
MAFFT, trimAl, and IQ-TREE workflow over those packaged inputs. Its JSON
metrics report:

- `artifact_count`
- `sequence_count`
- `sequence_type`
- `selected_model`
- `minimum_support`
- `maximum_support`
- `weakly_supported_clade_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/sequences.fasta`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/influenza-a-ha-reference-panel.aln`
- `workflow/influenza-a-ha-reference-panel.trimmed.aln`
- `workflow/influenza-a-ha-reference-panel.tree`
- `workflow/influenza-a-ha-reference-panel.model.tsv`
- `workflow/influenza-a-ha-reference-panel.support.tsv`
- `workflow/influenza-a-ha-reference-panel.log`
- `workflow/influenza-a-ha-reference-panel.manifest.json`
- `overview.md`

The packaged dataset carries raw unaligned viral nucleotide sequences and uses
the following governed inference controls:

- sequence type `dna`
- IQ-TREE seed `1`
- IQ-TREE threads `1`
- bootstrap replicates `1000`

This command requires MAFFT, trimAl, and IQ-TREE executables. Use
`--mafft-executable`, `--trimal-executable`, and `--iqtree-executable` when
they are not available on the default `PATH`.

`demo gnathostome-ortholog-protein-benchmark` is the governed packaged protein
sequence-to-tree surface. It materializes one shipped gnathostome ortholog
amino-acid FASTA panel into one output directory and reruns the owned MAFFT,
trimAl, and IQ-TREE workflow over those packaged protein inputs. Its JSON
metrics report:

- `artifact_count`
- `sequence_count`
- `sequence_type`
- `selected_model`
- `alignment_length`
- `trimmed_alignment_length`
- `minimum_support`
- `maximum_support`
- `weakly_supported_clade_count`
- `state_space`
- `model_selection_scope`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/sequences.fasta`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/molecular-assumptions.tsv`
- `workflow/gnathostome-ortholog-protein-benchmark.aln`
- `workflow/gnathostome-ortholog-protein-benchmark.trimmed.aln`
- `workflow/gnathostome-ortholog-protein-benchmark.tree`
- `workflow/gnathostome-ortholog-protein-benchmark.model.tsv`
- `workflow/gnathostome-ortholog-protein-benchmark.support.tsv`
- `workflow/gnathostome-ortholog-protein-benchmark.log`
- `workflow/gnathostome-ortholog-protein-benchmark.manifest.json`
- `overview.md`

The packaged dataset carries raw unaligned amino-acid sequences and uses the
following governed inference controls:

- sequence type `protein`
- IQ-TREE sequence keyword `AA`
- IQ-TREE seed `1`
- IQ-TREE threads `1`
- bootstrap replicates `1000`

This command writes one explicit molecular-assumption ledger because it is
meant to distinguish amino-acid inference from the repository's DNA demos. The
workflow starts from protein FASTA directly, searches protein models only, and
does not apply coding-DNA translation, codon-position partitioning, or
nucleotide-specific interpretation such as GC composition.

This command requires MAFFT, trimAl, and IQ-TREE executables. Use
`--mafft-executable`, `--trimal-executable`, and `--iqtree-executable` when
they are not available on the default `PATH`.

`demo pleistocene-bear-cytb-fragments` is the governed packaged ancient-DNA
sequence-to-tree surface. It materializes the shipped degraded bear
cytochrome b panel into one output directory and reruns the owned MAFFT,
trimAl, and IQ-TREE workflow with explicit missingness review over those
packaged inputs. Its JSON metrics report:

- `artifact_count`
- `sequence_count`
- `degraded_sequence_count`
- `selected_model`
- `minimum_support`
- `maximum_support`
- `removed_column_count`
- `cleaned_missing_data_fraction`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/sequences.fasta`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/missingness-effects.tsv`
- `workflow/pleistocene-bear-cytb-fragments.aln`
- `workflow/pleistocene-bear-cytb-fragments.trimmed.aln`
- `workflow/pleistocene-bear-cytb-fragments.cleaned.aln`
- `workflow/pleistocene-bear-cytb-fragments.tree`
- `workflow/pleistocene-bear-cytb-fragments.model.tsv`
- `workflow/pleistocene-bear-cytb-fragments.support.tsv`
- `overview.md`

The packaged dataset carries raw unaligned bear cytochrome b sequences and
uses the following governed review controls:

- sequence type `dna`
- site missingness threshold `0.15`
- sequence missingness threshold `0.15`
- IQ-TREE seed `1`
- IQ-TREE threads `1`
- bootstrap replicates `1000`

This command requires MAFFT, trimAl, and IQ-TREE executables. Use
`--mafft-executable`, `--trimal-executable`, and `--iqtree-executable` when
they are not available on the default `PATH`.

`demo catarrhine-mitogenome-five-locus-panel` is the governed packaged
multi-locus phylogenomics surface. It materializes the shipped catarrhine
mitochondrial panel into one output directory and reruns the owned
concatenation, occupancy, and partitioned IQ-TREE workflow over those packaged
inputs. Its JSON metrics report:

- `artifact_count`
- `taxon_count`
- `locus_count`
- `alignment_length`
- `partition_count`
- `selected_model`
- `minimum_support`
- `maximum_support`
- `weakly_supported_clade_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/taxa.csv`
- `dataset/loci/*.fasta`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/catarrhine-mitogenome-five-locus-panel.supermatrix.fasta`
- `workflow/catarrhine-mitogenome-five-locus-panel.partitions.txt`
- `workflow/occupancy-taxa.tsv`
- `workflow/occupancy-loci.tsv`
- `workflow/occupancy-matrix.tsv`
- `workflow/catarrhine-mitogenome-five-locus-panel.partition-summary.tsv`
- `workflow/catarrhine-mitogenome-five-locus-panel.model-candidates.tsv`
- `workflow/catarrhine-mitogenome-five-locus-panel.supported.tree`
- `workflow/catarrhine-mitogenome-five-locus-panel.support.tsv`
- `overview.md`

The packaged dataset carries one explicit aligned-locus contract:

- five mitochondrial coding loci
- one shared six-taxon identifier set across every locus
- partitioned IQ-TREE inference over the concatenated supermatrix

This command requires only IQ-TREE. The loci are already aligned, so the
governed demo focuses on multi-locus assembly and partitioned inference rather
than raw alignment generation.

`demo catarrhine-data-quality-stress-panel` is the governed packaged dirty-data
stress surface. It materializes the shipped catarrhine stress panel into one
output directory and reruns the owned audit-and-cleanup workflow over its raw
alignment, raw FASTA validation, coding-sequence, tree, and trait inputs. Its
JSON metrics report:

- `artifact_count`
- `raw_taxon_count`
- `cleaned_taxon_count`
- `duplicate_sequence_identifier_count`
- `illegal_character_count`
- `empty_sequence_count`
- `raw_sequence_length_outlier_count`
- `duplicate_trait_taxon_count`
- `missing_trait_value_count`
- `sequence_outlier_count`
- `tree_zero_length_branch_count`
- `tree_negative_branch_count`
- `tree_long_branch_outlier_count`
- `coding_frame_error_count`
- `coding_internal_stop_count`
- `raw_trait_missing_from_traits_count`
- `raw_trait_extra_taxon_count`
- `dropped_taxon_count`
- `repaired_branch_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/raw/alignment.fasta`
- `dataset/raw/sequence-input.fasta`
- `dataset/raw/coding-sequences.fasta`
- `dataset/raw/tree.nwk`
- `dataset/raw/traits.csv`
- `dataset/raw/traits-mismatch.csv`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/raw-sequence-findings.tsv`
- `workflow/raw-sequence-repair.tsv`
- `workflow/repaired-sequence-input.fasta`
- `workflow/repaired-sequence-validation.tsv`
- `workflow/coding-sequence-exclusions.tsv`
- `workflow/prepared-coding-sequences.fasta`
- `workflow/raw-trait-linkage.tsv`
- `workflow/trait-duplicates.tsv`
- `workflow/trait-missing-values.tsv`
- `workflow/sequence-outliers.tsv`
- `workflow/tree-issues.tsv`
- `workflow/repair-actions.tsv`
- `workflow/cleaned-traits.csv`
- `workflow/cleaned-alignment.fasta`
- `workflow/cleaned-tree.nwk`
- `workflow/cleaned-linkage.tsv`
- `workflow/cleaned-validation.tsv`
- `overview.md`

The packaged stress contract is explicit about scope:

- the raw alignment is already aligned and is stress-tested through composition review, not alignment generation
- the raw sequence-validation FASTA is intentionally dirty in duplicate identifiers, illegal characters, empty records, and length outliers
- the raw coding FASTA is intentionally dirty in frame consistency and premature stop codons
- the raw tree is intentionally dirty in zero and negative branch lengths plus one extreme long branch, not in syntax
- the raw trait table is intentionally dirty in duplicates and missingness
- the raw trait-mismatch table is intentionally wrong in taxon overlap and is kept as a failure-review surface
- the workflow resolves or excludes fixable inputs deterministically, records strict mismatch failure where repair would be dishonest, and writes one cleaned comparative subset instead of mutating raw inputs in place

`demo known-answer-reference-panel` is the governed packaged known-answer
simulation surface. It materializes the shipped deterministic simulation panel
into one output directory and reruns the owned recovery workflow over the
packaged true tree, simulated alignment, and simulated traits. Its JSON
metrics report:

- `artifact_count`
- `taxon_count`
- `sequence_length`
- `distance_method`
- `distance_model`
- `rooted_topology_equal`
- `same_unrooted_topology`
- `same_taxa_different_rooting`
- `robinson_foulds_distance`
- `parameter_row_count`
- `threshold_pass_count`
- `threshold_row_count`
- `continuous_internal_node_mean_absolute_error`
- `discrete_internal_node_accuracy`
- `host_internal_node_accuracy`
- `host_event_accuracy`
- `geographic_internal_node_accuracy`
- `geographic_event_accuracy`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/true-tree.nwk`
- `dataset/simulated-alignment.fasta`
- `dataset/continuous-traits.tsv`
- `dataset/ou-traits.tsv`
- `dataset/discrete-traits.tsv`
- `dataset/host-traits.tsv`
- `dataset/geographic-traits.tsv`
- `dataset/true-parameters.tsv`
- `dataset/true-continuous-nodes.tsv`
- `dataset/true-ou-nodes.tsv`
- `dataset/true-discrete-nodes.tsv`
- `dataset/true-host-nodes.tsv`
- `dataset/true-geographic-nodes.tsv`
- `dataset/true-host-switch-events.tsv`
- `dataset/true-geographic-transition-events.tsv`
- `dataset/recovery-thresholds.tsv`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/recovered-distance-tree.nwk`
- `workflow/tree-recovery.tsv`
- `workflow/parameter-recovery.tsv`
- `workflow/brownian-fit-summary.tsv`
- `workflow/ou-fit-summary.tsv`
- `workflow/continuous-ancestral-summary.tsv`
- `workflow/continuous-ancestral-uncertainty.tsv`
- `workflow/continuous-node-recovery.tsv`
- `workflow/discrete-ancestral-summary.tsv`
- `workflow/discrete-ancestral-probabilities.tsv`
- `workflow/discrete-node-recovery.tsv`
- `workflow/host-switch-summary.tsv`
- `workflow/host-state-nodes.tsv`
- `workflow/host-switch-branches.tsv`
- `workflow/host-node-recovery.tsv`
- `workflow/host-event-recovery.tsv`
- `workflow/geographic-ancestral-summary.tsv`
- `workflow/geographic-state-probabilities.tsv`
- `workflow/geographic-transition-summary.tsv`
- `workflow/geographic-node-recovery.tsv`
- `workflow/geographic-event-recovery.tsv`
- `workflow/recovery-threshold-evaluation.tsv`
- `overview.md`

The packaged truth contract is explicit rather than inferred from one recovery
score:

- Brownian and OU parameter recovery is measured against stored generating
  values
- discrete, host, and geographic internal-node recovery is measured against
  stored node truths
- host-switch and geographic-transition recovery is measured branch by branch
  against stored simulated events
- recovery pass and fail thresholds are declared in
  `dataset/recovery-thresholds.tsv` and evaluated in
  `workflow/recovery-threshold-evaluation.tsv`
- `workflow/discrete-ancestral-summary.tsv`
- `workflow/discrete-node-recovery.tsv`
- `overview.md`

The packaged simulation contract is explicit:

- tree model `birth-death`
- alignment model `jukes-cantor-like`
- continuous trait model `brownian-motion`
- discrete trait model `symmetric-discrete`
- recovery tree method `neighbor-joining`
- recovery distance model `p-distance`

This command does not require external executables because both the truth
surface and the governed recovery checks run entirely inside the owned runtime.

`demo continuous-mode-recovery-panel` is the governed packaged continuous
trait-model recovery surface. It materializes the shipped deterministic
simulation panel into one output directory and reruns the owned Brownian,
Ornstein-Uhlenbeck, and early-burst recovery workflow over one shared rooted
tree and four packaged simulation cases. Its JSON metrics report:

- `artifact_count`
- `taxon_count`
- `case_count`
- `selection_match_count`
- `parameter_pass_count`
- `parameter_row_count`
- `expected_warning_case_count`
- `expected_warning_present_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/reference-tree.nwk`
- `dataset/simulation-cases.tsv`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/recovery-summary.tsv`
- `workflow/parameter-recovery.tsv`
- `workflow/model-choice.tsv`
- `workflow/warning-review.tsv`
- `workflow/simulated-traits/*.tsv`
- `overview.md`

The packaged recovery contract is explicit:

- Brownian cases are judged on sigma-squared recovery and Brownian model choice.
- OU cases are judged on alpha, sigma-squared, optimum recovery, and OU model choice.
- Early-burst cases are judged on rate-change recovery and early-burst model choice.
- Weak OU cases are judged on warning transparency and Brownian-like model support rather than fake strong parameter certainty.

`simulate traits-early-burst` is the owned continuous-trait simulator for one
early-burst branch-rate case. It writes one tip-trait table and reports the
declared `rate_change` in JSON output so reviewers can tie downstream recovery
rows back to the generating parameter instead of inferring it indirectly.

`simulate traits-brownian` is the owned one-trait Brownian simulator. It writes
one tip-trait table, accepts one root state plus either `--sigma` or
`--sigma-squared`, and reports the resolved Brownian rate in JSON output so the
generated trait table keeps one explicit covariance-generating parameter
contract.

`simulate traits-brownian-correlated` is the owned multivariate Brownian
simulator. It writes one long-form replicate tip-trait table, accepts one
fixed-tree evolutionary covariance contract either directly through repeated
`--covariance-row` values or indirectly through repeated `--correlation-row`
values plus one `--trait-standard-deviation` per trait, and can also write one
summary ledger over root states, evolutionary covariance, tip distributions,
and tip covariances. Invalid covariance inputs fail explicitly instead of being
coerced into one fallback matrix.

`simulate tree-random` and `simulate tree-coalescent` are the owned governed
tree-simulation review surfaces for random rooted trees and coalescent trees.
They can each write one tree-set output plus one per-tree record ledger and one
envelope ledger through `--record-table-out` and `--envelope-table-out`. The
live `ape` parity lane now checks those envelope ledgers against governed
`ape::rtree` and `ape::rcoal` cases from
`shared_tree_simulation_fixture_catalog.json`, so simulation parity is tracked
as a machine-readable distribution review surface rather than one unstable
literal Newick target.

`demo rabies-cross-host-panel` is the governed packaged pathogen
host-switching surface. It materializes the shipped rabies nucleoprotein panel
into one output directory and reruns the owned host-switching workflow over
the packaged rooted tree and grouped host metadata. Its JSON metrics report:

- `artifact_count`
- `taxon_count`
- `workflow_trait`
- `observed_host_group_count`
- `analysis_constraint_mode`
- `root_host`
- `root_confidence`
- `host_switch_count`
- `certain_host_switch_count`
- `uncertain_host_switch_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/sequences.fasta`
- `dataset/tree.nwk`
- `dataset/hosts.csv`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/host-switch-summary.tsv`
- `workflow/host-state-nodes.tsv`
- `workflow/host-switch-branches.tsv`
- `workflow/host-switch-counts.tsv`
- `workflow/host-switch-fits.tsv`
- `workflow/host-switch-unsupported.tsv`
- `workflow/host-switch-exclusions.tsv`
- `overview.md`

The packaged dataset carries both exact `host_species` labels and one grouped
workflow trait:

- workflow trait `host_group`
- discrete ancestral model `ard`

This command does not require external inference executables because the
rooted rabies tree is packaged directly with the dataset.

`demo rabies-geographic-transition-panel` is the governed packaged pathogen
geography surface. It materializes the shipped rabies nucleoprotein panel into
one output directory and reruns the owned biogeography workflow over the
packaged rooted tree and grouped region metadata. Its JSON metrics report:

- `artifact_count`
- `taxon_count`
- `workflow_trait`
- `observed_region_group_count`
- `root_region`
- `root_region_probability`
- `changed_branch_count`
- `strongly_supported_transition_count`
- `migration_event_count`
- `strongly_supported_migration_event_count`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/sequences.fasta`
- `dataset/tree.nwk`
- `dataset/regions.csv`
- `dataset/expected/*`
- `workflow/workflow-summary.tsv`
- `workflow/geographic-state-summary.tsv`
- `workflow/geographic-region-probabilities.tsv`
- `workflow/geographic-transition-rates.tsv`
- `workflow/geographic-transition-events.tsv`
- `workflow/geographic-state-exclusions.tsv`
- `workflow/geographic-migration-summary.tsv`
- `workflow/geographic-migration-events.tsv`
- `workflow/geographic-migration-exclusions.tsv`
- `overview.md`

The packaged dataset carries both raw `country` provenance and one grouped
workflow trait:

- workflow trait `region_group`
- discrete ancestral model `ard`

This command does not require external inference executables because the
rooted rabies tree is packaged directly with the dataset.

`demo rabies-method-sensitivity-panel` is the governed packaged method-sensitivity
surface for the compact rabies nucleoprotein panel. It materializes the
packaged FASTA, metadata, and declared workflow-config matrix into one output
directory, reruns four alignment-and-trimming variants, compares IQ-TREE
against FastTree on each trimmed alignment, roots both engine trees on the
packaged outgroup, and writes one reviewer-facing bundle that separates rooted
preprocessing stability from unrooted engine-sensitive clade differences. Its
JSON metrics report:

- `artifact_count`
- `taxon_count`
- `variant_count`
- `parallel_workers`
- `execution_mode`
- `stable_clade_count`
- `changed_clade_count`
- `preprocessing_change_pair_count`
- `rooted_engine_change_variant_count`
- `serious_conflict_variant_count`
- `report_linked_artifact_count`
- `report_html_size_bytes`
- `report_linked_artifact_bytes`
- `report_total_output_bytes`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/workflow-config.json`
- `dataset/sequences.fasta`
- `dataset/metadata.csv`
- `dataset/expected/**`
- `workflow/workflow-summary.tsv`
- `workflow/variant-summary.tsv`
- `workflow/parallel-execution-summary.tsv`
- `workflow/rabies-method-sensitivity-panel.run.json`
- `workflow/preprocessing-rooted-comparisons.tsv`
- `workflow/stable-clades.tsv`
- `workflow/changed-clades.tsv`
- `workflow/method-conclusion-summary.tsv`
- `workflow/workflow-config.resolved.json`
- `workflow/rabies-method-sensitivity.manifest.json`
- `workflow/report-artifacts/rabies-method-sensitivity-report.manifest.json`
- `workflow/rabies-method-sensitivity-report.html`
- `workflow/parallel-logs/<variant-id>.log`
- `workflow/variants/<variant-id>/*.aln`
- `workflow/variants/<variant-id>/*.trimmed.aln`
- `workflow/variants/<variant-id>/fasttree.nwk`
- `workflow/variants/<variant-id>/iqtree-support.nwk`
- `workflow/variants/<variant-id>/rooted-fasttree.nwk`
- `workflow/variants/<variant-id>/rooted-iqtree-support.nwk`
- `workflow/variants/<variant-id>/rooting-summary.tsv`
- `workflow/variants/<variant-id>/rooted-engine-comparison.tsv`
- `workflow/variants/<variant-id>/unrooted-comparison.tsv`
- `workflow/variants/<variant-id>/unrooted-shared-clades.tsv`
- `workflow/variants/<variant-id>/unrooted-conflicting-clades.tsv`
- `workflow/variants/<variant-id>/unrooted-support-weighted-conflicts.tsv`
- `workflow/variants/<variant-id>/unrooted-conclusions.tsv`
- `workflow/variants/<variant-id>/unrooted-stability-summary.tsv`
- `overview.md`

The governed workflow now rejects concurrent reuse of the same `--out`
directory while one run is still active. Parallel execution stays safe because
each declared variant uses its own isolated output root inside that workflow
directory, and the raw `workflow/rabies-method-sensitivity-panel.run.json`
execution record keeps the worker count, execution mode, successful variants,
failed variants, and per-variant task logs auditable even when one isolated
variant task fails.

The packaged workflow matrix currently declares four durable variants:

- `auto-gap-threshold`
- `ginsi-gap-threshold`
- `auto-gappyout`
- `ginsi-gappyout`

The outgroup is fixed to `bat_chile_rv108`, and the workflow treats that
rooting choice as explicit evidence rather than an implicit side effect of the
engine comparison. This command depends on external `mafft`, `trimal`,
`iqtree2`, and `FastTree` executables because it reruns the real
sequence-to-alignment-to-tree path for each declared method combination.

The HTML report is intentionally summary-first. Large ledgers remain in the
linked TSV and JSON artifacts, while the report manifest records their
relative paths, checksums, and byte counts.

`tree-set report` now follows the same scaling contract. The HTML keeps
top-level uncertainty summaries in-page, writes large tables to a sibling
`<report>.artifacts/` directory, links those artifacts explicitly, and reports
`linked_artifact_count`, `html_size_bytes`, `linked_artifact_bytes`, and
`total_output_bytes` in JSON mode. Tree sets with `1,000+` trees switch to
`scaled-summary` mode and replace the most expensive supplemental sensitivity
passes with linked note artifacts instead of expanding them inline.

`demo rabies-cross-host-geography-panel` is the governed packaged integrated
pathogen workflow surface. It materializes the shipped rabies nucleoprotein
panel into one output directory and reruns the full owned sequence-to-tree,
host-switching, and biogeography workflow from raw sequences plus one combined
metadata table. Its JSON metrics report:

- `artifact_count`
- `sequence_count`
- `config_path`
- `biological_question`
- `short_answer`
- `host_trait`
- `geography_trait`
- `selected_model`
- `aligned_quality_score`
- `trimmed_quality_score`
- `minimum_support`
- `maximum_support`
- `root_host`
- `root_region`
- `host_switch_count`
- `migration_event_count`
- `clade_row_count`
- `bootstrap_tree_count`
- `timeout_seconds`
- `max_bootstrap_tree_count`
- `max_report_table_rows`
- `budget_warning_count`
- `comparative_formula`
- `comparative_selected_model`
- `reference_output_count`

The command writes:

- `dataset/README.md`
- `dataset/workflow-config.json`
- `dataset/sequences.fasta`
- `dataset/metadata.csv`
- `dataset/region-centroids.csv`
- `dataset/source-accessions.tsv`
- `dataset/expected/**`
- `workflow/workflow-summary.tsv`
- `workflow/input-validation.tsv`
- `workflow/alignment-quality.tsv`
- `workflow/alignment-sequence-ranking.tsv`
- `workflow/rabies-cross-host-geography-panel.aln`
- `workflow/rabies-cross-host-geography-panel.trimmed.aln`
- `workflow/rabies-cross-host-geography-panel.rooted.tree`
- `workflow/rabies-cross-host-geography-panel.rooting.tsv`
- `workflow/rabies-cross-host-geography-panel.model.tsv`
- `workflow/rabies-cross-host-geography-panel.support.tsv`
- `workflow/clade-table.tsv`
- `workflow/bootstrap-review/bootstrap-review.summary.tsv`
- `workflow/bootstrap-review/bootstrap-review.consensus.nwk`
- `workflow/bootstrap-review/bootstrap-review.clade-frequencies.tsv`
- `workflow/bootstrap-review/bootstrap-review.unstable-branches.tsv`
- `workflow/bootstrap-review/bootstrap-review.unstable-clades.tsv`
- `workflow/bootstrap-review/bootstrap-review.distance-matrix.tsv`
- `workflow/bootstrap-review/bootstrap-review.topology-clusters.tsv`
- `workflow/host-switch-summary.tsv`
- `workflow/host-state-nodes.tsv`
- `workflow/host-switch-branches.tsv`
- `workflow/host-switch-counts.tsv`
- `workflow/host-switch-fits.tsv`
- `workflow/host-switch-unsupported.tsv`
- `workflow/host-switch-exclusions.tsv`
- `workflow/biogeography/biogeography-report.html`
- `workflow/biogeography/ancestral-region-tree.svg`
- `workflow/biogeography/geographic-region-map.html`
- `workflow/biogeography/summary.tsv`
- `workflow/biogeography/region-counts.tsv`
- `workflow/biogeography/ancestral-regions.tsv`
- `workflow/biogeography/transition-matrix.tsv`

The packaged `dataset/workflow-config.json` is also the governed resource
budget surface for this workflow. In addition to the biological settings it
accepts:

- `iqtree_threads`
- `timeout_seconds`
- `max_bootstrap_tree_count`
- `max_report_table_rows`
- `memory_warning_threshold_bytes`

When the bootstrap review or integrated HTML report exceeds one of those
budgets, the runtime now either fails with a structured workflow-budget error
or records one explicit warning in the workflow summary and JSON metrics.
- `workflow/biogeography/event-table.tsv`
- `workflow/biogeography/map-markers.tsv`
- `workflow/biogeography/map-lines.tsv`
- `workflow/biogeography/exclusions.tsv`
- `workflow/comparative-traits.tsv`
- `workflow/comparative-tree.nwk`
- `workflow/comparative-tree-adjustments.tsv`
- `workflow/comparative/comparative-report.html`
- `workflow/comparative/comparative-summary.tsv`
- `workflow/comparative/coefficient-table.tsv`
- `workflow/comparative/residual-summary.tsv`
- `workflow/comparative/signal-summary.tsv`
- `workflow/comparative/model-comparison.tsv`
- `workflow/comparative/interpretation-table.tsv`
- `workflow/comparative/audit-table.tsv`
- `workflow/comparative/contrast-table.tsv`
- `workflow/comparative/model-matrix.tsv`
- `workflow/comparative/categorical-contrasts.tsv`
- `workflow/comparative/lambda-profile.tsv`
- `workflow/comparative/comparative.manifest.json`
- `workflow/rabies-cross-host-geography-report.html`
- `workflow/rabies-cross-host-geography.manifest.json`
- `overview.md`
- `rabies-cross-host-geography-overview.html`
- `rabies-cross-host-geography-package.manifest.json`

The packaged dataset carries grouped workflow traits for both downstream
biological surfaces:

- host workflow trait `host_group`
- geography workflow trait `region_group`
- comparative formula `region_longitude ~ host_group`
- discrete ancestral model `ard` for both state-evolution analyses
- explicit outgroup rooting on `bat_chile_rv108`

The package root is intentionally part of the public review contract.
`dataset/source-accessions.tsv` keeps accession provenance machine-readable,
the overview HTML states one biological question plus one short answer in
plain language, and the package manifest records that same question and answer
alongside config provenance, output checksums, and high-level workflow
metrics.

This command does require external `mafft`, `trimal`, and `iqtree2`
executables because it reruns the full raw-sequence inference path instead of
starting from a packaged rooted tree. Use `--config dataset/workflow-config.json`
against a packaged export when the exact shipped workflow settings should drive
the rerun.

`ancestral continuous` is the governed reconstruction surface for one numeric
trait on one rooted dichotomous tree. It estimates internal-node values under
the selected continuous model, reports 95% uncertainty intervals, and prunes
tips with missing or non-numeric trait values instead of hiding them. Its JSON
metrics report:
- `taxon_count`
- `estimate_count`
- `internal_node_count`
- `excluded_taxon_count`
- `unstable_node_count`
- `model`
- `estimator`
- `tree_is_ultrametric`
- `covariance_near_singular`
- `covariance_condition_number`
- `log_likelihood`
- `residual_sigma_squared`
- `optimizer_name`
- `optimizer_converged`
- `optimizer_iteration_count`
- `optimizer_function_evaluation_count`

The command supports `brownian` and `ou` reconstruction modes. The Brownian
path is aligned to the governed `ape::ace(method='pic')` reference surface with
explicit bounded tolerance rather than an undocumented local convention. The
optional `--estimator` flag makes the estimator surface explicit: `ace-pic`
preserves the governed `ape::ace(type='continuous', method='pic')` lane,
`anc-ml` switches to the governed live `phytools::anc.ML` lane with
Brownian log-likelihood, fitted sigma-squared, and optimizer diagnostics,
`fast-anc` switches to the governed live `phytools::fastAnc` lane, and
`generalized-least-squares` is reserved for the `ou` model.

When `--table-out` is supplied, `ancestral continuous` writes one flat node
ledger as CSV or TSV with both tips and internal nodes. When `--summary-out` is
supplied, it also writes one summary ledger. The summary row preserves:
- `trait`
- `taxon_column`
- `model`
- `estimator`
- `alpha`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `missing_tip_taxon_count`
- `non_numeric_tip_taxon_count`
- `internal_node_count`
- `unstable_node_count`
- `tree_is_ultrametric`
- `covariance_near_singular`
- `covariance_condition_number`
- `log_likelihood`
- `residual_sigma_squared`
- `optimizer_name`
- `optimizer_converged`
- `optimizer_iteration_count`
- `optimizer_function_evaluation_count`
- `log_likelihood`
- `residual_sigma_squared`
- `weak_support_node_count`
- `root_node`
- `root_estimate`
- `root_standard_error`
- `root_lower_95_interval`
- `root_upper_95_interval`
- `warning_count`

When `--uncertainty-out` is supplied, the command writes one internal-node
uncertainty ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `estimate`
- `standard_error`
- `lower_95_interval`
- `upper_95_interval`
- `uncertainty_width`
- `confidence`
- `interpretation`
- `unstable`

When `--exclusions-out` is supplied, the command writes one explicit excluded
tip ledger. Each row preserves:
- `taxon`
- `reason`

`ancestral discrete` is the governed reconstruction surface for one categorical
trait on one rooted dichotomous tree. It supports Fitch parsimony for a fast
set-based reconstruction and supports `equal-rates`, `symmetric`, and
`all-rates-different` likelihood models for Mk-style marginal ancestral
probabilities. Its JSON metrics report:
- `taxon_count`
- `estimate_count`
- `internal_node_count`
- `ambiguous_internal_node_count`
- `excluded_taxon_count`
- `state_count`
- `minimal_change_count`
- `parsimonious_root_state_count`
- `unstable_node_count`
- `comparison_node_count`
- `comparison_differing_node_count`
- `model`
- `root_prior_mode`
- `fixed_root_state`
- `log_likelihood`
- `parameter_count`
- `aic`
- `transition_rate_count`
- `phytools_rerooting_method_comparable`

For the likelihood models, the owned runtime fits an explicit Mk rate matrix,
reports node-level marginal probabilities, and can export one fitted directed
transition-rate ledger. The governed parity surface is checked against
`ape::ace` on ER, SYM, and ARD reference cases with explicit bounded
tolerances instead of a vague compatibility claim. Within that lane, the live
`ape::ace` discrete surface is governed explicitly for ER, SYM, and ARD,
while root-prior controls remain an owned Bijux policy surface because
`ape::ace` does not expose the same runtime root-prior interface. The owned
fit surface also warns when multi-parameter likelihood fits hit optimizer
bounds so weakly identified ARD and SYM reconstructions are reviewable instead
of looking falsely settled. The same report now also states whether the
requested run is comparable to live `phytools::rerootingMethod`: ER and SYM
with `--root-prior-mode equal` are governed rerooting-parity surfaces, while
Fitch, ordered-state, ARD, empirical-root-prior, and fixed-root-prior runs are
flagged explicitly as non-comparable.

For Fitch, the owned runtime now also reports the exact minimum parsimony
change count for the analyzed tree and the number of parsimonious root states.
That keeps the fast path reviewable instead of reducing it to only one chosen
state label per internal node.

When `--table-out` is supplied, `ancestral discrete` writes one flat node
ledger as CSV or TSV with both tips and internal nodes. When `--summary-out` is
supplied, it also writes one summary ledger. The summary row preserves:
- `trait`
- `taxon_column`
- `model`
- `state_ordering`
- `root_prior_mode`
- `fixed_root_state`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `internal_node_count`
- `ambiguous_internal_node_count`
- `unstable_node_count`
- `weak_support_node_count`
- `observed_state_count`
- `sparse_state_count`
- `minimal_change_count`
- `parsimonious_root_state_count`
- `root_node`
- `root_most_likely_state`
- `root_confidence`
- `phytools_rerooting_method_comparable`
- `log_likelihood`
- `parameter_count`
- `aic`
- `warning_count`

When `--probabilities-out` is supplied, the command writes one internal-node
marginal-probability ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `most_likely_state`
- `state_set`
- `state_probabilities`
- `confidence`
- `ambiguous`
- `unstable`
- `interpretation`

When `--transitions-out` is supplied, the command writes one directed fitted
transition-rate ledger for likelihood models. Each row preserves:
- `source_state`
- `target_state`
- `transition_allowed`
- `step_distance`
- `rate`

When `--comparison-out` is supplied, the command writes one direct node-wise
comparison ledger between the requested `--model` and `--compare-model`. Each
row preserves:
- `node`
- `descendant_taxa`
- `left_model`
- `right_model`
- `left_state`
- `right_state`
- `left_state_set`
- `right_state_set`
- `left_confidence`
- `right_confidence`
- `left_ambiguous`
- `right_ambiguous`
- `differs`
- `ambiguity_changed`

When `--exclusions-out` is supplied, the command writes one explicit excluded
tip ledger. Each row preserves:
- `taxon`
- `reason`

`ancestral discrete-reference` reruns the governed discrete ancestral
reference suite before any user dataset is interpreted. It validates
equal-rates, symmetric, and all-rates-different likelihood reconstructions
against checked-in `ape::ace` probability fixtures and then reruns the owned
root-prior, ambiguity, ordered-state, and irreversible-transition policy
surfaces on known examples. Its JSON metrics report:
- `case_count`
- `external_case_count`
- `all_passed`

`ancestral ordered-discrete` is the governed ordered-state comparison surface
for one discrete likelihood ancestral reconstruction. It fits the requested
likelihood model twice on the same tree: once with the supplied ordered state
vocabulary and once with the unrestricted unordered baseline. Its JSON metrics
report:
- `model`
- `ordered_state_count`
- `fit_count`
- `differing_node_count`
- `ambiguity_change_count`
- `restricted_transition_count`
- `preferred_ordering`

The command supports `equal-rates`, `symmetric`, and
`all-rates-different`. It requires `--ordered-states`, which defines the
durable ordered state vocabulary used to restrict transitions to adjacent
states only.

When `--summary-out` is supplied, `ancestral ordered-discrete` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `analyzed_taxon_count`
- `state_count`
- `ordered_log_likelihood`
- `unordered_log_likelihood`
- `ordered_parameter_count`
- `unordered_parameter_count`
- `ordered_aic`
- `unordered_aic`
- `delta_aic`
- `preferred_ordering`
- `differing_node_count`
- `ambiguity_change_count`
- `restricted_transition_count`
- `warning_count`

When `--fits-out` is supplied, the command writes one two-row fit ledger. Each
row preserves:
- `ordering_mode`
- `model`
- `state_ordering`
- `ordered_states`
- `analyzed_taxon_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `root_most_likely_state`
- `root_confidence`

When `--nodes-out` is supplied, the command writes one node-wise comparison
ledger. Each row preserves:
- `node`
- `descendant_taxa`
- `ordered_state`
- `unordered_state`
- `ordered_confidence`
- `unordered_confidence`
- `confidence_delta`
- `differs`
- `ambiguity_changed`

When `--transitions-out` is supplied, the command writes one directed
transition ledger. Each row preserves:
- `source_state`
- `target_state`
- `step_distance`
- `ordered_transition_allowed`
- `unordered_transition_allowed`
- `ordered_rate`
- `unordered_rate`

`ancestral irreversible-discrete` is the governed irreversible-state review
surface for one discrete likelihood ancestral reconstruction. It fits the
requested likelihood model twice on the same tree: once under an explicit
directed allowed-transition graph and once under the unrestricted baseline.
Its JSON metrics report:
- `model`
- `allowed_transition_count`
- `fit_count`
- `differing_node_count`
- `ambiguity_change_count`
- `forbidden_transition_count`
- `preferred_constraint`

The command supports `equal-rates`, `symmetric`, and
`all-rates-different`. It requires `--allowed-transitions`, which accepts a
comma-delimited directed graph such as `present->absent` or
`north->south,south->island`. Under `symmetric`, every allowed edge must be
bidirectional because the fitted rates are shared across both directions.

When `--summary-out` is supplied, `ancestral irreversible-discrete` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `analyzed_taxon_count`
- `constrained_log_likelihood`
- `unconstrained_log_likelihood`
- `likelihood_difference`
- `constrained_parameter_count`
- `unconstrained_parameter_count`
- `constrained_aic`
- `unconstrained_aic`
- `delta_aic`
- `preferred_constraint`
- `differing_node_count`
- `ambiguity_change_count`
- `forbidden_transition_count`
- `warning_count`

When `--fits-out` is supplied, the command writes one two-row fit ledger. Each
row preserves:
- `constraint_mode`
- `model`
- `analyzed_taxon_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `root_most_likely_state`
- `root_confidence`

When `--nodes-out` is supplied, the command writes one node-wise comparison
ledger. Each row preserves:
- `node`
- `descendant_taxa`
- `constrained_state`
- `unconstrained_state`
- `constrained_confidence`
- `unconstrained_confidence`
- `confidence_delta`
- `differs`
- `ambiguity_changed`

When `--transitions-out` is supplied, the command writes one directed
transition ledger. Each row preserves:
- `source_state`
- `target_state`
- `constrained_transition_allowed`
- `unconstrained_transition_allowed`
- `constrained_rate`
- `unconstrained_rate`

`biogeography model` is the governed ancestral-region review surface for one
taxon-region table on one rooted tree. It accepts ER, SYM, and ARD model
aliases and reuses the owned geographic-state engine to produce explicit
internal-node region probabilities, pairwise transition-rate rows, branchwise
event rows, and excluded-taxon rows. Its JSON metrics report:
- `model`
- `observed_region_count`
- `internal_node_count`
- `transition_rate_row_count`
- `changed_branch_count`
- `strongly_supported_transition_count`
- `excluded_taxon_count`

The command supports `er`, `sym`, and `ard`. `--allowed-regions` is optional
and defines an explicit region vocabulary when the metadata should be
restricted to named states only.

When `--summary-out` is supplied, `biogeography model` writes one overall
summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `likelihood_method`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observed_region_count`
- `internal_node_count`
- `ambiguous_internal_node_count`
- `changed_branch_count`
- `strongly_supported_transition_count`
- `transition_rate_row_count`
- `root_region`
- `root_region_probability`
- `warning_count`

When `--nodes-out` is supplied, the command writes one internal-node region
probability ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `most_likely_region`
- `region_probabilities`
- `confidence`
- `ambiguous`
- `is_root`

When `--rates-out` is supplied, the command writes one pairwise transition-rate
ledger. Each row preserves:
- `source_region`
- `target_region`
- `rate`
- `lower_95_interval`
- `upper_95_interval`
- `effective_transition_count`

When `--events-out` is supplied, the command writes one branchwise geographic
event ledger. Each row preserves:
- `parent_node`
- `child_node`
- `source_region`
- `target_region`
- `changed`
- `support`
- `strongly_supported`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_region`
- `normalized_region`
- `reason`
- `note`

`biogeography constrained` is the governed constrained-versus-unconstrained
geographic review surface for one taxon-region table on one rooted tree. It
accepts one explicit region adjacency matrix and fits one constrained and one
unconstrained likelihood geography model on the same analyzed region set. Its
JSON metrics report:
- `model`
- `allowed_transition_count`
- `forbidden_transition_count`
- `unsupported_transition_claim_count`
- `preferred_constraint`
- `excluded_taxon_count`

When `--summary-out` is supplied, `biogeography constrained` writes one overall
summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observed_region_count`
- `allowed_transition_count`
- `forbidden_transition_count`
- `constrained_log_likelihood`
- `unconstrained_log_likelihood`
- `likelihood_difference`
- `constrained_parameter_count`
- `unconstrained_parameter_count`
- `constrained_aic`
- `unconstrained_aic`
- `delta_aic`
- `preferred_constraint`
- `unsupported_transition_claim_count`
- `warning_count`

When `--fits-out` is supplied, the command writes one fit-comparison ledger.
Each row preserves:
- `constraint_mode`
- `model`
- `analyzed_taxon_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `root_region`
- `root_confidence`

When `--transitions-out` is supplied, the command writes one directed
transition-comparison ledger. Each row preserves:
- `source_region`
- `target_region`
- `transition_allowed`
- `unconstrained_rate`
- `constrained_rate`
- `rate_delta`

When `--unsupported-out` is supplied, the command writes one forbidden-claim
ledger. Each row preserves:
- `parent_node`
- `child_node`
- `descendant_taxa`
- `unconstrained_source_region`
- `unconstrained_target_region`
- `unconstrained_support`
- `constrained_source_region`
- `constrained_target_region`
- `constrained_support`
- `claim_resolved`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_region`
- `normalized_region`
- `reason`
- `note`

This constrained geography surface is intentionally explicit about scope. It is
an adjacency-constrained likelihood review over the owned discrete ancestral
runtime, not a full historical biogeography process model. It should be used
to test whether unconstrained geographic transition claims are compatible with
an explicit region-connectivity contract and how much fit is lost when those
forbidden transitions are removed.

`biogeography events` is the governed geographic movement-event review surface.
On one rooted tree it extracts only changed source-target branches from the
owned ancestral geography reconstruction and reports:
- `report_mode`
- `model`
- `event_count`
- `strongly_supported_event_count`
- `mean_event_support`
- `excluded_taxon_count`

When `--tree-set` is supplied, the same command switches into retained-tree
review mode over a posterior or bootstrap tree set and reports:
- `report_mode`
- `model`
- `kept_tree_count`
- `event_row_count`
- `event_summary_count`
- `topology_sensitive_event_count`
- `excluded_taxon_count`
- `warning_count`

When `--summary-out` is supplied on one tree, `biogeography events` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `likelihood_method`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `tree_depth`
- `event_count`
- `strongly_supported_event_count`
- `mean_event_support`
- `earliest_midpoint_depth`
- `latest_midpoint_depth`
- `warning_count`

When `--events-out` is supplied on one tree, the command writes one
branchwise event ledger. Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `branch_length`
- `parent_depth`
- `child_depth`
- `midpoint_depth`
- `source_region`
- `target_region`
- `support`
- `strongly_supported`
- `confidence_class`

When `--summary-out` is supplied with `--tree-set`, the command writes one
tree-set summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `total_tree_count`
- `burnin_tree_count`
- `kept_tree_count`
- `shared_tree_taxon_count`
- `analysis_taxon_count`
- `rooted_topology_count`
- `unrooted_topology_count`
- `event_row_count`
- `event_summary_count`
- `topology_sensitive_event_count`
- `low_support_event_count`
- `excluded_taxon_count`
- `warning_count`

When `--trees-out` is supplied with `--tree-set`, the command writes one
retained-tree ledger. Each row preserves:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `event_count`
- `strongly_supported_event_count`

When `--events-out` is supplied with `--tree-set`, the command writes one
per-tree event ledger. Each row preserves:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `branch_length`
- `parent_depth`
- `child_depth`
- `midpoint_depth`
- `source_region`
- `target_region`
- `support`
- `strongly_supported`
- `confidence_class`

When `--event-summaries-out` is supplied with `--tree-set`, the command writes
one comparable-event summary ledger across retained trees. Each row preserves:
- `branch_id`
- `child_descendant_taxa`
- `source_region`
- `target_region`
- `tree_presence_count`
- `tree_presence_fraction`
- `strongly_supported_tree_count`
- `strongly_supported_tree_fraction`
- `mean_support`
- `lower_95_midpoint_depth`
- `upper_95_midpoint_depth`
- `minimum_parent_depth`
- `maximum_child_depth`
- `stability_class`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_region`
- `normalized_region`
- `reason`
- `note`

This movement-event surface is intentionally explicit about temporal precision.
Its `midpoint_depth` value is a deterministic branch-midpoint estimate for
review, not a claim of exact stochastic event time inference.

`biogeography report` is the governed full geographic-evolution package surface
for one rooted tree, one taxon-region table, and one explicit centroid table.
It composes the owned ancestral-region model, the owned migration-event
surface, the owned ancestral-region tree renderer, and the owned region-map
surface into one durable handoff. Its JSON metrics report:
- `report_kind`
- `model`
- `output_dir`
- `artifact_count`
- `observed_region_count`
- `transition_rate_row_count`
- `event_count`
- `visible_map_line_count`

The command requires one centroid table positional input after the trait
arguments. `--region-column`, `--latitude-column`, and `--longitude-column`
resolve that centroid table. It supports `er`, `sym`, and `ard`.

When `--out-dir` is supplied, `biogeography report` writes this fixed package:
- `biogeography-report.html`
- `ancestral-region-tree.svg`
- `geographic-region-map.html`
- `summary.tsv`
- `region-counts.tsv`
- `ancestral-regions.tsv`
- `transition-matrix.tsv`
- `event-table.tsv`
- `map-markers.tsv`
- `map-lines.tsv`
- `exclusions.tsv`
- `biogeography-report.manifest.json`

`summary.tsv` keeps the owned ancestral-region model summary row.
`region-counts.tsv` keeps one observed region count row per analyzed region at
the tip surface after tree overlap and exclusion auditing.
`ancestral-regions.tsv` keeps one internal-node ancestral region probability
row. `transition-matrix.tsv` keeps the directed geographic transition-rate
ledger. `event-table.tsv` keeps the branchwise migration-event ledger.

`ancestral-region-tree.svg` is the governed tree figure with tip and internal
region calls rendered directly on the tree. `geographic-region-map.html` is the
self-contained region-transition map companion artifact. `map-markers.tsv` and
`map-lines.tsv` expose the marker and line ledgers that drive that map. The
combined `exclusions.tsv` ledger keeps both state-model and map-placement
exclusions in one durable review surface instead of forcing reviewers to merge
separate omission tables by hand.

`biogeography time-stratified` is the governed interval-specific geographic
transition review surface for one taxon-region table on one rooted tree with
positive branch lengths. It accepts ER, SYM, and ARD model aliases plus one or
more explicit `--time-bin LABEL:START:END` definitions. The workflow reuses the
owned ancestral-region reconstruction, allocates branch exposure and inferred
branch changes across the requested root-depth intervals, and reports:
- `model`
- `time_bin_count`
- `matrix_row_count`
- `changed_branch_count`
- `allocated_transition_weight_total`
- `excluded_taxon_count`
- `warning_count`

When `--summary-out` is supplied, `biogeography time-stratified` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `likelihood_method`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `tree_depth`
- `time_bin_count`
- `matrix_row_count`
- `changed_branch_count`
- `allocated_transition_weight_total`
- `warning_count`

When `--matrix-out` is supplied, the command writes one interval-specific
transition matrix ledger. Each row preserves:
- `time_bin_label`
- `start_depth`
- `end_depth`
- `source_region`
- `target_region`
- `source_exposure_length`
- `allocated_transition_weight`
- `time_stratified_rate`
- `global_rate`

When `--branches-out` is supplied, the command writes one branch-interval
allocation ledger. Each row preserves:
- `time_bin_label`
- `start_depth`
- `end_depth`
- `parent_node`
- `child_node`
- `parent_depth`
- `child_depth`
- `source_region`
- `target_region`
- `changed`
- `overlap_length`
- `allocated_transition_weight`
- `support`
- `strongly_supported`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_region`
- `normalized_region`
- `reason`
- `note`

This interval-specific surface is intentionally explicit about scope. It is a
deterministic branch-allocation review over the owned geographic-state
reconstruction, not a full time-inhomogeneous stochastic biogeographic process
fit. If the requested intervals do not cover the full tree depth, the command
reports that omission as a warning instead of silently claiming full temporal
coverage.

`biogeography chronology` is the governed dated-tree geographic chronology
review surface for one taxon-region table on one rooted ultrametric time tree.
It accepts ER, SYM, and ARD model aliases, verifies that the tree is
time-scaled, extracts node ages, maps inferred geographic transitions to
automatic equal-width age bins, and reports:
- `model`
- `tree_is_time_scaled`
- `root_age`
- `event_count`
- `time_bin_count`
- `high_uncertainty_bin_count`
- `excluded_taxon_count`

When `--summary-out` is supplied, `biogeography chronology` writes one overall
summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `likelihood_method`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `rooted`
- `branch_length_status`
- `tree_is_time_scaled`
- `tip_count`
- `node_age_row_count`
- `root_age`
- `event_count`
- `time_bin_count`
- `empty_time_bin_count`
- `high_uncertainty_bin_count`
- `warning_count`

When `--nodes-out` is supplied, the command writes one dated node ledger. Each
row preserves:
- `node`
- `node_name`
- `is_tip`
- `descendant_taxa`
- `branch_length`
- `depth_from_root`
- `age_before_present`
- `most_likely_region`
- `region_confidence`
- `ambiguous`
- `is_root`

When `--events-out` is supplied, the command writes one dated event ledger.
Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `source_region`
- `target_region`
- `branch_length`
- `parent_depth`
- `child_depth`
- `parent_age_before_present`
- `child_age_before_present`
- `midpoint_age_before_present`
- `time_bin_label`
- `support`
- `strongly_supported`
- `confidence_class`

When `--bins-out` is supplied, the command writes one time-bin chronology
ledger. Each row preserves:
- `time_bin_label`
- `start_age_before_present`
- `end_age_before_present`
- `event_count`
- `strongly_supported_event_count`
- `low_support_event_count`
- `support_weight_total`
- `mean_support`
- `support_uncertainty`
- `earliest_event_age_before_present`
- `latest_event_age_before_present`
- `dominant_transition`
- `transition_diversity`
- `uncertainty_class`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_region`
- `normalized_region`
- `reason`
- `note`

This dated-tree surface is intentionally explicit about scope. Its equal-width
age bins are reviewer-facing chronology bins over one owned reconstruction, not
a claim that the command fitted a fully time-varying stochastic biogeographic
process.

`biogeography sampling-bias` is the governed weighted-versus-unweighted
geographic review surface for one taxon-region table on one rooted tree. It
accepts ER, SYM, and ARD model aliases plus an optional explicit region-weight
table. When `--weights` is absent, the command applies automatic
inverse-frequency region weights. Its JSON metrics report:
- `model`
- `weighting_mode`
- `region_dominated`
- `dominant_region`
- `dominant_region_fraction`
- `root_region_changed`
- `changed_internal_node_count`
- `changed_transition_count`
- `excluded_taxon_count`

When `--summary-out` is supplied, `biogeography sampling-bias` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `weighting_mode`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observed_region_count`
- `region_dominated`
- `dominant_region`
- `dominant_region_fraction`
- `weighted_region_dominated`
- `weighted_dominant_region`
- `weighted_dominant_region_fraction`
- `root_region_unweighted`
- `root_region_weighted`
- `root_region_changed`
- `compared_internal_node_count`
- `changed_internal_node_count`
- `compared_transition_count`
- `changed_transition_count`
- `warning_count`

When `--regions-out` is supplied, the command writes one region-count and
weight ledger. Each row preserves:
- `region`
- `sample_count`
- `sample_fraction`
- `applied_weight`
- `weighted_sample_count`
- `weighted_sample_fraction`
- `dominant_unweighted`
- `dominant_weighted`

When `--nodes-out` is supplied, the command writes one weighted-versus-unweighted
internal-node ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `is_root`
- `unweighted_region`
- `weighted_region`
- `unweighted_confidence`
- `weighted_confidence`
- `confidence_delta`
- `changed`
- `unweighted_region_probabilities`
- `weighted_region_probabilities`

When `--transitions-out` is supplied, the command writes one weighted-versus-unweighted
branch transition ledger. Each row preserves:
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `unweighted_source_region`
- `unweighted_target_region`
- `weighted_source_region`
- `weighted_target_region`
- `unweighted_transition`
- `weighted_transition`
- `unweighted_changed`
- `weighted_changed`
- `changed_by_weighting`
- `unweighted_support`
- `weighted_support`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_region`
- `normalized_region`
- `reason`
- `note`

This weighted review surface is intentionally explicit about scope. It
reweights the owned deterministic geographic reconstruction so reviewers can
see how uneven sampling changes the conclusion surface. It is not presented as
one uniquely correct sampling model or as a substitute for richer generative
biogeographic process inference.

`host-association switches` is the governed host-switch review surface for one
host metadata table on one rooted parasite or pathogen tree. It reconstructs
internal host states, classifies branchwise host switches as certain or
uncertain, aggregates directed switch counts, and optionally compares a
constrained host-transition policy against the unconstrained fit. Its JSON
metrics report:
- `model`
- `analysis_constraint_mode`
- `observed_host_count`
- `host_switch_count`
- `certain_host_switch_count`
- `uncertain_host_switch_count`
- `preferred_constraint`
- `unsupported_switch_claim_count`
- `excluded_taxon_count`

The command supports `er`, `sym`, and `ard`. `--constraints` is optional and
accepts one CSV or TSV ledger with `source_host` and `target_host` columns,
plus an optional `transition_allowed` column. When present, the command fits
the same host-state model twice on the shared analyzed taxa: once
unconstrained and once under the explicit allowed-transition graph.

When `--summary-out` is supplied, `host-association switches` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `analysis_constraint_mode`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observed_host_count`
- `internal_node_count`
- `ambiguous_internal_node_count`
- `host_switch_count`
- `certain_host_switch_count`
- `uncertain_host_switch_count`
- `allowed_transition_count`
- `forbidden_transition_count`
- `constrained_log_likelihood`
- `unconstrained_log_likelihood`
- `constrained_aic`
- `unconstrained_aic`
- `preferred_constraint`
- `unsupported_switch_claim_count`
- `root_host`
- `root_confidence`
- `warning_count`

When `--nodes-out` is supplied, the command writes one internal-node host
probability ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `most_likely_host`
- `host_probabilities`
- `confidence`
- `ambiguous`
- `is_root`

When `--branches-out` is supplied, the command writes one branchwise
host-switch ledger. Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `branch_length`
- `parent_most_likely_host`
- `child_most_likely_host`
- `parent_host_set`
- `child_host_set`
- `overlapping_hosts`
- `changed`
- `transition`
- `certainty_class`
- `parent_confidence`
- `child_confidence`
- `transition_allowed`

When `--counts-out` is supplied, the command writes one directed switch-count
ledger. Each row preserves:
- `transition`
- `source_host`
- `target_host`
- `transition_allowed`
- `certain_switch_count`
- `uncertain_switch_count`
- `total_switch_count`

When `--fits-out` is supplied, the command writes one fit-comparison ledger.
Each row preserves:
- `constraint_mode`
- `model`
- `analyzed_taxon_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `root_host`
- `root_confidence`

When `--unsupported-out` is supplied, the command writes one forbidden-claim
ledger. Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `unconstrained_source_host`
- `unconstrained_target_host`
- `unconstrained_certainty_class`
- `constrained_source_host`
- `constrained_target_host`
- `constrained_certainty_class`
- `claim_resolved`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_host`
- `normalized_host`
- `reason`
- `note`

This host-association surface is intentionally explicit about scope. It is a
one-tree host-state evolution review over the owned discrete ancestral
runtime, not a full cophylogenetic reconciliation or host-parasite
transmission-history inference model.

`ecological-niche transitions` is the governed ecological niche transition
review surface for one ecological-state table on one rooted tree. It fits one
likelihood discrete transition model, reconstructs internal niche states,
counts branchwise niche changes, and ranks internal clades by concentrated
shift burden. Its JSON metrics report:
- `model`
- `observed_niche_count`
- `transition_rate_row_count`
- `changed_branch_count`
- `certain_transition_count`
- `uncertain_transition_count`
- `repeated_shift_clade_count`
- `excluded_taxon_count`

The command supports `er`, `sym`, and `ard`. It is intentionally likelihood
only so the transition-rate surface, likelihood, and AIC remain explicit
review artifacts rather than hidden assumptions.

When `--summary-out` is supplied, `ecological-niche transitions` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `internal_model`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observed_niche_count`
- `internal_node_count`
- `ambiguous_internal_node_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `transition_rate_row_count`
- `changed_branch_count`
- `certain_transition_count`
- `uncertain_transition_count`
- `strongly_supported_transition_count`
- `clade_shift_row_count`
- `repeated_shift_clade_count`
- `root_niche`
- `root_confidence`
- `warning_count`

When `--nodes-out` is supplied, the command writes one internal-node niche
probability ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `most_likely_niche`
- `niche_probabilities`
- `confidence`
- `ambiguous`
- `is_root`

When `--rates-out` is supplied, the command writes one fitted niche
transition-rate ledger. Each row preserves:
- `source_niche`
- `target_niche`
- `transition_allowed`
- `step_distance`
- `rate`

When `--branches-out` is supplied, the command writes one branchwise niche
transition ledger. Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `branch_length`
- `parent_most_likely_niche`
- `child_most_likely_niche`
- `parent_niche_set`
- `child_niche_set`
- `overlapping_niches`
- `changed`
- `transition`
- `certainty_class`
- `support`
- `strongly_supported`
- `parent_confidence`
- `child_confidence`

When `--counts-out` is supplied, the command writes one directed niche
transition-count ledger. Each row preserves:
- `transition`
- `source_niche`
- `target_niche`
- `certain_transition_count`
- `uncertain_transition_count`
- `total_transition_count`
- `strongly_supported_transition_count`

When `--clades-out` is supplied, the command writes one ranked internal-clade
shift ledger. Each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `descendant_taxon_count`
- `descendant_internal_node_count`
- `changed_branch_count`
- `certain_transition_count`
- `uncertain_transition_count`
- `strongly_supported_transition_count`
- `transition_diversity`
- `dominant_transition`
- `dominant_transition_count`
- `shift_burden_score`
- `contains_repeated_shifts`
- `rank`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_niche`
- `normalized_niche`
- `reason`
- `note`

This ecological-niche surface is intentionally explicit about scope. It is a
one-tree discrete niche-evolution review over the owned ancestral runtime, not
a full macroecological process model or a direct habitat-dependent
diversification inference surface.

`phylogeography coordinates` is the governed continuous-coordinate review
surface for one latitude/longitude table on one rooted tree. It reconstructs
ancestral coordinates under Brownian or OU continuous evolution, measures
branchwise great-circle displacement, flags jump outliers, and can render one
coordinate-space movement visualization. Its JSON metrics report:
- `model`
- `analyzed_taxon_count`
- `outlier_jump_count`
- `impossible_jump_count`
- `flagged_branch_count`
- `maximum_jump_km`
- `excluded_taxon_count`

The command supports `brownian` and `ou`. `--alpha` is accepted for the OU
path and remains explicit in the review summary. `--visualization-out` accepts
`.svg` or `.html` and produces a coordinate-space movement artifact rather
than a projected map.

When `--summary-out` is supplied, `phylogeography coordinates` writes one
overall summary ledger. The row preserves:
- `taxon_column`
- `latitude_column`
- `longitude_column`
- `model`
- `alpha`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `internal_node_count`
- `weak_support_node_count`
- `outlier_jump_count`
- `impossible_jump_count`
- `flagged_branch_count`
- `maximum_jump_km`
- `root_latitude`
- `root_longitude`
- `root_radial_standard_error_km`
- `warning_count`

When `--estimates-out` is supplied, the command writes one coordinate-estimate
ledger for tips and internal nodes. Each row preserves:
- `node`
- `node_name`
- `is_tip`
- `descendant_taxa`
- `latitude`
- `longitude`
- `latitude_standard_error`
- `longitude_standard_error`
- `radial_standard_error_km`
- `lower_95_latitude`
- `upper_95_latitude`
- `lower_95_longitude`
- `upper_95_longitude`
- `confidence`
- `unstable`
- `is_root`

When `--branches-out` is supplied, the command writes one branchwise movement
ledger. Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `branch_length`
- `parent_latitude`
- `parent_longitude`
- `child_latitude`
- `child_longitude`
- `great_circle_km`
- `branch_rate_km_per_unit`
- `support`
- `impossible_jump`
- `outlier_jump`
- `flag_codes`

When `--outliers-out` is supplied, the command writes one flagged movement
ledger. Each row preserves:
- `branch_id`
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `great_circle_km`
- `branch_rate_km_per_unit`
- `median_distance_km`
- `distance_threshold_km`
- `median_rate_km_per_unit`
- `rate_threshold_km`
- `impossible_jump`
- `outlier_jump`
- `flag_codes`

When `--exclusions-out` is supplied, the command writes one excluded-taxa
ledger. Each row preserves:
- `taxon`
- `raw_latitude`
- `raw_longitude`
- `reason`
- `note`

This phylogeography surface is intentionally explicit about scope. It is a
one-tree continuous-coordinate review over the owned continuous ancestral
runtime, not a projected map output, direct route reconstruction, or a full
spatial diffusion engine with historical cartography.

`phylogeography coordinates-map` is the governed map-rendering surface for one
continuous latitude/longitude reconstruction. It reuses the owned continuous
phylogeography review, projects tip and ancestral coordinates onto one fixed
world latitude/longitude extent, and can filter the visible branchwise
movement layer by midpoint depth. Its JSON metrics report:
- `map_mode`
- `model`
- `tip_marker_count`
- `internal_marker_count`
- `line_count`
- `visible_line_count`
- `time_filter_applied`
- `excluded_record_count`

The command supports `brownian` and `ou`. `--minimum-midpoint-depth` and
`--maximum-midpoint-depth` are optional visible-line filters over reconstructed
branch midpoints. `--html-out` writes one self-contained HTML map artifact.

When `--summary-out` is supplied, `phylogeography coordinates-map` writes one
overall map summary ledger. The row preserves:
- `mode`
- `model`
- `analyzed_taxon_count`
- `excluded_record_count`
- `tip_marker_count`
- `internal_marker_count`
- `root_marker_count`
- `line_count`
- `visible_line_count`
- `tree_depth`
- `time_filter_applied`
- `minimum_midpoint_depth_filter`
- `maximum_midpoint_depth_filter`
- `earliest_visible_midpoint_depth`
- `latest_visible_midpoint_depth`
- `warning_count`

When `--markers-out` is supplied, the command writes one geographic marker
ledger. Each row preserves:
- `marker_id`
- `label`
- `marker_kind`
- `latitude`
- `longitude`
- `state_label`
- `descendant_taxa`
- `confidence`
- `uncertainty_km`
- `is_tip`
- `is_root`
- `active_line_count`

When `--lines-out` is supplied, the command writes one geographic line ledger.
Each row preserves:
- `line_id`
- `line_kind`
- `source_label`
- `target_label`
- `source_latitude`
- `source_longitude`
- `target_latitude`
- `target_longitude`
- `child_descendant_taxa`
- `support`
- `midpoint_depth`
- `branch_length`
- `distance_km`
- `state_transition`
- `flag_codes`
- `visible`

When `--exclusions-out` is supplied, the command writes one map exclusion
ledger. Each row preserves:
- `subject_id`
- `subject_kind`
- `raw_left`
- `raw_right`
- `reason`
- `note`

`phylogeography regions-map` is the governed map-rendering surface for one
discrete ancestral geographic-region reconstruction plus one explicit region
centroid table. It reuses the owned biogeography reconstruction and movement
event surfaces, places tip and ancestral regions at their supplied centroids,
and renders directed transition lines on the same fixed world extent. Its JSON
metrics report the same map metrics as `phylogeography coordinates-map`.

The command requires `--trait` and `--centroids`. It supports `er`, `sym`,
and `ard`. `--region-column`, `--latitude-column`, and `--longitude-column`
resolve centroid-table columns. The same midpoint-depth filters control which
transition lines remain visible on the rendered map layer.

This map-rendering surface is intentionally explicit about scope. It does not
replace the underlying continuous or discrete geographic reconstruction.
Instead, it turns owned markers and owned movement or transition lines into a
reviewable HTML map without depending on a network tile service. Depth
filtering is a reviewer-facing visibility filter, not a new temporal model fit.

`ancestral confidence` is the governed ancestral-confidence review surface for
either one tree or one posterior/bootstrap tree set. It does not fit a new
ancestral model. Instead, it reuses the owned reconstruction surfaces and
turns their uncertainty into one ranked evidence ledger. Its JSON metrics
report:
- `kind`
- `source_kind`
- `model`
- `kept_tree_count`
- `confidence_row_count`
- `low_confidence_count`
- `unstable_count`
- `high_entropy_count`
- `top_uncertain_id`

The command requires `--kind continuous` or `--kind discrete`. On one tree it
reads `ancestral continuous` or `ancestral discrete`. With `--tree-set` it
reads `ancestral tree-set`, and `--burnin-fraction` applies only in that
tree-set mode.

When `--summary-out` is supplied, `ancestral confidence` writes one overall
summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `source_kind`
- `reconstruction_kind`
- `target_kind`
- `model`
- `state_ordering`
- `alpha`
- `analyzed_taxon_count`
- `kept_tree_count`
- `confidence_row_count`
- `low_confidence_count`
- `unstable_count`
- `high_entropy_count`
- `top_uncertain_id`
- `top_uncertain_label`
- `top_uncertain_score`
- `warning_count`

When `--confidence-out` is supplied, the command writes one ranked confidence
ledger. For one continuous tree, each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `estimate`
- `standard_error`
- `lower_95_interval`
- `upper_95_interval`
- `uncertainty_width`
- `relative_uncertainty`
- `confidence`
- `uncertainty_score`
- `uncertainty_rank`
- `confidence_class`
- `unstable`

For one discrete tree, each row preserves:
- `node`
- `node_name`
- `descendant_taxa`
- `most_likely_state`
- `state_set`
- `state_probabilities`
- `max_posterior_probability`
- `runner_up_probability`
- `probability_margin`
- `entropy`
- `normalized_entropy`
- `uncertainty_score`
- `uncertainty_rank`
- `confidence_class`
- `ambiguous`
- `unstable`

For a continuous tree set, each row preserves:
- `clade_id`
- `clade_taxa`
- `tree_presence_count`
- `tree_presence_fraction`
- `mean_confidence`
- `mean_standard_error`
- `empirical_interval_width`
- `normalized_empirical_interval_width`
- `unstable_tree_count`
- `unstable_tree_fraction`
- `instability_score`
- `uncertainty_score`
- `uncertainty_rank`
- `confidence_class`
- `stability_class`

For a discrete tree set, each row preserves:
- `clade_id`
- `clade_taxa`
- `tree_presence_count`
- `tree_presence_fraction`
- `dominant_state`
- `dominant_state_fraction`
- `unique_state_count`
- `ambiguous_tree_fraction`
- `unstable_tree_fraction`
- `state_distribution`
- `entropy`
- `normalized_entropy`
- `instability_score`
- `uncertainty_score`
- `uncertainty_rank`
- `confidence_class`
- `stability_class`

`ancestral root-sensitivity` is the governed root-assumption review surface
for one discrete likelihood ancestral reconstruction. It reruns the owned
likelihood path under an equal root prior, an empirical root prior derived
from the analyzed tip-state counts, and an optional user-supplied fixed-root
scenario. Its JSON metrics report:
- `model`
- `state_ordering`
- `analyzed_taxon_count`
- `assumption_count`
- `compared_node_count`
- `state_changed_node_count`
- `support_changed_node_count`
- `top_sensitive_node`
- `fixed_root_state`

The command supports `equal-rates`, `symmetric`, and
`all-rates-different`. It does not accept `fitch`, because root priors are a
likelihood-model concern rather than a Fitch set-propagation concern.
`--fixed-root-state` is optional. When it is absent, the surface compares only
equal and empirical root priors. When it is present, the command treats the
named state as a scenario assumption and adds that scenario to the comparison.

When `--summary-out` is supplied, `ancestral root-sensitivity` writes one
overall summary ledger. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `state_ordering`
- `analyzed_taxon_count`
- `assumption_count`
- `compared_node_count`
- `state_changed_node_count`
- `support_changed_node_count`
- `top_sensitive_node`
- `top_sensitive_score`
- `warning_count`

When `--assumptions-out` is supplied, the command writes one root-assumption
ledger. Each row preserves:
- `assumption_id`
- `root_prior_mode`
- `fixed_root_state`
- `root_prior_distribution`
- `root_most_likely_state`
- `root_confidence`
- `root_entropy`
- `unstable_node_count`
- `weak_support_node_count`

When `--nodes-out` is supplied, the command writes one node-wise comparison
ledger. Each row preserves:
- `node`
- `descendant_taxa`
- `assumption_states`
- `assumption_confidences`
- `assumption_entropies`
- `unique_state_count`
- `state_changed`
- `max_confidence_delta`
- `max_entropy_delta`
- `sensitivity_score`
- `sensitivity_rank`
- `stability_class`

`ancestral tree-set` is the governed reconstruction-stability surface for one
trait across a posterior or bootstrap tree set. It reruns ancestral
reconstruction on every retained tree, maps comparable internal clades by
descendant taxa, and reports how often each clade and ancestral conclusion
survive topology uncertainty. Its JSON metrics report:
- `kind`
- `model`
- `total_tree_count`
- `kept_tree_count`
- `rooted_topology_count`
- `unrooted_topology_count`
- `node_row_count`
- `clade_summary_count`
- `excluded_taxon_count`
- `unstable_clade_count`

The command requires `--kind continuous` or `--kind discrete`. Continuous mode
supports `brownian` and `ou`. Discrete mode supports `fitch`, `equal-rates`,
`symmetric`, and `all-rates-different`. `--burnin-fraction` removes the
requested leading fraction of trees before reconstruction, and every retained
tree keeps both its original one-based source index and its post-burnin index.

When `--summary-out` is supplied, `ancestral tree-set` writes one overall
summary ledger. In continuous mode the row preserves:
- `trait`
- `taxon_column`
- `model`
- `alpha`
- `total_tree_count`
- `burnin_tree_count`
- `kept_tree_count`
- `shared_tree_taxon_count`
- `analysis_taxon_count`
- `rooted_topology_count`
- `unrooted_topology_count`
- `clade_summary_count`
- `unstable_clade_count`
- `top_unstable_clade`
- `warning_count`

In discrete mode the same summary ledger preserves `state_ordering` instead of
`alpha`.

When `--trees-out` is supplied, the command writes one retained-tree ledger.
Each row preserves:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `internal_clade_count`

When `--nodes-out` is supplied, the command writes one per-tree internal-node
ledger. Continuous rows preserve:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `clade_id`
- `clade_taxa`
- `estimate`
- `standard_error`
- `lower_95_interval`
- `upper_95_interval`
- `confidence`
- `unstable`

Discrete rows preserve:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `clade_id`
- `clade_taxa`
- `most_likely_state`
- `state_set`
- `confidence`
- `ambiguous`
- `unstable`

When `--clades-out` is supplied, the command writes one comparable-clade
summary ledger. Continuous rows preserve:
- `clade_id`
- `clade_taxa`
- `tree_presence_count`
- `tree_presence_fraction`
- `mean_estimate`
- `median_estimate`
- `standard_deviation`
- `minimum_estimate`
- `maximum_estimate`
- `lower_95_empirical_estimate`
- `upper_95_empirical_estimate`
- `empirical_interval_width`
- `mean_standard_error`
- `unstable_tree_count`
- `unstable_tree_fraction`
- `instability_score`
- `stability_class`

Discrete rows preserve:
- `clade_id`
- `clade_taxa`
- `tree_presence_count`
- `tree_presence_fraction`
- `unique_state_count`
- `dominant_state`
- `dominant_state_tree_count`
- `dominant_state_fraction`
- `ambiguous_tree_count`
- `ambiguous_tree_fraction`
- `unstable_tree_count`
- `unstable_tree_fraction`
- `state_distribution`
- `instability_score`
- `stability_class`

When `--exclusions-out` is supplied, the command writes one explicit excluded
tip ledger with `taxon` and `reason`.

`ancestral transitions` is the governed categorical transition-counting surface
for one rooted tree or a retained tree set. It reruns the owned discrete
ancestral reconstruction path, converts each non-root branch into one explicit
parent-versus-child state comparison, and preserves whether each inferred
change is certain or uncertain. Its JSON metrics report:
- `tree_set`
- `model`
- `excluded_taxon_count`

For one-tree runs, the JSON metrics also report:
- `total_branch_count`
- `changed_branch_count`
- `certain_change_count`
- `uncertain_change_count`
- `transition_pair_count`

For tree-set runs, the JSON metrics also report:
- `total_tree_count`
- `kept_tree_count`
- `rooted_topology_count`
- `unrooted_topology_count`
- `transition_pair_count`
- `topology_sensitive_transition_pair_count`
- `uncertainty_sensitive_transition_pair_count`

The command supports `fitch`, `equal-rates`, `symmetric`, and
`all-rates-different`. Ordered-state transition counting requires a likelihood
model, because ordered-state support comes from the Mk likelihood path rather
than the Fitch set path. `--tree-set` switches the surface from one analyzed
tree to one retained posterior or bootstrap tree set, and `--burnin-fraction`
is only valid with `--tree-set`.

When `--summary-out` is supplied, `ancestral transitions` writes one overall
summary ledger. For one-tree runs the row preserves:
- `trait`
- `taxon_column`
- `model`
- `state_ordering`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `total_branch_count`
- `changed_branch_count`
- `certain_change_count`
- `uncertain_change_count`
- `transition_pair_count`
- `top_transition`
- `warning_count`

For tree-set runs the same summary ledger preserves:
- `trait`
- `taxon_column`
- `model`
- `state_ordering`
- `total_tree_count`
- `burnin_tree_count`
- `kept_tree_count`
- `shared_tree_taxon_count`
- `analysis_taxon_count`
- `rooted_topology_count`
- `unrooted_topology_count`
- `transition_pair_count`
- `topology_sensitive_transition_pair_count`
- `uncertainty_sensitive_transition_pair_count`
- `stable_transition_pair_count`
- `top_transition`
- `warning_count`

When `--branches-out` is supplied, the command writes one branch ledger. For
one-tree runs each row preserves:
- `parent_node`
- `child_node`
- `child_descendant_taxa`
- `branch_length`
- `parent_most_likely_state`
- `child_most_likely_state`
- `parent_state_set`
- `child_state_set`
- `overlapping_states`
- `changed`
- `certainty_class`
- `transition`

For tree-set runs the same branch ledger preserves those columns plus:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`

When `--counts-out` is supplied, the command writes one transition-pair count
ledger. For one-tree runs each row preserves:
- `transition`
- `source_state`
- `target_state`
- `certain_change_count`
- `uncertain_change_count`
- `total_change_count`

For tree-set runs each row preserves:
- `transition`
- `source_state`
- `target_state`
- `tree_presence_count`
- `tree_presence_fraction`
- `mean_certain_change_count`
- `mean_uncertain_change_count`
- `mean_total_change_count`
- `minimum_total_change_count`
- `maximum_total_change_count`
- `lower_95_empirical_total_change_count`
- `upper_95_empirical_total_change_count`
- `stability_class`

When `--trees-out` is supplied with `--tree-set`, the command writes one
retained-tree ledger. Each row preserves:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `branch_count`
- `changed_branch_count`
- `certain_change_count`
- `uncertain_change_count`
- `transition_pair_count`

When `--exclusions-out` is supplied, the command writes one explicit excluded
tip ledger. Each row preserves:
- `taxon`
- `reason`

`ancestral render` is the governed visualization surface for one ancestral
reconstruction. It accepts one continuous or discrete ancestral model input and
emits one reviewer-facing tree figure as SVG, PNG, or HTML depending on the
`--out` suffix. Its JSON metrics report:
- `tip_count`
- `format`
- `layout`
- `rendered_internal_annotation_count`
- `rendered_internal_pie_count`
- `rendered_branch_color_count`

The command requires `--kind continuous` or `--kind discrete`. Continuous mode
supports `brownian` and `ou`. Discrete mode supports `fitch`, `equal-rates`,
`symmetric`, and `all-rates-different`. The `--layout` surface remains
`cladogram`, `phylogram`, or `circular`.

For discrete visualization, `--discrete-node-style` chooses:
- `labels`: one internal text label per node
- `pies`: one marginal-state pie marker per node

For branch coloring, `--branch-coloring` chooses:
- `none`: keep branch lines neutral
- `state`: color descendant branches by inferred discrete state
- `regime`: color descendant branches by reconstructed continuous value regime

Continuous rendering rejects `state` branch coloring, and discrete rendering
rejects `regime` branch coloring, so the branch palette always matches the
underlying ancestral evidence type instead of silently coercing one mode into
another.

When `--out` ends in `.svg`, the command writes one standalone SVG figure. When
`--out` ends in `.png`, the command writes one raster PNG and also writes the
governed sibling SVG used to create that raster. When `--out` ends in `.html`,
the command writes one standalone HTML review page with the governed sibling
SVG embedded directly in the figure section.

`ancestral package` is the governed publication-bundle companion surface for
ancestral visualization. Its JSON metrics report:
- `output_dir`
- `artifact_count`

The package writes one figure bundle plus the node and uncertainty ledgers. The
visual artifacts are:
- `ancestral-figure.svg`
- `ancestral-figure.png`
- `ancestral-figure.html`

The tabular and descriptive artifacts remain:
- `node-states.tsv`
- `uncertainty.tsv`
- `legend.md`
- `model-description.md`
- `figure-caption.md`
- `figure-manifest.json`

The package uses the same owned visualization contract as `ancestral render`.
Continuous bundles color branches by reconstructed value regime. Discrete
bundles render marginal-state pies and color branches by inferred descendant
state, so the publication bundle preserves the richer visual evidence surface
without needing a separate styling step.

`ancestral report` is the governed full reconstruction-review surface for one
continuous or discrete ancestral analysis. It can still write one standalone
HTML report through `--out`, but `--out-dir` activates the complete package
surface that keeps the report, visualization files, node ledger, uncertainty
ledger, transition or branch-change ledgers, exclusion ledger, and manifest in
one directory. Its JSON metrics report:
- `report_kind`
- `reconstruction_kind`
- `output_dir`
- `artifact_count`
- `transition_count_row_count`

The command requires `--kind continuous` or `--kind discrete`. Continuous mode
supports `brownian` and `ou`. Discrete mode supports `fitch`, `equal-rates`,
`symmetric`, and `all-rates-different`. The same `--state-ordering`,
`--ordered-states`, `--compare-model`, `--compare-tree`, `--drop-taxa`, and
`--coding-map` surfaces remain available on the standalone and packaged paths.

When `--out-dir` is supplied, `ancestral report` writes this fixed package:
- `ancestral-report.html`
- `ancestral-figure.svg`
- `ancestral-figure.png`
- `ancestral-figure.html`
- `summary.tsv`
- `node-table.tsv`
- `uncertainty-table.tsv`
- `transition-counts.tsv`
- `transition-branches.tsv`
- `exclusions.tsv`
- `ancestral-report.manifest.json`

The HTML report embeds the governed SVG figure directly and turns the core
ancestral review surfaces into one durable handoff. `summary.tsv` keeps the
owned continuous or discrete summary row. `node-table.tsv` keeps the flattened
node-state ledger from the owned ancestral report surface.
`uncertainty-table.tsv` keeps either the continuous uncertainty ledger or the
discrete marginal-state probability ledger.

`transition-counts.tsv` and `transition-branches.tsv` preserve different but
explicit evidence depending on the reconstruction kind. For discrete traits,
they reuse the owned ancestral transition workflow, so each row is a directed
state-change count or branch record. For continuous traits, the same filenames
stay stable across the report package, but the rows become branch-delta review
ledgers with `increase`, `decrease`, and `stable` directions rather than
pretending that continuous values imply categorical state transitions.

When `--out` is supplied together with `--out-dir`, the command also writes one
copy of the packaged HTML report to the requested `--out` path and one sibling
SVG beside it. This keeps the complete package surface while still supporting
one explicitly named review artifact path.

`comparative pgls` is the governed regression surface for continuous trait
association under phylogenetic covariance. Its JSON metrics now report
`coefficient_count`, `confidence_interval_count`,
`residual_degrees_of_freedom`, `coefficient_inference_distribution`, and `aic`
so review tooling can distinguish a minimally identified model from one with
meaningful residual support and can compare model fit without scraping the
coefficient table.

`comparative covariance-audit` is the governed pre-fit review surface for
PGLS, Brownian trait evolution, and OU trait evolution. Its purpose is to
show whether the tree-trait overlap and induced covariance matrix are safe to
trust before coefficient or parameter interpretation starts.

Its JSON metrics report:
- `analysis`
- `covariance_model`
- `matrix_dimension`
- `matrix_rank`
- `condition_number`
- `fit_strategy`
- `singular`
- `near_singular`
- `matched_taxon_count`
- `analysis_taxon_count`
- `duplicate_tree_taxon_count`
- `duplicate_trait_taxon_count`
- `candidate_row_count`
- `blocker_count`
- `warning_count`

The command requires `--analysis pgls`, `--analysis brownian-trait`, or
`--analysis ou-trait`. PGLS accepts `--formula` or the
`--response` plus `--predictors` surface. Brownian and OU trait audits accept
`--trait`. `--lambda-value` accepts `estimate` or one numeric Pagel's lambda
for PGLS. `--alpha` accepts `estimate` or one positive numeric OU alpha for OU
trait audits.

The command reports:
- matched taxa
- tree taxa missing from the trait table
- extra trait-table taxa absent from the tree
- duplicate tree or trait taxa
- zero-length and negative branch counts
- minimum and maximum branch length
- whether the covariance is singular or near-singular
- whether the fitting path would proceed by `exact`, `regularization`,
  `pseudoinverse`, or `failure`

The candidate-level audit is explicit rather than narrative. When the audit
profiles estimated Pagel's lambda or OU alpha, each candidate row records its
matrix rank, raw condition number, stabilized fit condition number,
positive-definiteness before fitting, and the fit-strategy details used for
that candidate.

When `--summary-out` is supplied, `comparative covariance-audit` writes one
flat summary row as CSV or TSV. The row preserves:
- `analysis`
- `covariance_model`
- `analysis_label`
- `matrix_dimension`
- `matrix_rank`
- `condition_number`
- `fit_strategy`
- `singular`
- `near_singular`
- `tree_taxon_count`
- `trait_taxon_count`
- `matched_taxon_count`
- `analysis_taxon_count`
- `missing_from_traits_count`
- `extra_trait_taxon_count`
- `duplicate_tree_taxon_count`
- `duplicate_trait_taxon_count`
- `empty_trait_taxon_row_count`
- `zero_length_branch_count`
- `negative_branch_length_count`
- `minimum_branch_length`
- `maximum_branch_length`
- `blocker_count`
- `warning_count`

When `--candidates-out` is supplied, the command also writes one candidate
ledger with:
- `candidate_label`
- `parameter_name`
- `parameter_value`
- `matrix_dimension`
- `matrix_rank`
- `condition_number`
- `fit_condition_number`
- `positive_definite_before_fit`
- `singular`
- `near_singular`
- `fit_strategy`
- `fit_strategy_details`

When `--excluded-taxa-out` is supplied, the command writes one explicit
excluded-taxa ledger with:
- `taxon`
- `reason`
- `details`

The fit-strategy field is intentionally honest about the current runtime. The
governed Brownian, OU, and Pagel-lambda fitting paths currently stabilize
covariance inversion with diagonal epsilon regularization where needed, and the
audit reports that directly instead of implying an exact closed-form solve when
stabilization was required.

`comparative logistic` is the governed binary-response companion surface. It
preserves the same formula and predictor-encoding contract as `comparative
pgls`, but the fitted model is an explicit
`phylogenetic-working-correlation-gee` approximation rather than continuous
generalized least squares. It does not currently claim `ape::compar.gee`
parity and should not be treated as a drop-in `ape::compar.gee`
implementation. Its JSON metrics report:
- `taxon_count`
- `success_count`
- `failure_count`
- `coefficient_count`
- `fitted_row_count`
- `lambda_value`
- `approximation_method`
- `converged`
- `iteration_count`
- `binomial_log_likelihood`
- `separation_detected`
- `warning_count`
- `coefficient_inference_distribution`
- `method_excluded_reference_surfaces`

The command requires the response to be encoded as `0` and `1`. One-class
responses are rejected instead of silently producing degenerate output.

When `--coefficients-out` is supplied, `comparative logistic` writes one flat
coefficient ledger as CSV or TSV. Each row preserves:
- `response`
- `term`
- `estimate`
- `standard_error`
- `test_statistic`
- `p_value`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`
- `inference_distribution`
- `approximation_method`
- `lambda_value`
- `taxon_count`
- `success_count`
- `failure_count`
- `converged`
- `iteration_count`
- `binomial_log_likelihood`
- `separation_detected`

When `--fitted-out` is supplied, the command also writes one taxon-level
probability ledger as CSV or TSV. Each row preserves:
- `taxon`
- `observed_response`
- `fitted_probability`
- `linear_predictor`
- `residual`

When `--excluded-taxa-out` is supplied, the command writes one explicit
excluded-taxa ledger with:
- `taxon`
- `reason`
- `details`

Warnings are not hidden. If the fit requires information-matrix stabilization,
reaches the iteration limit, drives fitted probabilities to the `0/1`
boundary, or produces very large coefficients, the JSON output marks that as
`separation_detected` and preserves the warning messages directly.

`comparative correlated-traits` is the governed review surface for pairwise
trait-evolution coupling on one tree. It supports two analysis families:
- two numeric traits: `continuous-brownian-contrasts`
- two binary traits: `binary-joint-state`

Its JSON metrics report:
- `analysis_kind`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observation_row_count`
- `comparison_row_count`
- `association_measure_name`
- `association_measure_value`
- `evolutionary_covariance`
- `evolutionary_correlation`
- `better_model`
- `likelihood_ratio_p_value`
- `joint_state_count`
- `warning_count`

The command takes `--left-trait` and `--right-trait` instead of a regression
formula because this surface is not a predictor-response fit. `--analysis-kind`
defaults to `auto`, but may be forced to `continuous` or `binary`. For binary
review, `--binary-model` chooses the governed discrete transition surface from
`equal-rates`, `symmetric`, or `all-rates-different`.

The continuous path compares diagonal versus full Brownian contrast covariance.
The binary path compares separate discrete pseudo-likelihood fits against one
joint-state pseudo-likelihood fit over the observed `00`, `01`, `10`, and `11`
state combinations. The binary path is reviewer-facing and explicit about that
approximation boundary; it is not presented as a hidden full Pagel binary
correlation likelihood.

When `--summary-out` is supplied, `comparative correlated-traits` writes one
flat summary ledger as CSV or TSV. The row preserves:
- `analysis_kind`
- `left_trait`
- `right_trait`
- `taxon_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observation_row_count`
- `association_measure_name`
- `association_measure_value`
- `evolutionary_covariance`
- `evolutionary_correlation`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`
- `independent_parameter_count`
- `independent_log_likelihood`
- `independent_aic`
- `correlated_parameter_count`
- `correlated_log_likelihood`
- `correlated_aic`
- `better_model`
- `likelihood_ratio_statistic`
- `likelihood_ratio_degrees_of_freedom`
- `likelihood_ratio_p_value`
- `likelihood_ratio_p_value_method`
- `left_root_estimate`
- `right_root_estimate`
- `left_state_order`
- `right_state_order`
- `joint_state_count`
- `warning_count`

When `--comparison-out` is supplied, the command also writes one flat
model-comparison ledger with:
- `model_kind`
- `model_description`
- `parameter_count`
- `log_likelihood`
- `aic`
- `delta_aic`
- `selected`

When `--observations-out` is supplied, the command also writes one reviewer
evidence ledger. Each row preserves:
- `row_kind`
- `label`
- `taxon`
- `left_taxa`
- `right_taxa`
- `left_numeric_value`
- `right_numeric_value`
- `expected_variance`
- `left_state`
- `right_state`
- `joint_state`

When `--excluded-taxa-out` is supplied, the command also writes one explicit
excluded-taxa ledger with:
- `taxon`
- `reason`
- `missing_traits`

`comparative model-selection` is the governed review surface for comparing
competing comparative regression formulas on one shared response and one shared
complete-case taxon set. Its JSON metrics report:
- `response`
- `model_family`
- `model_count`
- `analysis_taxon_count`
- `excluded_taxon_count`
- `pairwise_comparison_count`
- `best_formula`
- `selected_criterion`
- `selected_log_likelihood`

The command requires at least two `--formula` values. Every candidate formula
must use the same response column. The model family is inferred from the shared
response after complete-case pruning:
- binary `0/1` response on the shared taxon set: phylogenetic logistic ranking
- otherwise: PGLS ranking

This surface is explicit about fairness. Model ranking only happens after one
common analysis set is fixed across all formulas, so a candidate does not gain
an artificial AIC advantage by silently dropping a different subset of taxa.

When `--ranking-out` is supplied, `comparative model-selection` writes one flat
ranking ledger as CSV or TSV. Each row preserves:
- `formula`
- `model_family`
- `parameter_count`
- `taxon_count`
- `phylogenetic_parameter_name`
- `phylogenetic_parameter_value`
- `phylogenetic_parameter_estimated`
- `log_likelihood`
- `aic`
- `aicc`
- `bic`
- `delta_aicc`
- `delta_bic`
- `akaike_weight`
- `rank`
- `selected`
- `encoded_columns`
- `warning_count`
- `separation_detected`

When `--pairwise-out` is supplied, the command also writes one pairwise
comparison ledger as CSV or TSV. Each row preserves:
- `left_formula`
- `right_formula`
- `comparison_kind`
- `preferred_formula`
- `left_rank`
- `right_rank`
- `left_parameter_count`
- `right_parameter_count`
- `delta_parameter_count`
- `left_log_likelihood`
- `right_log_likelihood`
- `left_aicc`
- `right_aicc`
- `left_bic`
- `right_bic`
- `likelihood_ratio_statistic`

`comparison_kind` is one of:
- `identical`
- `left_nested_in_right`
- `right_nested_in_left`
- `non_nested`

When `--excluded-taxa-out` is supplied, the command writes one explicit
shared-complete-case exclusion ledger with:
- `taxon`
- `reason`
- `missing_columns`

`comparative contrasts` is the governed review surface for phylogenetic
independent contrasts. Its JSON metrics report:
- `taxon_count`
- `contrast_count`
- `regression_row_count`
- `regression_slope`
- `regression_p_value`

Without `--predictor-trait`, the command returns one contrast report for the
requested trait. When `--predictor-trait` is supplied, the command also fits a
through-origin regression on the matched node-level contrasts and preserves
that under `data.regression`. The returned `data.contrast_report.input_audit`
also preserves the owned input-policy surface for the run:
- `tree_is_ultrametric`
- `minimum_root_to_tip_depth`
- `maximum_root_to_tip_depth`
- `ultrametric_policy`
- `missing_value_policy`
- `pruned_missing_value_taxa`
- `warnings`

When `--contrasts-out` is supplied, `comparative contrasts` writes one flat
contrast ledger as CSV or TSV. Each row preserves:
- `trait`
- `node_id`
- `node`
- `left_taxa`
- `right_taxa`
- `contrast`
- `expected_variance`
- `ancestral_value`
- `root_estimate`

When `--regression-out` is supplied together with `--predictor-trait`,
`comparative contrasts` also writes one regression-through-origin ledger as CSV
or TSV. Each row preserves:
- `response_trait`
- `predictor_trait`
- `node`
- `predictor_contrast`
- `response_contrast`
- `fitted_response_contrast`
- `residual`
- `leverage_fraction`
- `slope`
- `standard_error`
- `test_statistic`
- `p_value`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`
- `residual_sum_of_squares`
- `r_squared_through_origin`

The regression output is explicit rather than inferred. `--regression-out`
without `--predictor-trait` is rejected instead of silently writing nothing.
The owned Bijux surface now also has governed live `ape::pic` parity on
balanced rooted ultrametric, pectinate rooted non-ultrametric, and six-taxon
clean trait-table fixtures. Missing trait values remain an explicit owned
pruning policy rather than a literal live-`ape` parity lane, and negative
branch lengths are rejected as an invalid comparative-analysis boundary.

`comparative signal` is the governed review surface for one-trait phylogenetic
signal. Its JSON metrics report:
- `taxon_count`
- `blombergs_k`
- `pagels_lambda`
- `lambda_log_likelihood`
- `lambda_likelihood_ratio_statistic`
- `signal_p_value`
- `tree_is_ultrametric`
- `ultrametric_policy`
- `missing_value_policy`
- `pruned_missing_value_taxon_count`
- `signal_seed`
- `signal_null_k_minimum`
- `signal_null_k_mean`
- `signal_null_k_maximum`
- `lambda_likelihood_ratio_p_value`
- `lambda_optimizer_name`
- `lambda_optimizer_function_evaluation_count`
- `lambda_optimizer_hit_lower_boundary`
- `lambda_optimizer_hit_upper_boundary`
- `permutation_row_count`

The command preserves four distinct surfaces under `data`:
- `input_audit` for the rootedness, ultrametricity, pruning, and warning policy
- `blombergs_k` for the fitted K summary
- `pagels_lambda` for the fitted lambda summary
- `signal_test` for the permutation-based K test

The `signal_test` report now also preserves the seeded null-distribution
summary that the live `phytools::phylosig(method='K')` lane compares:
- `observed_k`
- `p_value`
- `permutations`
- `seed`
- `permuted_k_at_or_above_observed`
- `null_distribution_minimum`
- `null_distribution_mean`
- `null_distribution_maximum`

The `pagels_lambda` report now also preserves the fixed-lambda likelihood
context that the live `phytools::phylosig(method='lambda')` lane compares:
- `lambda_value`
- `log_likelihood`
- `null_log_likelihood`
- `brownian_log_likelihood`
- `likelihood_ratio_statistic`
- `likelihood_ratio_p_value`
- `p_value_method`
- `optimizer_diagnostics`
- `profile_rows`

The governed signal policy is explicit rather than implicit:
- rooted trees with branch lengths are accepted whether or not they are ultrametric
- ultrametric status is reported, not silently assumed
- overlapping missing trait values are pruned and reported under `input_audit`
- permutation rows are reproducible from `--seed`
- constant post-pruning trait vectors fail with `comparative_method_error`

When `--summary-out` is supplied, `comparative signal` writes one flat summary
ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_count`
- `blombergs_k`
- `blombergs_generalized_mean`
- `blombergs_observed_mean_square`
- `blombergs_phylogenetic_mean_square`
- `blombergs_expected_mean_square_ratio`
- `signal_permutation_p_value`
- `pagels_lambda`
- `lambda_log_likelihood`
- `lambda_null_log_likelihood`
- `lambda_brownian_log_likelihood`
- `lambda_likelihood_ratio_statistic`
- `lambda_likelihood_ratio_p_value`
- `lambda_p_value_method`
- `permutations`
- `permuted_k_at_or_above_observed`

When `--permutations-out` is supplied, `comparative signal` also writes one
permutation ledger as CSV or TSV. Each row preserves:
- `trait`
- `observed_k`
- `estimated_lambda`
- `permutations`
- `signal_permutation_p_value`
- `permutation_index`
- `permuted_k`
- `at_or_above_observed`

This surface exists so phylogenetic signal review does not collapse into one
scalar. Reviewers can inspect the fitted K and lambda values, the permutation
null distribution, and the explicit p-value contract without rerunning the
analysis manually, while still seeing whether the fit depended on pruning or a
non-ultrametric rooted tree.

`comparative discrete-mk` is the governed standalone discrete Mk fit surface
for one rooted tree and one categorical tip trait. Its JSON metrics report:
- `taxon_count`
- `model`
- `observed_state_count`
- `sparse_state_count`
- `pruned_missing_value_taxon_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `aicc`
- `optimizer_name`
- `optimizer_converged`
- `optimizer_iteration_count`
- `optimizer_function_evaluation_count`
- `optimizer_hit_lower_parameter_bound`
- `optimizer_hit_upper_parameter_bound`
- `overparameterized`
- `transition_rate_count`
- `baseline_model`
- `baseline_aic`
- `delta_aic`
- `preferred_model_by_aic`

The command preserves the full fit report under `data`, including:
- `input_audit`
- `transition_rate_rows`
- `optimizer_diagnostics`
- `baseline_comparison`

The governed policy is explicit rather than implicit:
- rooted trees with branch lengths are required
- overlapping missing trait values are pruned and reported under `input_audit`
- sparse states are surfaced explicitly instead of hidden
- optimizer non-convergence or boundary hits remain warnings, not silent success
- ER fits are the governed live `phytools::fitMk(model='ER')` parity surface
- unordered multistate SYM fits are the governed live `phytools::fitMk(model='SYM')` parity surface

When `--summary-out` is supplied, `comparative discrete-mk` writes one flat
summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `state_ordering`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `observed_state_count`
- `sparse_state_count`
- `log_likelihood`
- `parameter_count`
- `aic`
- `aicc`
- `optimizer_name`
- `optimizer_converged`
- `optimizer_iteration_count`
- `optimizer_function_evaluation_count`
- `optimizer_hit_lower_parameter_bound`
- `optimizer_hit_upper_parameter_bound`
- `overparameterized`
- `warning_count`

When `--rates-out` is supplied, the command also writes one directed
rate-matrix ledger as CSV or TSV. Each row preserves:
- `source_state`
- `target_state`
- `transition_allowed`
- `step_distance`
- `rate`

`comparative brownian` is the governed standalone Brownian trait-evolution
surface for one numeric trait on a rooted tree with branch lengths. Its JSON
metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `root_state`
- `sigma_squared`
- `log_likelihood`
- `aic`
- `aicc`

The command preserves the full summary under `data`, including:
- `analyzed_taxa`
- `excluded_taxa`
- `confidence_intervals`
- `residual_diagnostics`
- `readiness`

When `--summary-out` is supplied, `comparative brownian` writes one flat
summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `root_state`
- `root_state_lower_95`
- `root_state_upper_95`
- `sigma_squared`
- `sigma_squared_lower_95`
- `sigma_squared_upper_95`
- `log_likelihood`
- `aic`
- `aicc`
- `residual_variance`
- `max_abs_standardized_residual`
- `phylogenetic_residual_lambda`

When `--excluded-taxa-out` is supplied, `comparative brownian` also writes one
excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so Brownian trait-evolution review preserves both the fit
statistics and the taxon-pruning contract instead of reducing the workflow to a
single reported rate.

`comparative brownian-regimes` is the governed standalone multi-rate Brownian
surface for one numeric trait on a rooted tree with branch lengths plus a
user-supplied branch regime map. Its JSON metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `regime_count`
- `root_state`
- `log_likelihood`
- `aic`
- `aicc`
- `better_model`
- `likelihood_ratio_statistic`
- `likelihood_ratio_p_value`
- `identifiability_warning_count`
- `profile_row_count`

The command preserves the full summary under `data`, including:
- `analyzed_taxa`
- `excluded_taxa`
- `branch_rows`
- `regime_rows`
- `comparison_rows`
- `profile_rows`
- `identifiability_warnings`
- `residual_diagnostics`
- `readiness`

The branch regime map must assign every non-root branch exactly one regime. By
default the table uses:
- `branch_id`
- `regime`

Each `branch_id` is the normalized descendant-tip signature for the branch,
such as `A|B` or `A|B|C|D`. The optional `--branch-id-column` and
`--regime-column` flags let the same contract be read from differently named
columns without changing the underlying branch identity rule.

When `--summary-out` is supplied, `comparative brownian-regimes` writes one
flat summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `branch_id_column`
- `regime_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `regime_count`
- `root_state`
- `root_state_lower_95`
- `root_state_upper_95`
- `log_likelihood`
- `aic`
- `aicc`
- `better_model`
- `likelihood_ratio_statistic`
- `likelihood_ratio_degrees_of_freedom`
- `likelihood_ratio_p_value`
- `identifiability_warning_count`
- `residual_variance`
- `max_abs_standardized_residual`
- `phylogenetic_residual_lambda`

When `--rates-out` is supplied, `comparative brownian-regimes` also writes one
per-regime rate ledger as CSV or TSV. Each row preserves:
- `regime`
- `branch_count`
- `contributing_branch_count`
- `total_branch_length`
- `contributing_branch_length`
- `sigma_squared`
- `sigma_squared_lower_95`
- `sigma_squared_upper_95`
- `interval_method`

When `--comparison-out` is supplied, `comparative brownian-regimes` also writes
one single-rate versus multi-rate comparison ledger as CSV or TSV. The ledger
preserves:
- model-fit rows with `model`, `parameter_count`, `log_likelihood`, `aic`,
  `aicc`, `delta_aicc`, and `selected`
- one likelihood-ratio row with `comparison_id`, `left_model`, `right_model`,
  `statistic`, `degrees_of_freedom`, `p_value`, and `p_value_method`

When `--profile-out` is supplied, `comparative brownian-regimes` also writes
one conditional regime-rate profile as CSV or TSV. Each row preserves:
- `regime`
- `sigma_squared`
- `log_likelihood`
- `delta_log_likelihood`
- `in_support_interval`
- `selected`

When `--branches-out` is supplied, `comparative brownian-regimes` also writes
one normalized branch-assignment ledger as CSV or TSV. Each row preserves:
- `branch_id`
- `regime`
- `branch_length`
- `descendant_taxa`
- `analyzed_descendant_taxa`
- `contributes_to_analysis`

When `--excluded-taxa-out` is supplied, `comparative brownian-regimes` also
writes one excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so regime-aware Brownian review preserves the exact
branch-to-regime contract, the per-regime uncertainty surface, and the explicit
comparison against single-rate Brownian instead of collapsing into one claimed
rate shift.

`comparative regime-map` is the governed review surface for constructing or
validating branch regime assignments before a downstream regime-aware
comparative fit. It accepts exactly one source:
- `--table` plus `--trait` for discrete tip-state reconstruction
- `--regime-map` for a user-provided branch regime table

Its JSON metrics report:
- `source_kind`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `regime_count`
- `branch_count`
- `node_count`
- `ambiguous_branch_count`
- `rendered_internal_annotation_count`
- `rendered_categorical_trait_count`

The command preserves the full regime assignment under `data`, including:
- `observed_regimes`
- `branch_rows`
- `node_rows`
- `excluded_taxa`
- `warnings`

The normalized branch identity rule is explicit and shared with the
regime-aware Brownian workflow. Every non-root branch is identified by the
descendant-tip signature of the child node, such as `A|B` or `A|B|C|D`. When a
user-provided map is used, every non-root branch must appear exactly once.

When `--summary-out` is supplied, `comparative regime-map` writes one flat
summary ledger as CSV or TSV. The row preserves:
- `source_kind`
- `trait`
- `taxon_column`
- `reconstruction_model`
- `state_ordering`
- `ordered_states`
- `branch_id_column`
- `regime_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `branch_count`
- `regime_count`
- `ambiguous_branch_count`
- `node_count`

When `--branches-out` is supplied, `comparative regime-map` also writes one
normalized branch-regime ledger as CSV or TSV. Each row preserves:
- `branch_id`
- `child_node_name`
- `is_tip_branch`
- `branch_length`
- `regime`
- `candidate_regimes`
- `assignment_confidence`
- `ambiguous_assignment`
- `assignment_origin`
- `descendant_taxa`
- `analyzed_descendant_taxa`
- `contributes_to_analysis`

When `--nodes-out` is supplied, `comparative regime-map` also writes one
node-reconstruction ledger as CSV or TSV. Each row preserves:
- `node_id`
- `node_name`
- `is_tip`
- `descendant_taxa`
- `regime`
- `candidate_regimes`
- `assignment_confidence`
- `ambiguous_assignment`
- `state_probabilities`

This ledger is only populated when the source is `--table`. A user-provided
branch map does not infer hidden ancestral nodes and therefore keeps `node_rows`
empty.

When `--excluded-taxa-out` is supplied, `comparative regime-map` writes one
explicit excluded-taxa ledger for tip-state reconstruction. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states:
- `missing_from_state_table`
- `missing_state_value`
- `absent_from_tree`

When `--svg-out` is supplied, `comparative regime-map` renders one SVG tree
with the regime assignment overlaid on the analyzed tree. `--layout` accepts:
- `cladogram`
- `phylogram`
- `circular`

This surface exists so branch regimes can be reviewed as a first-class input
artifact before they are used to claim rate differences, clade shifts, or
ecological context in downstream comparative workflows.

`comparative ou` is the governed standalone OU trait-evolution surface for one
numeric trait on a rooted tree with branch lengths. Its JSON metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `alpha`
- `theta`
- `sigma_squared`
- `log_likelihood`
- `aic`
- `aicc`

The command preserves the full summary under `data`, including:
- `analyzed_taxa`
- `excluded_taxa`
- `confidence_intervals`
- `identifiability_warnings`
- `residual_diagnostics`
- `readiness`

When `--summary-out` is supplied, `comparative ou` writes one flat summary
ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `alpha`
- `alpha_lower_95`
- `alpha_upper_95`
- `theta`
- `theta_lower_95`
- `theta_upper_95`
- `sigma_squared`
- `sigma_squared_lower_95`
- `sigma_squared_upper_95`
- `log_likelihood`
- `aic`
- `aicc`
- `convergence_status`
- `identifiability_warning_count`
- `residual_variance`
- `max_abs_standardized_residual`
- `phylogenetic_residual_lambda`

When `--excluded-taxa-out` is supplied, `comparative ou` also writes one
excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so OU trait-evolution review preserves the fitted optimum,
pull strength, diffusion rate, and the pruning contract instead of reducing the
workflow to one preferred-alpha summary line.

`comparative early-burst` is the governed standalone early-burst trait-evolution
surface for one numeric trait on a rooted tree with branch lengths. Its JSON
metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `rate_change`
- `root_state`
- `sigma_squared`
- `log_likelihood`
- `aic`
- `aicc`
- `better_model`
- `identifiability_warning_count`
- `profile_row_count`

The command preserves the full summary under `data`, including:
- `analyzed_taxa`
- `excluded_taxa`
- `confidence_intervals`
- `comparison_rows`
- `likelihood_ratio_tests`
- `profile_rows`
- `identifiability_warnings`
- `residual_diagnostics`
- `readiness`

When `--summary-out` is supplied, `comparative early-burst` writes one flat
summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `rate_change`
- `rate_change_lower_95`
- `rate_change_upper_95`
- `root_state`
- `sigma_squared`
- `sigma_squared_lower_95`
- `sigma_squared_upper_95`
- `log_likelihood`
- `aic`
- `aicc`
- `better_model`
- `identifiability_warning_count`
- `residual_variance`
- `max_abs_standardized_residual`
- `phylogenetic_residual_lambda`

When `--excluded-taxa-out` is supplied, `comparative early-burst` also writes
one excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

When `--comparison-out` is supplied, `comparative early-burst` also writes one
combined comparison ledger as CSV or TSV. The ledger preserves:
- model-fit rows with `model`, `parameter_count`, `log_likelihood`, `aic`,
  `aicc`, `delta_aicc`, and `selected`
- likelihood-ratio rows with `comparison_id`, `left_mode`, `right_mode`,
  `statistic`, `degrees_of_freedom`, and `p_value`

When `--profile-out` is supplied, `comparative early-burst` also writes one
bounded rate-change likelihood profile as CSV or TSV. Each row preserves:
- `trait`
- `rate_change`
- `log_likelihood`
- `aic`
- `aicc`
- `delta_log_likelihood`
- `in_support_interval`
- `selected`

This surface exists so early-burst review preserves the bounded rate-change
profile, the explicit BM/OU comparison contract, and weak-identifiability
warnings instead of reducing the workflow to one optimistic point estimate.

`comparative rate-through-time` is the governed review surface for inspecting
whether reconstructed continuous-trait change is concentrated deeper or
shallower in one rooted tree. Its JSON metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `interval_count`
- `nonempty_interval_count`
- `tree_depth`
- `trend_direction`
- `earliest_interval_rate`
- `latest_interval_rate`
- `latest_to_earliest_rate_ratio`
- `weighted_rate_slope`
- `normalized_rate_slope`

The command requires a rooted tree with strictly positive branch lengths and at
least three analyzed taxa with numeric values for the requested trait. It uses
the Brownian continuous ancestral-state surface to reconstruct internal values,
then bins branches by root depth and allocates reconstructed squared change
across those depth intervals.

When `--summary-out` is supplied, `comparative rate-through-time` writes one
flat summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `interval_count`
- `nonempty_interval_count`
- `tree_depth`
- `ancestral_model`
- `earliest_interval_rate`
- `latest_interval_rate`
- `latest_to_earliest_rate_ratio`
- `weighted_rate_slope`
- `normalized_rate_slope`
- `trend_direction`
- `peak_interval_index`
- `trough_interval_index`
- `assumptions`
- `warnings`

When `--intervals-out` is supplied, the command also writes one depth-binned
interval ledger as CSV or TSV. Each row preserves:
- `interval_index`
- `start_depth`
- `end_depth`
- `midpoint_depth`
- `branch_length_in_interval`
- `change_sum`
- `estimated_rate`
- `branch_count`
- `is_empty`

When `--excluded-taxa-out` is supplied, the command also writes one explicit
excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so rate-through-time review preserves the interval ledger,
trend metrics, and pruning contract instead of reducing the workflow to one
visual impression from a plot or one unqualified scalar trend claim.

`comparative clade-traits` is the governed review surface for summarizing one
continuous or categorical trait across internal non-root clades in the analyzed
tree. Its JSON metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `trait_kind`
- `clade_count`
- `exceptional_clade_count`
- `top_exceptional_clade`
- `top_exceptionality_score`

The command analyzes one trait at a time. `--trait-kind auto` uses the table
schema inference, while `--trait-kind continuous` or `--trait-kind categorical`
can be used when a dirty column would otherwise infer the wrong family. Only
internal non-root clades meeting `--min-clade-size` are ranked.

When `--summary-out` is supplied, `comparative clade-traits` writes one flat
summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `trait_kind`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `minimum_clade_size`
- `clade_count`
- `exceptional_clade_count`
- `top_exceptional_clade`
- `top_exceptionality_score`
- `baseline_mean`
- `baseline_median`
- `baseline_minimum`
- `baseline_maximum`
- `baseline_range_width`
- `baseline_dominant_state`
- `baseline_dominant_state_fraction`
- `assumptions`
- `warnings`

When `--clades-out` is supplied, the command also writes one internal-clade
ledger as CSV or TSV. Each row preserves:
- `clade_id`
- `node_label`
- `trait_kind`
- `taxon_count`
- `taxa`
- `coverage_fraction`
- `mean`
- `median`
- `minimum`
- `maximum`
- `range_width`
- `mean_delta_from_global`
- `dominant_state`
- `dominant_state_count`
- `dominant_state_fraction`
- `dominant_state_enrichment`
- `distinct_state_count`
- `state_counts`
- `distribution_shift`
- `exceptionality_score`
- `exceptional`
- `rank`

When `--excluded-taxa-out` is supplied, the command also writes one explicit
excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so clade-level trait review preserves the ranking
heuristic, sample sizes, and pruning contract instead of reducing the question
to a visual impression or an undocumented spreadsheet sort.

`comparative trait-outliers` is the governed review surface for ranking
continuous-trait taxa by leave-one-taxon-out conditional phylogenetic residual
size. Its JSON metrics report:
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `selected_model`
- `outlier_count`
- `top_outlier_taxon`
- `top_abs_standardized_residual`

The command requires one rooted tree with complete branch lengths and one
numeric trait. Bijux fits standalone Brownian and OU continuous-trait models,
selects the better fit by AICc, then conditions each analyzed tip on all other
retained tips under that selected covariance surface. The taxon rank is thus a
real model-based conditional residual review, not a flat z-score over raw trait
values.

When `--summary-out` is supplied, `comparative trait-outliers` writes one flat
summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `tree_taxon_count`
- `analyzed_taxon_count`
- `excluded_taxon_count`
- `selected_model`
- `selected_mean_parameter`
- `selected_mean_value`
- `selected_alpha`
- `selected_sigma_squared`
- `brownian_aicc`
- `ou_aicc`
- `outlier_threshold`
- `outlier_count`
- `top_outlier_taxon`
- `top_abs_standardized_residual`

When `--outliers-out` is supplied, the command also writes one ranked taxon
ledger as CSV or TSV. Each row preserves:
- `taxon`
- `observed_value`
- `conditional_expected_value`
- `residual`
- `residual_direction`
- `conditional_variance`
- `conditional_standard_error`
- `standardized_residual`
- `abs_standardized_residual`
- `context_clade_id`
- `context_node_label`
- `context_taxon_count`
- `context_taxa`
- `context_mean`
- `sibling_context_id`
- `sibling_taxon_count`
- `sibling_taxa`
- `sibling_mean`
- `context_mean_shift`
- `outlier`
- `rank`

When `--excluded-taxa-out` is supplied, the command also writes one explicit
excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `missing_from_trait_table`
- `missing_trait_value`
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so taxon-level anomaly review preserves the selected model,
the conditional residual contract, and the local clade context instead of
reducing the question to an informal scan of raw values.

`comparative trait-imputation` is the governed review surface for imputing
missing continuous traits under a Brownian phylogenetic model. Its JSON metrics
report:
- `tree_taxon_count`
- `observed_taxon_count`
- `imputed_taxon_count`
- `excluded_taxon_count`
- `holdout_validation_status`
- `holdout_count`
- `holdout_mean_absolute_error`
- `holdout_interval_coverage`

The command requires one rooted tree with complete branch lengths and at least
three observed taxa with numeric values for the requested trait. Bijux fits the
Brownian mean and diffusion rate on the observed taxa, predicts every missing
tree taxon from the conditional Brownian distribution, and then validates that
prediction path by refitting after holding out each observed taxon in turn when
enough observed taxa remain.

When `--summary-out` is supplied, `comparative trait-imputation` writes one
flat summary ledger as CSV or TSV. The row preserves:
- `trait`
- `taxon_column`
- `model`
- `tree_taxon_count`
- `observed_taxon_count`
- `imputed_taxon_count`
- `excluded_taxon_count`
- `root_state`
- `sigma_squared`
- `log_likelihood`
- `aic`
- `aicc`
- `holdout_validation_status`
- `holdout_count`
- `holdout_mean_absolute_error`
- `holdout_root_mean_squared_error`
- `holdout_interval_coverage`

When `--imputations-out` is supplied, the command also writes one imputed-value
ledger as CSV or TSV. Each row preserves:
- `taxon`
- `missing_reason`
- `observed_support_taxon_count`
- `predicted_value`
- `conditional_variance`
- `conditional_standard_error`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`

When `--holdout-out` is supplied, the command also writes one
leave-one-observed-out validation ledger as CSV or TSV. Each row preserves:
- `taxon`
- `observed_value`
- `predicted_value`
- `residual`
- `absolute_error`
- `conditional_variance`
- `conditional_standard_error`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`
- `covered_by_95_confidence_interval`
- `observed_support_taxon_count`
- `rank`

When `--excluded-taxa-out` is supplied, the command also writes one explicit
excluded-taxa ledger as CSV or TSV. Each row preserves:
- `taxon`
- `reason`

The reasons are explicit reviewer-facing states rather than generic failures:
- `non_numeric_trait_value`
- `absent_from_tree`

This surface exists so missing-value prediction preserves the Brownian fit,
per-taxon uncertainty intervals, and holdout evidence instead of reducing the
question to one opaque filled-in spreadsheet column.

`comparative multivariate` is the governed review surface for fitting the same
comparative predictor set across multiple response traits on one shared
complete-case taxon set. Its JSON metrics report:
- `response_count`
- `predictor_count`
- `analysis_taxa`
- `excluded_taxa`
- `residual_covariance_response_count`
- `residual_covariance_matrix_rank`
- `residual_covariance_condition_number`
- `residual_covariance_singular`
- `residual_covariance_near_singular`
- `response_model_count`
- `coefficient_row_count`
- `residual_covariance_row_count`
- `residual_correlation_row_count`
- `residual_association_count`
- `warning_count`

The command preserves:
- `missing_value_policy` with the governed shared complete-case rule
- `numerical_tolerance` with the governed multivariate comparison tolerance
- `response_models` with one fitted PGLS result per requested response
- `response_model_rows` with one explicit fit-summary row per response
- `coefficient_rows` with one explicit coefficient row per response-term pair
- `covariance_rows` with one residual covariance row per ordered response pair
- `covariance_diagnostics` with residual covariance matrix rank, condition number, and singular-versus-near-singular state
- `correlation_rows` with one residual correlation row per ordered response pair
- `association_rows` with one residual association row per unique response pair
- `excluded_taxa` with one explicit complete-case exclusion row per dropped taxon
- `warnings` with explicit weak-sample-size and singular-covariance warnings when present

This surface is explicit about missing values. A taxon is analyzed only when
every requested response and every requested predictor term can be evaluated.
Predictor terms are interpreted with the comparative formula parser used by
PGLS, so categorical predictors, transformed numeric predictors, and explicit
interaction terms remain governed instead of being treated as raw columns.
Taxa that are missing from the trait table, missing from the tree, or missing
one required value are preserved in the excluded-taxa report instead of
disappearing into a generic row-count difference.

When `--response-models-out` is supplied, `comparative multivariate` writes one
per-response model ledger as CSV or TSV. Each row preserves:
- `response`
- `formula`
- `predictor_term_count`
- `encoded_term_count`
- `taxon_count`
- `lambda_value`
- `log_likelihood`
- `residual_variance`
- `r_squared`
- `residual_degrees_of_freedom`

When `--coefficients-out` is supplied, the command writes one per-response
coefficient ledger as CSV or TSV. Each row preserves:
- `response`
- `formula`
- `term`
- `estimate`
- `standard_error`
- `test_statistic`
- `p_value`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`
- `degrees_of_freedom`
- `inference_distribution`

When `--covariance-out` is supplied, `comparative multivariate` writes one
residual covariance ledger as CSV or TSV. Each row preserves:
- `left_response`
- `right_response`
- `pair_count`
- `is_diagonal`
- `covariance`
- `correlation`

When `--correlation-out` is supplied, the command writes one residual
correlation ledger as CSV or TSV. Each row preserves:
- `left_response`
- `right_response`
- `pair_count`
- `is_diagonal`
- `correlation`

When `--associations-out` is supplied, the command also writes one residual
trait-association ledger as CSV or TSV. Each row preserves:
- `left_response`
- `right_response`
- `pair_count`
- `covariance`
- `correlation`
- `test_statistic`
- `p_value`
- `lower_95_confidence_interval`
- `upper_95_confidence_interval`

When `--excluded-taxa-out` is supplied, the command writes one explicit
excluded-taxa ledger with:
- `taxon`
- `reason`
- `missing_columns`
- `blocking_responses`
- `details`

This surface exists so correlated-evolution review can inspect which traits
still move together after the fitted predictors are accounted for, while also
making the shared complete-case taxon policy auditable.

`comparative report` is the governed integrated comparative-analysis review
surface. It fits one comparative regression and then preserves the formula,
coefficient evidence, residual diagnostics, phylogenetic signal, model
comparison, audit trail, and independent contrasts in one reviewer-facing
package. Its JSON metrics report:
- `taxon_count`
- `selected_model`
- `audit_row_count`
- `excluded_taxa`
- `limitation_count`
- `coefficient_count`
- `package_output_count`

The command preserves:
- `snapshot` with the integrated comparative fit, signal, contrasts, and model comparison
- `influence` with taxon and predictor influence review
- `limitations` as explicit reviewer-facing cautions

When `--out` is supplied, `comparative report` writes one standalone HTML
report on the existing comparative-method surface.

When `--out-dir` is supplied, the command writes one full review package
directory with:
- `comparative-report.html`
- `comparative-summary.tsv`
- `coefficient-table.tsv`
- `residual-summary.tsv`
- `signal-summary.tsv`
- `model-comparison.tsv`
- `interpretation-table.tsv`
- `audit-table.tsv`
- `contrast-table.tsv`
- `comparative-report.manifest.json`

The package tables preserve these reviewer-facing contracts:
- `comparative-summary.tsv`: `response`, `formula`, `predictor_count`, `analysis_taxa`, `selected_model`, and fit/signal summary
- `coefficient-table.tsv`: one row per coefficient with estimate, standard error, test statistic, p-value, and 95% interval
- `residual-summary.tsv`: one row per residual-diagnostic surface with variance, leverage, outlier taxa, and warnings
- `signal-summary.tsv`: one row with Blomberg's K, Pagel's lambda, likelihood values, and contrast count
- `model-comparison.tsv`: one row per Brownian or OU candidate fit with AIC and AICc evidence
- `interpretation-table.tsv`: one row per reviewer-facing claim with supporting evidence and caution text
- `audit-table.tsv`: one row per audit surface with taxa used, excluded taxa, assumptions, and warnings
- `contrast-table.tsv`: one row per internal node with the standardized contrast and ancestral value

This surface exists so comparative conclusions can be reviewed from one durable
artifact bundle instead of being reassembled manually from separate regression,
signal, contrast, and diagnostics commands.

`comparative clade-residuals` is the governed review surface for asking whether
one fitted comparative model leaves concentrated residual burden in particular
subtrees. Its JSON metrics report:
- `model_family`
- `taxon_count`
- `clade_count`
- `residual_heavy_clade_count`
- `top_influential_clade`
- `standardized_residual_method`

The command fits one comparative model first and then aggregates residuals
across every internal non-root clade in the analyzed tree. The model family is
inferred from the fitted response:
- binary `0/1` response => phylogenetic logistic residual surface
- otherwise => PGLS residual surface

When `--taxa-out` is supplied, `comparative clade-residuals` writes one taxon
ledger as CSV or TSV. Each row preserves:
- `taxon`
- `observed_value`
- `fitted_value`
- `residual`
- `standardized_residual`

When `--clades-out` is supplied, the command also writes one internal-clade
aggregation ledger as CSV or TSV. Each row preserves:
- `clade_id`
- `node_label`
- `taxon_count`
- `taxa`
- `mean_residual`
- `mean_abs_residual`
- `mean_standardized_residual`
- `mean_abs_standardized_residual`
- `max_abs_standardized_residual`
- `residual_sum_of_squares`
- `residual_sum_of_squares_share`
- `positive_residual_taxa`
- `negative_residual_taxa`
- `influence_score`
- `residual_heavy`
- `rank`

This surface exists so residual review does not stop at single outlier taxa.
It makes subtree-level model misspecification explicit and ranks the clades
that carry the largest residual burden.

`comparative clade-stability` is the governed review surface for asking whether
one comparative conclusion depends on one major subtree. Its JSON metrics
report:
- `model_family`
- `baseline_taxon_count`
- `baseline_term_count`
- `candidate_clade_count`
- `blocked_clade_count`
- `coefficient_change_row_count`
- `top_influential_clade`
- `major_clade_fraction`
- `minimum_major_clade_size`

The command fits one baseline comparative model first, derives major internal
non-root clades from the baseline analyzed tree, removes each candidate clade,
and attempts the same refit on the retained taxa. The model family is inferred
from the fitted response:
- binary `0/1` response => phylogenetic logistic stability surface
- otherwise => PGLS stability surface

When `--clades-out` is supplied, `comparative clade-stability` writes one
leave-one-clade-out summary ledger as CSV or TSV. Each row preserves:
- `clade_id`
- `node_label`
- `dropped_taxon_count`
- `dropped_taxa`
- `retained_taxon_count`
- `fit_status`
- `blocker`
- `baseline_term_count`
- `coefficient_comparison_count`
- `missing_baseline_term_count`
- `missing_baseline_terms`
- `sign_changed_term_count`
- `significance_changed_term_count`
- `max_abs_delta_estimate`
- `max_abs_delta_p_value`
- `delta_log_likelihood`
- `influence_score`
- `rank`

When `--terms-out` is supplied, the command also writes one coefficient-delta
ledger as CSV or TSV. Each row preserves:
- `clade_id`
- `node_label`
- `term`
- `baseline_estimate`
- `dropped_estimate`
- `delta_estimate`
- `baseline_p_value`
- `dropped_p_value`
- `delta_p_value`
- `baseline_significant`
- `dropped_significant`
- `sign_changed`
- `significance_changed`

This surface exists so subtree dependence is explicit at the model level
rather than inferred indirectly from residuals or taxon-level exclusion
screens. Blocked rows are preserved on purpose: if removing one candidate
clade leaves too few taxa or collapses a binary response to one class, the
review surface records that failure instead of pretending the clade was
uninfluential.

`comparative posterior-pgls` is the governed review surface for propagating
tree-set uncertainty into one continuous-trait PGLS conclusion. Its JSON
metrics report:
- `total_tree_count`
- `burnin_tree_count`
- `kept_tree_count`
- `analysis_taxon_count`
- `rooted_topology_count`
- `unrooted_topology_count`
- `tree_fit_row_count`
- `coefficient_row_count`
- `coefficient_summary_count`
- `stable_supported_term_count`
- `direction_conflict_term_count`
- `lambda_mode`
- `significance_threshold`

The command reads one posterior or bootstrap tree set, optionally discards a
leading burn-in fraction, reduces retained trees to their shared taxa, fits the
same PGLS specification on every retained tree, and then summarizes the
coefficient distribution across those fits. This is intentionally a
continuous-trait PGLS surface, not a generic comparative bucket:
- the response must satisfy the normal PGLS input contract
- the same formula is reused on every retained tree
- lambda is either estimated independently per retained tree or fixed at one
  user-supplied value

When `--trees-out` is supplied, `comparative posterior-pgls` writes one
per-tree fit ledger as CSV or TSV. Each row preserves:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `unrooted_topology_id`
- `lambda_value`
- `log_likelihood`

When `--coefficients-out` is supplied, the command also writes one per-tree
coefficient ledger as CSV or TSV. Each row preserves:
- `source_tree_index`
- `post_burnin_index`
- `rooted_topology_id`
- `term`
- `estimate`
- `p_value`
- `significant`
- `direction`

When `--summary-out` is supplied, the command also writes one
coefficient-distribution summary ledger as CSV or TSV. Each row preserves:
- `term`
- `tree_fit_count`
- `positive_tree_count`
- `negative_tree_count`
- `zero_tree_count`
- `dominant_direction`
- `direction_consistency`
- `significant_tree_count`
- `significance_fraction`
- `conclusion_stability`
- `mean_estimate`
- `median_estimate`
- `standard_deviation`
- `minimum_estimate`
- `maximum_estimate`
- `lower_95_empirical_estimate`
- `upper_95_empirical_estimate`
- `mean_p_value`
- `median_p_value`
- `minimum_p_value`
- `maximum_p_value`

This surface exists so coefficient support can be reviewed against the full
tree-set uncertainty instead of being inferred from one summary tree. The
`conclusion_stability` field is deliberately reviewer-facing rather than
opaque:
- `stable_supported`: every retained tree supports the same directional effect
- `stable_unsupported`: every retained tree keeps the same direction but none
  cross the support threshold
- `mixed_support`: one direction is retained, but support calls differ across
  trees
- `direction_conflict`: positive and negative estimates both occur across the
  retained tree set

`comparative brownian-pgls` is the fixed-covariance companion surface when the
scientific question specifically assumes Brownian shared-path covariance rather
than an estimated Pagel lambda. Its JSON metrics report:
- `covariance_model`
- `lambda_value`
- `covariance_row_count`
- `tree_is_ultrametric`
- `minimum_root_to_tip_depth`
- `maximum_root_to_tip_depth`
- `raw_log_determinant`
- `positive_definite_before_stabilization`

The command does not estimate or optimize lambda. It fixes lambda at `1.0`,
audits the raw Brownian covariance before stabilization, and fails explicitly if
zero or negative branch lengths make that covariance invalid.

`comparative ou-pgls` is the analogous regression surface for stationary-root
OU covariance. It accepts either a fixed positive `--alpha` or `--alpha
estimate`. Its JSON metrics report:
- `alpha`
- `alpha_estimation_mode`
- `alpha_profile_point_count`
- `alpha_lower_95_confidence_interval`
- `alpha_upper_95_confidence_interval`
- `covariance_model`
- `covariance_row_count`
- `log_likelihood`
- `aic`

This surface exists so OU-style residual covariance can be reviewed directly
instead of being approximated through Brownian or Pagel-lambda summaries. The
reported `aic` uses the fitted regression coefficient count plus one extra
parameter only when alpha was estimated rather than fixed.

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

The same `comparative pgls` result now also preserves the fitted Pagel lambda
surface directly under `data.model.lambda_fit`. That report makes the covariance
choice auditable instead of treating lambda as one opaque scalar:
- `mode` distinguishes fixed-lambda review from estimated-lambda review
- `lambda_value` and `log_likelihood` identify the chosen optimum or fixed
  covariance strength
- `lower_95_confidence_interval` and `upper_95_confidence_interval` report the
  likelihood-ratio-supported interval when lambda was estimated
- `profile_rows` preserves the full bounded likelihood profile used for review

When `--lambda-profile-out` is supplied, `comparative pgls` also writes that
profile as CSV or TSV. Its JSON metrics then report
`lambda_estimation_mode`, `lambda_profile_point_count`,
`lambda_lower_95_confidence_interval`, and
`lambda_upper_95_confidence_interval`.

When `--model-matrix-out` is supplied, `comparative pgls` writes the encoded
design matrix as CSV or TSV. Its JSON metrics then also report
`intercept_included`, `model_matrix_row_count`, and
`model_matrix_column_count`, while `data.inputs.model_matrix` preserves the
encoded column names and one taxon-level row per analyzed observation. This
surface exists so reviewers can inspect the actual fitted predictors instead of
inferring the design matrix indirectly from coefficient names.

When `--categorical-contrasts-out` is supplied, `comparative pgls` also writes
one categorical-contrast ledger. Its JSON metrics then report
`categorical_contrast_predictor_count` and
`categorical_contrast_row_count`, while `data.categorical_contrasts.rows`
preserves one row per reported group level.

The contrast ledger is explicit about interpretation:
- `encoding_scheme` distinguishes treatment-coded `reference-level` rows from
  `full-indicator` rows used by intercept-free formulas
- `is_reference_level` marks the baseline group directly instead of forcing
  reviewers to infer it from a missing coefficient
- `coefficient_name` is blank only for a true baseline row
- `missing_category_taxa` keeps blank or absent category assignments visible for
  that predictor

This surface exists so categorical coefficients are read as group contrasts with
clear baselines, not as unlabeled free-floating numbers.

When `--interaction-coefficients-out` is supplied, `comparative pgls` also
writes one interaction-coefficient ledger. Its JSON metrics then report
`interaction_term_count` and `interaction_coefficient_row_count`, while
`data.interaction_coefficients.rows` preserves one row per fitted interaction
coefficient.

The interaction ledger is explicit about effect modification:
- `interaction_kind` distinguishes continuous-by-continuous,
  continuous-by-categorical, and categorical-by-categorical effects
- `component_columns` preserves the exact encoded columns that generated each
  fitted interaction coefficient
- `component_levels` shows which categorical levels participate in a given row
- `omitted_reference_levels` keeps treatment-coded baseline groups visible when
  an interaction term omits them from the coefficient table

This surface exists so interaction coefficients are interpreted as explicit
effect-modification terms instead of opaque colon-delimited names.

When `--covariance-out` is supplied, `comparative brownian-pgls` writes one
pairwise Brownian covariance ledger as CSV or TSV. Each row preserves:
- `left_taxon`
- `right_taxon`
- `is_diagonal`
- `shared_path_length`
- `left_root_depth`
- `right_root_depth`

The written ledger also repeats the tree-level covariance audit fields on every
row so a single extracted table still preserves whether the fitted tree was
ultrametric, the root-depth range, the branch-length range, and whether the raw
covariance was positive definite before stabilization.

When `--covariance-out` is supplied, `comparative ou-pgls` writes the analogous
pairwise OU covariance ledger as CSV or TSV. Each row preserves:
- `left_taxon`
- `right_taxon`
- `is_diagonal`
- `covariance_value`
- `shared_path_length`
- `left_root_depth`
- `right_root_depth`

When `--alpha-profile-out` is supplied, `comparative ou-pgls` also writes the
bounded alpha likelihood profile as CSV or TSV. Each row preserves:
- `alpha_estimation_mode`
- `alpha`
- `log_likelihood`
- `delta_log_likelihood`
- `within_95_confidence_interval`

The written ledgers therefore expose both the fitted covariance surface and the
alpha-selection surface directly instead of requiring reviewers to infer them
from one terminal regression summary.

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

Use `topology distance-reference` when the question is whether the runtime
still matches the governed tree-distance evidence surface instead of only one
pair of user trees. That command reruns checked rooted RF, unrooted RF,
normalized RF, and branch-score cases from the repository fixture set and
reports:

- `case_count`
- `external_case_count`
- `policy_case_count`
- `all_passed`

The governed cases cover binary trees, rooting-only disagreement, topology
disagreement, polytomies, star-tree collapse, and shared-taxa pruning. The
policy cases keep `--taxon-overlap-policy require-identical` explicit for both
RF and branch-score comparisons, so a future regression cannot silently turn
those mismatch rejections into pruning behavior.

Use `topology support-reference` when the question is whether branch-support
parsing still attaches IQ-TREE, FastTree, and posterior support values to the
correct clades instead of only parsing one tree ad hoc. That command reruns the
checked support fixture set and reports:

- `case_count`
- `reference_case_count`
- `policy_case_count`
- `all_passed`

The governed cases cover plain IQ-TREE UFBoot labels, compound SH-aLRT/UFBoot
labels, FastTree local support, and posterior clade frequencies from a checked
tree set. The policy cases keep two review guarantees explicit: rotated trees
with the same clades must still compare by clade rather than node order, and
bootstrap-versus-posterior comparison must flag topology mismatch when support
appears on clades that the other uncertainty surface does not contain.

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

`tree-set diversity` is the compact topology-dispersion surface for one
posterior or bootstrap tree set. It reads the tree set iteratively, buckets the
pairwise rooted RF signal into one frequency ledger, and reports rooted
topology count, pair count, runtime, peak memory, and skipped malformed-tree
count without requiring the full pairwise matrix as the primary review surface.
When `--out` is supplied, it writes one `.tsv` RF-distribution ledger with one
row per distinct rooted RF bucket.

The command writes a stable artifact bundle under `--out-dir`:
- `.summary.tsv` with one row of tree-count, runtime, peak memory, skipped malformed-tree count, topology-diversity, threshold, and consensus summary fields
- `.consensus.nwk` with the consensus topology labeled by bootstrap support percentages
- `.clade-frequencies.tsv` with one row per informative clade frequency
- `.unstable-branches.tsv` with one row per non-robust consensus branch
- `.unstable-clades.tsv` with the broader conflicting-clade ledger across the full replicate set
- `.rf-distribution.tsv` with one row per rooted RF bucket across distinct tree pairs
- `.distance-matrix.tsv` and `.topology-clusters.tsv` for direct topology-variation review

The unstable-branch contract is explicit. A consensus branch is omitted from the
unstable-branch ledger only when its replicate frequency reaches the robust
threshold and no conflicting alternative clade is present. That keeps a
bootstrap summary from overstating a majority-rule consensus as if every branch
were equally stable.

Both `tree-set diversity` and `tree-set bootstrap-summary` skip malformed
line-oriented Newick tree records instead of aborting the full review surface.
Their JSON payloads and summary ledgers record `runtime_seconds`,
`peak_memory_bytes`, and `skipped_malformed_tree_count` so large-tree review is
measurable instead of implicit.

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

The governed external execution adapters also share one explicit execution
control contract. `adapter align`, `trim`, `model-select`, `infer-ml`,
`bootstrap`, `sh-alrt`, `fasta-to-tree`, `consensus`, `infer-fast`,
`infer-large`, `compare-engines`, `mrbayes-run`, and `beast-run` accept:

- `--timeout-seconds` for a wall-clock execution budget
- `--resume` to reuse only one verified completed run
- `--incomplete-run-policy reject|clean` to stop on or remove partial outputs from a failed, timed-out, killed, or malformed-output earlier run

The JSON payloads now expose the applied timeout budget and the resolved resume
status, so automation can distinguish a fresh execution from a verified reuse.
If the executable itself cannot be resolved, the command fails before any
incomplete-run marker is written because no engine run started.

Success on these governed adapter commands now means more than "the process
exited and a file appeared." MAFFT and trimAl must emit non-empty valid
alignments, IQ-TREE must emit the required `.iqtree`, `.log`, tree, model,
and support artifacts for the selected workflow, FastTree must emit a valid
tree with parseable local-support annotations, and BEAST or MrBayes must emit
their full required posterior artifact sets. Missing required files, empty
required files, missing IQ-TREE model results, and missing required support
annotations surface stable structured error codes before manifests or
review-facing reports are written.

For coding nucleotide inputs, `adapter align --codon-aware` is the supported
alignment entrypoint. It excludes frame-broken sequences and sequences with
ambiguous or invalid codons plus sequences with premature stop codons, aligns a
translated amino-acid guide, and back-translates guide gaps as nucleotide
triplets so the resulting alignment stays codon-safe for downstream inference
steps.

This surface accepts:

- `--sequence-type dna|rna` when the nucleotide alphabet must be forced
- `--genetic-code` with an NCBI code id or codon-table name

The codon-aware workflow writes reviewer-facing sidecars in addition to the
final codon alignment:

- one translated amino-acid guide input FASTA
- one aligned amino-acid guide FASTA
- one exclusion ledger
- one codon summary ledger

`alignment coding` and `alignment translate` also accept `--genetic-code`, so
the standalone coding diagnostics and amino-acid translation surfaces can use
the same explicit codon table as the codon-aware alignment workflow.

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
`taxon_count`, `character_count`, `calibration_count`, `tip_date_count`,
`warning_count`, `starting_tree_source`, `beast_data_type`,
`substitution_model`, `clock_model`, `tree_prior`, `chain_length`, and
`log_every`.

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

`adapter beast-xml` is the governed XML evidence surface for one prepared BEAST
analysis file. It parses and validates the written XML directly, then reports
the assumed substitution model, clock model, tree prior, starting-tree source,
chain length, calibration count, tip-date count, and logger outputs. Its JSON
metrics expose `valid`, `issue_count`, `taxon_count`, `character_count`,
`calibration_count`, `tip_date_count`, `chain_length`, and `logger_count`.

`adapter beast-run` is the governed execution surface for one prepared BEAST
XML analysis. It runs BEAST, validates the posterior log and posterior tree
outputs, preserves the execution manifest beside the XML, and exposes
`warning_count`, `threads`, `seed`, `overwrite`, `resumed`, and
`timeout_seconds` in JSON metrics.

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

When a sibling BEAST XML is present, downstream Bayesian reviewer-text
surfaces reuse that XML so model assumptions and chain settings are stated
from the prepared analysis itself instead of from placeholder CLI defaults.

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
copied snippets. Its JSON metrics now also expose `resumed` and
`timeout_seconds`, matching the shared execution-control contract.

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

Those IQ-TREE commands are now also strict about bundle completeness. `adapter
model-select` fails if the best-fit model or candidate table cannot be parsed;
`adapter infer-ml` fails if the tree, report, log, or model result is missing
or empty; `adapter bootstrap` fails if the bootstrap tree set is empty or the
supported tree does not contain parseable support labels; `adapter sh-alrt`
fails if joint SH-aLRT/UFBoot labels are missing; and `adapter consensus`
fails if the consensus tree lacks parseable support values.

`adapter fasta-to-tree` is the governed raw-FASTA workflow surface above those
direct IQ-TREE adapters. It keeps `.aln`, `.trimmed.aln`, `.tree`, `.log`,
`.model.tsv`, `.support.tsv`, `.manifest.json`, and `.run.json` sidecars in
one review bundle, while leaving the step-specific engine artifacts under
`engine-artifacts/`. Its public validation corpus compares those reviewer-facing
artifacts semantically rather than byte-for-byte so stable scientific results
are not rejected because of harmless path or timestamp differences.

Its JSON payload now reports the workflow as `supported`, with explicit
real-engine validation basis, so downstream reviewers can distinguish this
validated external-engine lane from advisory, experimental, or parser-only
surfaces.

Its composite manifest now also records `stage_fingerprints` for raw-input
validation, alignment, trimming, model selection, inference, support, and the
final reviewer-facing report. Those fingerprints explain the resolved resume
state directly: a stage is reused only when the recorded inputs, config,
command, and detected engine version still match the current run, and changed
upstream fingerprints invalidate the downstream stages automatically.

The same manifest is now the entrypoint for `phylo bundle`, which exports one
portable workflow-result directory. That bundle keeps the copied workflow
manifest, extracted config, bundle-local rerun ledger, reviewer-facing HTML
report, copied inputs when still available, final workflow outputs, and
declared step-level engine artifacts together. `phylo validate-bundle` then
checks both checksum integrity and the required workflow entries before the
bundle is treated as a valid handoff.

`comparative logistic` intentionally reports a different trust state. It is
`experimental`, emits a warning in JSON output, and repeats the
`phylogenetic-working-correlation-gee` approximation method in both the direct
metrics and the method-tier payload, together with an explicit
`ape::compar.gee` non-claim. Bayesian report commands such as
`adapter mrbayes-report` and `adapter beast-calibration-report` report
`parser-only` because they summarize external posterior artifacts rather than
claiming that Bijux executed the inference itself.

`phylo run` is the governed one-command workflow-config surface above those
manifest tools. It takes one YAML or JSON config file, validates it before
engine preflight begins, executes the canonical `fasta-to-tree` workflow, and
then exports one validated result bundle automatically. Its config contract
currently supports:

- one input FASTA
- optional metadata and traits tables
- external engine executable choices
- alignment and trimming settings
- inference seed, threads, and bootstrap replicates
- output directory and optional bundle directory
- timeout and incomplete-run policy controls

Its JSON metrics report:

- `workflow`
- `selected_workflow_status`
- `metadata_present`
- `traits_present`
- `alignment_mode`
- `trimming_mode`
- `bootstrap_replicates`
- `iqtree_seed`
- `iqtree_threads`
- `timeout_seconds`
- `bundle_file_count`
- `bundle_validation_passed`

The exported bundle now includes the resolved workflow config plus copied
config-source, metadata, and traits files when they were supplied. Those
auxiliary files are recorded honestly as reviewer context and downstream
comparative inputs; the current `fasta-to-tree` execution path does not use
metadata or traits during tree building itself.

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
the recorded workflow stages. It also honors the same
`--incomplete-run-policy reject|clean` control as the smaller external adapter
surfaces, so stale partial FastTree outputs can be rejected or cleaned before a
fresh rerun.

Use `phylo preflight` before any external-engine workflow when you need a
single governed compatibility scan instead of discovering missing tools halfway
through a run. It inspects MAFFT, trimAl, IQ-TREE, FastTree, MrBayes, and
BEAST, records the resolved executable path and native version text for each
engine, and classifies each engine as `tested`, `untested`, `unsupported`, or
`missing` against the package's current support policy.

That same command also publishes workflow readiness for the governed
engine-backed workflows. Each workflow row names its required engines and is
classified as `ready`, `caution`, or `blocked`, with explicit notes when a
missing or untested engine drives the result. When `--workflow` is supplied,
the command fails early for blocked workflows instead of leaving the user to
discover the missing tool inside a longer adapter run.

Its JSON contract intentionally exposes both `selected_workflow_status` and
`overall_status`. The first answers whether the chosen workflow is runnable in
the current environment, while the second reflects the health of the broader
external-engine inventory across all supported workflows.

`phylo replay` is the governed rerun surface for any workflow manifest emitted
by the external-engine adapters and composite workflows above. It reads the
recorded workflow identifier, input checksums, structured config, resolved
commands, engine versions, seeds, runtime metadata, and output checksums from
the manifest, reruns the same workflow into a replay directory, and reports
whether the new outputs are scientifically equivalent under workflow-specific
comparison rules.

The replay contract is intentionally strict about provenance drift. Changed
inputs stop the replay before any engine rerun and surface the
`manifest_replay_input_changed` error. Engine-version drift is reported in the
structured output rather than treated as an automatic failure, because a newer
or older executable can still reproduce the same scientific result. Replay
equivalence is semantic rather than byte-for-byte: tree workflows compare
topology and support, alignment workflows compare aligned records, model
selection compares the chosen model, and reproducibility workflows compare the
governed classification outcome.

`phylo bundle` and `phylo validate-bundle` sit next to replay on that same
manifest contract. `phylo bundle` writes one reviewer-facing result directory
from a governed workflow manifest, while `phylo validate-bundle` fails if the
bundle is missing its required report, workflow outputs, or declared step
manifests even when the remaining copied files still hash correctly.

Across those governed engine-backed commands, JSON error payloads now carry one
scientific failure block instead of only a bare exception surface. The details
include `failure_reason`, `scientific_explanation`, `likely_causes`,
`actionable_fixes`, and `evidence`, so a blocked workflow can distinguish:

- invalid FASTA record problems such as duplicate identifiers, illegal characters, empty records, or length outliers
- trimming failures that removed every retained site versus workflows that never wrote the trimmed artifact
- tree-inference failures that never wrote a tree versus tree-like outputs that are present but unparsable
- comparative taxon-linkage mismatches listing tree taxa missing from the trait table and extra trait-table taxa
- BEAST and MrBayes parser failures naming the missing file, header, sampled row, or posterior-tree block section

`adapter compare-engines` is the governed side-by-side inference mode for one
aligned matrix. It runs IQ-TREE model selection, IQ-TREE ultrafast bootstrap
support inference, and FastTree approximate inference on the same input, then
emits the two inferred trees, an HTML comparison report, a flat comparison
table, a shared-clade ledger, a conflicting-clade ledger, a support-weighted
conflict ledger, a clade-conclusion ledger, a stability-summary ledger, a
taxon-influence ledger when shared-taxon pruning can rank conflict drivers, and
a manifest in one command. Its JSON summary exposes the selected model,
shared-taxon count, Robinson-Foulds distance, shared-clade count,
conflicting-clade count, stable-clade count, unstable-clade count,
engine-specific-clade count, high-support conflict count,
low-support disagreement count, serious-conflict count, the shared timeout
budget, and whether any governed engine step was resumed.

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
`--sequence-type` before the workflow can continue. The raw validation path now
scans FASTA inputs linearly instead of building one full record list just to
answer preflight questions, and the higher-level alignment-quality surface
reuses one loaded matrix with an explicit warning when near-duplicate pairwise
review is skipped above the governed large-alignment threshold.
