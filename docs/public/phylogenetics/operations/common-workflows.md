---
title: Common Workflows
audience: public
type: how-to
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-16
---

# Common Workflows

Typical public workflows include:

- validate and inspect a tree before downstream use
- validate core numerical outputs against established reference tools
- build a full reviewer-facing tree report package from one supplied tree
- trim and inspect an alignment before reporting
- run one command from raw FASTA to a supported inference bundle
- run one supported inference workflow from one YAML or JSON config file
- assemble aligned loci into a concatenated supermatrix
- audit a concatenated multi-locus matrix before inference
- run a comparative model and capture JSON plus report artifacts
- build a full comparative analysis review package from one tree and trait table
- reconstruct ancestral regions and inspect transition evidence
- generate a review-ready figure package for a tree
- benchmark large owned workloads under governed stress tiers

The public workflow contract is that important outputs should be inspectable
after the command finishes.

Serious workflow JSON outputs now also publish one explicit method tier:

- `supported` for surfaces backed by reference parity or real-engine validation
- `experimental` for approximate exploratory surfaces that emit a warning
- `advisory` for review packages that do not perform inference
- `parser-only` for summaries of external-engine artifacts that do not claim Bijux ran the inference

When present, inspect `method_tier`, `method_inference_mode`,
`method_validation_basis`, and `method_approximation` in the JSON metrics.

The major reviewer-facing TSV and JSON artifacts are also schema-governed now.
Downstream scripts can anchor on stable headers and top-level keys for the
canonical model, support, clade, branch, comparative-traits,
comparative-summary, event, and manifest outputs instead of scraping ad hoc
examples.

When the same governed workflow needs to run inside a notebook, script, or
pipeline instead of a shell command, use the stable Python surface under
`bijux_phylogenetics.api`. The workflow entry points intentionally mirror the
serious CLI categories rather than exposing only low-level helpers.

```python
from pathlib import Path

from bijux_phylogenetics.api import (
    run_alignment_workflow,
    run_comparative_model_workflow,
    run_tree_comparison_workflow,
)

alignment = run_alignment_workflow(
    Path("dataset/sequences.fasta"),
    Path("artifacts/alignment.fasta"),
)

comparison = run_tree_comparison_workflow(
    Path("dataset/tree-a.nwk"),
    Path("dataset/tree-b.nwk"),
)

comparative = run_comparative_model_workflow(
    Path("dataset/tree.nwk"),
    Path("dataset/traits.tsv"),
    response="response",
    predictors=["predictor_one"],
    lambda_value=1.0,
)

alignment.write_tsv(Path("artifacts/alignment-workflow.tsv"))
comparison.write_json(Path("artifacts/tree-comparison.json"))
comparative.write_tsv(Path("artifacts/comparative-model.tsv"))
```

Those Python functions return typed workflow result objects that delegate the
underlying CLI-grade reports and add stable JSON plus TSV serialization where
appropriate, so downstream code can inspect one stable runtime contract instead
of parsing command output.

When the goal is to prove that the repository's core numerical methods still
match established external tools on small governed fixtures, use `parity`.
This workflow keeps the fixtures CI-sized by default and writes reviewer-facing
summary and observation ledgers when requested.

```bash
bijux-phylogenetics parity \
  --summary-out artifacts/reference-parity-summary.tsv \
  --observations-out artifacts/reference-parity-observations.tsv \
  --json
```

Use `--extended` when the review needs the governed primate comparative cases
and the larger posterior-tree bundle in addition to the CI-sized default suite.

The core suite covers:

- Robinson-Foulds distance
- branch-score distance
- PGLS
- Pagel's lambda
- Brownian trait evolution
- OU trait evolution
- phylogenetic independent contrasts
- Blomberg's K
- posterior clade frequencies
- consensus tree generation

The PGLS portion of the suite now covers four distinct trust questions:

- fixed Brownian covariance on a simple numeric regression
- treatment-coded categorical predictor encoding against R model-matrix behavior
- treatment-coded interaction encoding and coefficients against R output
- one governed estimated-lambda primate regression against R `ape` plus `nlme`

The live `phytools` lane now adds a separate honest Brownian-only check for
those fixed-lambda regression shapes through `phytools::pgls.SEy`. Installed
`phytools 2.5.2` does not export a general `phytools::pgls` function, so the
governed live lane stays on `pgls.SEy` with `lambda = 1.0`, while the broader
exact PGLS contract for estimated lambda and full coefficient parity remains
the checked-in `ape` plus `nlme` suite.

The observation ledger records the expected failure classification, overlap
policy, and shared-versus-exclusive taxa so reviewer-visible mismatches can be
sorted quickly into topology, branch-length, missing-taxa-policy, numerical,
or model-assumption causes.

Each observation records the input fixtures, reference tool and version,
expected output, tolerance, tolerance rationale, observed output, and the
governed mismatch class.

When the goal is to run the same governed fixtures against a live local R
installation instead of only checked-in reference outputs, use the live `ape`
harness on the same command surface:

```bash
bijux-phylogenetics parity \
  --reference-source ape-live \
  --ape-rscript-executable Rscript \
  --summary-out artifacts/ape-parity-summary.tsv \
  --observations-out artifacts/ape-parity-observations.tsv \
  --json
```

This lane runs the checked-in R parity runner, loads `ape`, writes
machine-readable JSON, TSV, and normalized Newick outputs, compares them
against the Bijux runtime, and records the R version, `ape` version, Bijux
version, Bijux commit, function name, input fixture, tolerance, pass or fail
state, and mismatch reason for every governed case. If `Rscript` or `ape` is
unavailable, the case is recorded explicitly as skipped and one small
reproducible artifact bundle is written for review.

The live lane now uses two governed shared fixture catalogs. Tree parity cases
resolve durable fixture ids from `shared_tree_fixture_catalog.json`, while DNA
parity cases resolve durable fixture ids from
`shared_dna_alignment_fixture_catalog.json`. The lane now also resolves one
governed simulation catalog from
`shared_tree_simulation_fixture_catalog.json` for `ape::rtree` and
`ape::rcoal` envelope checks. The DNA portion of the lane covers
shared lowercase, ambiguity, gap, missing-data, identical-sequence,
high-divergence, invariant, one-variable-site, and coding-translation
fixtures across `ape::as.DNAbin`, `ape::base.freq`, `ape::seg.sites`,
`ape::dist.dna`, and `ape::trans`. The coding-translation portion now also
covers ambiguous-codon, frame-truncation, and alternate-genetic-code
fixtures. The unequal-length DNA fixture is now a
governed live `ape::dist.dna` failure case so Bijux and `ape` both prove the
same DNA-distance stop condition explicitly. The owned distance surface now
accepts the ape-compatible `raw`, `jc69`, `k80`, `f81`, and `tn93` aliases,
keeps `p-distance`, `jukes-cantor`, `kimura-2-parameter`,
`felsenstein-81`, and `tamura-nei-93` as the canonical internal labels,
reports saturated JC69, K80, F81, and TN93 pairs explicitly as undefined or
infinite instead of flattening them into one generic missing value, writes
one `--parameters-out` TSV ledger for model base frequencies and coefficients,
and can write one `--components-out` TSV ledger with pairwise mismatch,
transition, transversion, AG-transition, CT-transition, ambiguity, and
saturation fields alongside the distance matrix. TN93 warns explicitly when
the resolved alignment composition omits a nucleotide instead of silently
falling back to a simpler model.
The distance side of the same governed lane now also covers `ape::nj` over one
shared analytical three-taxon matrix plus four-taxon ultrametric and
non-ultrametric matrices. The owned Bijux side now builds neighbor-joining
trees with one in-repo deterministic NJ algorithm that validates zero-diagonal
and nonnegative matrix assumptions, produces branch lengths, and resolves tied
joins by stable taxon ordering rather than delegating the NJ method through
Biopython. When one review workflow already holds one loaded
`GeneticDistanceMatrix`, the same owned core is also reusable directly through
`build_distance_tree_from_genetic_distance_matrix(...)` instead of requiring a
restart from one path-based alignment or matrix file.
The comparative side of the same governed lane now also uses
`shared_trait_table_fixture_catalog.json` for `ape::pic` over balanced rooted
ultrametric, pectinate rooted non-ultrametric, and six-taxon clean trait
fixtures. The owned `comparative contrasts` surface writes one
`independent-contrasts.tsv` ledger with stable ape-style `node_id` values,
left-versus-right descendant partitions, standardized contrasts, and expected
variances, while its JSON report preserves one explicit input audit for
ultrametric reporting and missing-value pruning. Missing trait values remain an
owned pruning policy surface, and negative branch lengths remain an owned hard
stop instead of being pushed through as live `ape::pic` parity cases.
That same workflow boundary now states explicitly that `bionj` is out of scope
for this round. Reviewers therefore see one governed supported-method set
(`neighbor-joining`, `upgma`) and one explicit `ape::bionj` exclusion instead
of an ambiguous missing feature.
The same owned DNA surfaces now also share one DNAbin-compatible nucleotide
matrix instead of reparsing FASTA differently in each workflow. That matrix
preserves taxon order and alignment length, normalizes case, keeps gaps,
ambiguity codes, and explicit missing states literal, writes FASTA back
without nucleotide-state loss, and rejects unsupported symbols explicitly. The
same matrix is also reusable directly from Python through
`load_dna_bin_alignment(...)`,
`compute_alignment_base_frequency_report_from_dna_bin_alignment(...)`,
`compute_alignment_segregating_site_report_from_dna_bin_alignment(...)`, and
`compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(...)`. The
same matrix now also feeds
`inspect_coding_alignment_from_dna_bin_alignment(...)` and
`translate_coding_alignment_from_dna_bin_alignment(...)`, so one serious
nucleotide review workflow can carry one loaded matrix through composition,
segregating-site, nucleotide-distance, aligned coding diagnostics, and
aligned translation inspection.
The same `parity` command now also exposes one governed
`--reference-source phytools-live` lane for real `phytools` execution through
the checked-in R runner. That initial registry is intentionally narrow for
goal 201: it currently covers `phytools::phylosig(method='lambda')`,
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
`phytools::fastAnc`, and `phytools::anc.ML` on governed
twenty-four-taxon comparative fixtures drawn from the shared `phytools`
comparative fixture catalog. The live lambda lane now includes one
non-ultrametric strong-signal case and one ultrametric weak-signal case, so
reviewers can see both a boundary-adjacent high-signal fit and a near-zero
lambda fit validated against real `phytools` execution. The live K lane now
also includes seeded strong-signal and weak-signal permutation cases so
reviewers can compare K, permutation p-value, and null-distribution summary
behavior under one governed replicate count. The live `fitMk` lane now also
includes clean binary, clean multistate, and missing-value-pruned ER cases
plus clean multistate and missing-value-pruned multistate SYM cases. It now
also includes clean and missing-value-pruned binary ARD cases at full
rate-row parity plus clean and missing-value-pruned multistate ARD cases at
summary parity when the owned optimizer reports weakly identified boundary
rates, so reviewers can compare flat-root log-likelihood, AIC, AICc,
excluded taxa, ER-versus-SYM-versus-ARD model identity, and directed-rate
evidence honestly against real `phytools`. The live
`make.simmap` lane now also includes clean binary, clean multistate, and
missing-value-pruned binary ER cases; clean multistate and
missing-value-pruned multistate SYM cases; and binary plus
missing-value-pruned binary ARD cases at one governed seed and one governed
replicate count of 128 maps per case. Governed multistate ARD cases stay on
summary-envelope parity only when weakly identified boundary rates make
row-level transition summaries unstable across optimizers. It compares summary
envelopes only: excluded taxa, total-transition-count mean plus interval,
transition-count summary rows, and time-in-state summary rows. It does not
claim exact stochastic-history identity with `phytools`. On the owned Bijux
side, `discrete-evolution stochastic-map` now also carries fitted-model audit
fields, writes one flat branch-segment ledger, writes one per-state
time-summary ledger, and writes one per-branch state-occupancy ledger in
addition to the JSON collection and transition summary.
The same owned surface now also supports `discrete-evolution count-maps` over
saved stochastic-map collections, writing one per-replicate count matrix, one
aggregate transition matrix, one per-branch directional transition table, and
one flat event ledger. The live `countSimmap` lane now covers clean binary,
clean multistate, clean multistate SYM, and missing-value-pruned binary cases.
It compares total-transition envelopes plus directional transition-count rows,
including zero diagonal state pairs, without claiming exact stochastic-history
identity.
The same owned surface now also supports `discrete-evolution density-maps` over
saved stochastic-map collections, writing one branch-probability table, one
branch-level density envelope, one slice-level probability table at the
requested resolution, and one report-ready HTML or SVG artifact. The live
`densityMap` lane stays intentionally narrow: it currently covers governed
binary ER collections only, including one missing-value-pruned case. It
compares per-branch posterior probability summaries and branch-level
uncertainty against live `phytools::densityMap`, and does not claim
pixel-perfect plotting parity.
The live `describe.simmap` lane now also covers clean binary, clean multistate,
clean multistate SYM, and missing-value-pruned binary cases. It compares the
owned stochastic-map summary contract directly, including total changes,
transition rows, time-in-state rows, and per-branch state-occupancy rows.
The same owned surface now also supports `simulate history-discrete` for
fixed-tree discrete-history truth generation from one explicit rate matrix. It
writes one tip-state truth table, one node-state truth table, one branch
history table, one transition-event ledger, one branch-segment ledger, and
one parity summary table. The live `sim.history` lane now covers governed
binary and multistate no-change plus high-rate fixed-tree cases, and compares
total-transition summaries, transition-count rows, time-in-state rows, and
tip-state-frequency rows against real `phytools::sim.history` without
claiming exact cross-language history identity.
The same simulation surface now also supports one Brownian replicate-review
contract over tip distributions and tip covariances. One-off Brownian runs now
accept either `sigma` or explicit `sigma_squared`, keep the resolved Brownian
rate parameter in the report, and share the same covariance contract as the
replicate collection surface. The live `fastBM` lane now covers governed
low-variance, root-shift high-variance, and six-taxon fixed trees, and compares
summary envelopes plus tip-covariance rows against real `phytools::fastBM`
without claiming exact cross-language draws.
The same simulation surface now also supports one correlated Brownian
collection contract for two or more continuous traits on one fixed tree from
one explicit evolutionary covariance matrix. The owned CLI accepts either one
covariance matrix directly or one correlation matrix plus per-trait standard
deviations, keeps the generating covariance contract in the report, and can
write multivariate summary rows over root states, evolutionary covariance, tip
distributions, and tip covariances. The live `phytools::sim.corrs` lane now
covers governed low-correlation, negative-correlation root-shift, and
three-trait six-taxon fixed-tree cases, comparing summary envelopes,
tip-covariance rows, and tip-correlation rows against real `phytools::sim.corrs`
without claiming exact cross-language draws.
The same governed live `phytools` lane now also covers fixed-lambda Brownian
PGLS through `phytools::pgls.SEy` on one simple numeric regression plus one
categorical and one interaction-coded regression. That boundary is explicit:
installed `phytools 2.5.2` does not export a general `phytools::pgls`
function, so the live lane proves only `pgls.SEy` with `lambda = 1.0`, while
the broader exact PGLS contract for estimated lambda and full coefficient
parity remains the checked-in `ape` plus `nlme` suite.
The live
`rerootingMethod` lane now also covers governed ER binary, ER multistate, ER
missing-value-pruned, SYM multistate, and SYM missing-value-pruned cases. It
compares one node-probability table keyed by stable node signature and state
label against real `phytools` output. That governed claim stays narrow on
purpose: only ER or SYM with the equal root prior inherited from `fitMk` are
treated as live rerooting parity, while ARD, Fitch, ordered-state,
empirical-root-prior, and fixed-root-prior ancestral runs remain owned Bijux
review surfaces. The live
`fastAnc` lane now includes ultrametric strong-signal, ultrametric
weak-signal, non-ultrametric strong-signal, and missing-value pruning cases so
reviewers can compare node-signature estimates and standard errors directly
against real `phytools`. The live `anc.ML` lane now covers the same four
fixture shapes so reviewers can compare node-signature estimates, standard
errors, 95% intervals, Brownian log-likelihood, and fitted sigma-squared
directly against real `phytools`. It writes the same
summary-versus-observation TSV artifacts plus reproducible failure bundles as
the live `ape` lane.
The same baseline surface now also exposes
`simulate_random_tree(...)` and `simulate_coalescent_tree(...)` for one-tree
native simulation review, so callers that need one governed random or
coalescent draw plus its summary record do not have to route through one
batch-only simulation wrapper.
The owned `alignment composition` surface now also exposes
`--base-frequency-out` for one combined alignment-plus-sequence literal-state
frequency ledger that mirrors `ape::base.freq`. Lowercase, ambiguity-bearing,
missing-data, and all-gap-or-missing alignments all stay on the governed live
parity path, and the all-gap or missing edge case now warns explicitly instead
of fabricating canonical nucleotide content.
The owned `alignment segregating-sites` surface now also exposes
`--site-table-out` for one reviewer-facing segregating-site ledger with site
positions plus literal and ape-normalized state summaries. That same surface
follows live `ape::seg.sites` by normalizing leading and trailing gaps to `N`,
keeping explicit missing states from creating segregating sites by themselves,
and still treating incompatible ambiguity states or internal gaps as real
segregating-site evidence.
The owned `alignment translate` surface now also exposes
`--codon-validation-out` and `--excluded-sequences-out` for one amino-acid
translation review bundle. That aligned-translation surface now follows live
`ape::trans` by truncating trailing partial codons with an explicit warning,
while the stricter `prepare_coding_sequences_for_alignment` gate still owns
the serious workflow exclusion policy for frame errors, ambiguous codons, and
premature stop codons before codon-aware alignment.
The `ape::read.tree` portion of the same lane now validates structured clade
rows for rooted and unrooted trees, branch lengths, internal labels, support
labels, quoted labels, one governed multiple-tree Newick input, and one
governed malformed-Newick rejection case. Those cases now run through one
owned native Newick parser and writer on top of `PhyloTree`, including
location-aware parse failures, instead of depending on an external tree reader.
The owned tree-set side of the same workflow now also reads one native
`PhyloTree` per Newick record for `tree-set inspect`, `tree-set consensus`,
`tree-set support-map`, posterior diversity review, topology clustering, and
posterior tree-set comparison. Plain `.nwk` and plain `.trees` inputs are both
accepted when the content is one Newick record per tree. Strict consensus and
support workflows fail explicitly when the tree set does not share one exact
taxon set, while tolerant inspection workflows keep one reviewer-visible
malformed-record counter instead of silently dropping parse failures.
The `ape::root` portion now uses the
same shared tree catalog for single-tip outgroups, monophyletic multi-tip
outgroups, already-rooted trees, missing outgroups, and non-monophyletic
outgroups, and compares rooted clades plus branch lengths against live
`ape::root`. Bijux now rejects ambiguous non-monophyletic outgroups explicitly
instead of only warning about them. The `ape::unroot` portion now covers
rooted binary trees, post-outgroup-rooting trees, already-unrooted inputs, and
malformed parse failures, and it follows the same explicit branch-length
contract as live `ape::unroot`: the removed root-edge length is merged into
the retained sibling branch rather than folded into the expanded clade. The
owned tree-manipulation surface for rooting, unrooting, keep-tip pruning,
drop-tip pruning, clade extraction, MRCA lookup, and monophyly review now
runs directly on `PhyloTree` instead of splitting those operations across
native and Biopython-backed tree objects. The
`ape::consensus` portion now covers majority-rule and strict consensus over
governed conflicting and posterior-style tree sets, writes one normalized
consensus Newick plus one clade-frequency TSV ledger per case, and fails
explicitly when the tree set does not share one exact taxon set. The
`ape::prop.clades` portion now covers reference-tree clade support mapping
over duplicate, reordered, posterior-style, and mismatched shared tree sets.
Use `tree-set support-map` when one review tree needs one governed
`reference-tree-support.tsv` ledger keyed by descendant tip set. That surface
also preserves the real `ape` edge case explicitly: unsupported root-adjacent
splits stay unscored instead of being forced into zero-support rows. The
`ape::drop.tip` portion now covers rooted and unrooted taxon exclusion cases,
unknown excluded taxon names, and rootedness changes after pruning, while
keeping one explicit workflow safety rule: Bijux stops the run if pruning would
leave fewer than two retained taxa instead of carrying a one-tip tree forward.
The `ape::keep.tip` portion now covers valid keep-set pruning for rooted and
unrooted trees, selected-tip order differences, and rootedness changes after
pruning, while keeping two explicit workflow-side rules outside the live
parity subset: absent requested taxa are still reported for tree and trait
matching work, and fewer than two retained taxa still stop the run clearly.
The `ape::extract.clade` portion now covers rooted root-clade and internal-node
subtree extraction plus clear tip-node and out-of-bounds failures. Bijux also
keeps one owned selector next to that live parity lane: the same subtree can
be resolved by exact descendant taxa when a workflow has stable taxon identity
but not a durable ape-style node number.
The `ape::getMRCA` portion now covers stable internal-node identity for
two-tip, many-tip, full-tip-set, duplicate-tip, rooted-polytomy, and
already-rooted-outgroup cases. Bijux keeps one explicit workflow-side
difference next to that live lane: missing requested taxa fail clearly instead
of surfacing as a low-level parser-side condition.
The `ape::is.monophyletic` portion now covers rooted and unrooted monophyly
calls with explicit reroot policy, full-tip-set behavior, singleton and
mixed-missing requests, rooted-polytomy behavior, post-rooting behavior, and
all-missing reroot failures. Bijux also records the matched MRCA node and any
extra descendant taxa that make a direct clade fail monophyly, so the same
lane is usable for review instead of only for pass-fail checks.
The `ape::cophenetic.phylo` portion now covers rooted and unrooted
branch-length trees, compares one governed long-form tip-distance ledger
instead of only a printed matrix, and keeps the taxon order explicit in the
review payload. Bijux also exposes the same owned tip-distance surface
directly, and it rejects missing branch lengths unless the caller explicitly
opts into a unit-length fallback policy.
The `ape::dist.topo` portion now covers identical rooted trees, rooted
child-order rotations, one-conflict rooted pairs, rooted tree-versus-polytomy
pairs, one governed unrooted split conflict, and one governed 128-tip rooted
pair. It compares one explicit RF-style split ledger rather than only a
scalar distance, keeps rooted-versus-unrooted policy explicit per case, and
aligns directly with the owned `adapter compare --split-table-out` review
surface. The owned tree-distance, support-comparison, tree-set support, and
posterior clade-summary workflows now all consume one native clade-set core
for canonical split identities instead of keeping separate rooted-clade and
unrooted-split helper paths.
The `ape::vcv.phylo` portion now covers rooted ultrametric, rooted
non-ultrametric, unrooted branch-length, and singular zero-branch trees. It
compares one governed long-form Brownian shared-ancestry covariance ledger,
persists the compared covariance tables automatically when parity fails, and
keeps the taxon order explicit in the review payload. Bijux also exposes the
same owned Brownian covariance surface directly through
`summarize_brownian_covariance(...)`, and that runtime rejects missing or
negative branch lengths while reporting singular-versus-near-singular state on
the raw matrix instead of silently regularizing it away.
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
The `ape::gammaStat` portion now covers the same governed rooted ultrametric
fixture family, including internal-node-label and zero-internal-branch cases.
It compares one governed one-row `gamma-statistic.tsv` ledger, and the owned
Bijux surface `compute_diversification_gamma_statistic(...)` stays explicit
about two workflow boundaries instead of inheriting live `ape`'s looser edge
behavior: only fully bifurcating trees are accepted, and incomplete sampling
is surfaced as a warning rather than silently folded into the statistic.
The `ape::is.ultrametric` portion now covers exact ultrametric,
near-ultrametric, tight-tolerance near-ultrametric, and clearly
non-ultrametric trees. It compares one governed tip-depth diagnostic table,
and the owned Bijux surface `assess_tree_ultrametricity(...)` reports the
criterion name, criterion value, tolerance, maximum tip-depth deviation,
offending taxa, and a deterministic `ultrametric-diagnostics.tsv` ledger.
That same ape-style surface is now reused before rooted Brownian, OU, and
diversification workflows claim time-tree compatibility.
The
`ape::write.tree` portion now
roundtrips Bijux-written Newick through live `ape` for rooted, unrooted,
internal-label, support-label, quoted-label, and multiple-tree cases. Bijux
rejects unnamed tips, empty tree sets, and non-finite branch lengths before
that live comparison so malformed internal trees cannot be written silently.
Those tree checks are structural rather than text-based, so reordered
equivalent children pass while rootedness, tip-set, clade or split,
branch-length, and internal-label drift are reported explicitly.

When the goal is to check resource behavior across the repository's largest
owned workload families, use `benchmark stress-suite`. This workflow does not
invent a separate benchmark harness. It drives the owned large-alignment,
supermatrix, tree-set, comparative, and table-generation surfaces and records
their runtime, peak memory, input size, and output row counts in one JSON
report.

```bash
bijux-phylogenetics benchmark stress-suite \
  --tier small \
  --json
```

The `small` tier is intended for routine checks. The `heavy` tier raises the
same workload families to `1,000+` sequences, `1,000+` sampled trees, and
large reviewer-facing table sizes. Each observation reports:

- input size in bytes
- sequence count when applicable
- alignment length when applicable
- tree count when applicable
- taxon count when applicable
- locus count when applicable
- runtime
- peak memory
- memory observation kind
- output row count

This benchmark surface is intentionally reviewer-facing rather than a hidden CI
only lane. Users can run it directly and inspect the same governed workload
classes that back the repository stress tiers.

When the goal is to produce one honest pre-release artifact from actual test
and workflow outputs, use `report release-truth`. This workflow consumes pytest
JUnit XML reports for the full test lane and the real-engine lane, reruns the
owned workflow-validation, release-gate, parity, and stress-suite checks, and
summarizes supported workflows, experimental workflows, flagship datasets, and
known limitations in one HTML report plus one machine manifest.

```bash
bijux-phylogenetics report release-truth \
  --test-report artifacts/pytest/full-suite.xml \
  --real-engine-test-report artifacts/pytest/real-engine.xml \
  --out artifacts/release-truth-report.html \
  --json
```

When the goal is to turn one supplied or inferred tree into a review-ready
artifact bundle instead of only a diagnostic HTML page, use
`report tree-package`. This workflow keeps the tree itself as the single input
surface and materializes one coherent report directory with an image plus TSV
ledgers.

```bash
bijux-phylogenetics report tree-package example_tree.nwk \
  --out-dir artifacts/tree-report-package \
  --json
```

The package writes:

- one HTML report with the rendered tree image embedded
- one SVG tree image rendered on the owned phylogram surface
- one branch-support ledger
- one clade ledger
- one branch-statistics ledger
- one machine-readable manifest

This workflow is intentionally different from the older `report tree` audit
surface. `report tree` remains the lightweight structural and forensic HTML
diagnostic. `report tree-package` is the richer reviewer-facing output used
when users need the image, support table, clade table, and branch-length
statistics together in one durable bundle. Its JSON and HTML outputs classify
the package as `advisory`, because it audits one supplied tree rather than
inferring a new one.

When the goal is to start from a real packaged mammal dataset rather than
assemble a tree and trait table by hand, use `demo primate-comparative`. This
workflow materializes the packaged primate dataset, which is the first
repository-owned mammal comparative dataset shipped with the runtime, and then
recomputes the governed reference outputs for PGLS, Brownian motion, OU,
phylogenetic signal, and ancestral reconstruction.

```bash
bijux-phylogenetics demo primate-comparative \
  --out artifacts/primate-comparative \
  --json
```

The `dataset/` directory contains the tree, trait table, and checked reference
outputs that ship with the package. The `workflow/` directory contains a fresh
rerun of the same comparative surfaces over those packaged inputs:

- PGLS with `longevity ~ social_group_size`
- Brownian trait evolution for `longevity`
- OU trait evolution for `longevity`
- phylogenetic signal for `longevity`
- continuous ancestral reconstruction for `longevity`
- discrete ancestral reconstruction for `mating_system`

This keeps the real mammal dataset usable in two directions. Users can run the
full packaged workflow immediately, or they can point individual comparative
commands at `dataset/tree.nwk` and `dataset/traits.csv` with
`--taxon-column species`.

When the goal is to start from a real packaged bird dataset focused on trait
evolution and clade-pattern review, use `demo avian-reproductive-traits`. This
workflow materializes the packaged avian reproductive dataset, which ships with
one rooted 94-taxon bird tree and one cleaned reproductive trait table, and
then recomputes the governed workflow bundle for regression, trait evolution,
ancestral reconstruction, and clade summaries.

```bash
bijux-phylogenetics demo avian-reproductive-traits \
  --out artifacts/avian-reproductive-traits \
  --json
```

The packaged bird workflow writes the same dataset and workflow split as the
primate surface, but uses bird-specific review choices:

- PGLS with `testes_mass ~ body_mass`
- Brownian trait evolution for `testes_mass`
- OU trait evolution for `testes_mass`
- phylogenetic signal for `testes_mass`
- continuous ancestral reconstruction for `testes_mass`
- discrete ancestral reconstruction for `mating_system`
- clade-specific trait summaries for `mating_system`

This makes the packaged bird dataset usable for both trait-evolution examples
and clade-pattern review without requiring users to assemble a bird tree and
trait matrix by hand.

When the goal is to start from a real packaged plant dataset rather than an
animal comparative example, use `demo central-european-seashore-flora`. This
workflow materializes the packaged Central European plant dataset, which ships
with a rooted saltwater-and-seashore flora subset from a published phylogeny
and trait matrix, and then recomputes the governed workflow bundle for
regression, trait evolution, ancestral reconstruction, and clade summaries.

```bash
bijux-phylogenetics demo central-european-seashore-flora \
  --out artifacts/central-european-seashore-flora \
  --json
```

The packaged plant workflow keeps the source subset rule explicit and writes
the same dataset/workflow split as the other packaged comparative surfaces,
using the following governed review choices:

- PGLS with `seed_mass ~ plant_height`
- Brownian trait evolution for `seed_mass`
- OU trait evolution for `seed_mass`
- phylogenetic signal for `seed_mass`
- continuous ancestral reconstruction for `seed_mass`
- discrete ancestral reconstruction for `lifeform`
- clade-specific trait summaries for `lifeform`

This gives users a real non-animal comparative surface without requiring them
to assemble or subset the published flora tree and trait matrix by hand.

When the goal is to start from a real packaged viral FASTA panel rather than a
checked example alignment, use `demo influenza-a-ha-reference-panel`. This
workflow materializes the packaged Influenza A hemagglutinin dataset and reruns
the owned raw-sequence workflow for alignment, trimming, maximum-likelihood
inference, and bootstrap support review.

```bash
bijux-phylogenetics demo influenza-a-ha-reference-panel \
  --out artifacts/influenza-a-ha-reference-panel \
  --json
```

This packaged viral workflow requires executable access to MAFFT, trimAl, and
IQ-TREE. The governed run keeps the source accession panel explicit and uses
the following deterministic inference controls:

- sequence type `dna`
- MAFFT alignment on the raw FASTA panel
- trimAl trimming on the aligned FASTA
- IQ-TREE model selection and maximum-likelihood inference
- IQ-TREE bootstrap support with `1000` replicates
- IQ-TREE random seed `1`
- IQ-TREE threads `1`

This gives users a real viral sequence-to-tree surface without requiring them
to assemble accession panels or wire the external engine chain by hand.

When the goal is to exercise the same owned sequence-to-tree runtime under
short degraded inputs and explicit missing-data burden, use
`demo pleistocene-bear-cytb-fragments`. This workflow materializes the
packaged ancient-DNA-style bear cytochrome b panel and reruns the owned raw
FASTA workflow with explicit missingness-review outputs.

```bash
bijux-phylogenetics demo pleistocene-bear-cytb-fragments \
  --out artifacts/pleistocene-bear-cytb-fragments \
  --json
```

This packaged ancient-DNA-style workflow also requires executable access to
MAFFT, trimAl, and IQ-TREE. The dataset keeps its provenance and degradation
rule explicit by pairing three modern bear cytochrome b references with two
real cave-bear cytochrome b fragments shortened and interrupted by internal
`N` blocks. The governed run makes missingness effects reviewable through:

- one raw unaligned FASTA panel
- one aligned output
- one trimAl-trimmed alignment
- one missingness-cleaned alignment with explicit cleanup thresholds
- one missingness-effects ledger per sequence
- one supported tree, selected-model table, and support table

This gives users a real degraded-sequence stress surface without requiring
them to assemble ancient and modern accession panels or invent their own
missingness reporting convention.

When the goal is to review a real multi-locus assembly and partitioned
inference surface instead of a single-gene example, use
`demo catarrhine-mitogenome-five-locus-panel`. This workflow materializes the
packaged catarrhine mitochondrial panel, which ships with five aligned coding
genes over the same six taxa, and then reruns the owned concatenation,
occupancy, model-selection, and partitioned-support surfaces.

```bash
bijux-phylogenetics demo catarrhine-mitogenome-five-locus-panel \
  --out artifacts/catarrhine-mitogenome-five-locus-panel \
  --json
```

The packaged multi-locus workflow writes the same dataset/workflow split as
the other demos, but focuses on reviewer-facing supermatrix assembly and
partitioned inference:

- concatenated supermatrix output
- explicit locus partition file
- per-taxon, per-locus, and matrix occupancy ledgers
- partition summary ledger
- model-candidate ledger from partitioned model selection
- supported tree and branch-support ledger from partitioned bootstrap review

The packaged locus set is intentionally small and explicit: five
mitochondrial coding genes extracted from stable RefSeq mitochondrial genomes.
That keeps the surface honest about scale while still exercising the full
multi-locus workflow contract.

When the goal is to inspect how the package handles messy comparative inputs
instead of assuming clean data, use
`demo catarrhine-data-quality-stress-panel`. This workflow materializes the
packaged catarrhine stress panel and reruns the owned audit-and-cleanup
workflow over its raw alignment, tree, and traits surface.

```bash
bijux-phylogenetics demo catarrhine-data-quality-stress-panel \
  --out artifacts/catarrhine-data-quality-stress-panel \
  --json
```

The packaged stress workflow keeps the messiness explicit rather than hiding
it in preprocessing:

- one raw aligned FASTA with a deliberate sequence composition outlier
- one raw FASTA validation surface with duplicate identifiers, illegal characters, empty sequences, and explicit length outliers
- one raw coding FASTA surface with one frame error and one premature stop codon
- one raw tree with a zero-length branch, one negative branch, and one extreme terminal branch
- one raw trait table with a duplicate taxon row and missing values
- one raw trait-linkage mismatch table with one missing tree taxon and one extra trait taxon
- raw-sequence findings and repair ledgers
- coding-sequence exclusion ledger
- raw trait-linkage mismatch ledger
- duplicate-trait and missing-trait ledgers
- sequence-outlier and tree-issue ledgers
- explicit repair-actions ledger
- one repaired raw FASTA subset, one prepared coding-sequence subset, and one cleaned alignment, cleaned tree, and cleaned trait table
- cleaned linkage and cleaned validation ledgers for the resolved subset

This gives users a governed dirty-data review surface that proves how
`bijux-phylogenetics` identifies and handles realistic pathologies without
pretending the raw inputs were already analysis-ready.

When the goal is to check recovery behavior against stored truth instead of
against only other inferred outputs, use `demo known-answer-reference-panel`.
This workflow materializes the packaged deterministic simulation panel and
reruns the owned recovery surfaces over its true tree, simulated alignment,
and simulated traits.

```bash
bijux-phylogenetics demo known-answer-reference-panel \
  --out artifacts/known-answer-reference-panel \
  --json
```

The packaged known-answer workflow keeps the truth contract explicit:

- one true birth-death tree
- one simulated DNA alignment generated on that tree
- one simulated Brownian continuous trait
- one simulated OU continuous trait
- one simulated three-state discrete trait
- one simulated host-association trait with stored branch-change truth
- one simulated geographic trait with stored branch-change truth
- one true-parameter ledger with the exact simulation seeds and values
- one true continuous-node ledger and one true OU-node ledger
- one true discrete-node ledger, one true host-node ledger, and one true geographic-node ledger
- one host-switch event ledger, one geographic transition-event ledger, and one declared recovery-threshold ledger

The governed recovery bundle then exposes how well the owned runtime recovers
those known answers:

- neighbor-joining distance-tree recovery against the true tree
- Brownian parameter recovery on the true tree
- OU parameter recovery on the true tree
- continuous ancestral reconstruction against stored node truth
- discrete ancestral reconstruction against stored node truth
- host-switch and geographic-state recovery against stored node and branch-event truth
- explicit pass or fail evaluation against declared thresholds instead of only a narrative claim

This gives users one durable internal reference surface for checking that
topology and ancestral-state recovery remain interpretable when the answers are
actually known in advance.

When the goal is to check whether Brownian, OU, and early-burst trait models
recover known generating parameters instead of only fitting one observed trait
table, use `demo continuous-mode-recovery-panel`. This workflow materializes
the packaged deterministic recovery panel and reruns the owned simulation plus
refit bundle over one shared rooted tree and four governed scenarios.

```bash
bijux-phylogenetics demo continuous-mode-recovery-panel \
  --out artifacts/continuous-mode-recovery-panel \
  --json
```

The packaged recovery workflow keeps the trust contract explicit:

- one shared rooted reference tree
- one Brownian recovery case judged on sigma-squared recovery
- one strong OU recovery case judged on alpha, sigma-squared, and optimum recovery
- one strong early-burst recovery case judged on rate-change recovery
- one weak OU case judged on identifiability warnings and Brownian-like support

The governed recovery bundle then exposes:

- one case summary ledger with expected-versus-selected best model
- one parameter-recovery ledger with absolute and relative error plus declared tolerances
- one BM-versus-OU-versus-early-burst model-choice ledger
- one warning ledger for OU and early-burst identifiability review
- one simulated trait table per case so users can rerun the same recovery checks directly

When the goal is to generate a standalone early-burst trait table for custom
downstream review, use `simulate traits-early-burst`.

```bash
bijux-phylogenetics simulate traits-early-burst tree.nwk \
  --root-state 1.0 \
  --sigma 0.5 \
  --rate-change 4.0 \
  --out artifacts/simulated-early-burst.tsv \
  --json
```

This simulator is the owned generating surface used by the packaged recovery
panel. It keeps the declared `rate_change` explicit in JSON output so recovery
reports can be tied back to known truth without reconstructing provenance from
notes or filenames.

When the goal is to review host-state evolution on a real pathogen panel, use
`demo rabies-cross-host-panel`. This workflow materializes the packaged rabies
dataset, which ships with a rooted nucleoprotein tree, raw sequences, and host
metadata, and then reruns the owned host-switching review surface over the
grouped host trait.

```bash
bijux-phylogenetics demo rabies-cross-host-panel \
  --out artifacts/rabies-cross-host-panel \
  --json
```

The packaged pathogen workflow writes the same dataset/workflow split as the
other demos, but focuses on ancestral host reconstruction and branchwise
transition evidence:

- host-switch summary ledger
- internal-node host-state probability ledger
- branchwise host-switch ledger
- directed host-switch count ledger
- host-transition fit ledger
- unsupported-claim and exclusion ledgers

The packaged metadata carries both exact `host_species` labels and the grouped
`host_group` trait used by the governed workflow. That grouping is explicit and
intentional: it keeps the compact rabies panel interpretable for host-switch
review instead of overstating species-level certainty from a tiny demo panel.

When the goal is to review geographic state evolution on a real pathogen
panel, use `demo rabies-geographic-transition-panel`. This workflow
materializes the packaged rabies geography dataset, which ships with a rooted
nucleoprotein tree, raw sequences, and grouped macroregion metadata, and then
reruns the owned biogeography review surfaces over that grouped region trait.

```bash
bijux-phylogenetics demo rabies-geographic-transition-panel \
  --out artifacts/rabies-geographic-transition-panel \
  --json
```

The packaged geography workflow writes the same dataset/workflow split as the
other demos, but focuses on ancestral region reconstruction and branchwise
movement evidence:

- geographic-state summary ledger
- internal-node region probability ledger
- geographic transition-rate ledger
- branchwise geographic transition-event ledger
- geographic migration summary and event ledgers
- geographic exclusion ledgers

The packaged metadata carries both raw `country` provenance and the grouped
`region_group` trait used by the governed workflow. That grouping is explicit
and intentional: it keeps the compact rabies panel interpretable for geography
review instead of overstating locality-level movement certainty from a tiny
demo panel.

When the goal is to test whether your biological conclusion survives
reasonable preprocessing and engine choices on the same compact pathogen
panel, use `demo rabies-method-sensitivity-panel`. This workflow reruns a
declared four-variant matrix over the packaged rabies nucleoprotein FASTA,
compares `auto` versus `ginsi` alignment, compares trimAl gap-threshold versus
`gappyout` trimming, and then compares rooted IQ-TREE versus rooted FastTree
topologies for each variant. The packaged workflow config now declares
`parallel_workers`, so those independent variant roots can be executed in
parallel safely without reusing the same output prefix.

```bash
bijux-phylogenetics demo rabies-method-sensitivity-panel \
  --out artifacts/rabies-method-sensitivity-panel \
  --mafft-executable mafft \
  --trimal-executable trimal \
  --iqtree-executable iqtree2 \
  --fasttree-executable FastTree \
  --json
```

The resulting bundle is designed for reviewer-facing method-sensitivity
inspection rather than only raw tree production:

- one workflow summary with stable-clade, changed-clade, rooted preprocessing, and rooted engine counts
- one variant summary with alignment length, trimmed length, selected model, support range, and rooted engine RF distance
- one parallel-execution summary plus one workflow manifest that record the worker count, execution mode, and per-variant isolated task outputs
- one raw workflow execution record at `workflow/rabies-method-sensitivity-panel.run.json` that records run status, worker count, successful variants, failed variants, and per-variant task logs before the reviewer bundle is even inspected
- one rooted preprocessing comparison table across every declared variant pair
- one stable-clade ledger and one changed-clade ledger aggregated across variants
- one method-conclusion ledger that states which claims remained stable and which remained engine-sensitive
- one report manifest that records the reviewer HTML links, checksums, and byte counts for the governed workflow ledgers
- one per-variant task log so concurrent variant runs remain inspectable without mixing batch execution output
- one per-variant evidence package with aligned FASTA, trimmed FASTA, unrooted comparison tables, clade ledgers, support-weighted conflicts, rooted trees, and rooting review

The reviewer HTML is intentionally compact. It surfaces top-level counts and
conclusion text, then links the large TSV and JSON artifacts instead of
embedding those tables directly. The JSON metrics expose the same report-size
surface through `report_linked_artifact_count`, `report_html_size_bytes`,
`report_linked_artifact_bytes`, and `report_total_output_bytes`.
If a second caller tries to reuse the same workflow output root while one run
is still active, the command now fails explicitly instead of interleaving task
logs and variant outputs.

When the goal is to review posterior or bootstrap uncertainty at larger scale,
`tree-set report` now uses the same summary-first rule. Large tables are
written into one sibling `<report>.artifacts/` directory as TSV or JSON, the
HTML report links those artifacts directly, and the CLI reports the resulting
HTML bytes, linked-artifact bytes, and total output bytes. For inputs with
`1,000+` trees the command switches to `scaled-summary` mode and replaces the
highest-cost supplemental sensitivity passes with linked note artifacts so the
report remains reviewable instead of attempting an unbounded inline expansion.

This workflow is intentionally narrower than the flagship integrated rabies
demo. It does not fit comparative models or reconstruct host and geography
states. Its job is to answer a different scientific question: whether the
rooted high-level rabies conclusion is robust to declared alignment, trimming,
and engine choices, and where unrooted internal structure remains sensitive to
those choices.

When the goal is to prove one full sequence-to-result biological workflow on a
real compact pathogen panel, use `demo rabies-cross-host-geography-panel`.
This workflow starts from raw rabies nucleoprotein sequences plus one combined
host-and-geography metadata table, reruns alignment, trimming,
maximum-likelihood inference, bootstrap support estimation, explicit outgroup
rooting, bootstrap topology review, clade extraction, host-switching review,
geographic transition review, migration-event extraction, one comparative
model, and one final HTML handoff.

```bash
bijux-phylogenetics demo rabies-cross-host-geography-panel \
  --out artifacts/rabies-cross-host-geography-panel \
  --config packages/bijux-phylogenetics/src/bijux_phylogenetics/resources/datasets/pathogens/rabies_cross_host_geography_panel/workflow-config.json \
  --mafft-executable mafft \
  --trimal-executable trimal \
  --iqtree-executable iqtree2 \
  --json
```

The packaged integrated workflow writes the same dataset/workflow split as the
other demos, but it keeps the full scientific chain together:

- one accession-provenance ledger at `dataset/source-accessions.tsv`
- FASTA validation, sequence-type detection, and alignment-quality review ledgers
- aligned FASTA and trimmed FASTA
- maximum-likelihood rooted tree with bootstrap support ledger
- bootstrap topology review with consensus, clade-frequency, instability, and distance tables
- rooted clade table annotated with host and region metadata
- explicit rooting-evidence ledger
- host-switch summary, node, branch, count, fit, unsupported-claim, and exclusion ledgers
- full biogeography package with ancestral-region tree, transition matrix, migration events, and self-contained map
- comparative trait table, comparative-ready tree, branch-adjustment ledger, and comparative report package
- conclusion-stability ledgers for key clades, support values, ancestral states, and comparative coefficients
- integrated HTML report and machine-readable manifest

Its `workflow-config.json` is also the place to govern execution budgets
honestly. That config now carries the real runtime controls for this workflow:

- `iqtree_threads` for the IQ-TREE execution lane
- `timeout_seconds` for the engine-backed sequence-to-tree steps
- `max_bootstrap_tree_count` for bootstrap-summary inputs
- `max_report_table_rows` for reviewer-facing HTML tables
- `memory_warning_threshold_bytes` for measured bootstrap-review memory

When one of those limits is crossed, the workflow now raises one structured
budget failure or records one explicit warning in `workflow-summary.tsv`
instead of silently continuing.

At the package root it also writes one reviewer-facing overview HTML page and
one package manifest. Those top-level artifacts state the biological question,
one short answer, the exact reproduction config path and checksum, and the key
workflow metrics without forcing a reviewer to open the deeper integrated
report first.

The conclusion-stability surface combines two uncertainty lanes that already
exist elsewhere in the repository:

- rooted bootstrap topology review
- governed rabies method-sensitivity variants across alignment, trimming, and engine choices

Its reviewer-facing HTML report separates `stable`, `weak`, and `unstable`
conclusions directly and the JSON metrics expose
`conclusion_stable_count`, `conclusion_weak_count`, and
`conclusion_unstable_count`.

The packaged metadata carries both raw `host_species` and `country` provenance
plus the grouped workflow traits `host_group` and `region_group`. Those grouped
traits are explicit and intentional: they keep one small real rabies panel
interpretable for end-to-end host and geography review without overstating
species-level or locality-level resolution. The comparative layer derives one
regional longitude trait from the shipped centroid table and records any
nonpositive branch-length repair needed before comparative fitting. This
integrated workflow does depend on external `mafft`, `trimal`, and `iqtree2`
executables because it reruns the real sequence-to-tree path instead of
relying on a pre-rooted tree.

When the goal is to fit a phylogenetic regression rather than only measure
signal, use `comparative pgls`. The command inspects the requested response and
predictors, fits generalized least squares on the phylogenetic covariance, and
reports coefficient estimates with standard errors, Student-t test statistics,
p-values, explicit 95% confidence intervals, and AIC.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size brain_mass_g \
  --taxon-column species \
  --json
```

When `--lambda-value estimate` is used, Bijux first chooses the covariance
strength that best fits the data and then reports coefficient uncertainty for
that fitted covariance. When an external workflow already fixed one
phylogenetic covariance, pass that exact numeric lambda instead. That keeps
coefficient and p-value review tied to one covariance assumption rather than
silently mixing model-selection differences with coefficient-level inference.

When the goal is to inspect the covariance assumptions themselves before
trusting a comparative fit, use `comparative covariance-audit`. This workflow
does not fit coefficients first and explain problems later. It audits the tree,
trait table, taxon overlap, branch lengths, and candidate covariance matrices
up front for PGLS, Brownian trait models, or OU trait models.

```bash
bijux-phylogenetics comparative covariance-audit \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --analysis pgls \
  --formula "longevity ~ social_group_size + mating_system" \
  --summary-out artifacts/primates.covariance-audit-summary.tsv \
  --candidates-out artifacts/primates.covariance-audit-candidates.tsv \
  --excluded-taxa-out artifacts/primates.covariance-audit-excluded.tsv \
  --json
```

The JSON result reports:

- matrix dimension and rank
- matched taxa and analysis-taxon counts
- taxa missing from the trait table or absent from the tree
- duplicate tree or trait taxa
- zero-length and negative branch counts
- covariance singularity and near-singularity
- raw condition number
- whether fitting would proceed by exact solve, regularization, pseudoinverse, or failure

For `--analysis pgls`, use `--formula` or `--response` plus `--predictors`.
For `--analysis brownian-trait` or `--analysis ou-trait`, use `--trait`.
Estimated Pagel's lambda and estimated OU alpha produce one candidate ledger
row per audited covariance candidate so reviewers can see whether only one part
of the profile is unstable.

The ledger outputs keep the audit reviewable:

- `--summary-out` writes one overall audit row
- `--candidates-out` writes one row per lambda or alpha candidate, or one row for fixed Brownian covariance
- `--excluded-taxa-out` writes the explicit taxon-pruning reasons used before covariance construction

This workflow is intentionally honest about stabilization. The current
governed comparative fitting paths may regularize covariance inversion with a
small diagonal epsilon, and the audit reports that as `regularization` instead
of hiding it behind a successful coefficient table.

When the response is binary rather than continuous, use `comparative logistic`.
This workflow keeps the comparative formula surface, but fits a binary
working-correlation approximation instead of reusing continuous-trait PGLS
output. It requires the response to be encoded explicitly as `0` and `1`,
supports the same predictor encoding rules as `comparative pgls`, and reports
Wald-normal coefficient uncertainty plus explicit separation-risk warnings.
It does not currently claim `ape::compar.gee` parity, so it should be used as
an exploratory approximate surface rather than a drop-in `ape::compar.gee`
replacement.

```bash
bijux-phylogenetics comparative logistic \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response sociality_present \
  --predictors brain_mass_g habitat \
  --taxon-column species \
  --lambda-value 1.0 \
  --coefficients-out artifacts/primates.logistic-coefficients.tsv \
  --fitted-out artifacts/primates.logistic-fitted.tsv \
  --excluded-taxa-out artifacts/primates.logistic-excluded.tsv \
  --json
```

The coefficient ledger keeps one row per fitted term with the estimate,
standard error, Wald test statistic, p-value, and 95% confidence interval. The
fitted ledger keeps one row per analyzed taxon with the observed binary
response, fitted probability, linear predictor, and raw residual. The
excluded-taxa ledger keeps the same explicit pruning contract as the other
comparative workflows. If the fitted probabilities collapse toward `0` or `1`,
or the working information matrix becomes singular enough to require
stabilization, the JSON result marks that as separation risk instead of
pretending the approximation is as stable as an ordinary continuous-trait fit.

When the goal is to compare competing comparative hypotheses rather than fit
just one formula, use `comparative model-selection`. This workflow keeps one
shared complete-case taxon set across every candidate formula, auto-detects
whether the shared response should use continuous-trait PGLS or the binary
working-correlation logistic surface, and then ranks the candidate formulas by
information criteria on that one fixed analysis set.

```bash
bijux-phylogenetics comparative model-selection \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "sociality_present ~ brain_mass_g" \
  --formula "sociality_present ~ habitat" \
  --formula "sociality_present ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value 1.0 \
  --ranking-out artifacts/primates.model-ranking.tsv \
  --pairwise-out artifacts/primates.model-pairwise.tsv \
  --excluded-taxa-out artifacts/primates.model-excluded.tsv \
  --json
```

The ranking ledger keeps one row per candidate with log-likelihood, AIC, AICc,
BIC, delta values, Akaike weight, selected-model status, and the encoded model
columns actually used in the fit. The pairwise ledger makes the comparison
contract explicit by marking whether each candidate pair is identical, nested,
or non-nested, and by preserving a likelihood-ratio statistic only where a
nested comparison is real. The excluded-taxa ledger records the shared
complete-case rule directly so reviewers can see which taxa were dropped before
any candidate formula was ranked.

When the goal is to ask whether two traits changed in coupled ways across one
tree, use `comparative correlated-traits`. This workflow is pairwise
trait-evolution review rather than another regression wrapper. For two numeric
traits it evaluates correlated versus independent Brownian contrast covariance.
For two binary traits it evaluates independent versus joint-state discrete
evolution on the shared analyzed taxon set.

```bash
bijux-phylogenetics comparative correlated-traits \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --left-trait body_mass_g \
  --right-trait home_range_km2 \
  --taxon-column species \
  --summary-out artifacts/primates.correlated-traits-summary.tsv \
  --comparison-out artifacts/primates.correlated-traits-comparison.tsv \
  --observations-out artifacts/primates.correlated-traits-observations.tsv \
  --excluded-taxa-out artifacts/primates.correlated-traits-excluded.tsv \
  --json
```

The summary ledger keeps the selected analysis kind, the association measure,
the fitted covariance or correlation surface, and the independent-versus-
correlated likelihood comparison. The comparison ledger keeps the two model
rows directly so reviewers can inspect parameter counts, log-likelihoods, AIC,
delta AIC, and selection outcome without recomputing them from JSON. The
observation ledger keeps the evidence rows that drove the fit: one
contrast-per-node ledger for continuous traits or one tip-state ledger for
binary traits. The excluded-taxa ledger preserves the shared filtering reasons
explicitly, including missing trait values and taxa present in only one input.

For binary traits, the surface is intentionally explicit about scope. Bijux
fits a joint-state discrete pseudo-likelihood review surface rather than
claiming a full native reimplementation of an external Pagel binary-correlation
fit. That approximation boundary is preserved in the warnings and JSON output
instead of being hidden behind the same wording as the continuous-trait path.

When the goal is to reconstruct one continuous ancestral trait directly, use
`ancestral continuous`. This workflow estimates internal-node values for one
numeric trait on the analyzed rooted tree, reports 95% uncertainty intervals
for internal nodes, and keeps explicit ledgers for dropped missing or
non-numeric tips instead of silently excluding them.

```bash
bijux-phylogenetics ancestral continuous \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --model brownian \
  --estimator fast-anc \
  --table-out artifacts/primates.ancestral-nodes.tsv \
  --summary-out artifacts/primates.ancestral-summary.tsv \
  --uncertainty-out artifacts/primates.ancestral-uncertainty.tsv \
  --exclusions-out artifacts/primates.ancestral-excluded.tsv \
  --json
```

The node ledger keeps one row per analyzed tip or internal node so reviewers
can inspect the reconstructed value surface directly. The summary ledger keeps
one row with the analyzed taxon count, excluded taxon count, internal node
count, unstable node count, the Brownian ultrametric verdict, covariance
conditioning fields, one GLS likelihood summary, and the root estimate with
its 95% interval. The uncertainty ledger keeps one row per internal node with
the standard error, interval width, confidence score, and interpretation label
so broad intervals cannot be hidden behind only point estimates. The
excluded-taxa ledger keeps one row per dropped tip with an explicit reason
such as `missing_trait_value` or `non_numeric_trait_value`.

The governed live `ape` lane for this surface is `ape::ace(type='continuous',
method='pic', CI=TRUE)` over balanced, pectinate, six-taxon, and pruned
missing-value fixtures. Bijux therefore uses one explicit closed-form
Brownian parity target while still surfacing extra owned diagnostics that live
`ape::ace(..., method='pic')` does not emit. When `--estimator fast-anc` is
selected instead, the Brownian runtime switches to the governed live
`phytools::fastAnc` lane and keeps stable node-signature rows plus standard
errors for parity review on the shared comparative fixture corpus. When
`--estimator anc-ml` is selected instead, the Brownian runtime switches to
the governed live `phytools::anc.ML` lane and adds fitted sigma-squared,
Brownian log-likelihood, explicit optimizer diagnostics, stable
node-signature rows, standard errors, and 95% intervals under the same
missing-value pruning policy.

For Python-native work, the same continuous ancestral runtime is also
available directly through
`reconstruct_continuous_ancestral_states_from_dataset(...)` once one
`AncestralContinuousDataset` has already been loaded.

When the goal is to reconstruct one categorical ancestral trait directly, use
`ancestral discrete`. This workflow supports a fast Fitch path for parsimony
review and supports likelihood reconstructions under `equal-rates`,
`symmetric`, or `all-rates-different` Mk models when reviewers need marginal
node probabilities instead of only one state label. When the baseline model is
Fitch, the workflow also reports the minimum parsimony change count, the number
of ambiguous internal nodes, and an optional direct node-by-node comparison
against one likelihood model. For likelihood models, the same surface also
reports the fitted log-likelihood, AIC, one directed transition-rate ledger,
and an owned root-prior policy through `--root-prior-mode equal|empirical|fixed`
plus `--fixed-root-state`.

```bash
bijux-phylogenetics ancestral discrete \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --taxon-column species \
  --model equal-rates \
  --root-prior-mode empirical \
  --table-out artifacts/primates.ancestral-discrete-nodes.tsv \
  --summary-out artifacts/primates.ancestral-discrete-summary.tsv \
  --probabilities-out artifacts/primates.ancestral-discrete-probabilities.tsv \
  --transitions-out artifacts/primates.ancestral-discrete-transitions.tsv \
  --exclusions-out artifacts/primates.ancestral-discrete-excluded.tsv \
  --json
```

The node ledger keeps one row per analyzed tip or internal node with the
reported state label and state set. The summary ledger keeps one row with the
analyzed taxon count, excluded taxon count, internal node count, ambiguous
internal node count, unstable node count, observed state count, sparse state
count, the root state with its confidence, the root-prior mode, and the fitted
likelihood summary. Under Fitch, the same summary ledger also preserves the
minimum parsimony change count and the number of parsimonious root states so
the fast path still exposes its evidence burden explicitly. The probabilities
ledger keeps one row per internal node with the full marginal probability
vector so reviewers can distinguish a narrow state assignment from a weakly
resolved one. For likelihood models, the transition ledger keeps one directed
row per fitted state-to-state rate so reviewers can inspect the ER, SYM, or
ARD fit directly instead of inferring it from node calls alone. When
`--comparison-out` is supplied, the command also writes one direct node-wise
comparison ledger between the baseline model and the requested comparison
model. The excluded-taxa ledger keeps one row per dropped tip with an explicit
reason such as `missing_discrete_trait_state`. The same summary ledger now also
states whether the run is comparable to live `phytools::rerootingMethod`: ER
and SYM with the equal root prior are governed rerooting-parity surfaces,
while ARD, Fitch, ordered-state, empirical-root-prior, and fixed-root-prior
runs are flagged explicitly as non-comparable.

For Python-native work, the same discrete ancestral runtime is also available
directly through `reconstruct_discrete_ancestral_states_from_dataset(...)`
once one `AncestralDiscreteDataset` has already been loaded.

When the goal is to verify that the owned discrete ancestral likelihood
surface still matches governed external references and known policy examples,
use `ancestral discrete-reference`. This workflow reruns the checked-in
`ape::ace` ER, SYM, and ARD reference cases and then reruns the owned
root-prior, ambiguity, ordered-state, and irreversible-loss review surfaces
on packaged known examples.

```bash
bijux-phylogenetics ancestral discrete-reference \
  --json
```

The JSON report keeps the trust contract compact and explicit:

- one governed case count across the whole suite
- one external-case count for the `ape::ace` parity subset
- one all-passed flag for the full discrete ancestral reference lane

The governed live `ape::ace` discrete lane is intentionally explicit rather
than broad. ER, SYM, and ARD parity now cover governed balanced, pectinate,
six-taxon, and pruned missing-value cases, and the owned fit review warns when
multi-parameter likelihood surfaces hit optimizer bounds. Root-prior controls
are still validated as owned Bijux policy because `ape::ace` does not expose
the same runtime root-prior interface.

When the goal is to decide whether a discrete trait should be reconstructed as
ordered rather than unordered, use `ancestral ordered-discrete`. This workflow
fits the same discrete likelihood model twice on one tree: once with the
supplied ordered state vocabulary and once with the unrestricted unordered
baseline. It then writes the fit comparison, node-wise ancestral differences,
and directed transition restrictions so reviewers can see exactly what the
ordered assumption changed.

```bash
bijux-phylogenetics ancestral ordered-discrete \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait habitat_stage \
  --taxon-column species \
  --model equal-rates \
  --ordered-states absent,partial,complete \
  --summary-out artifacts/primates.ordered-discrete-summary.tsv \
  --fits-out artifacts/primates.ordered-discrete-fits.tsv \
  --nodes-out artifacts/primates.ordered-discrete-nodes.tsv \
  --transitions-out artifacts/primates.ordered-discrete-transitions.tsv \
  --json
```

The summary ledger keeps one row with the ordered and unordered likelihoods,
AIC values, delta AIC, preferred ordering, differing node count, ambiguity
change count, and restricted transition count. The fit ledger keeps one row
for the ordered fit and one for the unordered baseline, including the ordered
state vocabulary, parameter count, root state, and root confidence. The node
ledger keeps one row per internal node with the ordered state, unordered
state, confidence delta, and ambiguity change flag. The transition ledger
keeps one row per directed state change so reviewers can see which transitions
the ordered model forbids entirely and which adjacent transitions remain
allowed.

When the goal is to model irreversible gains or losses, use
`ancestral irreversible-discrete`. This workflow fits one constrained discrete
likelihood model under an explicit directed transition graph and compares it to
the unconstrained baseline of the same model family. It is the owned review
surface for one-way losses, forbidden reversals, and similar trait histories
where some transitions should be impossible.

```bash
bijux-phylogenetics ancestral irreversible-discrete \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait trait_state \
  --taxon-column species \
  --model all-rates-different \
  --allowed-transitions present->absent \
  --summary-out artifacts/primates.irreversible-summary.tsv \
  --fits-out artifacts/primates.irreversible-fits.tsv \
  --nodes-out artifacts/primates.irreversible-nodes.tsv \
  --transitions-out artifacts/primates.irreversible-transitions.tsv \
  --json
```

The summary ledger keeps one row with the constrained and unconstrained
likelihoods, likelihood difference, AIC difference, preferred constraint,
differing node count, ambiguity change count, and forbidden transition count.
The fit ledger keeps one row for the constrained fit and one for the
unconstrained baseline, including root state and root confidence. The node
ledger keeps one row per internal node with the constrained state,
unconstrained state, confidence delta, and ambiguity change flag. The
transition ledger keeps one row per directed state change so reviewers can see
which transitions the constraint forbids and what rate the unconstrained model
would otherwise assign to them.

When the goal is to reconstruct ancestral geographic regions rather than only a
generic discrete trait, use `biogeography model`. This workflow accepts one
taxon-region table, fits one ER, SYM, or ARD geographic transition model on
the rooted tree, and writes reviewer-facing ledgers for internal-node region
probabilities, pairwise transition rates, branchwise geographic events, and
explicit excluded taxa.

```bash
bijux-phylogenetics biogeography model \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait region \
  --taxon-column species \
  --model ard \
  --summary-out artifacts/primates.biogeography-summary.tsv \
  --nodes-out artifacts/primates.biogeography-nodes.tsv \
  --rates-out artifacts/primates.biogeography-rates.tsv \
  --events-out artifacts/primates.biogeography-events.tsv \
  --exclusions-out artifacts/primates.biogeography-excluded.tsv \
  --json
```

The summary ledger keeps one row with the ER/SYM/ARD choice, analyzed taxon
count, observed region count, internal-node count, changed-branch count,
strongly supported transition count, root region, and root-region support. The
node ledger keeps one row per internal node with the most likely ancestral
region and the full region-probability vector. The rate ledger keeps one row
per directed source-target region pair with the fitted rate and uncertainty
interval. The event ledger keeps one row per branch with the inferred source
and target regions plus branch-level support. The excluded-taxa ledger keeps
raw region rows that were dropped because the taxon was absent from the tree or
the region coding was unusable.

The current biogeography surface is explicit about method scope. It reuses the
owned deterministic geographic-state engine already in the repository. The
one-tree `biogeography model` workflow should therefore be interpreted as
model-conditioned geographic state evidence rather than as direct proof of
dispersal timing or mechanism.

When the goal is to forbid impossible region changes and compare that
constraint against the unconstrained fit, use `biogeography constrained`. This
workflow accepts one rooted tree, one taxon-region table, and one explicit
region adjacency matrix. It reuses the owned likelihood discrete ancestral
runtime to fit one constrained and one unconstrained geographic model on the
same analyzed taxa, then writes reviewer-facing ledgers for fit comparison,
pairwise transition-rate comparison, forbidden unconstrained branch claims, and
excluded taxa.

```bash
bijux-phylogenetics biogeography constrained \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  artifacts/region-adjacency.tsv \
  --trait region \
  --taxon-column species \
  --model ard \
  --summary-out artifacts/primates.biogeography-constrained-summary.tsv \
  --fits-out artifacts/primates.biogeography-constrained-fits.tsv \
  --transitions-out artifacts/primates.biogeography-constrained-transitions.tsv \
  --unsupported-out artifacts/primates.biogeography-constrained-unsupported.tsv \
  --exclusions-out artifacts/primates.biogeography-constrained-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, observed region
count, allowed and forbidden transition counts, constrained and unconstrained
likelihoods, AIC difference, preferred constraint regime, unsupported
transition-claim count, and warning count. The fit ledger keeps one row per
constraint regime with root region and fit evidence. The transition ledger
keeps one row per directed region pair with the allowed flag plus constrained
and unconstrained rates. The unsupported-claims ledger keeps one row per
unconstrained branchwise transition claim that violates the supplied adjacency
matrix, including the constrained replacement claim on the same branch. The
excluded-taxa ledger keeps dropped or unsupported raw region rows.

The constrained geography surface is explicit about scope. It is an
adjacency-constrained likelihood review over the owned discrete ancestral
runtime, not a richer time-stratified dispersal-extinction-cladogenesis model.
Its value is to show where unconstrained geographic transition claims are not
compatible with the supplied region-connectivity contract and whether the
constrained fit remains competitive.

When the goal is to list inferred dispersal or migration events directly, use
`biogeography events`. This workflow accepts one rooted tree and one
taxon-region table, or the same inputs plus `--tree-set` when the tree input is
a posterior or bootstrap tree set. It reuses the owned ancestral geographic
reconstruction, keeps only changed source-target branches as event rows, and
writes reviewer-facing ledgers for branch identity, root-depth interval,
midpoint depth estimate, support, and excluded taxa. In tree-set mode it also
writes retained-tree ledgers and comparable event summaries across trees.

```bash
bijux-phylogenetics biogeography events \
  artifacts/primates.posterior.nwk \
  artifacts/primates.csv \
  --trait region \
  --taxon-column species \
  --model ard \
  --tree-set \
  --burnin-fraction 0.25 \
  --summary-out artifacts/primates.biogeography-events-summary.tsv \
  --trees-out artifacts/primates.biogeography-events-trees.tsv \
  --events-out artifacts/primates.biogeography-events.tsv \
  --event-summaries-out artifacts/primates.biogeography-event-summaries.tsv \
  --exclusions-out artifacts/primates.biogeography-events-excluded.tsv \
  --json
```

On one tree, the summary ledger keeps one row with the analyzed taxon count,
tree depth, event count, strongly supported event count, mean event support,
and earliest and latest event midpoint depths. The event ledger keeps one row
per inferred geographic movement event with the source region, target region,
branch identity, descendant taxa, branch length, parent and child depths,
midpoint depth, support, and confidence class. The excluded-taxa ledger keeps
dropped raw region rows.

In tree-set mode, the retained-tree ledger keeps one row per retained tree with
rooted and unrooted topology identifiers plus tree-level event counts. The
event ledger keeps one row per inferred event per retained tree. The
event-summary ledger groups comparable source-target events by branch identity
across retained trees and reports presence fractions, strongly supported tree
fractions, mean support, empirical midpoint-depth bounds, and stability class.

The migration-event surface is explicit about what its time column means. It
does not claim an exact stochastic event time on a branch. Instead, it reports
the branch depth interval and a deterministic midpoint-depth estimate so
reviewers can place likely movement events on the tree without overstating
temporal precision.

When the goal is to deliver one complete geographic evolution review package
rather than separate region, event, and map commands, use `biogeography
report`. This workflow accepts one rooted tree, one taxon-region table, and
one explicit centroid table. It reruns the owned ancestral-region model, the
owned migration-event extraction surface, the owned ancestral-region tree
renderer, and the owned region-map renderer into one governed package.

```bash
bijux-phylogenetics biogeography report \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait region \
  artifacts/region-centroids.tsv \
  --taxon-column species \
  --model ard \
  --out-dir artifacts/primates-biogeography-report \
  --json
```

The package writes:
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

That bundle makes the core evidence surfaces explicit in one place: observed
region counts at the tips, ancestral-region probabilities at internal nodes,
the fitted directed transition matrix, branchwise migration-event rows, the
ancestral-region tree figure, and the self-contained region-transition map.
The HTML report embeds the tree figure directly and links the same package to
the map artifact, so reviewers can move from rates and counts to geographic
placement without reassembling outputs by hand.

When the goal is to compare inferred geographic transition pressure across
explicit historical intervals, use `biogeography time-stratified`. This
workflow accepts one rooted tree with positive branch lengths, one taxon-region
table, and one or more explicit root-depth bins in `LABEL:START:END` form. It
reuses the governed ancestral-region reconstruction, allocates branch exposure
and inferred branch changes across the requested intervals, and writes
reviewer-facing ledgers for interval-specific transition matrices, branch-bin
overlaps, and excluded taxa.

```bash
bijux-phylogenetics biogeography time-stratified \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait region \
  --taxon-column species \
  --model ard \
  --time-bin deep-history:0.0:12.0 \
  --time-bin crown-radiation:12.0:24.0 \
  --summary-out artifacts/primates.biogeography-time-summary.tsv \
  --matrix-out artifacts/primates.biogeography-time-matrix.tsv \
  --branches-out artifacts/primates.biogeography-time-branches.tsv \
  --exclusions-out artifacts/primates.biogeography-time-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, tree depth,
time-bin count, changed-branch count, total allocated transition weight, and
warning count. The matrix ledger keeps one row per directed source-target
region pair inside each requested interval, including source exposure length,
allocated transition weight, interval-specific rate, and the corresponding
global one-tree rate for comparison. The branch ledger keeps one row per
branch-bin overlap so reviewers can see exactly which reconstructed changes and
which source-region exposure lengths were assigned to each interval. The
excluded-taxa ledger keeps the same dropped region rows as the one-tree model
workflow.

The time-stratified surface is also explicit about method scope. It is a
deterministic interval-allocation review surface layered on the owned
geographic-state reconstruction, not a full time-inhomogeneous stochastic
biogeographic process fit. When requested intervals do not cover the full tree
depth, the workflow reports that boundary as a warning instead of silently
pretending those uncovered branch segments were modeled.

When the goal is to review one dated biogeographic tree directly instead of
declaring explicit historical bins by hand, use `biogeography chronology`.
This workflow accepts one rooted ultrametric time tree and one taxon-region
table, verifies that the tree is time-scaled, extracts node ages, maps
inferred geographic transitions to automatic equal-width age bins, and reports
which bins are weakly supported. This keeps the dated-tree chronology surface
separate from `biogeography time-stratified`, which remains the explicit
interval-definition workflow.

```bash
bijux-phylogenetics biogeography chronology \
  artifacts/primates.time-tree.nwk \
  artifacts/primates.csv \
  --trait region \
  --taxon-column species \
  --model ard \
  --time-bin-count 4 \
  --summary-out artifacts/primates.biogeography-chronology-summary.tsv \
  --nodes-out artifacts/primates.biogeography-node-ages.tsv \
  --events-out artifacts/primates.biogeography-dated-events.tsv \
  --bins-out artifacts/primates.biogeography-time-bins.tsv \
  --exclusions-out artifacts/primates.biogeography-chronology-excluded.tsv \
  --json
```

The summary ledger keeps one row with the model, analyzed taxon count,
time-scaled-tree audit result, root age, node-age row count, dated event
count, time-bin count, empty-bin count, high-uncertainty-bin count, and
warning count. The node ledger keeps one row per reviewed node with branch
depth, age before present, most likely region, confidence, and root flag. The
event ledger keeps one row per inferred geographic transition with parent and
child depths, parent and child ages before present, midpoint event age,
automatic time-bin label, and support class.

The time-bin ledger is the chronology review surface. Each row preserves the
bin bounds, event totals, strongly supported versus low-support event counts,
mean support, support uncertainty, earliest and latest dated event in the bin,
dominant transition, transition diversity, and one reviewer-facing uncertainty
class. This is intentionally a dated-tree review over one owned reconstruction,
not a full time-varying stochastic biogeographic process fit. Equal-width age
bins are explicit reviewer-facing chronology bins rather than hidden model
parameters.

When the goal is to test whether ancestral geographic conclusions are being
driven by uneven region sampling instead of the phylogenetic signal itself, use
`biogeography sampling-bias`. This workflow accepts one rooted tree, one
taxon-region table, and an optional region-weight table. It reports raw region
sample counts, applies either automatic inverse-frequency weights or explicit
user weights, reruns the owned deterministic geographic reconstruction under
that weighting scheme, and compares the weighted and unweighted conclusions
directly.

```bash
bijux-phylogenetics biogeography sampling-bias \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait region \
  --taxon-column species \
  --model ard \
  --weights artifacts/region-weights.tsv \
  --summary-out artifacts/primates.biogeography-sampling-summary.tsv \
  --regions-out artifacts/primates.biogeography-region-counts.tsv \
  --nodes-out artifacts/primates.biogeography-weighted-nodes.tsv \
  --transitions-out artifacts/primates.biogeography-weighted-transitions.tsv \
  --exclusions-out artifacts/primates.biogeography-sampling-excluded.tsv \
  --json
```

The summary ledger keeps one row with the weighting mode, observed-region
count, dominant-region fraction before and after weighting, whether the root
region changes, how many internal nodes change, how many branch transitions
change, and the warning burden. The region-count ledger is the main sampling
audit surface. Each row preserves the raw sample count, raw sample fraction,
applied weight, weighted sample count, weighted sample fraction, and whether
the region dominates before or after weighting.

The node ledger compares weighted and unweighted ancestral region probabilities
for each internal node. The transition ledger compares weighted and unweighted
branchwise transitions, including whether the correction changes the inferred
source or target region on that branch. This surface is intentionally explicit
about scope: it is a reviewer-facing weighted reanalysis over the owned
deterministic geographic reconstruction, not a full hierarchical sampling model
or a hidden claim that the weighting scheme is uniquely correct.

When the goal is to reconstruct host association evolution and count inferred
host shifts on one pathogen or parasite tree, use `host-association switches`.
This workflow accepts one rooted tree, one host metadata table, and
optionally one explicit host-transition constraint ledger. It reconstructs
internal host states, classifies branchwise host shifts as certain or
uncertain, aggregates directed host-switch counts, and compares constrained
versus unconstrained host-transition fits when a constraint ledger is
supplied.

```bash
bijux-phylogenetics host-association switches \
  artifacts/pathogens.nwk \
  artifacts/pathogens-hosts.tsv \
  --trait host \
  --taxon-column taxon \
  --model er \
  --constraints artifacts/host-transition-constraints.tsv \
  --summary-out artifacts/pathogens.host-switch-summary.tsv \
  --nodes-out artifacts/pathogens.host-switch-nodes.tsv \
  --branches-out artifacts/pathogens.host-switch-branches.tsv \
  --counts-out artifacts/pathogens.host-switch-counts.tsv \
  --fits-out artifacts/pathogens.host-switch-fits.tsv \
  --unsupported-out artifacts/pathogens.host-switch-unsupported.tsv \
  --exclusions-out artifacts/pathogens.host-switch-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, observed host
count, internal-node ambiguity burden, certain and uncertain switch counts,
allowed and forbidden transition counts, constrained and unconstrained fit
evidence, preferred constraint regime, unsupported unconstrained claim count,
and root-host support. The node ledger keeps one row per internal node with
the most likely host and the full host-probability vector. The branch ledger
keeps one row per analyzed branch with the parent and child host sets,
certainty class, and allowed-transition flag. The count ledger aggregates
directed host-switch counts across branches. When `--constraints` is supplied,
the fit ledger keeps one row for the unconstrained fit and one for the
constrained fit, while the unsupported-claims ledger keeps one row for each
unconstrained branchwise host-switch claim that violates the supplied
transition policy. The excluded-taxa ledger keeps dropped raw host rows.

The host-association surface is explicit about scope. It is a one-tree
host-state evolution review over the owned discrete ancestral runtime, not a
full cophylogenetic reconciliation or a transmission-history inference model.
Its value is to make inferred host-switch evidence, uncertainty, and explicit
transition constraints reviewable on the analyzed tree.

When the goal is to reconstruct ecological niche evolution and identify which
clades carry concentrated niche-shift burden, use
`ecological-niche transitions`. This workflow accepts one rooted tree and one
ecological-state table, fits one ER, SYM, or ARD discrete transition model,
reconstructs ancestral niche states, counts branchwise niche changes, and
ranks internal clades by repeated or concentrated transition burden.

```bash
bijux-phylogenetics ecological-niche transitions \
  artifacts/lineages.nwk \
  artifacts/lineages-niche.tsv \
  --trait niche \
  --taxon-column taxon \
  --model er \
  --summary-out artifacts/lineages.niche-summary.tsv \
  --nodes-out artifacts/lineages.niche-nodes.tsv \
  --rates-out artifacts/lineages.niche-rates.tsv \
  --branches-out artifacts/lineages.niche-branches.tsv \
  --counts-out artifacts/lineages.niche-counts.tsv \
  --clades-out artifacts/lineages.niche-clades.tsv \
  --exclusions-out artifacts/lineages.niche-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, observed niche
count, internal-node ambiguity burden, fitted log-likelihood and AIC,
transition-rate row count, certain and uncertain niche-transition counts, and
the number of internal clades carrying repeated shifts. The node ledger keeps
one row per internal node with the most likely niche and full niche
probability vector. The rate ledger keeps one row per directed niche pair with
the fitted transition rate. The branch ledger keeps one row per analyzed
branch with the inferred parent and child niches, certainty class, and support
score. The count ledger aggregates directed niche-transition counts across
branches. The clade ledger ranks internal clades by concentrated shift burden,
including repeated-shift flags and dominant transition labels. The
excluded-taxa ledger keeps dropped raw niche rows.

The ecological-niche surface is explicit about scope. It is a one-tree
discrete niche-evolution review over the owned ancestral runtime, not a full
macroecological process model. Its value is to make niche-transition evidence,
uncertainty, and clade-localized shift concentration directly reviewable on
the analyzed tree.

When the goal is to reconstruct continuous spatial movement directly from
latitude and longitude traits on one tree, use `phylogeography coordinates`.
This workflow accepts one rooted tree plus one taxon-keyed coordinate table,
reconstructs ancestral coordinates under Brownian or OU continuous evolution,
reviews branchwise great-circle displacement, flags jump outliers, and renders
one coordinate-space movement visualization as SVG or HTML.

```bash
bijux-phylogenetics phylogeography coordinates \
  artifacts/lineages.nwk \
  artifacts/lineages-coordinates.tsv \
  --latitude-column latitude \
  --longitude-column longitude \
  --taxon-column taxon \
  --model brownian \
  --summary-out artifacts/lineages.phylogeography-summary.tsv \
  --estimates-out artifacts/lineages.phylogeography-estimates.tsv \
  --branches-out artifacts/lineages.phylogeography-branches.tsv \
  --outliers-out artifacts/lineages.phylogeography-outliers.tsv \
  --exclusions-out artifacts/lineages.phylogeography-excluded.tsv \
  --visualization-out artifacts/lineages.phylogeography-movement.svg \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, internal-node
count, weak-support-node count, outlier and impossible jump counts, maximum
branch displacement, and root-coordinate uncertainty. The estimate ledger
keeps one row per tip or internal node with latitude, longitude, paired
uncertainty, and radial coordinate uncertainty in kilometers. The branch
ledger keeps one row per branch with parent and child coordinates,
great-circle displacement, branchwise displacement rate, and explicit jump
flags. The outlier ledger keeps only flagged branches with the robust distance
and rate thresholds used for review. The excluded-taxa ledger keeps dropped
rows for missing, non-numeric, out-of-range, or out-of-tree coordinates. The
visualization output is intentionally a coordinate-space movement review, not a
projected base map.

The phylogeography surface is explicit about scope. It is a one-tree
continuous-coordinate review over the owned continuous ancestral runtime, not
a map product or a direct historical route inference model. Its value is to
make ancestral coordinates, uncertainty, extreme branch displacements, and a
reviewable movement trace visible before the later map-specific workflow.

When the goal is to render a real HTML map from geographic reconstruction
outputs, use `phylogeography coordinates-map` for continuous latitude and
longitude reconstructions or `phylogeography regions-map` for discrete region
transitions with an explicit region-centroid table. These workflows keep the
owned reconstruction surfaces intact, then project tips, ancestral states, and
movement or transition lines onto one fixed world latitude/longitude extent.

```bash
bijux-phylogenetics phylogeography coordinates-map \
  artifacts/lineages.nwk \
  artifacts/lineages-coordinates.tsv \
  --latitude-column latitude \
  --longitude-column longitude \
  --minimum-midpoint-depth 2.0 \
  --summary-out artifacts/lineages.phylogeography-map-summary.tsv \
  --markers-out artifacts/lineages.phylogeography-map-markers.tsv \
  --lines-out artifacts/lineages.phylogeography-map-lines.tsv \
  --exclusions-out artifacts/lineages.phylogeography-map-exclusions.tsv \
  --html-out artifacts/lineages.phylogeography-map.html \
  --json
```

```bash
bijux-phylogenetics phylogeography regions-map \
  artifacts/lineages.nwk \
  artifacts/lineages-regions.tsv \
  --trait region \
  --centroids artifacts/region-centroids.tsv \
  --model ard \
  --maximum-midpoint-depth 0.15 \
  --summary-out artifacts/lineages.region-map-summary.tsv \
  --markers-out artifacts/lineages.region-map-markers.tsv \
  --lines-out artifacts/lineages.region-map-lines.tsv \
  --exclusions-out artifacts/lineages.region-map-exclusions.tsv \
  --html-out artifacts/lineages.region-map.html \
  --json
```

The summary ledger keeps one row with the map mode, model, marker counts,
line counts, visible line counts after optional midpoint-depth filtering, tree
depth, and warning count. The marker ledger keeps one row per tip or ancestral
location with latitude, longitude, confidence, and active-line context. The
line ledger keeps one row per movement or transition path with source and
target coordinates, midpoint depth, support, path distance, and one explicit
`visible` flag. The exclusion ledger keeps dropped taxa, missing centroids, or
other records that could not be placed on the rendered map. The HTML output is
self-contained and intentionally uses one fixed world extent instead of a
network-dependent tile service.

When branch lengths represent dated depth, the midpoint-depth filters make one
time-restricted visible line layer without refitting the reconstruction. This
is explicit filtering over already reconstructed branches or events, not the
full dated-tree interval analysis that belongs to the later time-aware
biogeography workflow.

When the goal is to rank which ancestral nodes or comparable clades remain
weakly resolved, use `ancestral confidence`. This workflow reads the owned
ancestral reconstruction surfaces and turns their uncertainty into one ranked
confidence ledger instead of forcing reviewers to merge multiple ancestral
tables manually. On one tree, the continuous path ranks internal nodes by
relative interval width and confidence score, while the discrete path ranks
internal nodes by maximum marginal probability, probability margin, and node
entropy. Across a tree set, the same surface ranks comparable clades by
topology presence, within-tree uncertainty, and cross-tree state or value
dispersion.

```bash
bijux-phylogenetics ancestral confidence \
  artifacts/primates.posterior.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --kind discrete \
  --model equal-rates \
  --taxon-column species \
  --tree-set \
  --burnin-fraction 0.25 \
  --summary-out artifacts/primates.ancestral-confidence-summary.tsv \
  --confidence-out artifacts/primates.ancestral-confidence.tsv \
  --json
```

The summary ledger keeps one row with the source kind, reconstruction kind,
analyzed taxon count, retained tree count when relevant, the ranked-row count,
the low-confidence count, the unstable count, the high-entropy count, and the
top uncertain node or clade. The confidence ledger is the main review surface.
On one continuous tree it keeps one ranked internal-node row with the
reconstructed estimate, 95% interval, relative uncertainty, uncertainty score,
and confidence class. On one discrete tree it keeps one ranked internal-node
row with the most likely state, full marginal probability vector, maximum
posterior probability, runner-up probability, entropy, normalized entropy,
uncertainty score, and confidence class.

On a tree set, the confidence ledger shifts from internal nodes to comparable
clades. Continuous rows preserve clade presence, mean within-tree confidence,
empirical interval width across retained trees, unstable-tree fraction, and
final uncertainty rank. Discrete rows preserve clade presence, dominant-state
fraction, state-distribution entropy, ambiguous-tree fraction,
unstable-tree fraction, and final uncertainty rank. This keeps topology
uncertainty and within-tree state uncertainty visible in the same governed
surface.

When the goal is to decide whether one discrete ancestral conclusion depends
on the root assumption rather than the likelihood model itself, use
`ancestral root-sensitivity`. This workflow reruns the owned discrete
likelihood reconstruction path under an equal root prior, an empirical root
prior derived from the observed tip-state counts, and an optional fixed-root
scenario. It then compares every internal node directly across those root
assumptions so reviewers can see whether the state changes, whether only the
support changes, and which nodes remain stable.

```bash
bijux-phylogenetics ancestral root-sensitivity \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --taxon-column species \
  --model equal-rates \
  --fixed-root-state island \
  --summary-out artifacts/primates.root-sensitivity-summary.tsv \
  --assumptions-out artifacts/primates.root-assumptions.tsv \
  --nodes-out artifacts/primates.root-sensitive-nodes.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, assumption
count, compared node count, state-changed node count, support-changed node
count, the top sensitive node, and warning count. The assumption ledger keeps
one row per root assumption with the explicit prior distribution, the inferred
root state, its confidence, its entropy, and the unstable or weak-support node
counts under that assumption. The node ledger is the main review surface: it
keeps one row per internal node with the per-assumption states, per-assumption
confidences, per-assumption entropies, the maximum confidence and entropy
deltas, the final sensitivity score, the rank, and the stability class
(`stable`, `root_sensitive_support`, or `root_sensitive_state`).

When the goal is to see whether ancestral conclusions survive topology
uncertainty instead of relying on one summary tree, use `ancestral tree-set`.
This workflow reruns ancestral reconstruction across every retained
posterior/bootstrap tree, maps comparable internal clades by descendant taxon
set, and writes both per-tree node ledgers and cross-tree clade summaries. The
surface supports continuous reconstruction under `brownian` or `ou` and
discrete reconstruction under `fitch`, `equal-rates`, `symmetric`, or
`all-rates-different`.

```bash
bijux-phylogenetics ancestral tree-set \
  artifacts/primates.posterior.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --kind discrete \
  --model equal-rates \
  --taxon-column species \
  --burnin-fraction 0.25 \
  --summary-out artifacts/primates.ancestral-tree-set-summary.tsv \
  --trees-out artifacts/primates.ancestral-tree-set-trees.tsv \
  --nodes-out artifacts/primates.ancestral-tree-set-nodes.tsv \
  --clades-out artifacts/primates.ancestral-tree-set-clades.tsv \
  --exclusions-out artifacts/primates.ancestral-tree-set-excluded.tsv \
  --json
```

The tree ledger keeps one row per retained tree with the original tree index,
post-burnin index, and rooted/unrooted topology identifiers so reviewers can
trace any unstable ancestral call back to the exact sampled topology. The node
ledger keeps one row per internal node per retained tree. For continuous
traits, each row preserves the reconstructed estimate, standard error, 95%
interval, confidence, and unstable flag. For discrete traits, each row
preserves the most likely state, candidate state set, confidence, ambiguity,
and unstable flag.

The clade ledger is the cross-tree review surface. It keys every comparable
clade by descendant taxa, reports how often that clade is present across the
retained trees, and summarizes whether its ancestral conclusion is stable or
topology-sensitive. Continuous clade rows preserve empirical value ranges and
mean uncertainty; discrete clade rows preserve dominant-state fractions and
state distributions. The summary ledger keeps one row with total tree counts,
retained tree counts, topology diversity, analyzed taxon counts, and unstable
clade counts. The excluded-taxa ledger keeps one row per dropped tip with an
explicit reason, so tree-set uncertainty is never mixed with silent trait-table
loss.

When the goal is to count inferred categorical state changes instead of only
report one ancestral state per node, use `ancestral transitions`. This
workflow reruns the owned discrete ancestral reconstruction surface, converts
each non-root branch into one explicit parent-versus-child state comparison,
and separates certain from uncertain changes instead of collapsing every
ambiguous branch into the same transition count. The same surface also supports
posterior or bootstrap tree sets through `--tree-set`, where it repeats the
counting path across every retained tree and summarizes which transition pairs
are stable versus topology-sensitive.

```bash
bijux-phylogenetics ancestral transitions \
  artifacts/primates.posterior.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --taxon-column species \
  --model equal-rates \
  --tree-set \
  --burnin-fraction 0.25 \
  --summary-out artifacts/primates.transition-summary.tsv \
  --trees-out artifacts/primates.transition-trees.tsv \
  --branches-out artifacts/primates.transition-branches.tsv \
  --counts-out artifacts/primates.transition-counts.tsv \
  --exclusions-out artifacts/primates.transition-excluded.tsv \
  --json
```

On one tree, the branch ledger keeps one row per non-root branch with the
parent state set, child state set, overlap, changed flag, and certainty class.
That makes it possible to distinguish a branch that certainly changed from one
that still overlaps across candidate states. The count ledger then collapses
those branch rows into one transition pair row per `source->target` direction
with separate certain and uncertain change totals. On a tree set, the tree
ledger keeps one retained-tree row with source and post-burnin indices plus
topology identifiers, the branch ledger keeps one branch row per retained
tree, and the count ledger adds tree-presence and stability summaries so
reviewers can see whether one reported transition depends on one sampled
topology or persists across the retained tree distribution.

When the goal is to inspect reconstructed histories directly on the tree, use
`ancestral render`. This workflow turns one continuous or discrete ancestral
reconstruction into a reviewer-facing figure and now supports three durable
export paths from the same owned visualization surface:
- SVG for publication or editing
- PNG for slideware and static sharing
- HTML for inline figure review with the governed visualization manifest

```bash
bijux-phylogenetics ancestral render \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --kind discrete \
  --model equal-rates \
  --discrete-node-style pies \
  --branch-coloring state \
  --layout phylogram \
  --out artifacts/primates.ancestral-visualization.html \
  --json
```

For continuous traits, the workflow keeps numeric internal-node labels and can
color branches by the descendant reconstructed value regime. For discrete
traits, the workflow supports either state labels or marginal-state pie charts
at internal nodes and can color branches by the inferred descendant state. The
HTML path writes the governed SVG as a sibling artifact and embeds that same
figure directly into the review page. The PNG path also writes the sibling SVG
so the raster export still preserves one editable source figure.

When the goal is to hand off a complete reviewer or publication bundle rather
than one standalone figure, use `ancestral package`. The package now writes the
same node-state and uncertainty ledgers as before, but it also emits:
- `ancestral-figure.svg`
- `ancestral-figure.png`
- `ancestral-figure.html`

That bundle keeps one rich visual surface alongside the node tables, legend,
caption, and manifest instead of forcing reviewers to regenerate alternate
formats from the SVG by hand.

When the goal is to deliver one complete ancestral reconstruction review
package instead of a figure-only bundle, use `ancestral report --out-dir`.
This workflow turns one continuous or discrete ancestral reconstruction into
one governed directory that keeps the reconstructed node ledger, the
uncertainty ledger, the transition or branch-change review ledgers, one
embedded HTML report, and the same SVG, PNG, and HTML tree visualization
surfaces in one handoff.

```bash
bijux-phylogenetics ancestral report \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --kind discrete \
  --model equal-rates \
  --taxon-column species \
  --out-dir artifacts/primates-ancestral-report \
  --json
```

The package always writes:
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

For discrete traits, the transition ledgers preserve directed ancestral state
changes using the owned transition-counting surface. For continuous traits, the
same filenames remain stable, but the transition ledgers become branch-change
review ledgers over reconstructed value deltas instead of pretending that a
continuous trajectory has categorical state transitions. The HTML report embeds
the governed SVG figure directly and keeps the summary, uncertainty, and
transition review surfaces visible in one place so reviewers do not need to
reassemble them manually.

When the goal is to review phylogenetic independent contrasts directly, use
`comparative contrasts`. The base workflow computes one standardized contrast
row per internal node for one numeric trait and can optionally fit one
regression-through-origin over matched predictor and response contrasts.

```bash
bijux-phylogenetics comparative contrasts \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --predictor-trait brain_mass_g \
  --taxon-column species \
  --contrasts-out artifacts/primates.independent-contrasts.tsv \
  --regression-out artifacts/primates.independent-contrast-regression.tsv \
  --json
```

The contrast ledger keeps one row per internal node with the left and right
descendant taxa, the standardized contrast, the expected variance, the local
ancestral value, and the shared root estimate for the analyzed trait. When
`--predictor-trait` is supplied, the regression ledger keeps one row per
matched node with the predictor contrast, response contrast, fitted
through-origin response contrast, residual, and leverage fraction, while the
JSON metrics report the fitted slope and p-value explicitly.

For Python-native comparative work, the same Brownian comparative core is also
available without restarting from file paths: use
`summarize_brownian_covariance_from_tree(...)` once you already hold one
`PhyloTree`, and use
`compute_phylogenetic_independent_contrasts_from_dataset(...)` once you already
hold one `ComparativeDataset`.

When the goal is to measure whether one numeric trait shows phylogenetic
structure rather than fit a regression, use `comparative signal`. That
workflow reports Blomberg's K, Pagel's lambda, a permutation p-value for the
observed K value, and a likelihood-ratio-style p-value for the fitted lambda
against the zero-signal lambda boundary. It also records whether the rooted
tree is ultrametric, keeps overlapping missing trait values as an explicit
pruning decision instead of a hidden cleanup step, reproduces permutation rows
from the supplied seed, and rejects constant post-pruning trait vectors with a
typed comparative-method error.

```bash
bijux-phylogenetics comparative signal \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --permutations 199 \
  --summary-out artifacts/primates.signal-summary.tsv \
  --permutations-out artifacts/primates.signal-permutations.tsv \
  --json
```

The summary ledger keeps one row with the fitted K, lambda, log-likelihood
context, both p-values, and the permutation exceedance count. The permutation
ledger keeps one row per shuffled trait realization so reviewers can see the
null K distribution directly instead of only one final exceedance count.
The JSON and summary surfaces now also keep Pagel-lambda optimizer diagnostics,
including the bounded grid-search name, function-evaluation count, lower versus
upper boundary hits, and one explicit fine-grid likelihood profile, so
reviewers can tell when the fitted lambda is an interior optimum versus one
boundary-supported solution.

When the goal is to review one categorical tip trait under a direct Mk
likelihood fit instead of one ancestral-state report, use
`comparative discrete-mk`. That workflow fits one ER, SYM, or ARD discrete Mk
surface on one rooted tree plus one categorical trait column, keeps the
missing-value pruning policy explicit, writes one flat fit summary plus one
directed rate-matrix ledger when requested, and exposes the same optimizer
diagnostics that the owned discrete Mk runtime uses internally. Its JSON
metrics now also expose overparameterization status plus ER baseline
comparison fields so model-comparison review does not require reopening the
TSV summary. The governed live `phytools` claim for this surface is explicit
and narrow: ER fits track `phytools::fitMk(model='ER')` on governed binary,
multistate, and missing-value-pruned fixtures, and unordered multistate SYM
fits track `phytools::fitMk(model='SYM')` on governed clean and
missing-value-pruned multistate fixtures. ARD fits now track
`phytools::fitMk(model='ARD')` on governed binary and missing-value-pruned
binary fixtures at full rate-row parity, while the governed multistate ARD
claim stays at summary parity when the optimizer reports weakly identified
boundary rates.

When the goal is to fit a standalone continuous-trait evolution model rather
than a regression, use `comparative brownian`. This workflow keeps the
Brownian motion surface explicit by reporting the fitted root state,
evolutionary rate `sigma²`, log-likelihood, AIC, and AICc, while preserving the
taxa pruned before fitting.

```bash
bijux-phylogenetics comparative brownian \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.brownian-summary.tsv \
  --excluded-taxa-out artifacts/primates.brownian-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted root state, `sigma²`, log-likelihood, AIC, AICc, and residual
diagnostic summary fields. The excluded-taxa ledger keeps one row per taxon
that was absent from the tree, absent from the trait table, missing the target
trait value, or pruned because the trait value was non-numeric. This makes the
Brownian fit auditable as a tree-plus-trait workflow rather than a detached
scalar estimate.

When the goal is to test whether named branches or clades evolve under
different Brownian rates, use `comparative brownian-regimes`. This workflow
fits one rate per user-supplied branch regime, compares that fit against the
single-rate Brownian baseline on the same analyzed taxon set, and exposes the
rate-profile uncertainty directly instead of reducing the outcome to one model
winner.

The regime map is a separate tabular input keyed by deterministic branch
identifiers. By default the workflow expects a `branch_id` column and a
`regime` column. Each `branch_id` is the normalized descendant-tip signature
for one non-root branch, such as `A|B` or `A|B|C|D`. Every non-root branch in
the tree must be assigned exactly one regime.

```bash
bijux-phylogenetics comparative brownian-regimes \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  artifacts/primates.branch-regimes.tsv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.brownian-regimes-summary.tsv \
  --rates-out artifacts/primates.brownian-regimes-rates.tsv \
  --comparison-out artifacts/primates.brownian-regimes-comparison.tsv \
  --profile-out artifacts/primates.brownian-regimes-profile.tsv \
  --branches-out artifacts/primates.brownian-regimes-branches.tsv \
  --excluded-taxa-out artifacts/primates.brownian-regimes-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted root state, multi-rate log-likelihood, AIC, AICc, the preferred
model between single-rate Brownian and regime-specific Brownian, and the
likelihood-ratio context for that comparison. The regime-rate ledger keeps one
row per regime with the fitted `sigma²`, conditional 95% profile interval, and
the amount of contributing branch length actually retained after taxon
filtering. The comparison ledger keeps both model-fit rows and the
Brownian-vs-regime likelihood-ratio row. The profile ledger keeps one row per
evaluated regime-rate candidate so weak identifiability can be reviewed
directly. The branch ledger keeps the normalized branch-to-regime assignment
that was actually used, including filtered descendant taxa, so the rate fit is
auditable end to end.

When the goal is to build or validate the branch regime map itself, use
`comparative regime-map`. This workflow either reconstructs branch regimes from
a discrete tip-state table or normalizes a user-provided branch regime map to
the repository's deterministic branch identity contract. It exists so
regime-aware downstream fits can start from one explicit review surface instead
of hiding the regime assignment step inside another model.

```bash
bijux-phylogenetics comparative regime-map \
  artifacts/primates.nwk \
  --table artifacts/primates.geography.tsv \
  --trait region \
  --taxon-column species \
  --summary-out artifacts/primates.regime-map-summary.tsv \
  --branches-out artifacts/primates.regime-map-branches.tsv \
  --nodes-out artifacts/primates.regime-map-nodes.tsv \
  --excluded-taxa-out artifacts/primates.regime-map-excluded.tsv \
  --svg-out artifacts/primates.regime-map.svg \
  --json
```

When `--table` is used, the workflow reconstructs one regime assignment per
node and then projects those assignments onto every non-root branch. The node
ledger preserves ambiguous internal states directly instead of flattening them
into one claimed historical narrative, while the branch ledger keeps the
normalized branch ids, chosen regime labels, candidate regimes, and filtered
descendant taxa that downstream comparative workflows would actually consume.
The SVG output renders the same regime assignment on the tree so reviewers can
inspect geography or ecology mapping visually instead of reading only a TSV.

When a branch regime map already exists, replace `--table ... --trait ...` with
`--regime-map artifacts/primates.branch-regimes.tsv`. In that mode the
workflow validates that every non-root branch is present exactly once, keeps
the normalized branch-to-regime ledger, and can still render the review SVG.
The governing branch identity is the descendant-tip signature for each
non-root branch, such as `A|B` or `A|B|C|D`, so user-provided maps and
downstream regime-aware models share one stable branch naming surface.

When the goal is to test whether a continuous trait is better explained by
constrained evolution toward an optimum, use `comparative ou`. This workflow
fits the stationary-root Ornstein-Uhlenbeck surface directly and reports the
fitted pull strength `alpha`, optimum `theta`, diffusion rate `sigma²`,
log-likelihood, AIC, and AICc, while preserving explicit pruning and
identifiability warnings.

```bash
bijux-phylogenetics comparative ou \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.ou-summary.tsv \
  --excluded-taxa-out artifacts/primates.ou-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted `alpha`, `theta`, `sigma²`, log-likelihood, AIC, AICc, and
residual diagnostic summary fields. The same row also preserves the
convergence-status label and the count of OU identifiability warnings. The
excluded-taxa ledger keeps one row per taxon that was absent from the tree,
absent from the trait table, missing the target trait value, or pruned because
the trait value was non-numeric. This keeps OU review grounded in both the fit
statistics and the data-retention contract.

When the goal is to test whether trait evolution decelerates through time in an
early-burst or adaptive-radiation style pattern, use `comparative
early-burst`. This workflow fits the bounded rate-change surface directly,
compares the retained fit against Brownian and OU on the same pruned taxon set,
and flags weak identifiability when the likelihood surface stays broad or
slides onto the search boundary.

```bash
bijux-phylogenetics comparative early-burst \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --summary-out artifacts/primates.early-burst-summary.tsv \
  --excluded-taxa-out artifacts/primates.early-burst-excluded.tsv \
  --comparison-out artifacts/primates.early-burst-comparison.tsv \
  --profile-out artifacts/primates.early-burst-profile.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, fitted `rate_change`, root state, `sigma²`, log-likelihood, AIC, AICc,
the selected best model among Brownian/OU/early-burst, and the count of
identifiability warnings. The excluded-taxa ledger keeps one row per taxon that
was absent from the tree, absent from the trait table, missing the target trait
value, or pruned because the trait value was non-numeric. The comparison ledger
keeps both model-fit rows and likelihood-ratio rows so reviewers can see
whether early-burst is actually preferred over Brownian or OU. The profile
ledger keeps one fixed `rate_change` row per bounded likelihood evaluation so
weak identifiability can be reviewed directly rather than inferred from one
point estimate.

When the goal is to inspect whether observed trait change is concentrated
deeper or shallower in the tree, use `comparative rate-through-time`. This
workflow reconstructs continuous ancestral states on the analyzed tree, bins
branch depth into explicit time intervals, and summarizes how much
branch-length-normalized reconstructed change falls into each interval. The
summary then classifies the pattern as `acceleration`, `slowdown`, `stable`, or
`insufficient_data` instead of leaving reviewers to infer that direction from a
raw table.

```bash
bijux-phylogenetics comparative rate-through-time \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait longevity \
  --taxon-column species \
  --interval-count 6 \
  --summary-out artifacts/primates.rate-through-time-summary.tsv \
  --intervals-out artifacts/primates.rate-through-time-intervals.tsv \
  --excluded-taxa-out artifacts/primates.rate-through-time-excluded.tsv \
  --json
```

The summary ledger keeps one row with the analyzed taxon count, excluded taxon
count, retained interval count, total tree depth, earliest and latest interval
rate estimates, the latest-to-earliest rate ratio, the weighted depth trend
slope, the normalized slope, and the classified trend direction. The interval
ledger keeps one row per depth bin with the bin bounds, contributing branch
length, contributing reconstructed change, and estimated rate so acceleration
or slowdown can be reviewed directly rather than inferred from only one final
label. The excluded-taxa ledger keeps one row per taxon that was absent from
the tree, absent from the trait table, missing the target trait value, or
pruned because the trait value was non-numeric.

When the goal is to ask whether particular internal clades have unusual trait
distributions, use `comparative clade-traits`. This workflow prunes the tree to
the analyzed taxa for one requested trait, summarizes every internal non-root
clade above a chosen size threshold, and ranks clades by reviewer-facing
exceptionality against the analyzed global trait distribution.

```bash
bijux-phylogenetics comparative clade-traits \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait habitat \
  --taxon-column species \
  --trait-kind categorical \
  --summary-out artifacts/primates.clade-traits-summary.tsv \
  --clades-out artifacts/primates.clade-traits.tsv \
  --excluded-taxa-out artifacts/primates.clade-traits-excluded.tsv \
  --json
```

For continuous traits the clade ledger keeps mean, median, minimum, maximum,
range width, and the mean shift from the analyzed global mean. For categorical
traits it keeps dominant state, dominant-state fraction, state counts, and the
distribution shift from the analyzed global state frequencies. In both cases
the rank is explicit and reviewer-facing rather than hidden behind a plot or an
informal scan of the tree.

The exceptionality score is intentionally descriptive. For continuous traits it
is a weighted standardized mean shift; for categorical traits it is a weighted
total-variation shift from the analyzed global state distribution. That makes
clade ranking auditable without overstating it as a formal one-size-fits-all
hypothesis test.

When the goal is to identify taxa whose continuous trait values are unusual
given the fitted phylogenetic covariance, use `comparative trait-outliers`.
This workflow fits Brownian and OU continuous-trait baselines, selects the
better standalone model by AICc, then scores each analyzed taxon by its
leave-one-taxon-out conditional residual under the selected covariance model.
The ranked taxon ledger also keeps local clade context so reviewers can see
whether an outlier departs from a nearby sister lineage or from a broader
cluster.

```bash
bijux-phylogenetics comparative trait-outliers \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait body_mass \
  --taxon-column species \
  --summary-out artifacts/primates.trait-outlier-summary.tsv \
  --outliers-out artifacts/primates.trait-outliers.tsv \
  --excluded-taxa-out artifacts/primates.trait-outlier-excluded.tsv \
  --json
```

The summary ledger keeps one row with the selected model, the selected mean
parameter name and value, the selected `sigma²`, the Brownian and OU AICc
values, the outlier threshold, the outlier count, and the top-ranked taxon.
The ranked taxon ledger keeps one row per analyzed tip with observed value,
conditional expected value, conditional variance, standardized residual, local
context clade, sibling context, and a stable reviewer-facing rank. The
excluded-taxa ledger keeps one row per taxon that was absent from the tree,
absent from the trait table, missing the target trait value, or pruned because
the trait value was non-numeric.

When the goal is to estimate missing continuous trait values from observed tips
and the phylogeny, use `comparative trait-imputation`. This workflow fits a
Brownian trait-evolution baseline on the observed taxa, imputes every missing
tree taxon from the conditional Brownian distribution, and then validates the
same prediction contract by leave-one-observed-out holdout runs.

```bash
bijux-phylogenetics comparative trait-imputation \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --trait body_mass \
  --taxon-column species \
  --summary-out artifacts/primates.trait-imputation-summary.tsv \
  --imputations-out artifacts/primates.trait-imputations.tsv \
  --holdout-out artifacts/primates.trait-imputation-holdout.tsv \
  --excluded-taxa-out artifacts/primates.trait-imputation-excluded.tsv \
  --json
```

The summary ledger keeps one row with the fitted Brownian root state,
`sigma²`, likelihood criteria, the count of imputed taxa, and the holdout
validation metrics. The imputation ledger keeps one row per missing tree taxon
with the predicted value, conditional variance, and 95% uncertainty interval.
The holdout ledger keeps one row per observed taxon with the leave-one-out
prediction error and interval coverage so the imputation contract can be
reviewed directly rather than assumed. The excluded-taxa ledger keeps one row
per taxon that could not enter the Brownian workflow because the trait value
was non-numeric or the row referred to a taxon absent from the tree.

When the goal is to fit the same comparative regression across several response
traits and then inspect how those fitted traits still co-vary, use
`comparative multivariate`. This workflow keeps one shared complete-case taxon
set across every requested response and predictor term, fits one comparative
regression per response on that exact taxon set, and then reports per-response
coefficients, per-response model summaries, residual covariance, residual
correlation, and residual trait-trait association explicitly. Predictor terms
are interpreted with the same comparative formula parser used by PGLS, so
categorical predictors, transformed numeric predictors such as `log(body_mass)`,
and explicit interaction terms such as `body_mass:habitat` remain governed
instead of falling back to raw-column guessing.

```bash
bijux-phylogenetics comparative multivariate \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --responses longevity range_size \
  --predictors brain_mass_g habitat brain_mass_g:habitat \
  --taxon-column species \
  --response-models-out artifacts/primates.multivariate-models.tsv \
  --coefficients-out artifacts/primates.multivariate-coefficients.tsv \
  --covariance-out artifacts/primates.multivariate-covariance.tsv \
  --correlation-out artifacts/primates.multivariate-correlation.tsv \
  --associations-out artifacts/primates.multivariate-associations.tsv \
  --excluded-taxa-out artifacts/primates.multivariate-excluded.tsv \
  --json
```

The missing-value policy is explicit: one taxon is retained only when every
requested response and every predictor term can be evaluated on that taxon. The
response-model ledger keeps one row per response with the fitted formula,
encoded-term count, taxon count, lambda value, log-likelihood, residual
variance, and residual degrees of freedom. The coefficient ledger keeps one row
per response-term coefficient with standard errors, test statistics, p-values,
and 95% intervals. The covariance ledger keeps one response-pair row with
residual covariance, residual correlation, pair count, and diagonal status. The
correlation ledger keeps the residual correlation matrix explicitly, separated
from covariance magnitudes. The association ledger keeps one unique
response-pair row with the same covariance and correlation plus a correlation
test statistic, p-value, and Fisher-style interval. The excluded-taxa ledger
makes the shared complete-case rule explicit by recording which taxa were
dropped, which responses they blocked, and which columns or terms failed.

The JSON report also preserves reviewer-facing warnings when the shared fit has
weak residual degrees of freedom or when the residual covariance matrix is
singular within the governed multivariate numerical tolerance. It also reports
residual covariance matrix rank, condition number, and singular-versus-near-
singular state directly, so reviewers do not need to infer matrix stability
only from the warning text or the raw covariance ledger.

When the goal is to hand a reviewer one durable comparative bundle rather than
separate regression, signal, contrast, and diagnostics outputs, use
`comparative report`. This workflow fits the integrated comparative surface and
materializes one HTML report plus flat TSV ledgers for the model formula,
coefficient table, residual summary, phylogenetic signal, model comparison,
biological interpretation, audit trail, and independent contrasts.

```bash
bijux-phylogenetics comparative report \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "longevity ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value estimate \
  --out-dir artifacts/primates.comparative-report \
  --json
```

The package directory keeps the comparative decision chain reviewable in one
place:

- `comparative-report.html` with integrated reviewer-facing narrative and tables
- `comparative-summary.tsv` with formula, selected model, signal, and fit summary
- `coefficient-table.tsv` with one explicit PGLS coefficient row
- `residual-summary.tsv` with one row per fitted residual-diagnostic surface
- `signal-summary.tsv` with Blomberg's K, Pagel's lambda, and contrast counts
- `model-comparison.tsv` with the Brownian-versus-OU fit comparison
- `interpretation-table.tsv` with explicit claim, evidence, and caution rows
- `audit-table.tsv` with taxa used, excluded taxa, assumptions, and warnings
- `contrast-table.tsv` with one independent-contrast row per internal node
- `comparative-report.manifest.json` with checksums and output inventory

This surface is intentionally different from a generic HTML export. It keeps
the regression formula, coefficient uncertainty, residual diagnostics,
phylogenetic signal, and model comparison together so biological interpretation
can be reviewed against the exact evidence rather than reconstructed later from
separate commands.

When the goal is to detect whether one fitted comparative model is leaving
systematic residual structure inside particular subtrees, use
`comparative clade-residuals`. This workflow keeps the fitted taxon-level
residuals from one comparative model and then aggregates those residuals across
every internal non-root clade in the analyzed tree.

```bash
bijux-phylogenetics comparative clade-residuals \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "longevity ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value estimate \
  --taxa-out artifacts/primates.residual-taxa.tsv \
  --clades-out artifacts/primates.residual-clades.tsv \
  --json
```

The taxon ledger keeps one row per analyzed taxon with the observed value,
fitted value, raw residual, and standardized residual used for clade
aggregation. The clade ledger keeps one row per internal clade with its member
taxa, residual averages, residual sum-of-squares share, influence score,
residual-heavy flag, and rank. That makes it possible to distinguish one
isolated outlier taxon from a whole subtree that is carrying consistent model
misspecification.

When the question is not residual burden after one fit but whether the fitted
comparative conclusion itself depends on one major subtree, use
`comparative clade-stability`. This workflow fits one comparative model on the
baseline analyzed taxon set, derives major internal non-root clades from that
baseline tree, removes each candidate clade in turn, and refits the same model
on the retained taxa.

```bash
bijux-phylogenetics comparative clade-stability \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula "longevity ~ brain_mass_g + habitat" \
  --taxon-column species \
  --lambda-value estimate \
  --clades-out artifacts/primates.clade-stability.tsv \
  --terms-out artifacts/primates.clade-stability-terms.tsv \
  --json
```

The clade ledger keeps one row per candidate removal with the dropped clade,
retained taxon count, fit status, blocked-refit reason when one exists,
coefficient comparison count, sign-change count, significance-change count,
largest coefficient and p-value shifts, delta log-likelihood, and the final
influence rank. The term ledger then keeps one row per comparable coefficient
for each successful clade removal, including baseline and refit estimates,
baseline and refit p-values, and explicit sign-change and
significance-change flags.

This surface matters when a biological interpretation looks fragile. A large
or influential clade can drive coefficient direction, apparent significance,
or fit quality without leaving obviously extreme single-taxon residuals. The
workflow keeps those subtree dependencies explicit and preserves blocked rows
instead of silently skipping removals that would collapse the remaining fit.

When the question is not whether one clade drives the fit but whether topology
uncertainty across a posterior or bootstrap tree set changes the comparative
conclusion, use `comparative posterior-pgls`. This workflow keeps one
continuous-trait PGLS formula fixed, applies it to every retained tree in a
tree set, and then summarizes how the coefficients and support calls vary
across those trees.

```bash
bijux-phylogenetics comparative posterior-pgls \
  artifacts/primates.posterior.trees \
  artifacts/primates.csv \
  --formula "longevity ~ social_group_size" \
  --taxon-column species \
  --lambda-value estimate \
  --burnin-fraction 0.25 \
  --significance-threshold 0.05 \
  --trees-out artifacts/primates.posterior-pgls-trees.tsv \
  --coefficients-out artifacts/primates.posterior-pgls-coefficients.tsv \
  --summary-out artifacts/primates.posterior-pgls-summary.tsv \
  --json
```

The per-tree ledger keeps one row per retained tree with its post-burn-in
position, topology identifiers, fitted lambda, and log-likelihood. The
coefficient ledger then keeps one row per coefficient per retained tree,
including estimate, p-value, direction, and whether the coefficient met the
chosen support threshold on that tree. The summary ledger collapses those
coefficient rows into reviewer-facing distributions with empirical estimate
ranges, direction consistency, support fractions, and a conclusion-stability
classification such as `stable_supported`, `stable_unsupported`,
`mixed_support`, or `direction_conflict`.

This surface matters when one MCC tree or one consensus tree would hide the
fact that coefficient support is topology-sensitive. It lets users propagate
tree uncertainty into the comparative conclusion directly instead of treating
the posterior tree set and the trait model as separate review steps.

When the covariance assumption itself must stay fixed to Brownian shared branch
lengths, use `comparative brownian-pgls`. This keeps the regression surface
separate from Pagel-lambda fitting and writes an explicit covariance ledger
when requested.

```bash
bijux-phylogenetics comparative brownian-pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size \
  --taxon-column species \
  --covariance-out artifacts/primates.brownian-covariance.tsv \
  --json
```

That workflow fits PGLS under the raw Brownian shared-path covariance implied
by the rooted tree, checks that the unstabilized covariance is positive
definite, and records one pairwise row per taxon pair with shared path length
and root-to-tip depths. Rooted ultrametric and rooted non-ultrametric trees are
both supported as long as the Brownian covariance stays valid. If zero or
negative branch lengths make the covariance invalid, the workflow fails
explicitly instead of silently regularizing the tree into a different
scientific assumption.

When the residual structure is expected to reflect pull toward an optimum rather
than unconstrained Brownian diffusion, use `comparative ou-pgls`. The workflow
supports either a fixed positive `--alpha` or a governed `estimate` mode that
profiles alpha across the bounded search grid already used by the OU trait-model
surface.

```bash
bijux-phylogenetics comparative ou-pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size \
  --taxon-column species \
  --alpha estimate \
  --covariance-out artifacts/primates.ou-covariance.tsv \
  --alpha-profile-out artifacts/primates.ou-alpha-profile.tsv \
  --json
```

That workflow fits the comparative regression under stationary-root OU
covariance, records the fitted log-likelihood and AIC, and writes two optional
review ledgers. The covariance ledger keeps one pairwise row per taxon pair
with the implied OU covariance, shared path length, and root-depth context. The
alpha profile ledger keeps one row per candidate alpha value, the
log-likelihood drop from the best-supported fit, and whether the row stays
inside the likelihood-ratio-supported 95% interval. Fixed-alpha review and
estimated-alpha review therefore stay explicit instead of being folded into one
opaque regression number.

When the main question is how strongly the regression residuals prefer
phylogenetic covariance, write the governed lambda profile in the same run.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size \
  --taxon-column species \
  --lambda-value estimate \
  --lambda-profile-out artifacts/primates.lambda-profile.tsv \
  --json
```

The profile ledger keeps one row per candidate lambda value across the bounded
search surface, records the log-likelihood drop from the best fit, and marks
which rows stay inside the likelihood-ratio-supported 95% confidence interval.

When the main question is whether the requested biological formula expands into
the encoded predictors you actually expect, use formula syntax directly and
write the governed model matrix. Intercept-free formulas use the standard
comparative spellings `0 + ...` or `... - 1`.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula 'longevity ~ 0 + social_group_size * habitat' \
  --taxon-column species \
  --model-matrix-out artifacts/primates.model-matrix.tsv \
  --json
```

The written matrix keeps one row per analyzed taxon, the response value, and
one encoded column per fitted predictor term. That makes it possible to review
continuous predictors, categorical indicator columns, interaction columns, and
intercept inclusion before interpreting the fitted coefficients.

When categorical predictors are present and the question is how each biological
group was encoded and interpreted, write the governed categorical-contrast
ledger in the same run.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --response longevity \
  --predictors social_group_size habitat \
  --taxon-column species \
  --categorical-contrasts-out artifacts/primates.categorical-contrasts.tsv \
  --json
```

That ledger keeps one row for the baseline group when treatment coding is used
and one row for every estimated non-baseline group coefficient. It also records
missing-category taxa explicitly, so a dropped or blank category value does not
disappear into a generic exclusion count.

When interaction terms are present and the question is how effect modification
was encoded and estimated, write the governed interaction-coefficient ledger in
the same run.

```bash
bijux-phylogenetics comparative pgls \
  artifacts/primates.nwk \
  artifacts/primates.csv \
  --formula 'longevity ~ social_group_size * habitat' \
  --taxon-column species \
  --interaction-coefficients-out artifacts/primates.interaction-coefficients.tsv \
  --json
```

That ledger keeps one row per fitted interaction coefficient, records whether
the interaction is continuous-by-continuous, continuous-by-categorical, or
categorical-by-categorical, and preserves any omitted treatment-coded baseline
levels so the effect-modification interpretation stays explicit.

When a tree is already inferred and you need to root it on one known outgroup
taxon or one expected outgroup clade, use `topology root-outgroup`. The
command writes the rooted tree and can also emit a one-row TSV report that
records which requested taxa were found, which requested taxa were absent,
whether the matched outgroup is monophyletic in the input tree, which extra
taxa fall inside the matched outgroup MRCA when it is not monophyletic, and
which taxa end up isolated on the rooted outgroup side.

```bash
bijux-phylogenetics topology root-outgroup \
  artifacts/mammals.unrooted.nwk \
  --taxa Ornithorhynchus_anatinus Tachyglossus_aculeatus \
  --out artifacts/mammals.rooted.nwk \
  --report-out artifacts/mammals.rooting.tsv \
  --json
```

If the requested outgroup taxa are missing, the TSV and JSON report that
explicitly instead of silently dropping them. If the requested outgroup taxa
are not monophyletic in the input tree, the workflow still records the rooted
tree but also reports the MRCA spillover taxa and a warning that the root does
not cleanly isolate the requested outgroup as one coherent clade.

When no explicit outgroup is available, use `topology reroot-midpoint` for an
exploratory rooted tree. The command writes the rerooted tree and can also
emit a one-row TSV report that records the anchor tip pair that defined the
selected midpoint path, the total tip-to-tip path length, the midpoint
distance from the anchor tip, whether the midpoint landed on an original node
or within an original branch, and which taxa ended up on each side of the new
root.

```bash
bijux-phylogenetics topology reroot-midpoint \
  artifacts/mammals.unrooted.nwk \
  --out artifacts/mammals.midpoint-rooted.nwk \
  --report-out artifacts/mammals.midpoint-rooting.tsv \
  --json
```

The midpoint report also records whether the input tree was suitable for
straightforward midpoint interpretation. Trees that are not strictly
bifurcating are still rerooted, but the report marks them as exploratory and
adds an explicit warning so downstream review does not overclaim the result.

When the goal is to inspect every clade directly instead of transforming the
tree, use `topology clades`. That command writes one row per node-derived
clade, including tips, internal clades, and the root. Each row keeps the
member taxa, any parsed support label, the incoming branch length, root depth,
descendant tip depths, and `node_age` when the tree is branch-length complete
and ultrametric.

```bash
bijux-phylogenetics topology clades \
  artifacts/mammals.supported.nwk \
  --metadata artifacts/mammals.metadata.tsv \
  --metadata-column species \
  --metadata-column location \
  --out artifacts/mammals.clades.tsv \
  --json
```

When `--metadata` is supplied, Bijux treats the table as a taxon-keyed
metadata or trait table and flattens the requested columns into per-clade
review fields such as `species_values`, `species_distinct_values`, and
`species_missing_taxa`. That keeps trait inspection tied directly to clade
membership instead of forcing a separate join step.

When the goal is to inspect tree shape directly rather than individual clades,
use `topology shape`. That command writes one review row with balance and shape
metrics including Sackin imbalance, Colless imbalance where defined, cherry
count, tree height in edges, branch-length tree height where available, and the
governed `balanced`, `skewed`, or `ladderized` summary.

```bash
bijux-phylogenetics topology shape \
  artifacts/mammals.supported.nwk \
  --out artifacts/mammals.shape.tsv \
  --json
```

The JSON output also preserves whether the tree is star-like, comb-like, or
unusually imbalanced. That keeps obvious ladderization or star-topology risks
visible without forcing users to infer them from raw node depths alone.

When the goal is to inspect branch-length patterns directly, use
`topology branch-lengths`. That command writes one row per non-root branch and
summarizes minimum, maximum, mean, and median branch length alongside explicit
zero-length, negative-length, and outlier flags.

```bash
bijux-phylogenetics topology branch-lengths \
  artifacts/mammals.supported.nwk \
  --out artifacts/mammals.branch-lengths.tsv \
  --json
```

This surface matters when the question is not just whether branch lengths exist
but whether they contain odd zeros, extreme long branches, or other scale
distortions that can mislead downstream interpretation.

When you already have two inferred trees and need an explicit topology
distance, use `compare`. The command now supports rooted and unrooted
Robinson-Foulds review directly, and it records which taxa were shared versus
present on only one side before the distance is computed.

```bash
bijux-phylogenetics compare \
  artifacts/mammals.iqtree.nwk \
  artifacts/mammals.fasttree.nwk \
  --rf-mode unrooted \
  --json
```

Use rooted RF when root placement is part of the scientific claim. Use
unrooted RF when the question is only whether the same splits were recovered
regardless of root placement. By default the command prunes both trees to
their shared taxa before computing RF distance, which is appropriate for
reviewing partially overlapping outputs. If taxon-set mismatch itself should
fail the comparison, add `--taxon-overlap-policy require-identical`.

```bash
bijux-phylogenetics compare \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --rf-mode rooted \
  --taxon-overlap-policy require-identical \
  --json
```

When the first question is whether two trees can be reduced to the same taxon
set safely before any deeper comparison, use `compare prune`. That command
writes the two pruned trees plus a governed pruning review bundle in one output
directory.

```bash
bijux-phylogenetics compare prune \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --out artifacts/mammals.shared-taxa \
  --json
```

The pruning bundle keeps the evidence explicit:
- `left-shared.nwk` and `right-shared.nwk` are the retained shared-taxon trees
- `shared-taxa-pruning.tsv` records one row per input tree with retained and
  removed taxa, branch-length audit fields, and information-loss fields
- `shared-taxa-removed.tsv` records one row per removed taxon with the side and
  removal reason
- `shared-taxa-comparison.tsv` compares the retained trees directly so the
  post-pruning topology review is durable instead of implicit

The JSON payload also preserves the full left and right pruning audits plus a
`post_pruning_comparison` report. That makes it possible to inspect branch
length preservation, removed taxa, and retained-tree topology without
reconstructing the review from separate commands.

When branch lengths matter as well as topology, use
`compare branch-lengths`. That surface preserves the per-split length table
for shared rooted clades and also computes Felsenstein branch-score distance on
the union of unrooted splits. Missing splits count as zero-length branches for
the score calculation, which makes topology disagreement contribute directly to
the final distance instead of disappearing from the branch-length review.

```bash
bijux-phylogenetics compare branch-lengths \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --json
```

By default the command prunes both trees to their shared taxa before computing
branch-score distance. If taxon-set mismatch itself should fail the review, add
`--taxon-overlap-policy require-identical`. When a matched split is present but
lacks a branch length on one side, the per-split ledger records that missing
length explicitly and the branch-score summary becomes unavailable instead of
silently treating the missing value as zero.

When the goal is to check whether the runtime still agrees with the governed
tree-distance reference corpus rather than compare one ad hoc pair of trees,
use `topology distance-reference`.

```bash
bijux-phylogenetics topology distance-reference --json
```

This surface reruns the checked hard cases for rooted RF, unrooted RF,
normalized RF, and branch-score distance. The fixture set includes binary
disagreement, rooting-only disagreement, polytomies, star-tree collapse, and
shared-taxa pruning, and it keeps the
`--taxon-overlap-policy require-identical` rejection path explicit for both RF
and branch-score comparisons.

When the goal is to verify that support values still stay attached to the right
branches across parser and topology surfaces, use `topology support-reference`.

```bash
bijux-phylogenetics topology support-reference --json
```

This surface reruns checked support cases for plain IQ-TREE UFBoot labels,
compound SH-aLRT/UFBoot labels, FastTree local support, and posterior clade
frequencies. It also keeps two policy checks explicit:

- support comparison must stay clade-mapped when two trees differ only by node
  ordering
- bootstrap-versus-posterior review must report topology mismatch explicitly
  when one support surface contains clades the other does not

When the main question is whether topology conflicts are serious or only weakly
supported, use `compare support`. That surface combines clade presence with the
support values parsed from each tree and writes one support-aware conflict
ledger with `--out`.

```bash
bijux-phylogenetics compare support \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --out artifacts/mammals.support-conflicts.tsv \
  --json
```

The JSON and TSV outputs separate three situations explicitly:
- shared clades, where both trees carry the same clade and Bijux reports the
  normalized support delta
- high-support conflicts, where a clade present in only one tree still carries
  normalized support of at least `0.9`
- low-support disagreements, where a conflicting clade is present with support
  below `0.7`

Conflicting clades between `0.7` and `0.9` are preserved as
`moderate_support_disagreement` so the report does not flatten moderate support
into either strong conflict or weak noise. If the present-side tree did not
carry a parseable support label, the conflict row is marked
`support_unavailable`.

When the question is which shared taxa are driving disagreement, use
`compare influence`. That workflow performs a leave-one-taxon-out comparison
for every shared taxon, recomputes the topology and support conflict surfaces
after each exclusion, and ranks taxa by how much the disagreement surface
changes.

```bash
bijux-phylogenetics compare influence \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.candidate.nwk \
  --out artifacts/mammals.taxon-influence.tsv \
  --json
```

The written ledger keeps one row per excluded taxon and preserves:
- whether the baseline rooted conflict disappeared or persisted after exclusion
- rooted and unrooted Robinson-Foulds deltas after the taxon was removed
- support-disagreement, conflicting-clade, and high-support-conflict deltas
- a transparent influence score and resulting rank

This review surface is useful when one taxon appears unstable, misplaced, or
poorly aligned and the practical question is whether disagreement is global or
concentrated around that single tip. The ranking is deliberately heuristic: it
adds the normalized topology shift to the absolute support-surface count shifts
so reviewers can see why one taxon outranks another instead of treating the
rank as an opaque diagnostic.

When the main question is which clades agree or conflict across several trees,
use `compare clades`. The command accepts two required tree paths plus
additional trees through repeated `--tree` flags, computes clade overlap on the
shared taxon set, preserves support values where the underlying tree labels
carry them, and can write one flat clade-overlap table with `--out`.

```bash
bijux-phylogenetics compare clades \
  artifacts/mammals.reference.nwk \
  artifacts/mammals.bootstrap.nwk \
  --tree artifacts/mammals.fasttree.nwk \
  --out artifacts/mammals.clade-overlap.tsv \
  --json
```

The JSON summary separates clades present in every tree from clades that are
conflicting or tree-specific. The written table then gives one row per
clade-per-tree observation, including whether that clade is present in the
given tree and which support value was observed when the tree format exposed
one.

When the goal is not overlap comparison but a direct ledger of every clade from
every sampled tree, use `tree-set clades`. That command writes one row per
clade observation per tree, preserving the tree index alongside the same clade
membership, support, branch-length, depth, and optional metadata summaries
used by the single-tree surface.

```bash
bijux-phylogenetics tree-set clades \
  artifacts/mammals.posterior.nwk \
  --metadata artifacts/mammals.metadata.tsv \
  --metadata-column species \
  --out artifacts/mammals.posterior-clades.tsv \
  --json
```

This surface is for reviewer-readable extraction, not just summary counts. It
is especially useful before consensus-building because it keeps minority and
tree-specific clades visible instead of collapsing them into one consensus or
frequency table immediately.

When the question is how shape varies across many sampled trees rather than how
individual clades vary, use `tree-set shape`. That command writes one row per
tree with the same balance and height metrics as `topology shape`, and its JSON
aggregate counts how many sampled trees are balanced, skewed, ladderized,
star-like, or comb-like.

```bash
bijux-phylogenetics tree-set shape \
  artifacts/mammals.posterior.nwk \
  --out artifacts/mammals.posterior-shape.tsv \
  --json
```

This surface is for tree-set review before or alongside consensus-building. It
keeps the distribution of imbalance and ladderization explicit instead of
reducing a posterior or bootstrap set to one representative tree immediately.

When the question is how branch-length distributions vary across many sampled
trees, use `tree-set branch-lengths`. That command writes one row per
non-root branch per tree and preserves the source tree index so long branches,
zero branches, and missing lengths can be traced back to the individual sample
that produced them.

```bash
bijux-phylogenetics tree-set branch-lengths \
  artifacts/mammals.posterior.nwk \
  --out artifacts/mammals.posterior-branch-lengths.tsv \
  --json
```

The JSON aggregate keeps the set-wide minimum, maximum, median, zero-length
count, negative-length count, and long-outlier count explicit. That makes it
possible to review whether one sampled tree is distorting the branch-length
distribution instead of assuming the full set is numerically homogeneous.

When the input is a bootstrap replicate tree file and the goal is one governed
review bundle rather than a chain of separate tree-set commands, use
`tree-set bootstrap-summary`. The command reads the bootstrap trees directly,
builds a consensus tree, computes clade frequencies and topology diversity, and
exports a dedicated unstable-branch ledger for consensus branches that are not
yet robust.

```bash
bijux-phylogenetics tree-set bootstrap-summary \
  artifacts/mammals.bootstrap.ufboot \
  --out-dir artifacts/mammals-bootstrap-review \
  --prefix mammals-bootstrap \
  --json
```

The written artifact set includes:
- `mammals-bootstrap.summary.tsv`
- `mammals-bootstrap.consensus.nwk`
- `mammals-bootstrap.clade-frequencies.tsv`
- `mammals-bootstrap.unstable-branches.tsv`
- `mammals-bootstrap.unstable-clades.tsv`
- `mammals-bootstrap.rf-distribution.tsv`
- `mammals-bootstrap.distance-matrix.tsv`
- `mammals-bootstrap.topology-clusters.tsv`

This bootstrap-specific bundle matters because majority-rule consensus alone can
hide whether its retained branches are only weakly recovered or compete with
clear alternative clades across the replicate set. The unstable-branch ledger
keeps that distinction explicit.

For large posterior or bootstrap tree sets where the full pairwise matrix is no
longer the right first review surface, use `tree-set diversity` to emit the
compact rooted RF distribution directly. The JSON payload and the written
summary ledgers from both `tree-set diversity` and `tree-set bootstrap-summary`
record `runtime_seconds`, `peak_memory_bytes`, and
`skipped_malformed_tree_count`, so 1,000-tree review stays explicit about cost
and malformed-record loss instead of hiding them behind a successful exit code.

```bash
bijux-phylogenetics tree-set diversity \
  artifacts/mammals.bootstrap.ufboot \
  --out artifacts/mammals-bootstrap.rf-distribution.tsv \
  --json
```

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

Use `adapter mrbayes-prepare` when you need a Bayesian NEXUS input file for
one aligned matrix and, optionally, one same-datatype partition file. The
command writes a ready-to-run MrBayes file with the data block, charset
definitions, partition declaration, model commands, and MCMC settings in one
step.

```bash
bijux-phylogenetics adapter mrbayes-prepare \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --partitions artifacts/mixed-locus-supermatrix.partitions.txt \
  --out artifacts/multilocus-bayesian.nex \
  --model gtr \
  --rates gamma \
  --ngen 200000 \
  --samplefreq 100 \
  --printfreq 100 \
  --json
```

When partitions are supplied, the preparation surface validates that the
coordinates fit the alignment and that any declared partition datatypes still
match the inferred alignment alphabet. The generated MrBayes block uses named
charsets plus one active partition declaration so the resulting file is
accepted by MrBayes as a partitioned Bayesian analysis input instead of as a
flat unpartitioned matrix.

Those direct adapter runs now preserve IQ-TREE's own `.iqtree` and `.log`
artifacts, plus the model-selection sidecar, a generated
`.model-candidates.tsv`, and `.treefile`, `.ufboot`, or `.contree` where the
invoked step produces them. The emitted JSON and manifests also include the
parsed `selected_model`, `selected_criterion`, `candidate_model_count`,
`best_model_aic`, `best_model_aicc`, `best_model_bic`, `log_likelihood`, and
support-value counts so reviewers can verify the ML result against structured
fields instead of manually scraping engine text files.

The governed external execution adapters now share one explicit restart and
failure-control contract. `adapter align`, `trim`, `model-select`, `infer-ml`,
`bootstrap`, `sh-alrt`, `fasta-to-tree`, `consensus`, `infer-fast`,
`infer-large`, `compare-engines`, `mrbayes-run`, and `beast-run` all accept:

- `--timeout-seconds` to stop one engine execution after a wall-clock budget
- `--resume` to reuse one completed manifest only when the recorded command, input checksums, and output checksums still match
- `--incomplete-run-policy reject|clean` to either stop on stale incomplete outputs or remove them before a fresh rerun

That contract is intentionally strict. A resume is allowed only for verified
completed runs. Failed, timed-out, killed, or malformed-output runs leave an
explicit incomplete-run marker, and `--incomplete-run-policy clean` is the
governed way to discard that partial state before rerunning. A missing
executable stops before any incomplete-run marker is written because no engine
run started.

Workflow success is also strict about output completeness. Non-empty files are
not enough by themselves: alignment workflows must emit valid retained-site
matrices, IQ-TREE workflows must emit the required tree, log, report, model,
and support artifacts for the chosen mode, FastTree must emit parseable local
support labels, and BEAST or MrBayes must emit their full posterior artifact
sets. Missing required files, empty required files, missing model results, and
missing support annotations stop the workflow with stable structured error
codes before manifests or reviewer-facing reports are written.

Use `adapter beast-prepare` when you need a real BEAST2 XML template from one
aligned matrix plus optional dating metadata. The command writes a BEAST2-style
XML file with an explicit alignment block, a starting tree, a strict or
uncorrelated lognormal clock, a Yule or birth-death tree prior, MCMC loggers,
and one MRCA prior per validated calibration target.

```bash
bijux-phylogenetics adapter beast-prepare \
  artifacts/mixed-locus-supermatrix.aln.fasta \
  --tree artifacts/mixed-locus-guide.nwk \
  --calibrations artifacts/fossil-calibrations.tsv \
  --tip-dates artifacts/tip-dates.tsv \
  --out artifacts/multilocus-beast.xml \
  --clock-model relaxed-lognormal \
  --tree-prior birth-death \
  --chain-length 2000000 \
  --log-every 1000 \
  --json
```

When `--calibrations` or `--tip-dates` are supplied, `--tree` is required so
the prepared XML can use the same taxa and named clades that were validated
during preparation. If the alignment is nucleotide or RNA, the template uses
an HKY site model; if the alignment is protein, it uses a JTT site model. The
written XML also names the default BEAST run logs as `STEM.$(seed).log` and
`STEM.$(seed).trees` so later execution produces a predictable artifact set
beside the XML.

Use `adapter beast-xml` when you need a governed review surface over one
prepared BEAST XML before execution. The command parses the XML, checks that
the alignment block, site model, clock model, tree prior, MCMC run, and
posterior loggers are all present, and reports the assumed substitution model,
clock model, tree prior, starting-tree source, chain length, calibration
count, and dated-tip count directly.

```bash
bijux-phylogenetics adapter beast-xml \
  artifacts/multilocus-beast.xml \
  --json
```

That JSON payload is the durable preparation evidence surface. It turns the
XML into a reviewable assumptions contract instead of leaving prior choice,
clock choice, chain settings, and dated-tree inputs buried inside engine
markup.

Calibration handling is explicit rather than silent. If a calibration table
already provides both lower and upper bounds, the template preserves those hard
bounds as a BEAST uniform prior. If a calibration provides only a lower bound,
the template emits an explicit offset density above that minimum bound and
records a preparation warning in JSON so reviewers can see that the prior shape
was translated for template generation instead of copied as a literal hard
uniform interval.

Tip-dated workflows also surface one tree-prior boundary explicitly. If you
ask for the standard `birth-death` prior together with tip dates, the template
still writes valid XML, but the JSON warnings mark that combination as
exploratory because BEAST's own validator reports that the standard birth-death
prior is not serial-sampling aware.

Use `adapter beast-run` after preparation when you want governed BEAST
execution with explicit timeout, resume, and incomplete-run handling.

```bash
bijux-phylogenetics adapter beast-run \
  artifacts/multilocus-beast.xml \
  --threads 1 \
  --seed 7 \
  --timeout-seconds 1800 \
  --resume \
  --incomplete-run-policy clean \
  --json
```

That workflow writes the posterior log and posterior tree file beside the XML
using the governed `STEM.SEED.log` and `STEM.SEED.trees` names. The JSON
summary reports the thread count, seed, overwrite mode, resume status, timeout
budget, and warning count so BEAST execution can be reviewed as a structured
run instead of as an opaque shell invocation.

Use `adapter beast-log` once a BEAST run has produced a posterior log and you
need a governed review surface instead of manual spreadsheet work. The command
parses native BEAST log files, accepts BEAST's own `Sample` state header, and
can discard an explicit burn-in fraction before computing summary statistics
and effective sample sizes.

```bash
bijux-phylogenetics adapter beast-log \
  artifacts/multilocus-beast.1.log \
  --burnin-fraction 0.1 \
  --summary-out artifacts/multilocus-beast.log-summary.tsv \
  --json
```

The JSON output keeps the raw parsed log plus a summary block that separates
posterior, likelihood, prior, clock, and tree-related parameters into explicit
lists. The optional `--summary-out` table writes one row per sampled parameter
with the post-burn-in mean, median, sample standard deviation, 95% HPD
interval, min/max range, first-half and second-half means, standardized drift,
and effective sample size so reviewers can audit BEAST log behaviour without
re-parsing the engine file themselves.

When a sibling preparation XML is available, the downstream Bayesian methods
and diagnostics surfaces reuse that XML to state the actual BEAST substitution
model, clock model, tree prior, calibration count, tip-date count, and chain
settings explicitly in generated reviewer text.

Use `adapter beast-parameters` when the main goal is posterior parameter
diagnostics rather than raw log parsing.

```bash
bijux-phylogenetics adapter beast-parameters \
  artifacts/multilocus-beast.1.log \
  --burnin-fraction 0.1 \
  --summary-out artifacts/multilocus-beast.parameter-summary.tsv \
  --json
```

That command emits one burn-in-aware JSON report and, when requested, one TSV
table covering posterior mean, median, sample standard deviation, 95% HPD
interval, effective sample size, and retained state window for every parameter
that survived the requested burn-in cut.

Use `adapter beast-convergence` when you want the same burn-in handling applied
to convergence warnings directly.

```bash
bijux-phylogenetics adapter beast-convergence \
  artifacts/multilocus-beast.1.log \
  --burnin-fraction 0.1 \
  --ess-threshold 200 \
  --mean-shift-threshold 0.5 \
  --json
```

That command emits structured warnings for low ESS and mean drift after the
requested burn-in has been discarded, so convergence review stays aligned with
the same retained posterior window described in the summary table.

Use `adapter beast-burnin-sensitivity` when you want to compare posterior
parameter estimates and clade probabilities across the governed burn-in review
fractions `5%`, `10%`, `25%`, and `50%`, or across an explicit custom set.

```bash
bijux-phylogenetics adapter beast-burnin-sensitivity \
  artifacts/multilocus-beast.1.trees \
  --log artifacts/multilocus-beast.1.log \
  --slice-out artifacts/multilocus-beast.burnin-slices.tsv \
  --parameter-out artifacts/multilocus-beast.burnin-parameters.tsv \
  --clade-out artifacts/multilocus-beast.burnin-clades.tsv \
  --json
```

That workflow writes one per-fraction slice ledger plus separate parameter and
clade comparison ledgers. Parameter instability is reported when the tested
95% HPD intervals fail to share a common overlap. Clade instability is
reported when a clade posterior probability crosses the majority-rule
threshold across the tested burn-in fractions.

Use `adapter beast-trees` once a BEAST run has produced a posterior tree file
and you need a governed tree-sample surface rather than ad hoc NEXUS handling.
The command parses native `.trees` files, keeps the sampled `STATE_*`
generations, applies an explicit burn-in fraction, extracts clade frequencies,
and can write the retained trees back out as normalized Newick for downstream
MCC and consensus workflows.

```bash
bijux-phylogenetics adapter beast-trees \
  artifacts/multilocus-beast.1.trees \
  --burnin-fraction 0.1 \
  --tree-set-out artifacts/multilocus-beast.postburnin.nwk \
  --json
```

The JSON summary reports the total sampled tree count, the discarded burn-in
count, the retained tree count, rooted-tree count, sampled state count, and
the extracted clade table. The optional `--tree-set-out` file is a normalized
Newick tree set that can be handed directly to the existing tree-set
consensus, topology-diversity, and MCC review surfaces without re-scraping the
original BEAST NEXUS container.

Use `adapter beast-subsample` when the posterior tree file is already too large
for the next review step, but you still want a governed retained subset instead
of ad hoc copying. The command supports either evenly spaced thinning or a
seeded random subset and preserves the native `STATE_*` labels in the retained
sample ledger.

```bash
bijux-phylogenetics adapter beast-subsample \
  artifacts/multilocus-beast.1.trees \
  --method evenly-spaced \
  --burnin-fraction 0.1 \
  --thinning-interval 5 \
  --tree-set-out artifacts/multilocus-beast.subsample.nwk \
  --sample-table-out artifacts/multilocus-beast.subsample.tsv \
  --json
```

For random retained subsets, switch to `--method random`, pass
`--sample-count`, and set `--seed` explicitly when you need the same retained
trees again later. The retained tree-set file stays normalized Newick, while
the TSV ledger records the retained source index, post-burn-in index, tree
name, sampled state, and rooted flag for each kept tree.

Use `adapter beast-consensus` when you want a governed majority-rule summary
tree from BEAST posterior samples instead of chaining generic tree-set commands
by hand.

```bash
bijux-phylogenetics adapter beast-consensus \
  artifacts/multilocus-beast.1.trees \
  --burnin-fraction 0.1 \
  --out artifacts/multilocus-beast.consensus.nwk \
  --tree-set-out artifacts/multilocus-beast.postburnin.nwk \
  --clade-table-out artifacts/multilocus-beast.clades.tsv \
  --json
```

That command applies burn-in filtering once, writes the retained posterior tree
set, computes informative clade frequencies, and writes a consensus tree whose
internal labels are posterior clade probabilities on the `0..1` scale. The
clade-frequency ledger preserves every informative retained clade, not only the
majority clades that appear in the consensus topology, so reviewers can see
which alternative groupings remained credible after burn-in removal.

Use `adapter beast-diversity` when you want a governed uncertainty view over
posterior topology dispersion rather than only a single consensus summary.

```bash
bijux-phylogenetics adapter beast-diversity \
  artifacts/multilocus-beast.1.trees \
  --burnin-fraction 0.1 \
  --tree-set-out artifacts/multilocus-beast.postburnin.nwk \
  --distance-out artifacts/multilocus-beast.distances.tsv \
  --topology-out artifacts/multilocus-beast.topologies.tsv \
  --unstable-clade-out artifacts/multilocus-beast.unstable-clades.tsv \
  --json
```

That workflow keeps the retained posterior tree set, writes the full pairwise
RF distance ledger, clusters retained trees by rooted topology, and exports an
unstable-clade ledger for non-unanimous groupings. The JSON summary then adds
the retained tree count, number of unique rooted topologies, dominant topology
frequency, effective topology count, RF pair count, and unstable-clade count
so reviewers can judge whether the posterior is tightly concentrated or broadly
dispersed before treating one consensus tree as sufficient.

Use `adapter mrbayes-run` after preparation when you want the governed runtime
to execute MrBayes and preserve the native posterior artifacts for later
inspection instead of leaving the engine outputs implicit.

```bash
bijux-phylogenetics adapter mrbayes-run \
  artifacts/multilocus-bayesian.nex \
  --timeout-seconds 1800 \
  --resume \
  --incomplete-run-policy clean \
  --json
```

The run keeps the sampled posterior trees (`.run1.t`), parameter traces
(`.run1.p`), MCMC diagnostics (`.mcmc`), and consensus tree (`.con.tre`) in
the same output directory. Those files are then consumable through the parser
surfaces instead of through manual text scraping:

- `adapter mrbayes-traces` for tabular parameter traces from `.run1.p`
- `adapter mrbayes-parameters` for burn-in-aware posterior mean, median, SD, 95% HPD, and ESS summaries from `.run1.p`
- `adapter mrbayes-trees` for sampled posterior trees and generation tags from `.run1.t`
- `adapter mrbayes-subsample` for governed retained posterior subsets with generation metadata from `.run1.t`
- `adapter mrbayes-mcmc` for acceptance-rate and split-frequency diagnostics from `.mcmc`
- `adapter mrbayes-consensus` for the annotated consensus topology and posterior-probability range from `.con.tre`

That separation matters because the consensus-tree file uses MrBayes-specific
annotation syntax that common generic NEXUS readers may not accept directly.
The governed parser strips the native inline annotations only after recording
their posterior-probability values, so reviewers keep both a clean tree object
and the support summary that MrBayes actually emitted.

When you want posterior parameter summaries directly, use the dedicated
diagnostics surface:

```bash
bijux-phylogenetics adapter mrbayes-parameters \
  artifacts/multilocus-bayesian.nex.run1.p \
  --burnin-fraction 0.25 \
  --summary-out artifacts/multilocus-bayesian.parameter-summary.tsv \
  --json
```

That workflow writes one row per retained parameter with posterior mean,
median, sample standard deviation, 95% HPD interval, effective sample size,
first-half versus second-half mean split, and the retained generation window.

Use `adapter mrbayes-burnin-sensitivity` when you want the same cross-fraction
review for MrBayes posterior trees and traces.

```bash
bijux-phylogenetics adapter mrbayes-burnin-sensitivity \
  artifacts/multilocus-bayesian.nex.run1.t \
  --traces artifacts/multilocus-bayesian.nex.run1.p \
  --slice-out artifacts/multilocus-bayesian.burnin-slices.tsv \
  --parameter-out artifacts/multilocus-bayesian.burnin-parameters.tsv \
  --clade-out artifacts/multilocus-bayesian.burnin-clades.tsv \
  --json
```

This workflow tests the same governed default fractions `5%`, `10%`, `25%`,
and `50%` unless you pass an explicit custom set. It reports parameter
instability from non-overlapping 95% HPD intervals and clade instability from
posterior probabilities that move across the majority-rule threshold.

Use `adapter mrbayes-subsample` when a downstream audit only needs a retained
posterior subset instead of the full `.run1.t` file.

```bash
bijux-phylogenetics adapter mrbayes-subsample \
  artifacts/multilocus-bayesian.nex.run1.t \
  --method random \
  --burnin-fraction 0.25 \
  --sample-count 200 \
  --seed 7 \
  --tree-set-out artifacts/multilocus-bayesian.subsample.nwk \
  --sample-table-out artifacts/multilocus-bayesian.subsample.tsv \
  --json
```

That workflow supports both evenly spaced thinning and seeded random
subsampling. The retained tree-set file is normalized Newick, and the sample
ledger keeps the retained source index, post-burn-in index, tree name,
sampled generation, and rooted flag so the subset remains traceable back to the
native MrBayes posterior output.

For ultrafast bootstrap review specifically, `adapter bootstrap` now writes
three reviewer-facing TSV artifacts alongside the native IQ-TREE files:

- `PREFIX.support.tsv` for every supported internal branch
- `PREFIX.low-support.tsv` for the subset of branches below the governed weak-support threshold
- `PREFIX.support-histogram.tsv` for the support distribution buckets used in reports

Those artifacts map support values back onto explicit descendant-taxon clades,
flag low-support branches directly, and preserve the same `lt50`, `50to69`,
`70to89`, and `ge90` buckets that appear in the manifest and HTML workflow
report.

That bootstrap surface is now a real validation gate as well as an artifact
writer. If IQ-TREE returns a treefile and `.ufboot` set but the supported tree
does not actually contain parseable support labels, the workflow fails instead
of writing a superficially complete review bundle.

Use `adapter infer-fast` when you need a rapid approximate tree for one aligned
DNA or protein matrix and you want reviewer-facing local-support evidence
instead of a bare Newick file.

```bash
bijux-phylogenetics adapter infer-fast \
  aligned-matrix.fasta \
  --out artifacts/mammals.fasttree.nwk \
  --sequence-type protein \
  --json
```

That workflow runs FastTree directly, writes the inferred tree to the requested
output path, and emits three sidecar TSV artifacts next to it:

- `mammals.fasttree.support.tsv` for every internal branch with a parsable local-support label
- `mammals.fasttree.low-support.tsv` for the subset of branches below the governed weak-support threshold
- `mammals.fasttree.support-histogram.tsv` for the local-support distribution buckets used in reports

The structured manifest and JSON also expose the FastTree approximation
contract explicitly: the method is approximately maximum-likelihood, the native
support labels are SH-like local-support proportions, and the governed support
scale is `0..1` rather than bootstrap percentages.

Use `adapter infer-large` when the matrix is already aligned and you need a
large-alignment inference path that avoids copying the matrix through multiple
Python-side structures before FastTree runs.

```bash
bijux-phylogenetics adapter infer-large \
  aligned-matrix.fasta \
  --out-dir artifacts/large-alignment \
  --prefix mammals \
  --sequence-type protein \
  --timeout-seconds 600 \
  --resume \
  --incomplete-run-policy clean \
  --json
```

That workflow performs a streamed preflight scan of the aligned FASTA, runs
FastTree in place on the original matrix, and writes these review outputs:

- `mammals.tree`
- `mammals.support.tsv`
- `mammals.low-support.tsv`
- `mammals.support-histogram.tsv`
- `mammals.resources.tsv`
- `mammals.log`
- `mammals.manifest.json`

The streamed preflight records sequence count, alignment width, total site
cells, and inferred sequence type without materializing the full matrix as an
in-memory Python alignment object. The resource ledger separates preflight
allocation observations from sampled FastTree process RSS so large runs expose
both wall time and memory pressure directly. When `--resume` is used, the
workflow reuses a completed manifest only if the input checksum, command, and
recorded outputs still agree; otherwise it reruns the inference step.

Use `adapter compare-engines` when you need a governed side-by-side comparison
between IQ-TREE and FastTree on the same aligned matrix.

```bash
bijux-phylogenetics adapter compare-engines \
  aligned-matrix.fasta \
  --out-dir artifacts/engine-comparison \
  --prefix mammals \
  --sequence-type dna \
  --bootstrap-replicates 1000 \
  --timeout-seconds 1800 \
  --resume \
  --incomplete-run-policy clean \
  --json
```

That workflow runs IQ-TREE model selection, IQ-TREE ultrafast bootstrap
support inference, and FastTree approximate inference on the same alignment.
It then writes these user-facing outputs:

- `mammals.fasttree.nwk`
- `mammals.iqtree-support.nwk`
- `mammals.comparison.html`
- `mammals.comparison.tsv`
- `mammals.shared-clades.tsv`
- `mammals.conflicting-clades.tsv`
- `mammals.support-weighted-conflicts.tsv`
- `mammals.conclusions.tsv`
- `mammals.stability-summary.tsv`
- `mammals.taxon-influence.tsv`
- `mammals.manifest.json`

Use `phylo preflight` before any governed external-engine workflow when the
environment itself is uncertain.

```bash
bijux-phylogenetics phylo preflight \
  --workflow fasta-to-tree \
  --json
```

That command inspects MAFFT, trimAl, IQ-TREE, FastTree, MrBayes, and BEAST
before the workflow starts. The report includes the resolved executable path,
detected version text, engine support class, and a workflow-readiness table so
you can see exactly which workflows are runnable in the current environment.

If you install an engine outside `PATH`, pass its explicit
`--mafft-executable`, `--trimal-executable`, `--iqtree-executable`,
`--fasttree-executable`, `--mrbayes-executable`, or `--beast-executable`
override here first. When the selected workflow is blocked, the command exits
early instead of letting a longer alignment or Bayesian run fail partway
through.

Every governed external-engine workflow in this section writes a
`*.manifest.json` file with input checksums, structured workflow config,
resolved engine commands, detected engine versions, seed values, runtime, and
output checksums. Treat that manifest as the durable provenance record for
review, reruns, and governed downstream reports.

Use `phylo run` when you want one governed workflow to start from a single
config file instead of a long command line.

```bash
bijux-phylogenetics phylo run workflow-config.yaml --json
```

The current one-config surface targets the canonical `fasta-to-tree` workflow.
One config file can declare:

- the input FASTA
- optional metadata and traits tables
- engine executable choices
- alignment and trimming settings
- inference seed, threads, and bootstrap replicates
- output directory and optional bundle directory
- timeout and incomplete-run controls

Invalid config files fail before engine preflight or alignment starts. A valid
run executes the same governed workflow as `adapter fasta-to-tree`, then
exports and validates one complete result bundle automatically. That bundle now
includes the resolved workflow config plus copied config-source, metadata, and
traits files when they were supplied.

When one of these governed workflows fails and `--json` is enabled, the error
surface now explains the biological blocker directly. The structured payload
includes `failure_reason`, `scientific_explanation`, `likely_causes`,
`actionable_fixes`, and `evidence`, so the same run can say whether the input
FASTA contains duplicate or empty records, whether trimming removed every
retained site, whether inference never produced a tree or produced an
unparsable tree artifact, whether the tree and trait table disagree on taxon
membership, or whether a BEAST or MrBayes artifact is missing its expected
file, header, sampled row, or posterior-tree block.

Use `phylo replay` when you need to rerun one of those manifests and verify the
result against the original workflow output.

```bash
bijux-phylogenetics phylo replay \
  artifacts/fasta-to-tree/mammals.manifest.json \
  --out-dir artifacts/fasta-to-tree-replay \
  --json
```

Replay refuses to continue when the recorded input checksums no longer match
the current inputs. It also reports engine-version drift when the current
executable version differs from the one captured in the manifest, while still
checking whether the replayed outputs remained scientifically equivalent.
Equivalence is workflow-specific rather than byte-for-byte: trees are compared
by topology and support, alignments by aligned records, and model-selection
workflows by the selected substitution model.

Use `phylo bundle` when you need one portable workflow handoff directory for a
reviewer or collaborator.

```bash
bijux-phylogenetics phylo bundle \
  artifacts/fasta-to-tree/mammals.manifest.json \
  --out-dir artifacts/fasta-to-tree-bundle \
  --json
bijux-phylogenetics phylo validate-bundle \
  artifacts/fasta-to-tree-bundle \
  --json
```

The exported bundle keeps the copied workflow manifest, extracted config,
bundle-local rerun ledger, reviewer-facing HTML report, final workflow outputs,
and declared step-level engine artifacts in one directory. Validation is not
just a checksum pass: `phylo validate-bundle` also checks the required workflow
entries for completeness, so a missing report, tree, or step manifest fails the
handoff before it is shared.

The shared-clade ledger preserves both engines' support values for clades that
appear in both trees. The conflicting-clade ledger separates clades that appear
in only one tree from shared clades whose normalized support fractions diverge
enough to merit review. The support-weighted conflict ledger ranks disagreements
by normalized support burden so weak disagreements stay distinct from serious
high-support conflicts. The conclusion ledger classifies each reviewed clade as
stable, unstable, or engine-specific, and the stability summary reduces that
evidence to reviewer-facing counts plus the top conflict-driver taxa when they
can be detected by shared-taxon pruning.

The normalization rule is explicit and limited: FastTree SH-like local support
and IQ-TREE UFBoot support are shown together as fractions only for side-by-side
review, not as proof that the two support methods are interchangeable.

Use `adapter reproducibility` when you need to test whether repeated
bootstrap-supported IQ-TREE inference stays deterministic under fixed settings.

```bash
bijux-phylogenetics adapter reproducibility \
  aligned-matrix.fasta \
  --out-dir artifacts/inference-reproducibility \
  --prefix mammals \
  --sequence-type dna \
  --bootstrap-replicates 1000 \
  --repeats 3 \
  --json
```

That workflow runs model selection once to choose a fixed model, then reruns
the same supported IQ-TREE inference settings multiple times on the same
alignment. It writes these review outputs:

- `mammals.runs.tsv`
- `mammals.comparisons.tsv`
- `mammals.support-deltas.tsv`
- `mammals.manifest.json`

The run ledger records one manifest and supported tree per rerun. The
comparison ledger classifies each rerun relative to the baseline as
`deterministic`, `equivalent`, or `unstable` after checking topology,
log-likelihood, and clade support values. The support-delta ledger preserves
the per-clade support shifts as normalized fractions so support drift can be
reviewed directly instead of being inferred from summary text alone.

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
  --genetic-code 1 \
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

The codon-aware alignment workflow validates each coding sequence under one
explicit genetic code before MAFFT runs. It excludes frame-broken sequences,
sequences with ambiguous or invalid codons, and sequences with premature stop
codons. It then aligns a translated amino-acid guide and projects guide gaps
back as nucleotide triplets, so the final alignment length stays divisible by
three and codon boundaries are retained. The workflow also writes:

- the translated amino-acid guide input
- the aligned amino-acid guide
- a TSV exclusion ledger with invalid-codon and stop-codon counts
- a TSV codon summary ledger covering every input sequence

Use `alignment coding --genetic-code ...` to inspect one nucleotide dataset
under the same codon table before alignment, and use
`alignment translate --genetic-code ... --codon-validation-out <table.tsv>
--excluded-sequences-out <table.tsv>` when you need the amino-acid translation
itself as a reviewable artifact.

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
  --timeout-seconds 1800 \
  --resume \
  --incomplete-run-policy clean \
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
- `mammals.run.json`

It also retains step-specific engine artifacts under
`artifacts/fasta-to-tree/engine-artifacts/mammals/` for auditability.

The JSON payload for this workflow reports a `supported` method tier with its
real-engine validation basis. That tier is intended to separate this validated
workflow from approximate exploratory commands such as `comparative logistic`,
which now report `experimental` and emit a warning naming their approximation,
and from Bayesian report commands such as `adapter mrbayes-report`, which are
classified as `parser-only`.

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

`--resume` is stage-aware for this composite workflow. Raw-input validation,
alignment, trimming, model selection, inference, support, and final reporting
each contribute one deterministic fingerprint to `mammals.manifest.json`. A
downstream stage is reused only when its recorded inputs, config, command, and
detected engine version still match; changing the raw FASTA invalidates the
downstream tree-building stages, while changing bootstrap replicates invalidates
the support and report stages without forcing a fresh alignment.

When you need to send this workflow to another scientist, export from that same
`mammals.manifest.json` with `phylo bundle`. The resulting bundle includes the
copied input FASTA when it is still available, the workflow config, the
reviewer-facing `.aln`, `.trimmed.aln`, `.tree`, `.log`, `.model.tsv`,
`.support.tsv`, `.run.json` outputs, the copied step manifests, the native
engine artifacts, and one bundle-local rerun ledger plus HTML summary report.

Use `alignment sequence-type` when you need the raw FASTA type decision before
any engine is invoked. It reports compatible types, the selected default, the
confidence level, and mixed or invalid blocking signals.

Use `alignment validate-input` when you need one broader report of duplicate
identifiers, illegal sequence characters, empty records, raw-sequence length
outliers, and sequence-type detection before any engine is invoked.

That raw validation path now streams FASTA inputs linearly, so it remains
useful on thousand-sequence inputs before alignment. The downstream
`alignment quality` surface also reuses one loaded alignment instead of
re-reading the same matrix for each sub-diagnostic, and it emits an explicit
warning if near-duplicate pairwise review is skipped above the governed
large-alignment threshold.

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

That governed corpus now covers at least three real DNA inputs and three real
protein inputs. The validation surface compares alignment, tree, support,
model, manifest, and run artifacts semantically instead of byte-for-byte so
harmless path and timestamp variation does not masquerade as scientific drift.

Those checks pin reviewer-facing output bundles for:

- `gnathostome-ortholog-proteins`
- `gnathostome-ortholog-coding-sequences`
- `strnog-enog411bqtj-proteins`

They verify that the workflow emits the aligned matrix, trimmed matrix, tree,
log, model table, and support table with stable names on real DNA and protein
inputs.

When you want one public packaged protein benchmark rather than the broader
test-only regression corpus, use
`demo gnathostome-ortholog-protein-benchmark`. It materializes one shipped
amino-acid FASTA panel plus governed expected outputs, reruns MAFFT, trimAl,
and IQ-TREE on those packaged protein inputs, and writes one explicit
`molecular-assumptions.tsv` ledger so reviewers can see that the workflow runs
with `-st AA`, searches protein models only, does not translate coding DNA,
and does not depend on nucleotide-specific assumptions such as codon position
or GC interpretation.

For coding DNA, prefer the codon-aware alignment workflow above and then run
the downstream adapter steps explicitly. The current `adapter fasta-to-tree`
path is still a generic trimAl-based pipeline rather than the codon-preserving
entrypoint.
