---
title: Common Workflows
audience: public
type: how-to
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-12
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
  --json
```

The run keeps the sampled posterior trees (`.run1.t`), parameter traces
(`.run1.p`), MCMC diagnostics (`.mcmc`), and consensus tree (`.con.tre`) in
the same output directory. Those files are then consumable through the parser
surfaces instead of through manual text scraping:

- `adapter mrbayes-traces` for tabular parameter traces from `.run1.p`
- `adapter mrbayes-parameters` for burn-in-aware posterior mean, median, SD, 95% HPD, and ESS summaries from `.run1.p`
- `adapter mrbayes-trees` for sampled posterior trees and generation tags from `.run1.t`
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

For ultrafast bootstrap review specifically, `adapter bootstrap` now writes
three reviewer-facing TSV artifacts alongside the native IQ-TREE files:

- `PREFIX.support.tsv` for every supported internal branch
- `PREFIX.low-support.tsv` for the subset of branches below the governed weak-support threshold
- `PREFIX.support-histogram.tsv` for the support distribution buckets used in reports

Those artifacts map support values back onto explicit descendant-taxon clades,
flag low-support branches directly, and preserve the same `lt50`, `50to69`,
`70to89`, and `ge90` buckets that appear in the manifest and HTML workflow
report.

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
- `mammals.manifest.json`

The shared-clade ledger preserves both engines' support values for clades that
appear in both trees. The conflicting-clade ledger separates clades that appear
in only one tree from shared clades whose normalized support fractions diverge
enough to merit review. The normalization rule is explicit and limited:
FastTree SH-like local support and IQ-TREE UFBoot support are shown together as
fractions only for side-by-side review, not as proof that the two support
methods are interchangeable.

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
