---
title: Repository Overview
audience: public
type: handbook
status: active
owner: bijux-phylogenetics
last_reviewed: 2026-04-29
---

# Repository Overview

`bijux-phylogenetics` provides a governed Python surface for tree validation,
inspection, comparison, metadata linkage, alignment trimming, coding-sequence
translation, explicit rooting transforms, comparative trait analysis,
ancestral-state reconstruction, external engine orchestration, deterministic
tree rendering, dataset crosswalk and completeness auditing, publication
figure packaging, evidence bundles, and HTML report generation.

The repository intentionally does not reimplement inference engines. Its
current product surface is the reproducible orchestration and evidence layer
around trees, alignments, and trait tables.

## Current CLI Surface

- `bijux-phylogenetics validate tree.nwk`
- `bijux-phylogenetics inspect tree.nwk`
- `bijux-phylogenetics compare left.nwk right.nwk`
- `bijux-phylogenetics annotate tree.nwk --metadata traits.tsv`
- `bijux-phylogenetics render tree.nwk --layout phylogram --support-labels --metadata samples.tsv --label-column species --metadata-strip-columns location --traits traits.tsv --categorical-column habitat --continuous-column height_cm --heatmap-columns height_cm,status --package-dir artifacts/tree-figure --out artifacts/tree.svg`
- `bijux-phylogenetics bundle run/ --out evidence-pack/`
- `bijux-phylogenetics report --tree tree.nwk --alignment alignment.fasta --traits traits.tsv --metadata samples.tsv --out phylo-report.html`
- `bijux-phylogenetics alignment trim alignment.fasta --out trimmed.fasta --sequence-missingness-threshold 0.4`
- `bijux-phylogenetics alignment coding coding.fasta --json`
- `bijux-phylogenetics alignment translate coding.fasta --out translated.fasta`
- `bijux-phylogenetics report dataset --tree tree.nwk --metadata samples.tsv --traits traits.tsv --alignment alignment.fasta --tip-dates tip-dates.tsv --calibrations calibrations.tsv --out artifacts/dataset-report.html --json`
- `bijux-phylogenetics comparative readiness tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative contrasts tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative signal tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative brownian tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative compare-models tree.nwk traits.tsv --trait height_cm --json`
- `bijux-phylogenetics comparative pgls tree.nwk traits.tsv --response height_cm --predictors body_mass log_range --json`
- `bijux-phylogenetics comparative pgls tree.nwk traits.tsv --formula "height_cm ~ body_mass * habitat" --json`
- `bijux-phylogenetics comparative multiple-testing tree.nwk traits.tsv --responses height_cm range_km --predictors body_mass log_range --json`
- `bijux-phylogenetics comparative report tree.nwk traits.tsv --formula "height_cm ~ body_mass + habitat" --out artifacts/comparative-report.html --json`
- `bijux-phylogenetics ancestral continuous tree.nwk traits.tsv --trait height_cm --model brownian --json`
- `bijux-phylogenetics ancestral discrete tree.nwk traits.tsv --trait habitat --model symmetric --state-ordering ordered --ordered-states low,medium,high --json`
- `bijux-phylogenetics ancestral compare tree.nwk traits.tsv --trait height_cm --left-model brownian --right-model ou --json`
- `bijux-phylogenetics ancestral sensitivity tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --json`
- `bijux-phylogenetics ancestral report tree.nwk traits.tsv --trait height_cm --kind continuous --compare-model ou --compare-tree tree-alt.nwk --out artifacts/ancestral-report.html`
- `bijux-phylogenetics ancestral package tree.nwk traits.tsv --trait habitat --kind discrete --model symmetric --state-ordering ordered --ordered-states low,medium,high --out-dir artifacts/ancestral-package --json`
- `bijux-phylogenetics adapter align unaligned.fasta --out aligned.fasta --json`
- `bijux-phylogenetics adapter model-select alignment.fasta --out-dir artifacts/model-select --prefix mammals --json`
- `bijux-phylogenetics adapter infer-ml alignment.fasta --out-dir artifacts/ml --model GTR+G --prefix mammals --json`
- `bijux-phylogenetics adapter bootstrap alignment.fasta --out-dir artifacts/bootstrap --model GTR+G --replicates 1000 --prefix mammals --json`
- `bijux-phylogenetics adapter consensus artifacts/bootstrap/mammals.ufboot --out-dir artifacts/consensus --prefix mammals --json`
- `bijux-phylogenetics adapter infer-fast alignment.fasta --out artifacts/fasttree.nwk --json`
- `bijux-phylogenetics adapter compare --fast-tree artifacts/fasttree.nwk --ml-tree artifacts/ml/mammals.treefile --out artifacts/engine-comparison.html --json`
- `bijux-phylogenetics adapter mrbayes-prepare alignment.fasta --out artifacts/mrbayes/analysis.nex --ngen 20000 --samplefreq 100 --json`
- `bijux-phylogenetics adapter mrbayes-run artifacts/mrbayes/analysis.nex --resume --json`
- `bijux-phylogenetics adapter mrbayes-summarize artifacts/mrbayes/analysis.run1.t --burnin-fraction 0.25 --json`
- `bijux-phylogenetics adapter mrbayes-traces artifacts/mrbayes/analysis.run1.p --json`
- `bijux-phylogenetics adapter mrbayes-ess artifacts/mrbayes/analysis.run1.p --json`
- `bijux-phylogenetics adapter mrbayes-convergence artifacts/mrbayes/analysis.run1.p --ess-threshold 200 --json`
- `bijux-phylogenetics adapter mrbayes-report artifacts/mrbayes/analysis.run1.t --traces artifacts/mrbayes/analysis.run1.p --out artifacts/mrbayes/posterior-report.html --json`
- `bijux-phylogenetics adapter beast-prepare alignment.fasta --out artifacts/beast/analysis.xml --tree tree.nwk --calibrations calibrations.tsv --tip-dates tip-dates.tsv --clock-model relaxed-lognormal --tree-prior birth-death --json`
- `bijux-phylogenetics adapter beast-calibrations tree.nwk calibrations.tsv --json`
- `bijux-phylogenetics adapter beast-tip-dates tree.nwk tip-dates.tsv --alignment alignment.fasta --json`
- `bijux-phylogenetics adapter beast-log artifacts/beast/run.log --json`
- `bijux-phylogenetics adapter beast-convergence artifacts/beast/run.log --ess-threshold 200 --json`
- `bijux-phylogenetics adapter beast-calibration-report tree.nwk calibrations.tsv --tip-dates tip-dates.tsv --alignment alignment.fasta --out artifacts/beast/calibration-audit.html --json`
- `bijux-phylogenetics adapter bayesian-evidence --out-dir artifacts/bayesian-bundle --inputs alignment.fasta calibrations.tsv tip-dates.tsv --configs artifacts/beast/analysis.xml --trees tree.nwk --logs artifacts/beast/run.log --diagnostics diagnostics.json --reports artifacts/beast/calibration-audit.html --json`
- `bijux-phylogenetics adapter report artifacts/mrbayes/analysis.manifest.json --out artifacts/mrbayes/inference-report.html --json`
- `bijux-phylogenetics alignment identity-matrix alignment.fasta --out identity.tsv`
- `bijux-phylogenetics topology root-outgroup tree.nwk --taxa OutgroupA OutgroupB --out rooted.nwk`
- `bijux-phylogenetics topology reroot-midpoint tree.nwk --out midpoint-rooted.nwk`
- `bijux-phylogenetics inspect tree-with-support.nwk --json`

## Tree Diagnostics Highlights

- internal-node child counts are reported explicitly for every internal node
- missing branch lengths are separated into internal-branch and terminal-branch diagnostics
- singleton internal nodes are detected instead of being silently treated as ordinary branching structure
- long and short nonzero branch outliers are reported as concrete affected nodes
- numeric internal labels are classified as likely support values, textual internal labels as likely clade names
- support-like labels are checked for suspicious ranges and mixed probability-versus-percentage scales

## Tree Figure Highlights

- rectangular cladograms, branch-length phylograms, and circular tree layouts are supported from one render surface
- phylogram output includes a branch-length scale bar instead of leaving absolute lengths visually ambiguous
- internal support values can be rendered directly on branch junctions for review-ready trees
- tip annotations can combine categorical markers, continuous bars, metadata strips, and aligned heatmap columns
- named clades can be collapsed into figure summaries without mutating the input tree
- figure bundles can be emitted as a reusable package containing `figure.svg`, a JSON manifest, a caption draft, and aligned tip-annotation tables

## Comparative Analysis Highlights

- rooted trees and numeric traits can be checked explicitly for comparative readiness before modeling
- phylogenetic independent contrasts are available as one deterministic internal-node table per trait
- trait signal can be summarized with Blomberg's K, Pagel's lambda, and a permutation-based significance surface
- standalone Brownian and OU fits now report likelihoods, confidence intervals, identifiability warnings, residual structure, and leave-one-taxon-out sensitivity
- phylogenetic generalized least-squares now supports formula-audited categorical and interaction terms with explicit encoded-column tracking
- comparative workflows can emit multiple-testing summaries, integrated audit tables, influence reports, alternative-tree comparisons, pruning comparisons, and reviewer-facing limitations text

## Dataset Audit Highlights

- dataset reports now include explicit taxon crosswalk rows across tree, alignment, metadata, traits, dates, geography, and calibration targets
- completeness matrices expose which taxa are present on each surface instead of leaving reviewers to infer omissions from separate tables
- exclusion tables now attach exact causes, affected analysis families, pruning step counts, and ordering-drift warnings to dataset review artifacts
- named readiness levels summarize whether the dataset is inspection-ready, inference-ready, comparative-ready, time-tree-ready, or publication-ready

## Ancestral-State Highlights

- continuous ancestral-state reconstruction supports Brownian and OU-style trait models over a rooted pruned analysis tree
- discrete ancestral-state reconstruction supports Fitch parsimony plus likelihood-style ER, SYM, and ARD models with explicit ordered versus unordered state assumptions
- uncertainty is surfaced directly through continuous confidence intervals, low-confidence node warnings, discrete state-probability tables, and reviewer-facing downstream-risk summaries
- ancestral workflows can now emit sensitivity summaries, supplement-style reports, and publication-ready figure bundles in addition to annotated trees and deterministic node tables

## External Engine Highlights

- MAFFT-, trimAl-, IQ-TREE-, and FastTree-style workflows are exposed through one governed adapter surface instead of ad hoc shell snippets
- every engine workflow records the resolved executable, exact command vector, version output, stdout, stderr, and extracted warning lines in a deterministic manifest
- model selection emits a stable selected-model artifact, ML inference emits validated Newick trees, and bootstrap workflows retain both support trees and bootstrap tree sets
- fast approximate and ML trees can be compared through the same topology, clade, support, and branch-length report surface already used for native tree comparison

## Bayesian Workflow Highlights

- aligned FASTA inputs can now be turned into deterministic MrBayes NEXUS analyses with explicit model and MCMC settings
- MrBayes posterior runs emit posterior trees, parameter traces, resumable manifests, and HTML inference workflow reports through the same adapter layer as other external engines
- posterior tree sets can be burn-in filtered and summarized into consensus trees, rooted-topology counts, and shared-taxon summaries
- parameter traces can be parsed directly and converted into per-parameter effective sample size summaries without leaving the native Python/runtime surface
- posterior trace diagnostics now flag low ESS and unstable mean drift explicitly, and posterior HTML reports package those warnings beside consensus and clade-support summaries
- BEAST-style XML preparation now validates fossil calibrations, tip dates, and impossible age constraints before emitting a time-tree configuration
- BEAST logs can be parsed and checked for ESS and drift problems through the same deterministic trace-diagnostics surface used for MrBayes
- calibration audit reports and full Bayesian evidence packages now bundle dated-tree assumptions, diagnostics, logs, configs, and rendered reviewer-facing artifacts together
