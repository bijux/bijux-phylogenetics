---
title: APE Correspondence
audience: public
type: reference
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# APE Correspondence

APE is the principal live reference for tree structure, DNA utilities, and
several foundational comparative operations. The governed registry currently
contains 180 cases across 35 Bijux function shards. Each shard identifies one
owned Python surface and the APE function used as its reference.

## Governed Surface

| Territory | APE references | Compared behavior |
| --- | --- | --- |
| tree I/O | `read.tree`, `write.tree` | taxa, rootedness, branch lengths, clades, tree-set cardinality and round trips |
| topology operations | `root`, `unroot`, `drop.tip`, `keep.tip`, `extract.clade`, `getMRCA`, `is.monophyletic` | resulting clades, stable node correspondence, retained taxa and length policy |
| tree sets and support | `consensus`, `prop.clades` | strict/majority clades, frequencies, reference-tree support mapping |
| tree geometry | `cophenetic.phylo`, `dist.topo`, `vcv.phylo`, `node.depth.edgelength`, `branching.times`, `is.ultrametric`, `gammaStat` | matrices, topology distance, covariance, node depth, times and scalar diagnostics |
| DNA operations | `as.DNAbin`, `base.freq`, `seg.sites`, `dist.dna`, `trans` | state normalization, frequencies, variable sites, distance models and translation |
| foundational methods | `nj`, `pic`, `ace` | neighbor-joining clades, independent contrasts, continuous estimates and discrete probabilities |
| simulation | `rtree`, `rcoal` | distribution-summary envelopes and structural invariants, not identical random draws |

The `dist.dna` registry covers raw, JC69, K80, F81, F84, T92, TN93, pairwise
deletion, ambiguity, gap, and variance-related conditions where the two
interfaces expose comparable semantics. The `ace` registry separates
continuous and discrete reconstruction because their outputs and root-prior
boundaries differ.

## Representative Function Correspondence

| APE reference | Bijux registry surface | Compared identity |
| --- | --- | --- |
| `read.tree`, `write.tree` | `load_tree`, `write_newick`, tree-set load/write variants | taxa, rootedness, clades, lengths, tree count, round trip |
| `root`, `unroot` | `root_tree_on_outgroup`, `unroot_tree` | outgroup, resulting rootedness, clades, length policy |
| `drop.tip`, `keep.tip` | `drop_tree_taxa`, `prune_tree_to_requested_taxa` | requested, retained and excluded taxa, resulting topology |
| `getMRCA`, `extract.clade`, `is.monophyletic` | `find_tree_mrca`, `extract_tree_clade_by_node_id`, `assess_tree_monophyly` | descendant-taxon identity rather than incidental node number |
| `cophenetic.phylo`, `dist.topo` | `compute_tree_tip_distance_matrix`, `compare_topology_distance` | ordered taxa, path lengths, rooted/unrooted distance convention |
| `vcv.phylo`, `pic` | `summarize_brownian_covariance`, `compute_phylogenetic_independent_contrasts` | covariance/contrast order, branch-length scale, numerical rule |
| `dist.dna`, `nj` | `compute_pairwise_genetic_distance_matrix`, `build_tree_from_imported_distance_matrix` | alphabet and deletion policy, pair counts, distance model, clades |
| `ace` | continuous and discrete ancestral reconstruction surfaces | model/state order, node keys, estimates or probabilities |
| `rtree`, `rcoal` | `simulate_random_trees`, `simulate_coalescent_trees` | replication population and structural/distribution summaries |

This table describes registered correspondence families, not every argument
accepted by the R functions. The case registry is authoritative for fixtures,
parameters, normalization, observations, and expected availability.

## How Correspondence Is Established

1. A checked case selects a Bijux function, APE function, and shared fixture.
2. The R runner and Bijux surface execute independently.
3. Taxon ordering, node identity, matrices, states, and support scales are
   normalized by case-specific rules.
4. Exact, tolerance-based, or structural observations are evaluated.
5. Summary rows, observation rows, and reproducible mismatch artifacts retain
   the case status and mismatch reason.

Newick strings are not compared as raw text when topology is the scientific
object. Random trees are not expected to be identical across random-number
generators. APE node numbers are translated through stable structural identity
where direct numbering would be brittle.

## Run The Live Lane

```bash
RUN_LIVE_PARITY=1 bijux-phylogenetics parity \
  --reference-source ape-live \
  --ape-shard compare_topology_distance \
  --summary-out artifacts/parity/ape-summary.tsv \
  --observations-out artifacts/parity/ape-observations.tsv \
  --ape-failure-root artifacts/parity/ape-failures \
  --json
```

Use `--ape-case` for one registered fixture or repeat `--ape-shard` to combine
owned surfaces. The live lane requires `Rscript`, APE, and
`RUN_LIVE_PARITY=1`. An unavailable reference produces a visible skip or
failure record; it does not fall back to a checked fixture while retaining the
`ape-live` label.

## Programmatic Surface

| Public locator | Purpose |
| --- | --- |
| `bijux_phylogenetics.parity.ape:list_ape_parity_cases` | enumerate the 180 eligible case contracts |
| `bijux_phylogenetics.parity.ape:list_ape_parity_shards` | enumerate 35 function-owned case groups |
| `bijux_phylogenetics.parity.ape:run_ape_parity_cases` | execute selected live cases and retain every terminal status |
| `bijux_phylogenetics.parity.ape:write_ape_parity_summary_table` | persist shard and case totals |
| `bijux_phylogenetics.parity.ape:write_ape_parity_observation_table` | persist field-, row-, topology-, and tolerance-level observations |

Use discovery output to construct the requested denominator before execution.
Case and shard filters that resolve to no eligible entries are selection
errors, not empty successful runs.

## Interpret The Case Status

An APE case passes only when the selected R and Bijux executions complete and
every registered observation rule holds. A failed case retains the specific
execution, field, row, topology, or tolerance mismatch. A skipped case records
that the live reference could not execute; it is excluded from neither the
selected denominator nor the report.

Quote the shard and case count with any result. “APE parity passed” is too broad
when the execution selected one of 35 shards or only a subset of its fixtures.
Case status does not create an Evidence Book verdict unless a governed bundle
later binds the observation to a scientific claim.

## Read One Observation Row

An observation row identifies the case and shard, APE and Bijux functions,
fixture, compared field or structural key, reference and observed values,
normalization, rule, tolerance where applicable, outcome, and mismatch reason.
Matrix and vector rows also retain label identity; topology rows retain the
taxon/rooting policy; stochastic rows retain the replication rule.

The summary table is derived from these rows and terminal case states. If a
summary total cannot be reconstructed from the selected cases and observations,
stop at that integrity failure rather than quoting the total.

## Accept A Structural Comparison

Before accepting an APE case, confirm that the observation rule matches the
scientific object. Taxon order may be normalized for a labelled matrix but not
discarded from its identity. Child order may be irrelevant to a clade but root
placement is not. Node numbers may be translated through descendant taxa but
must still resolve one-to-one. A rejected input can correspond with APE only
when the case explicitly classifies the behavioral difference.

The observation table, not the final Newick rendering, owns these decisions.

## Deliberate Non-Claims

- A matching case does not make Bijux a drop-in replacement for all of APE.
- Strict input rejection in Bijux may intentionally differ from permissive APE
  behavior and must be classified rather than hidden.
- Simulation correspondence covers declared distributions and invariants, not
  seed-for-seed identity across languages.
- Plotting APIs and rendered-figure equivalence are outside the live analytical
  registry unless a case explicitly says otherwise.
