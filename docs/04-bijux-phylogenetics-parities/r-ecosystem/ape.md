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

## How Correspondence Is Established

1. A checked case selects a Bijux function, APE function, and shared fixture.
2. The R runner and Bijux surface execute independently.
3. Taxon ordering, node identity, matrices, states, and support scales are
   normalized by case-specific rules.
4. Exact, tolerance-based, or structural observations are evaluated.
5. Summary rows, observation rows, and reproducible mismatch artifacts retain
   the verdict.

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

## Deliberate Non-Claims

- A matching case does not make Bijux a drop-in replacement for all of APE.
- Strict input rejection in Bijux may intentionally differ from permissive APE
  behavior and must be classified rather than hidden.
- Simulation correspondence covers declared distributions and invariants, not
  seed-for-seed identity across languages.
- Plotting APIs and rendered-figure equivalence are outside the live analytical
  registry unless a case explicitly says otherwise.
