# Primate Longevity Signal

This study is the governed PCM1 record for primate data preparation, tree
reconciliation, phylogenetic signal, and continuous ancestral reconstruction.
It establishes the prepared primate population reused by later studies and
keeps every population-changing operation visible before model output is
interpreted.

## Current Evidence State

| Record | Scientific surface | Verdict |
| --- | --- | --- |
| `evidence-001` | aggregate lecture reconstruction | `matched_with_tolerance` |
| `evidence-002` | workbook loading | `matched` |
| `evidence-003` | type repair | `matched` |
| `evidence-004` | missing-data accounting | `matched` |
| `evidence-005` | duplicate-species aggregation | `matched` |
| `evidence-006` | source and trimmed tree import | `matched` |
| `evidence-007` | rootedness, topology, and node-label diagnostics | `matched` |
| `evidence-008` | tree–trait correspondence | `matched` |
| `evidence-009` | governed processed exports | `matched` |

These values come from each bundle's `manifest.json`. The aggregate verdict is
not a claim that every lecture block is numerically comparable: its manifest
separately counts verified, tolerance-matched, plot-only, artifact-only,
seeded-input-only, and workflow-only blocks.

## Scientific Dependency Chain

```mermaid
flowchart LR
    source["raw workbook and tree"] --> repair["types, missingness,<br/>duplicate species"]
    repair --> population["reconciled 75-taxon<br/>analysis population"]
    population --> signal["Pagel lambda and<br/>zero-signal test"]
    population --> ancestral["continuous ancestral<br/>state reconstruction"]
    signal --> aggregate["evidence-001"]
    ancestral --> aggregate
```

A change to type coercion, row exclusion, duplicate aggregation, tree pruning,
or species mapping creates a different analytical population. Matching final
taxon counts do not establish that two populations are identical.

## Record Layout

| Path | Authority |
| --- | --- |
| `datasets/` | source copies and governed primate derivatives used by the bundles |
| `provenance/lund-course-sources.json` | origin and role of lecture sources, datasets, and repository adapters |
| `reference/` | lecture source plus executable R and Bijux reference paths |
| `evidence-001/` | aggregate claim, comparison records, block payloads, and review summaries |
| `evidence-002/` through `evidence-009/` | focused preparation and export contracts |

Within a bundle, begin with `manifest.json`, then read
`inputs.manifest.json`, `claims.json`, `checks.json`, and the machine-readable
files under `results/`. `report.md` supplies context after those records have
established identity and status.

## What The Aggregate Comparison Retains

The aggregate bundle records the prepared 75-taxon population, longevity
summary, estimated Pagel lambda, the lambda-zero likelihood-ratio comparison,
continuous ancestral estimates, branch-change summaries, and the R/Bijux
comparison rule for each scalar. The checked result inventory also retains
block-level classifications so rendering and workflow coverage cannot be
mistaken for numerical equivalence.

Use the descendant-taxon signature to identify an internal node. A lecture
node label is not durable enough on its own because parsing, pruning, or
rerooting can change internal numbering while preserving the clade.

## Claim Boundary

PCM1 supports the registered preparation, tree correspondence, processed
export, signal, and ancestral comparisons for this study family. It does not
establish:

- plot-for-plot equivalence with the R graphics stack;
- accuracy across unrelated trees, traits, or taxonomic populations;
- interchangeability of APE, GEIGER, phytools, and Bijux beyond the recorded
  observations;
- causal interpretation of phylogenetic signal or reconstructed change.

The aggregate bundle explicitly lists rendered-figure equivalence as out of
scope. Preserve that limitation whenever the result is summarized.

## Reviewer Route

1. Confirm source and derivative identities in `datasets/` and `provenance/`.
2. Establish the transformation chain with `evidence-002` through
   `evidence-009`.
3. Inspect the R and Bijux result JSON and the scalar parity table in
   `evidence-001/results/`.
4. Read the aggregate manifest verdict, block status counts, limitations, and
   freshness basis.
5. Use the reviewer summary and public
   [PCM1 guide](../../../docs/03-bijux-phylogenetics-evidence-book/pcm1-primate-longevity-signal/index.md)
   as interpretation aids, not substitutes for the governed records.
