---
title: Data And Structure
audience: public
type: index
status: canonical
owner: bijux-phylogenetics-docs
last_reviewed: 2026-07-21
---

# Data And Structure

Phylogenetic methods depend on a shared structural vocabulary: taxa, trees,
branch lengths, alignments, characters, partitions, and support annotations.
This domain defines how those objects enter the runtime and which invariants
must hold before a scientific method can consume them.

```mermaid
flowchart LR
    files["Newick · Nexus · FASTA<br/>tables · partitions"]
    parse["Parse without<br/>semantic guessing"]
    identity["Normalize and reconcile<br/>taxon identity"]
    validate["Topology · rootedness<br/>length · state checks"]
    objects["Typed trees, alignments<br/>and datasets"]

    files --> parse --> identity --> validate --> objects
```

Start with [Trees, Alignments, And Taxa](trees-alignments-and-taxa.md) for
rootedness, branch-length meaning, tree sets, alignment partitions, support
labels, and taxon reconciliation.

## Boundary Rules

- Parsing preserves source meaning; it does not silently repair scientific
  ambiguity.
- Taxon reconciliation emits a ledger of matches, exclusions, and conflicts.
- Rooting and branch-length assumptions belong to the consuming method.
- Tree-set identity and ordering remain stable across summaries.
- Support values retain their producing method and scale.
