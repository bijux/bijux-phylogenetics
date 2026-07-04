# Test Fixtures

This directory is the governed fixture corpus for the runtime package. It holds
the checked inputs, reference outputs, and catalog metadata that keep tests
reviewable, reproducible, and stable across refactors.

## Fixture Surfaces

- `trees/` holds single-tree and tree-set inputs used for parsing, validation,
  topology, and reporting coverage.
- `alignments/` holds aligned sequence fixtures for alignment, composition,
  translation, and distance review.
- `concatenation/` holds per-locus alignment inputs used to govern
  supermatrix-assembly behavior.
- `metadata/` holds trait tables, distance matrices, simulation catalogs, and
  other structured review inputs.
- `fasta_to_tree/real/` holds biologically real raw FASTA inputs for the
  end-to-end tree-inference workflow.
- `expected/` holds checked deterministic outputs and reviewer-facing golden
  bundles.
- `engine_outputs/` holds parser-stress artifacts for external-engine outputs.

## Why Catalogs Exist

Most governed fixture families are indexed through shared catalogs under
`metadata/`. Those catalogs assign durable fixture ids, define expected
behavior, and let validation and parity lanes resolve the same reviewed inputs
without hardcoded one-off paths.

The main catalogs cover:

- tree fixtures and tree-set fixtures for parsing, rooting, and topology review
- DNA alignments for composition, distance, translation, and invalid-input
  boundaries
- trait tables for comparative, ancestral, and contrast workflows
- distance matrices for neighbor-joining validation
- rooted random-tree and coalescent simulation envelopes
- shared `phytools` and `geiger` comparative parity panels
- the checked BEAST posterior corpus used when a live `beast` executable is not
  present

## Reviewer-Facing Goldens

The `expected/` surface holds more than byte snapshots. It records the stable
reader-facing artifacts promised by product workflows, including checked
workflow bundles, report assets, and benchmark-reference outputs. These files
are the evidence surface used to detect regressions in owned outputs.

## Usage Rule

Tests should resolve fixture paths through the local `fixture(...)` helper or a
catalog-backed loader rather than assuming a permanently flat directory
layout.
