# Test Fixtures

This fixture corpus is organized by durable surface:

- `trees/` for tree-format inputs
- `alignments/` for aligned FASTA inputs
- `metadata/` for metadata and trait tables
- `expected/` for checked-in expected outputs and goldens

The runtime now also uses these fixtures as a Level 1 reference-validation
corpus for tree, taxonomy, alignment, dataset, figure, and report regression
checks.

Tests should resolve files through the local `fixture(...)` helper rather than
assuming a flat directory layout.
