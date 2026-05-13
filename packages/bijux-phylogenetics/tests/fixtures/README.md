# Test Fixtures

This fixture corpus is organized by durable surface:

- `trees/` for tree-format inputs
- `alignments/` for aligned FASTA inputs
- `concatenation/` for aligned per-locus FASTA inputs used to assemble checked supermatrices
- `metadata/` for metadata and trait tables
- `fasta_to_tree/real/` for biologically real raw FASTA inputs used to exercise the canonical end-to-end workflow
- `expected/` for checked-in expected outputs and goldens

The runtime now also uses these fixtures as a Level 1 reference-validation
corpus for tree, taxonomy, alignment, dataset, figure, and report regression
checks.

The `expected/` directory now also carries benchmark-corpus regression
snapshots used to pin stable dataset summaries across releases.

The checked-in fixture surfaces are additionally grouped into clean, broken,
and messy benchmark corpora through `bijux_phylogenetics.validation_corpus`
so validation, warning-rich behavior, and failure signatures can be audited
from library code as well as tests.

Tests should resolve files through the local `fixture(...)` helper rather than
assuming a flat directory layout.
