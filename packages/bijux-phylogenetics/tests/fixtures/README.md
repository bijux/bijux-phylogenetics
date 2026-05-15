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

The tree fixture corpus now also has one governed shared catalog at
`metadata/shared_tree_fixture_catalog.json`. That catalog assigns durable
fixture ids, records parse and validation expectations, and marks the
structural features that matter for cross-tool parity, including balanced,
pectinate, star, polytomy, rooted, unrooted, ultrametric, non-ultrametric,
zero-branch, long-branch, internal-label, branch-support, quoted-label, and
malformed-Newick cases. The live `ape` parity harness resolves its tree inputs
through those fixture ids instead of hardcoding ad hoc file paths. That same
catalog now also owns the governed outgroup-rooting fixtures used to compare
Bijux and `ape::root` on one-tip outgroups, monophyletic multi-tip outgroups,
already-rooted trees, missing outgroups, and non-monophyletic outgroup
failures.

The tree fixture corpus now also has one governed shared tree-set catalog at
`metadata/shared_tree_set_fixture_catalog.json`. That catalog assigns durable
fixture ids for multiple-tree Newick inputs, records expected tree counts and
shared taxon sets, and gives the live `ape::read.tree` parity lane one owned
surface for multiple-tree structure checks instead of one-off fixture paths.

The DNA alignment fixture corpus now also has one governed shared catalog at
`metadata/shared_dna_alignment_fixture_catalog.json`. That catalog assigns
durable fixture ids, records load and translation expectations, and marks the
distance, composition, missingness, lowercase, identical-sequence,
high-divergence, unequal-length, valid-coding, frame-error, internal-stop, and
terminal-stop cases that matter for cross-tool DNA parity. The live `ape`
parity harness resolves its DNA inputs through those fixture ids for
`base.freq`, `dist.dna`, and `trans` instead of relying on loose path lists or
one-off inline sequences. The unequal-length fixture now also serves as a
governed `ape::dist.dna` failure case, so ragged DNA input is checked as an
explicit parity boundary rather than only as local validation diagnostics.
The same governed DNA catalog now also anchors both raw and JC69
`ape::dist.dna` parity, including the high-divergence fixture where JC69 must
distinguish undefined pairs from pairs that diverge toward infinity.

The metadata fixture corpus now also has one governed shared catalog at
`metadata/shared_trait_table_fixture_catalog.json`. That catalog assigns
durable fixture ids for continuous traits, binary and multistate discrete
traits, missing-value cases, extra-versus-missing taxon mismatches, duplicate
taxon negatives, constant-trait negatives, categorical predictors, and
misordered taxon rows. The owned discrete-reference validator and the live
`ape::ace` ancestral review tests now resolve their small governed trait
tables through those fixture ids instead of hardcoding one-off table paths.

The `expected/` directory now also carries benchmark-corpus regression
snapshots used to pin stable dataset summaries across releases.

The checked-in fixture surfaces are additionally grouped into clean, broken,
and messy benchmark corpora through `bijux_phylogenetics.validation_corpus`
so validation, warning-rich behavior, and failure signatures can be audited
from library code as well as tests.

The `metadata/beast2_strict_yule_posterior.xml`, `.log`, and `.trees` files
are the governed real-artifact BEAST corpus used by the Bayesian validation
matrix whenever a live `beast` executable is not present locally. They are not
toy parser strings; they are the checked evidence surface for real BEAST output
structure.

Tests should resolve files through the local `fixture(...)` helper rather than
assuming a flat directory layout.
