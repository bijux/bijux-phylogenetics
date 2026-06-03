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
The same governed DNA catalog now also anchors raw, JC69, K80, F81, and TN93
`ape::dist.dna` parity, including the high-divergence fixture where corrected
DNA distances must distinguish undefined pairs from pairs that diverge toward
infinity and the unequal-composition fixtures where F81 and TN93 must expose
their estimated base-frequency assumptions instead of silently falling back to
simpler models.
The catalog also includes one all-gap-or-missing alignment fixture so
`ape::base.freq` parity and the owned `alignment composition --base-frequency-out`
surface can prove their explicit edge-case policy without fabricating A/C/G/T
content from a biologically empty alignment.
That same governed DNA catalog now also includes invariant and one-variable-site
fixtures for `ape::seg.sites` parity and the owned
`alignment segregating-sites --site-table-out` surface, so segregating-site
review stays tied to durable fixture ids instead of ad hoc inline alignments.
The coding portion of the same DNA catalog now also includes ambiguous-codon
and alternate-genetic-code fixtures for `ape::trans` parity, plus the
frame-truncation fixture where aligned translation must drop trailing partial
codons with an explicit warning instead of silently mutating the review path.
The catalog now also includes one invalid-symbol alignment fixture for the
owned DNAbin-compatible nucleotide matrix boundary, so Bijux can reject
unsupported nucleotide states explicitly while the governed live
`ape::as.DNAbin` lane still proves basic structure parity on valid lowercase,
gap-bearing, and ambiguity-bearing inputs.

The metadata fixture corpus now also has one governed shared catalog at
`metadata/shared_trait_table_fixture_catalog.json`. That catalog assigns
durable fixture ids for continuous traits, binary and multistate discrete
traits, missing-value cases, extra-versus-missing taxon mismatches, duplicate
taxon negatives, constant-trait negatives, categorical predictors, and
misordered taxon rows. The owned discrete-reference validator and the live
`ape::ace` ancestral review tests now resolve their small governed trait
tables through those fixture ids instead of hardcoding one-off table paths.
That same trait-table catalog now also includes governed phylogenetic
independent-contrast fixtures for balanced, pectinate, six-taxon, and
missing-value comparative cases, so the owned `comparative contrasts` surface
and live `ape::pic` parity lane resolve the same durable trait identities.
The same catalog now also carries governed continuous ancestral Brownian cases
for balanced, pectinate, six-taxon, and pruned missing-value `ape::ace`
review, so `ancestral continuous` and the live continuous-`ace` parity lane
resolve the same durable trait identities instead of one-off table paths.
That same shared trait-table catalog now also carries governed discrete ER,
SYM, and ARD `ape::ace` cases for balanced binary, balanced multistate,
pectinate multistate, six-taxon, and pruned missing-value review, so
`ancestral discrete` and the live discrete-`ace` parity lane resolve the same
durable trait identities instead of one-off table paths.
The metadata fixture corpus now also has one governed shared catalog at
`metadata/shared_distance_matrix_fixture_catalog.json`. That catalog assigns
durable fixture ids to analytical, ultrametric, and non-ultrametric distance
matrices for owned neighbor-joining validation plus live `ape::nj` parity, so
distance-tree trust does not depend on loose file names or incidental matrix
paths.
The metadata fixture corpus now also has one governed shared simulation catalog
at `metadata/shared_tree_simulation_fixture_catalog.json`. That catalog assigns
durable fixture ids to rooted random-tree and coalescent simulation envelopes,
records replicate count, tip count, seed, branch-length model, and reference
function, and gives the live `ape::rtree` and `ape::rcoal` parity lane one
owned distribution-review surface instead of ad hoc seeded script fragments.
The same metadata corpus now also has one governed shared `phytools`
comparative catalog at
`metadata/shared_phytools_comparative_fixture_catalog.json`. That catalog pairs
shared tree fixtures and shared trait-table fixtures into one comparative
surface with durable fixture ids for twenty-four-taxon and one-hundred-
twenty-eight-taxon ultrametric signal cases, rooted non-ultrametric cases,
binary and multistate discrete traits, missing-value pruning, constant-trait
negatives, tree-versus-table mismatch, branch-length edge cases, and explicit
known-truth simulation metadata. The live `phytools` harness now resolves its
governed comparative cases through those ids instead of hardcoded ad hoc tree
and trait paths.
The same metadata corpus now also has one governed shared `geiger`
continuous-trait catalog at
`metadata/shared_geiger_continuous_fixture_catalog.json`. That catalog resolves
continuous `fitContinuous` review through durable fixture ids that point to
shared trait-table and tree fixtures, rather than letting the live `geiger`
lane hardcode recovery-panel paths. It covers twenty-four-taxon and one-
hundred-twenty-eight-taxon ultrametric surfaces, a rooted non-ultrametric
negative-control tree, Brownian, OU, and early-burst known-truth traits,
white-noise low-signal traits, missing-value pruning, constant-trait blockers,
one explicit outlier surface, one trend proxy for future model expansion, and
one governed per-taxon standard-error review surface retained for future
`fitContinuous(SE=...)` parity work.
The same metadata corpus now also has one governed shared `geiger`
discrete-trait catalog at
`metadata/shared_geiger_discrete_fixture_catalog.json`. That catalog resolves
future `fitDiscrete` review through durable fixture ids that point to shared
trait-table and tree fixtures instead of letting live parity lanes hardcode ad
hoc state tables. It covers ER binary known truth, SYM three-state and
four-state known truth, ARD binary and four-state known truth, one sparse
six-state overparameterization surface, missing-state pruning, one constant
negative, and one tree-versus-table mismatch panel, and it records explicit
transition-matrix metadata wherever the fixture reflects a governed simulation
surface.

The `expected/` directory now also carries benchmark-corpus regression
snapshots used to pin stable dataset summaries across releases.

The checked-in fixture surfaces are additionally grouped into clean, broken,
and messy benchmark corpora through `bijux_phylogenetics.validation`
so validation, warning-rich behavior, and failure signatures can be audited
from library code as well as tests.

The `engine_outputs/` fixture surface is the governed parser-stress corpus for
external tool outputs. It intentionally carries normal, warning-heavy,
truncated, malformed, and version-variant BEAST, MrBayes, and IQ-TREE parser
artifacts so external-output hardening stays tied to durable files instead of
ad hoc inline strings.

The `metadata/beast2_strict_yule_posterior.xml`, `.log`, and `.trees` files
are the governed real-artifact BEAST corpus used by the Bayesian validation
matrix whenever a live `beast` executable is not present locally. They are not
toy parser strings; they are the checked evidence surface for real BEAST output
structure.

That BEAST corpus is also cataloged at
`metadata/shared_beast_posterior_fixture_catalog.json`, including the governed
consensus tree, maximum clade credibility tree, burn-in counts, and posterior
summary references that anchor BEAST validation, publication, and reviewer
evidence surfaces to one durable fixture owner.

Tests should resolve files through the local `fixture(...)` helper rather than
assuming a flat directory layout.
