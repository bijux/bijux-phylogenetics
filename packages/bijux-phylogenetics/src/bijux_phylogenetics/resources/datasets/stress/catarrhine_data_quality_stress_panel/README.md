# Catarrhine Data Quality Stress Panel

This packaged dataset is a governed stress surface for dirty phylogenetic inputs.

The packaged Python surface lives under
`bijux_phylogenetics.datasets.data_quality_stress` and is intentionally split
by ownership: `panel` exposes dataset metadata, `traits` owns permissive raw
trait parsing and duplicate or missing-value review, `cleanup` owns sequence,
tree, and comparative-subset repair policy, `export` copies the packaged raw
panel, `bundle` writes the reviewer-facing ledgers, and `demo` materializes
the full review bundle.

It combines three deliberately different provenance layers:

- a real catarrhine mitochondrial supermatrix derived from the packaged multi-locus catarrhine panel
- a derived rooted tree whose topology follows the same catarrhine panel but whose branch lengths intentionally retain one zero-length branch, one negative branch, and one extreme long branch
- a synthetic trait table keyed to the same taxa and intentionally containing one duplicate taxon row plus missing values
- a separate raw FASTA input keyed to the same taxa and intentionally containing one duplicate identifier, one illegal character, one empty record, and explicit length outliers
- a separate raw coding FASTA input containing one frame error and one premature stop codon
- a separate raw trait-linkage table containing one tree taxon missing from traits and one extra trait taxon absent from the tree

The goal is not to pretend this is a clean scientific reference dataset.
The goal is to expose realistic data-quality failure modes in one durable public demo surface so reviewers can inspect how `bijux-phylogenetics` identifies and handles them.

Expected stress conditions:

- duplicate taxa in the raw traits table
- duplicate sequence identifiers in the raw FASTA validation input
- illegal FASTA characters in the raw FASTA validation input
- empty FASTA records in the raw FASTA validation input
- raw-sequence length outliers in the raw FASTA validation input
- one coding sequence with a frame error
- one coding sequence with a premature stop codon
- missing trait values in both required and non-required columns
- one raw tree-versus-traits mismatch surface with one missing tree taxon and one extra trait taxon
- one alignment composition outlier taxon
- one terminal branch-length outlier taxon
- one zero-length branch and one negative branch requiring repair before comparative reuse

The governed workflow writes explicit ledgers for discovery and handling, then materializes one cleaned comparative subset rather than hiding the exclusions.
