# Catarrhine Data Quality Stress Panel

This packaged dataset is a governed stress surface for dirty phylogenetic inputs.

It combines three deliberately different provenance layers:

- a real catarrhine mitochondrial supermatrix derived from the packaged multi-locus catarrhine panel
- a derived rooted tree whose topology follows the same catarrhine panel but whose branch lengths intentionally retain one zero-length branch and one extreme long branch
- a synthetic trait table keyed to the same taxa and intentionally containing one duplicate taxon row plus missing values

The goal is not to pretend this is a clean scientific reference dataset.
The goal is to expose realistic data-quality failure modes in one durable public demo surface so reviewers can inspect how `bijux-phylogenetics` identifies and handles them.

Expected stress conditions:

- duplicate taxa in the raw traits table
- missing trait values in both required and non-required columns
- one alignment composition outlier taxon
- one terminal branch-length outlier taxon
- one zero-length branch requiring repair before comparative reuse

The governed workflow writes explicit ledgers for discovery and handling, then materializes one cleaned comparative subset rather than hiding the exclusions.
