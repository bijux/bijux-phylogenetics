# Real FASTA-to-Tree Reference Inputs

These fixtures hold biologically real raw sequence inputs for end-to-end
tree-inference validation with installed engines.

They let the canonical `adapter fasta-to-tree` workflow run against reviewed
real data rather than only synthetic toy examples.

- Each FASTA file was derived from the source alignment recorded in
  `provenance.json` by removing alignment gap characters.
- Fixture names describe biological scope or stable ortholog identity, not
  temporary example numbering.
