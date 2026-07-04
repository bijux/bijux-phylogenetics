# Real FASTA-to-Tree Expected Bundles

Each directory in this surface holds the checked output bundle for one
biologically real raw FASTA dataset run through the canonical
`adapter fasta-to-tree` workflow.

## Governed Output Set

The checked bundle is intentionally limited to the durable outputs promised by
the workflow contract:

- `.aln`
- `.trimmed.aln`
- `.tree`
- `.log`
- `.methods-summary.md`
- `.model.tsv`
- `.support.tsv`
- `.manifest.json`
- `.run.json`

Engine-specific scratch files and other incidental intermediates remain runtime
artifacts and are not treated as governed expected outputs here.

## Comparison Policy

Real-engine validation compares these bundles semantically rather than
byte-for-byte. That keeps the review lane stable across harmless path,
timestamp, and environment variation while still governing the reader-facing
workflow contract.

The comparison surfaces are typed deliberately:

- trees must preserve compatible topology, tip identity, and support meaning
- tabular ledgers must preserve schema, row identity, and tolerated numeric
  values
- alignments must preserve identifier-to-sequence content
- manifests and run records must preserve structured output ownership and hash
  coverage
