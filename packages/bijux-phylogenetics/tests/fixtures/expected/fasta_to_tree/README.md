# Real FASTA-to-tree expected bundles

Each directory in this fixture surface contains the checked reviewer-facing output bundle for one real raw FASTA dataset run through the canonical `adapter fasta-to-tree` workflow.

The checked files are limited to the durable outputs promised by the workflow contract:

- `.aln`
- `.trimmed.aln`
- `.tree`
- `.log`
- `.model.tsv`
- `.support.tsv`
- `.manifest.json`
- `.run.json`

Engine-specific intermediates remain runtime artifacts and are not tracked in this governed corpus.

Real-engine validation compares these bundles semantically rather than byte-for-byte.
That keeps the scientific checks stable across harmless path and timestamp variation while still verifying the reviewer-facing workflow contract.
