# Real FASTA-to-tree expected bundles

Each directory in this fixture surface contains the checked reviewer-facing output bundle for one real raw FASTA dataset run through the canonical `adapter fasta-to-tree` workflow.

The checked files are limited to the durable outputs promised by the workflow contract:

- `.aln`
- `.trimmed.aln`
- `.tree`
- `.log`
- `.model.tsv`
- `.support.tsv`

Engine-specific intermediates remain runtime artifacts and are not tracked in this golden corpus.
