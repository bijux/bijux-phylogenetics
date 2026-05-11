# Real FASTA-to-tree reference inputs

These fixtures carry biologically real sequence records derived from the repository-local trimAl reference corpus. They are stored here as raw FASTA inputs so the canonical `adapter fasta-to-tree` workflow can be exercised end to end with real engines and checked outputs.

Each FASTA file was produced by removing alignment gap characters from the named source alignment listed in `provenance.json`.
The fixture names describe the biological scope or stable external ortholog identifier rather than the temporary source example number.
