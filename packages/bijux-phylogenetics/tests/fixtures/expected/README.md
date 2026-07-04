# Expected Outputs

This directory holds the checked reviewer-facing outputs used as deterministic
reference bundles for runtime tests.

Use these fixtures when a test needs to confirm a stable public artifact rather
than only an internal metric or intermediate structure.

- `beast2_strict_yule_posterior/` governs BEAST posterior reference summaries,
  including the checked 10% burn-in consensus and maximum-clade-credibility
  trees.
- `concatenation/` governs expected multi-locus supermatrix outputs.
- `fasta_to_tree/` governs end-to-end expected bundles for the canonical raw
  FASTA to inferred tree workflow.
