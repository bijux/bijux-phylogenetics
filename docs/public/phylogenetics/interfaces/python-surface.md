---
title: Python Surface
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-15
---

# Python Surface

The canonical Python import surface lives under `bijux_phylogenetics`.

Use the canonical package name when you need the durable runtime API. The
`phylogenetic` package exists as a compatibility alias, not as a second
independent runtime.

The owned native tree API now lives on `bijux_phylogenetics.PhyloTree`,
`bijux_phylogenetics.TreeNode`, and `bijux_phylogenetics.TaxonLabel`. That
surface is the single in-memory tree contract for native traversal, topology
transforms, branch-length review, comparative covariance, ancestral
reconstruction, simulation, and canonical Newick conversion. Stable node IDs,
parent-child links, node metadata, edge metadata, validation, deep-copy
behavior, native Newick parsing and writing, multi-tree Newick loading, and
location-aware Newick parse failures are part of that runtime promise.
Outgroup rooting, unrooting, keep-tip pruning, drop-tip pruning, clade
extraction, MRCA lookup, and monophyly review are also part of that same
owned tree-manipulation core, so those surfaces no longer depend on an
external tree object model for normal runtime behavior.
The same native ownership boundary now covers canonical rooted-clade
extraction, canonical unrooted-split extraction, Robinson-Foulds metrics, and
clade-support matching, so tree distance, topology comparison, tree-set
support, posterior clade frequencies, and live `ape::dist.topo` parity all
read one shared split identity contract.
The same owned runtime now also loads Newick tree sets directly into
`PhyloTree` records for consensus building, clade-frequency summaries,
reference-tree support mapping, topology clustering, and posterior tree-set
comparison. Strict consensus and support surfaces validate one exact taxon set
across the whole tree set, while tolerant inspection surfaces keep one
explicit malformed-record counter instead of failing silently.
The same native tree runtime now also owns direct baseline simulation entry
points through `simulate_random_tree(...)` and `simulate_coalescent_tree(...)`
when one caller needs one governed random or coalescent tree plus its summary
record without going through one batch simulation wrapper.

The owned native DNA-alignment API now also lives on one
`bijux_phylogenetics.DnaBinAlignment` runtime loaded through
`bijux_phylogenetics.load_dna_bin_alignment(...)`. That matrix preserves taxon
order and alignment length, normalizes case, keeps gaps, ambiguity codes, and
explicit missing states literal, and rejects unsupported symbols explicitly.
The same runtime object now feeds direct Python nucleotide review surfaces for
literal-state composition, ape-style segregating-site detection, and
nucleotide distances through
`compute_alignment_base_frequency_report_from_dna_bin_alignment(...)`,
`compute_alignment_segregating_site_report_from_dna_bin_alignment(...)`, and
`compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(...)`. It now
also feeds aligned coding diagnostics and aligned translation through
`inspect_coding_alignment_from_dna_bin_alignment(...)` and
`translate_coding_alignment_from_dna_bin_alignment(...)` instead of forcing
those workflows to reparse FASTA independently.
The same owned distance runtime now also exposes
`build_distance_tree_from_genetic_distance_matrix(...)`, so one loaded
`GeneticDistanceMatrix` can recover one governed NJ or UPGMA tree directly
without restarting from path-based input loaders.

The owned comparative runtime now also exposes one direct Brownian covariance
surface on an in-memory tree through
`summarize_brownian_covariance_from_tree(...)` and one direct PIC surface on a
loaded comparative dataset through
`compute_phylogenetic_independent_contrasts_from_dataset(...)`. Once a caller
already holds one `PhyloTree` or `ComparativeDataset`, covariance review and
independent-contrast analysis no longer need to restart from path-based
loading wrappers.
The same comparative runtime now also exposes one direct discrete Mk fit
surface on a loaded ancestral discrete dataset through
`fit_discrete_mk_model_from_dataset(...)`, so ER, SYM, and ARD review no
longer need to restart from tree and trait file paths once one
`AncestralDiscreteDataset` already exists in memory. That same owned surface
now underlies the governed live `phytools::fitMk(model='ER')` and
`phytools::fitMk(model='SYM')` parity lanes for the validated ER and
unordered-multistate SYM cases, and it now also underlies the governed live
`phytools::fitMk(model='ARD')` lane for one rate-row-parity binary surface
plus one summary-parity multistate surface when the optimizer flags weakly
identified boundary rates.

The owned ancestral runtime now also exposes direct dataset-backed
reconstruction surfaces through
`reconstruct_continuous_ancestral_states_from_dataset(...)` and
`reconstruct_discrete_ancestral_states_from_dataset(...)`. Once a caller
already holds one `AncestralContinuousDataset` or `AncestralDiscreteDataset`,
continuous and discrete ancestral reconstruction no longer need to restart
from path-based loading wrappers.

For end-to-end external-engine orchestration, the public engine surface includes
`bijux_phylogenetics.run_fasta_to_tree_workflow(...)`. That workflow owns the
raw-FASTA to aligned matrix, trimmed matrix, selected-model table, supported
tree, support-summary table, log, and manifest contract used by the CLI.

The public promise is ownership clarity: imports should resolve to one runtime
family, not to competing public surfaces with different meanings.
