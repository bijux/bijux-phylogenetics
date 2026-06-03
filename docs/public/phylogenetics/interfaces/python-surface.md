---
title: Python Surface
audience: public
type: reference
status: active
owner: bijux-phylogenetics-docs
last_reviewed: 2026-05-16
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
The same owned comparative runtime also underlies the governed live
`phytools::pgls.SEy` lane for fixed-lambda Brownian covariance regression over
simple numeric, categorical, and interaction-coded fixtures. That governed
claim is narrower than the full owned PGLS API on purpose: installed
`phytools 2.5.2` does not export a general `phytools::pgls` surface, so live
parity stays on `pgls.SEy` with `lambda = 1.0`, while estimated-lambda and
broader exact regression parity remain covered by the checked-in `ape` plus
`nlme` reference suite.

The owned ancestral runtime now also exposes direct dataset-backed
reconstruction surfaces through
`reconstruct_continuous_ancestral_states_from_dataset(...)` and
`reconstruct_discrete_ancestral_states_from_dataset(...)`. Once a caller
already holds one `AncestralContinuousDataset` or `AncestralDiscreteDataset`,
continuous and discrete ancestral reconstruction no longer need to restart
from path-based loading wrappers. The discrete ancestral report now also carries
one explicit `rerooting_method_compatibility` contract so Python callers can
tell whether one ER or SYM reconstruction with the equal root prior is
comparable to the governed live `phytools::rerootingMethod` lane. Fitch,
ordered-state, ARD, empirical-root-prior, and fixed-root-prior runs are marked
explicitly as non-comparable instead of being left to inference.
The same owned discrete-evolution runtime now also exposes seeded stochastic
character mapping through `simulate_discrete_stochastic_maps(...)` and
`simulate_discrete_stochastic_maps_from_fit_report(...)` plus review-friendly
writers for the resulting collection JSON, transition summary, flat
branch-segment ledger, per-state time-in-state ledger, and per-branch
state-occupancy ledger. That surface fits
one discrete CTMC, carries one explicit fitted-model audit with model
identity, parameter count, log-likelihood, AIC, AICc, baseline-model
comparison, optimizer convergence, and weak-fit warnings, conditions one
sampled history on the observed tip states, reports failed branch-history
draws explicitly, and now underlies the governed live
`phytools::make.simmap(model='ER')`, `phytools::make.simmap(model='SYM')`, and
`phytools::make.simmap(model='ARD')` summary-parity lanes without claiming
exact stochastic-history identity across languages. Governed multistate ARD
cases stay on summary-envelope parity only when weakly identified boundary
rates make row-level transition summaries unstable across optimizers.
The same owned runtime also exposes
`bijux_phylogenetics.count_discrete_stochastic_map_transitions(...)` plus
writers for one per-replicate transition-count matrix, one aggregate
transition matrix, one per-branch directional transition table, and one flat
event ledger. That surface now underlies the governed live
`phytools::countSimmap` lane, which compares total-transition envelopes and
directional transition-count rows over selected seeded map collections,
including zero diagonal state pairs, without claiming exact stochastic-history
identity.
The same owned stochastic-map runtime also exposes
`bijux_phylogenetics.summarize_discrete_stochastic_map_density(...)`,
`bijux_phylogenetics.render_stochastic_map_density_artifact(...)`, and writers
for one branch-probability table, one branch-level density envelope, and one
slice-level probability table. Binary collections can summarize one default
focal state directly; multistate collections keep branch probability summaries
for every state, and require one explicit focal state before generating
density-slice rows or one report-ready artifact. That surface now underlies
the governed live `phytools::densityMap` lane for selected binary ER
collections, including one missing-value-pruned case, and compares
branch-level posterior probability summaries plus branch uncertainty without
claiming pixel-perfect plot parity.
The same owned summary contract also underlies the governed live
`phytools::describe.simmap` lane, which compares total-change summary,
transition-count rows, time-in-state rows, and per-branch state-occupancy rows
over selected seeded map collections.
The same owned simulation surface also exposes
`bijux_phylogenetics.simulate_brownian_trait_collection(...)` plus writers for
one replicate trait ledger and one Brownian summary ledger over tip
distributions and tip covariances. The one-trait Brownian entrypoint
`bijux_phylogenetics.simulate_brownian_traits(...)` now accepts either
`sigma` or explicit `sigma_squared`, keeps the resolved Brownian rate
parameter in the returned report, and shares the same fixed-tree covariance
contract as the collection surface. That runtime now underlies the governed
live `phytools::fastBM` lane over selected low-variance, root-shift
high-variance, and six-taxon fixed-tree cases, comparing distribution
summaries and tip-covariance rows without claiming exact cross-language
draw identity.
The same owned simulation surface also exposes
`bijux_phylogenetics.simulate_correlated_brownian_trait_collection(...)` plus
writers for one long-form replicate tip-trait ledger and one multivariate
summary ledger over root states, evolutionary covariance, tip distributions,
and tip covariances. It accepts either one explicit covariance matrix or one
correlation matrix plus per-trait standard deviations, rejects invalid
non-positive-definite covariance inputs explicitly, and keeps the generating
parameter matrix in the returned report. That runtime now underlies the
governed live `phytools::sim.corrs` lane over selected low-correlation,
negative-correlation root-shift, and three-trait fixed-tree cases, comparing
distribution summaries, tip-covariance rows, and tip-correlation rows without
claiming exact cross-language draw identity.
The same owned simulation surface also exposes
`bijux_phylogenetics.simulate_discrete_histories(...)` plus writers for one
tip-state truth table, one node-state truth table, one branch-history truth
table, one transition-event ledger, one branch-segment ledger, and one parity
summary table. That surface simulates discrete histories on one fixed tree
from an explicit rate matrix, supports binary and multistate states, supports
fixed or probabilistic root states, and keeps one seeded truth surface for
downstream recovery tests. It now underlies the governed live
`phytools::sim.history` lane over selected no-change and high-rate fixed-tree
cases, comparing total-transition summaries, transition-count rows,
time-in-state rows, and tip-state-frequency rows without claiming exact
cross-language history identity.

For end-to-end external-engine orchestration, the public engine surface includes
`bijux_phylogenetics.run_fasta_to_tree_workflow(...)`. That workflow owns the
raw-FASTA to aligned matrix, trimmed matrix, selected-model table, supported
tree, support-summary table, log, and manifest contract used by the CLI.

The public promise is ownership clarity: imports should resolve to one runtime
family, not to competing public surfaces with different meanings.
