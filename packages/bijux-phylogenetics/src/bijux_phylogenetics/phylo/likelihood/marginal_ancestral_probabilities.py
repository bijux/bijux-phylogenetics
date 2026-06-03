from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    normalize_unambiguous_dna_records,
    one_hot_dna_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    MarginalAncestralSequenceProbabilityReport,
    MarginalAncestralStateProbabilityRow,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.posteriors import (
    compute_marginal_state_posteriors,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def validate_nucleotide_marginal_ancestral_probability_model(model_name: str) -> str:
    return validate_selected_nucleotide_likelihood_model(model_name)


def evaluate_nucleotide_marginal_ancestral_probabilities(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> MarginalAncestralSequenceProbabilityReport:
    """Evaluate internal-node marginal ancestral state probabilities per site."""
    normalized_model_name = validate_nucleotide_marginal_ancestral_probability_model(
        model_name
    )
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name=normalized_model_name.upper(),
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=(
            f"{normalized_model_name.upper()} marginal ancestral probability evaluation"
        ),
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return _evaluate_selected_nucleotide_marginal_probabilities_from_patterns(
        tree,
        compressed_patterns,
        model_name=specification.model_name,
        root_prior=specification.root_prior,
        parameter_values=specification.parameter_values,
        transition_matrix_for_child=lambda child: (
            specification.transition_matrix_for_branch_length(
                max(float(child.branch_length or 0.0), 0.0)
            )
        ),
    )


def evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> MarginalAncestralSequenceProbabilityReport:
    """Evaluate one fixed-topology nucleotide posterior report from file paths."""
    return evaluate_nucleotide_marginal_ancestral_probabilities(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def _evaluate_selected_nucleotide_marginal_probabilities_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
    root_prior: numpy.ndarray,
    parameter_values: dict[str, float],
    transition_matrix_for_child,
) -> MarginalAncestralSequenceProbabilityReport:
    validate_explicit_branch_lengths(tree, model_name=model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=model_name,
    )

    internal_nodes = list(tree.iter_internal_nodes(order="preorder"))
    posterior_rows: list[MarginalAncestralStateProbabilityRow] = []
    for pattern in compressed_patterns.patterns:
        states_by_taxon = dict(
            zip(compressed_patterns.taxon_order, pattern.states, strict=True)
        )
        current_states_by_taxon = states_by_taxon
        posterior_pass = compute_marginal_state_posteriors(
            tree,
            state_count=len(DNA_STATE_ORDER),
            leaf_likelihood=lambda node, current_states_by_taxon=current_states_by_taxon: (
                one_hot_dna_leaf_vector(
                    current_states_by_taxon,
                    model_name=model_name,
                    node_name=node.name,
                )
            ),
            transition_matrix_for_child=transition_matrix_for_child,
            root_prior=root_prior,
        )
        for site_position in pattern.site_positions:
            for node in internal_nodes:
                node_probabilities = posterior_pass.posterior_for_node(node)
                for state, posterior_probability in zip(
                    DNA_STATE_ORDER,
                    node_probabilities,
                    strict=True,
                ):
                    posterior_rows.append(
                        MarginalAncestralStateProbabilityRow(
                            node_id=node.node_id or "",
                            node_name=node.name,
                            descendant_taxa=node.descendant_taxa,
                            pattern_id=pattern.pattern_id,
                            site_position=site_position,
                            state=state,
                            posterior_probability=float(posterior_probability),
                        )
                    )
    return MarginalAncestralSequenceProbabilityReport(
        model_name=model_name,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        internal_node_count=len(internal_nodes),
        compression_used=True,
        expansion_policy="expanded-internal-node-site-state-rows",
        tree_newick=dumps_newick(tree),
        parameter_values=parameter_values,
        posterior_rows=posterior_rows,
    )
