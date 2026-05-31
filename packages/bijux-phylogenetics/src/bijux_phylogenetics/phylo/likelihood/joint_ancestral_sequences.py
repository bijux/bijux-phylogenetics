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
from bijux_phylogenetics.phylo.likelihood.joint_states import (
    compute_joint_state_assignment,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    JointAncestralSequenceRecord,
    JointAncestralSequenceReport,
    JointAncestralStateAssignmentRow,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def validate_nucleotide_joint_ancestral_sequence_model(model_name: str) -> str:
    return validate_selected_nucleotide_likelihood_model(model_name)


def reconstruct_nucleotide_joint_ancestral_sequences(
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
) -> JointAncestralSequenceReport:
    """Reconstruct one globally optimal internal-node state assignment per site."""
    normalized_model_name = validate_nucleotide_joint_ancestral_sequence_model(
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
        owner_name=f"{normalized_model_name.upper()} joint ancestral reconstruction",
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return _reconstruct_selected_nucleotide_joint_sequences_from_patterns(
        tree,
        compressed_patterns,
        model_name=specification.model_name,
        parameter_values=specification.parameter_values,
        root_prior=specification.root_prior,
        transition_matrix_for_child=lambda child: (
            specification.transition_matrix_for_branch_length(
                max(float(child.branch_length or 0.0), 0.0)
            )
        ),
    )


def reconstruct_nucleotide_joint_ancestral_sequences_from_alignment(
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
) -> JointAncestralSequenceReport:
    """Reconstruct one joint ancestral sequence report from file paths."""
    return reconstruct_nucleotide_joint_ancestral_sequences(
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


def _reconstruct_selected_nucleotide_joint_sequences_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
    parameter_values: dict[str, float],
    root_prior: numpy.ndarray,
    transition_matrix_for_child,
) -> JointAncestralSequenceReport:
    validate_explicit_branch_lengths(tree, model_name=model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=model_name,
    )

    internal_nodes = list(tree.iter_internal_nodes(order="preorder"))
    assignment_rows: list[JointAncestralStateAssignmentRow] = []
    sequence_by_node_id = {node.node_id or "": [] for node in internal_nodes}
    for pattern in compressed_patterns.patterns:
        states_by_taxon = dict(
            zip(compressed_patterns.taxon_order, pattern.states, strict=True)
        )
        current_states_by_taxon = states_by_taxon
        joint_assignment = compute_joint_state_assignment(
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
        assigned_state_by_node_id = {
            node.node_id or "": DNA_STATE_ORDER[
                joint_assignment.assigned_state_index_for_node(node)
            ]
            for node in internal_nodes
        }
        for site_position in pattern.site_positions:
            for node in internal_nodes:
                node_id = node.node_id or ""
                state = assigned_state_by_node_id[node_id]
                sequence_by_node_id[node_id].append(state)
                assignment_rows.append(
                    JointAncestralStateAssignmentRow(
                        node_id=node_id,
                        node_name=node.name,
                        descendant_taxa=node.descendant_taxa,
                        pattern_id=pattern.pattern_id,
                        site_position=site_position,
                        state=state,
                    )
                )
    sequence_records = [
        JointAncestralSequenceRecord(
            node_id=node.node_id or "",
            node_name=node.name,
            descendant_taxa=node.descendant_taxa,
            sequence="".join(sequence_by_node_id[node.node_id or ""]),
        )
        for node in internal_nodes
    ]
    return JointAncestralSequenceReport(
        model_name=model_name,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        internal_node_count=len(internal_nodes),
        compression_used=True,
        expansion_policy="expanded-internal-node-site-joint-state-rows",
        tree_newick=dumps_newick(tree),
        parameter_values=parameter_values,
        sequence_records=sequence_records,
        assignment_rows=assignment_rows,
    )
