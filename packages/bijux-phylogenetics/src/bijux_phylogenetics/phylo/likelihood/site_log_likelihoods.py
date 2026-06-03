from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    evaluate_fixed_topology_dna_site_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    normalize_dna_likelihood_records,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FixedTopologySiteLogLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
    validate_site_log_likelihood_reconstruction,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def validate_nucleotide_site_log_likelihood_model(model_name: str) -> str:
    return validate_selected_nucleotide_likelihood_model(model_name)


def evaluate_nucleotide_site_log_likelihoods(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    observation_policy: str = "reject",
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
) -> FixedTopologySiteLogLikelihoodReport:
    """Evaluate one fixed-topology nucleotide likelihood with expanded site rows."""
    normalized_model_name = validate_nucleotide_site_log_likelihood_model(model_name)
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name=normalized_model_name.upper(),
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} site log likelihood export",
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return evaluate_selected_dna_site_log_likelihoods_from_patterns(
        tree,
        compressed_patterns,
        model_name=specification.model_name,
        state_count=specification.state_count,
        observation_policy=specification.observation_policy,
        root_prior=specification.root_prior,
        parameter_values=specification.parameter_values,
        transition_matrix_for_child=lambda child: (
            specification.transition_matrix_for_branch_length(
                max(float(child.branch_length or 0.0), 0.0)
            )
        ),
    )


def evaluate_nucleotide_site_log_likelihoods_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    observation_policy: str = "reject",
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
) -> FixedTopologySiteLogLikelihoodReport:
    """Evaluate one fixed-topology nucleotide likelihood from paths with site rows."""
    return evaluate_nucleotide_site_log_likelihoods(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def evaluate_selected_dna_site_log_likelihoods_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
    state_count: int,
    observation_policy: str,
    root_prior: numpy.ndarray,
    parameter_values: dict[str, float],
    transition_matrix_for_child,
) -> FixedTopologySiteLogLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name=model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=model_name,
    )
    site_log_likelihoods, total_log_likelihood = (
        expanded_site_log_likelihood_rows_from_patterns(
            compressed_patterns,
            site_log_likelihood=lambda states: (
                evaluate_fixed_topology_dna_site_log_likelihood(
                    tree,
                    states,
                    taxon_order=compressed_patterns.taxon_order,
                    model_name=model_name,
                    observation_policy=observation_policy,
                    root_prior=root_prior,
                    transition_matrix_for_child=transition_matrix_for_child,
                )
            ),
        )
    )
    validate_site_log_likelihood_reconstruction(
        site_log_likelihoods,
        expected_total_log_likelihood=total_log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name=f"{model_name} site log likelihood export",
    )
    return FixedTopologySiteLogLikelihoodReport(
        model_name=model_name,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        expansion_policy="expanded-site-rows",
        tree_newick=dumps_newick(tree),
        state_count=state_count,
        observation_policy=observation_policy,
        parameter_values=parameter_values,
        log_likelihood=total_log_likelihood,
        site_log_likelihoods=site_log_likelihoods,
    )
