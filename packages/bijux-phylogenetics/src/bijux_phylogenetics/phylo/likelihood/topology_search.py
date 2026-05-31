from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from statistics import median

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import normalize_unambiguous_dna_records
from bijux_phylogenetics.phylo.likelihood.fixed_topology_branch_lengths import (
    BranchReoptimizationResult,
    evaluate_selected_nucleotide_log_likelihood_from_patterns,
    optimize_selected_nucleotide_branch_length_subset,
    optimize_selected_nucleotide_branch_lengths,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    SelectedNucleotideLikelihoodSpecification,
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.substitution_parameters import (
    optimize_nucleotide_substitution_parameters,
)
from bijux_phylogenetics.phylo.likelihood.validation import validate_explicit_branch_lengths
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id, split_sort_key
from bijux_phylogenetics.phylo.topology.tree import descendant_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES = frozenset({"coordinate-branch-lengths"})


@dataclass(frozen=True, slots=True)
class ResolvedNucleotideTopologySearchSurface:
    """Resolved nucleotide likelihood surface used during topology search."""

    model_name: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    specification: SelectedNucleotideLikelihoodSpecification


def resolve_nucleotide_topology_search_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    """Resolve one topology-search tree input and preserve its path when present."""
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_tree(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def resolve_nucleotide_topology_search_records(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    """Resolve one topology-search alignment input and preserve its path when present."""
    resolved_alignment_path = records if isinstance(records, Path) else None
    resolved_records = (
        load_fasta_alignment(records) if isinstance(records, Path) else list(records)
    )
    return resolved_records, resolved_alignment_path


def validate_nucleotide_topology_search_tree(
    tree: PhyloTree,
    *,
    workflow_name: str,
) -> None:
    """Require one structurally valid strictly binary rooted tree."""
    errors = tree.validation_errors()
    if errors:
        raise ValueError(
            f"{workflow_name} requires a structurally valid tree: {'; '.join(errors)}"
        )
    if len(tree.root.children) != 2:
        raise ValueError(f"{workflow_name} requires a rooted binary tree")
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ValueError(
            f"{workflow_name} requires a strictly binary tree: {invalid_internal_nodes}"
        )


def initialize_generated_nucleotide_topology_search_tree(
    tree: PhyloTree,
    *,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
) -> PhyloTree:
    """Seed missing and out-of-bound branch lengths on one generated topology candidate."""
    working_tree = tree.copy().refresh()
    existing_branch_lengths = [
        child.branch_length
        for _parent, child in working_tree.iter_edges()
        if child.branch_length is not None
    ]
    default_branch_length = (
        float(median(existing_branch_lengths))
        if existing_branch_lengths
        else lower_branch_length_bound
    )
    seeded_default_branch_length = min(
        upper_branch_length_bound,
        max(lower_branch_length_bound, default_branch_length),
    )
    for _parent, child in working_tree.iter_edges():
        branch_length = child.branch_length
        if branch_length is None:
            child.branch_length = seeded_default_branch_length
            continue
        child.branch_length = min(
            upper_branch_length_bound,
            max(lower_branch_length_bound, float(branch_length)),
        )
    return working_tree.refresh()


def validate_branch_reoptimization_policy(policy: str) -> str:
    normalized_policy = policy.strip().lower()
    if normalized_policy not in _SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES:
        raise ValueError(
            "branch_reoptimization_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES))
        )
    return normalized_policy


def resolve_nucleotide_topology_search_surface(
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
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
) -> ResolvedNucleotideTopologySearchSurface:
    """Resolve one selected nucleotide likelihood surface for topology search."""
    normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name=f"{normalized_model_name.upper()} likelihood NNI search",
    )

    if normalized_model_name == "jc69":
        specification = resolve_selected_nucleotide_likelihood_specification(
            normalized_records,
            model_name=normalized_model_name,
            owner_name="JC69 likelihood NNI search",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideTopologySearchSurface(
            model_name=specification.model_name,
            substitution_parameter_policy="fixed-from-model",
            substitution_parameter_values=specification.parameter_values,
            substitution_parameter_warnings=[],
            specification=specification,
        )

    if normalized_model_name == "f81" and base_frequencies is None:
        specification = resolve_selected_nucleotide_likelihood_specification(
            normalized_records,
            model_name=normalized_model_name,
            owner_name="F81 likelihood NNI search",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideTopologySearchSurface(
            model_name=specification.model_name,
            substitution_parameter_policy="estimated-from-alignment",
            substitution_parameter_values=specification.parameter_values,
            substitution_parameter_warnings=[],
            specification=specification,
        )

    parameters_are_fully_provided = (
        (normalized_model_name == "k80" and kappa is not None)
        or (
            normalized_model_name == "f81"
            and base_frequencies is not None
        )
        or (
            normalized_model_name == "hky85"
            and kappa is not None
        )
        or (
            normalized_model_name == "gtr"
            and exchangeabilities is not None
        )
    )
    if parameters_are_fully_provided:
        specification = resolve_selected_nucleotide_likelihood_specification(
            normalized_records,
            model_name=normalized_model_name,
            owner_name=f"{normalized_model_name.upper()} likelihood NNI search",
            kappa=kappa,
            base_frequencies=base_frequencies,
            exchangeabilities=exchangeabilities,
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideTopologySearchSurface(
            model_name=specification.model_name,
            substitution_parameter_policy="provided",
            substitution_parameter_values=specification.parameter_values,
            substitution_parameter_warnings=[],
            specification=specification,
        )

    optimization_report = optimize_nucleotide_substitution_parameters(
        tree,
        normalized_records,
        model_name=normalized_model_name,
        base_frequencies=base_frequencies,
        initial_kappa=1.0 if normalized_model_name in {"k80", "hky85"} else None,
        initial_exchangeabilities=exchangeabilities,
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} likelihood NNI search",
        kappa=_resolved_kappa(optimization_report),
        base_frequencies=_resolved_base_frequencies(optimization_report),
        exchangeabilities=_resolved_exchangeabilities(optimization_report),
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return ResolvedNucleotideTopologySearchSurface(
        model_name=specification.model_name,
        substitution_parameter_policy="optimized-on-start-topology",
        substitution_parameter_values=specification.parameter_values,
        substitution_parameter_warnings=list(optimization_report.warnings),
        specification=specification,
    )


def normalize_nucleotide_topology_search_records(
    records: list[AlignmentRecord],
    *,
    owner_name: str,
) -> tuple[list[AlignmentRecord], CompressedAlignmentSitePatterns]:
    """Normalize one nucleotide alignment and precompute compressed site patterns."""
    normalized_records = normalize_unambiguous_dna_records(records, model_name=owner_name)
    return normalized_records, compress_alignment_site_patterns_from_records(
        normalized_records
    )


def reoptimize_nucleotide_topology_tree(
    tree: PhyloTree,
    *,
    compressed_patterns: CompressedAlignmentSitePatterns,
    resolved_surface: ResolvedNucleotideTopologySearchSurface,
    branch_reoptimization_policy: str,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> BranchReoptimizationResult:
    """Reoptimize one topology-search tree under one resolved likelihood surface."""
    validated_branch_reoptimization_policy = validate_branch_reoptimization_policy(
        branch_reoptimization_policy
    )
    if validated_branch_reoptimization_policy not in _SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES:
        raise ValueError(
            "branch_reoptimization_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES))
        )
    return optimize_selected_nucleotide_branch_lengths(
        tree,
        compressed_patterns,
        specification=resolved_surface.specification,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def reoptimize_nucleotide_topology_tree_branch_subset(
    tree: PhyloTree,
    *,
    compressed_patterns: CompressedAlignmentSitePatterns,
    resolved_surface: ResolvedNucleotideTopologySearchSurface,
    optimized_branch_ids: list[str],
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> BranchReoptimizationResult:
    """Reoptimize one declared subset of branch lengths on one topology-search tree."""
    return optimize_selected_nucleotide_branch_length_subset(
        tree,
        compressed_patterns,
        specification=resolved_surface.specification,
        optimized_branch_ids=optimized_branch_ids,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def resolve_reoptimized_branch_clade_ids(
    tree: PhyloTree,
    optimized_branch_ids: list[str],
) -> list[str]:
    """Render one deterministic descendant-clade ledger for reoptimized branch identifiers."""
    signatures = [
        frozenset(descendant_taxa(tree.node_by_id(branch_id)))
        for branch_id in optimized_branch_ids
    ]
    return [
        canonical_clade_id(signature)
        for signature in sorted(signatures, key=split_sort_key)
    ]


def prefer_higher_likelihood(
    left_score: float,
    left_newick: str,
    right_score: float,
    right_newick: str,
) -> bool:
    """Prefer higher likelihoods, then deterministic Newick order, across search moves."""
    if left_score > right_score and not math.isclose(left_score, right_score):
        return True
    if right_score > left_score and not math.isclose(left_score, right_score):
        return False
    return left_newick < right_newick


def _resolved_base_frequencies(
    optimization_report,
) -> dict[str, float] | numpy.ndarray | None:
    if optimization_report.base_frequency_a is None:
        return None
    return {
        "A": optimization_report.base_frequency_a,
        "C": optimization_report.base_frequency_c or 0.0,
        "G": optimization_report.base_frequency_g or 0.0,
        "T": optimization_report.base_frequency_t or 0.0,
    }


def _resolved_kappa(optimization_report) -> float | None:
    row_by_name = {
        row.parameter_name: row.optimized_value
        for row in optimization_report.parameter_rows
    }
    return row_by_name.get("kappa")


def _resolved_exchangeabilities(
    optimization_report,
) -> dict[str, float] | None:
    if optimization_report.model_name != "GTR":
        return None
    row_by_name = {
        row.parameter_name: row.optimized_value
        for row in optimization_report.parameter_rows
    }
    return {
        "AC": optimization_report.fixed_parameter_values.get("AC", 1.0),
        "AG": row_by_name["AG"],
        "AT": row_by_name["AT"],
        "CG": row_by_name["CG"],
        "CT": row_by_name["CT"],
        "GT": row_by_name["GT"],
    }
