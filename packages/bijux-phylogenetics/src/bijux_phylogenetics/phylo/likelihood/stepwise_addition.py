from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology import (
    PhyloTree,
    StepwiseAdditionCandidateScore,
    StepwiseAdditionTraceRow,
    StepwiseAdditionTreeReport,
    TreeNode,
    apply_stepwise_addition_candidate,
    iter_stepwise_addition_edge_candidates,
    summarize_stepwise_addition_tree,
    validate_stepwise_addition_taxa,
)

from .topology_search import (
    normalize_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
)
from .fixed_topology_branch_lengths import optimize_selected_nucleotide_branch_lengths

_SUPPORTED_LIKELIHOOD_STEPWISE_ADDITION_MODELS = frozenset({"jc69"})


def validate_likelihood_stepwise_addition_model(model_name: str) -> str:
    """Validate the native likelihood model family supported for greedy starts."""
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name not in _SUPPORTED_LIKELIHOOD_STEPWISE_ADDITION_MODELS:
        raise ValueError(
            "likelihood stepwise addition model must be one of "
            + ", ".join(sorted(_SUPPORTED_LIKELIHOOD_STEPWISE_ADDITION_MODELS))
        )
    return normalized_model_name


def build_likelihood_stepwise_addition_tree(
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str = "jc69",
    insertion_order: list[str] | None = None,
    default_branch_length: float = 0.1,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> tuple[PhyloTree, StepwiseAdditionTreeReport]:
    """Build one rooted starting tree by maximizing optimized likelihood per insertion."""
    normalized_model_name = validate_likelihood_stepwise_addition_model(model_name)
    if default_branch_length <= 0.0:
        raise ValueError("default_branch_length must be strictly positive")
    if lower_branch_length_bound < 0.0:
        raise ValueError("lower_branch_length_bound must be nonnegative")
    if upper_branch_length_bound <= lower_branch_length_bound:
        raise ValueError("upper_branch_length_bound must be greater than lower bound")
    if max_coordinate_passes < 1:
        raise ValueError("max_coordinate_passes must be at least one")

    resolved_records, _resolved_alignment_path = resolve_nucleotide_topology_search_records(
        records
    )
    ordered_taxa = _resolve_likelihood_stepwise_insertion_order(
        resolved_records,
        insertion_order=insertion_order,
    )
    current_tree = _build_initial_likelihood_stepwise_tree(
        ordered_taxa[:2],
        default_branch_length=default_branch_length,
    )
    initial_records = _restrict_alignment_records_to_taxa(
        resolved_records,
        taxa=set(ordered_taxa[:2]),
    )
    normalized_records, compressed_patterns = normalize_nucleotide_topology_search_records(
        initial_records,
        owner_name="likelihood stepwise addition",
    )
    resolved_surface = resolve_nucleotide_topology_search_surface(
        current_tree,
        normalized_records,
        model_name=normalized_model_name,
    )
    current_tree = _initialized_likelihood_stepwise_tree(
        current_tree,
        default_branch_length=default_branch_length,
        branch_length_floor=lower_branch_length_bound,
    )
    initial_optimization = optimize_selected_nucleotide_branch_lengths(
        current_tree,
        compressed_patterns,
        specification=resolved_surface.specification,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    current_tree = initial_optimization.optimized_tree
    current_score = initial_optimization.log_likelihood
    trace_rows: list[StepwiseAdditionTraceRow] = []

    for step_index, taxon in enumerate(ordered_taxa[2:], start=1):
        step_records = _restrict_alignment_records_to_taxa(
            resolved_records,
            taxa=set(ordered_taxa[: step_index + 2]),
        )
        _normalized_step_records, step_patterns = normalize_nucleotide_topology_search_records(
            step_records,
            owner_name="likelihood stepwise addition",
        )
        tested_edge_rows: list[StepwiseAdditionCandidateScore] = []
        best_branch_id: str | None = None
        best_descendant_taxa: list[str] | None = None
        best_score: float | None = None
        best_tree: PhyloTree | None = None
        best_tree_newick: str | None = None
        for candidate in iter_stepwise_addition_edge_candidates(current_tree):
            candidate_tree = apply_stepwise_addition_candidate(current_tree, candidate, taxon)
            candidate_tree = _initialized_likelihood_stepwise_tree(
                candidate_tree,
                default_branch_length=default_branch_length,
                branch_length_floor=lower_branch_length_bound,
            )
            optimization_result = optimize_selected_nucleotide_branch_lengths(
                candidate_tree,
                step_patterns,
                specification=resolved_surface.specification,
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
                improvement_tolerance=improvement_tolerance,
                max_coordinate_passes=max_coordinate_passes,
            )
            candidate_score = optimization_result.log_likelihood
            candidate_tree_newick = optimization_result.optimized_tree.to_newick()
            tested_edge_rows.append(
                StepwiseAdditionCandidateScore(
                    branch_id=candidate.branch_id,
                    descendant_taxa=list(candidate.descendant_taxa),
                    score=candidate_score,
                    candidate_tree_newick=candidate_tree_newick,
                )
            )
            if (
                best_score is None
                or candidate_score > best_score
                or (
                    candidate_score == best_score
                    and candidate_tree_newick < (best_tree_newick or "")
                )
            ):
                best_branch_id = candidate.branch_id
                best_descendant_taxa = list(candidate.descendant_taxa)
                best_score = candidate_score
                best_tree = optimization_result.optimized_tree
                best_tree_newick = candidate_tree_newick
        if (
            best_branch_id is None
            or best_descendant_taxa is None
            or best_score is None
            or best_tree is None
            or best_tree_newick is None
        ):
            raise AssertionError(
                "likelihood stepwise addition must evaluate at least one insertion edge"
            )
        current_tree = best_tree
        current_score = best_score
        trace_rows.append(
            StepwiseAdditionTraceRow(
                step_index=step_index,
                taxon=taxon,
                inserted_taxa=ordered_taxa[: step_index + 2],
                tested_edge_rows=tested_edge_rows,
                best_edge_id=best_branch_id,
                best_edge_descendant_taxa=best_descendant_taxa,
                best_score=best_score,
                selected_tree_newick=best_tree_newick,
            )
        )

    report = summarize_stepwise_addition_tree(
        current_tree,
        insertion_order=ordered_taxa,
        objective_name=f"likelihood-{normalized_model_name}",
        objective_direction="maximize",
        final_score=current_score,
        trace_rows=trace_rows,
    )
    return current_tree, report


def build_likelihood_stepwise_addition_tree_from_alignment(
    alignment_path: Path,
    *,
    model_name: str = "jc69",
    insertion_order: list[str] | None = None,
    default_branch_length: float = 0.1,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> tuple[PhyloTree, StepwiseAdditionTreeReport]:
    """Build one likelihood stepwise-addition tree from one FASTA alignment path."""
    return build_likelihood_stepwise_addition_tree(
        alignment_path,
        model_name=model_name,
        insertion_order=insertion_order,
        default_branch_length=default_branch_length,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def _build_initial_likelihood_stepwise_tree(
    starting_taxa: list[str],
    *,
    default_branch_length: float,
) -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(name=starting_taxa[0], branch_length=default_branch_length),
                TreeNode(name=starting_taxa[1], branch_length=default_branch_length),
            ]
        ),
        rooted=True,
    ).refresh()


def _initialized_likelihood_stepwise_tree(
    tree: PhyloTree,
    *,
    default_branch_length: float,
    branch_length_floor: float,
) -> PhyloTree:
    working_tree = tree.copy().refresh()
    for _parent, child in working_tree.iter_edges():
        branch_length = child.branch_length
        if branch_length is None:
            child.branch_length = max(default_branch_length, branch_length_floor)
            continue
        child.branch_length = max(float(branch_length), branch_length_floor)
    return working_tree.refresh()


def _resolve_likelihood_stepwise_insertion_order(
    records: list[AlignmentRecord],
    *,
    insertion_order: list[str] | None,
) -> list[str]:
    record_order = validate_stepwise_addition_taxa([record.identifier for record in records])
    if insertion_order is None:
        return record_order
    validated_insertion_order = validate_stepwise_addition_taxa(insertion_order)
    record_taxa = set(record_order)
    insertion_taxa = set(validated_insertion_order)
    missing_taxa = sorted(record_taxa - insertion_taxa)
    unexpected_taxa = sorted(insertion_taxa - record_taxa)
    if missing_taxa or unexpected_taxa:
        raise ValueError(
            "likelihood stepwise addition insertion_order must match alignment taxa exactly"
        )
    return list(validated_insertion_order)


def _restrict_alignment_records_to_taxa(
    records: list[AlignmentRecord],
    *,
    taxa: set[str],
) -> list[AlignmentRecord]:
    return [record for record in records if record.identifier in taxa]
