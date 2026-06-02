from __future__ import annotations

import csv
from dataclasses import asdict
import json
import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    loads_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import normalize_unambiguous_dna_records
from bijux_phylogenetics.phylo.likelihood.models import (
    LikelihoodPlacementAlternativeRow,
    LikelihoodPlacementQuerySummary,
    LikelihoodPlacementReport,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    BoundedLikelihoodSearchResult,
    run_bounded_coordinate_likelihood_search,
    run_bounded_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.site_log_likelihoods import (
    evaluate_selected_dna_site_log_likelihoods_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
    InvalidBranchLengthError,
)

_PLACEMENT_LIKELIHOOD_EQUALITY_TOLERANCE = 1e-12
# Stable report field name; Bandit misclassifies the word "pass" here.
_OPTIMIZATION_PASS_COUNT_FIELD = "optimization_pass_count"  # nosec B105


class _PlacementEvaluationError(RuntimeError):
    """Internal sentinel used to keep placement search payloads structured."""


def place_queries_by_likelihood(
    reference_tree: PhyloTree,
    reference_records: list[AlignmentRecord],
    query_records: list[AlignmentRecord],
    *,
    model: str = "jc69",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    lower_pendant_length_bound: float = 1e-6,
    upper_pendant_length_bound: float = 5.0,
    max_coordinate_passes: int = 12,
) -> LikelihoodPlacementReport:
    """Place each query sequence on every reference-tree edge by nucleotide likelihood."""
    normalized_model = model.strip().lower()
    if normalized_model != "jc69":
        raise ValueError("likelihood placement currently supports only the JC69 model")
    if lower_pendant_length_bound <= 0.0:
        raise InvalidBranchLengthError(
            "likelihood placement lower pendant length bound must be positive"
        )
    if upper_pendant_length_bound <= lower_pendant_length_bound:
        raise InvalidBranchLengthError(
            "likelihood placement pendant length bounds must be strictly increasing"
        )
    if max_coordinate_passes < 1:
        raise ValueError(
            "likelihood placement max_coordinate_passes must be at least one"
        )

    working_tree = reference_tree.copy().refresh()
    validate_explicit_branch_lengths(working_tree, model_name="JC69 placement")
    normalized_reference_records = normalize_unambiguous_dna_records(
        list(reference_records),
        model_name="JC69 placement",
    )
    normalized_query_records = normalize_unambiguous_dna_records(
        list(query_records),
        model_name="JC69 placement",
    )
    _validate_query_alignment_inputs(
        working_tree,
        normalized_reference_records,
        normalized_query_records,
    )

    reference_patterns = compress_alignment_site_patterns_from_records(
        normalized_reference_records
    )
    validate_tree_taxa_against_patterns(
        working_tree,
        reference_patterns,
        model_name="JC69 placement reference tree",
    )
    reference_taxa = reference_patterns.taxon_order
    reference_tree_newick = dumps_newick(working_tree)
    edge_rows = list(_iter_reference_edges(working_tree))
    alternative_rows: list[LikelihoodPlacementAlternativeRow] = []
    query_summaries: list[LikelihoodPlacementQuerySummary] = []
    total_function_evaluation_count = 0

    for query_record in normalized_query_records:
        placement_rows, summary_row = _place_single_query(
            working_tree,
            normalized_reference_records,
            query_record,
            edge_rows=edge_rows,
            model_name=normalized_model,
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            lower_pendant_length_bound=lower_pendant_length_bound,
            upper_pendant_length_bound=upper_pendant_length_bound,
            max_coordinate_passes=max_coordinate_passes,
        )
        alternative_rows.extend(placement_rows)
        query_summaries.append(summary_row)
        total_function_evaluation_count += sum(
            row.function_evaluation_count for row in placement_rows
        )

    return LikelihoodPlacementReport(
        model_name="JC69",
        reference_tree_newick=reference_tree_newick,
        reference_taxa=reference_taxa,
        edge_count=len(edge_rows),
        query_count=len(normalized_query_records),
        site_count=reference_patterns.alignment_length,
        lower_pendant_length_bound=lower_pendant_length_bound,
        upper_pendant_length_bound=upper_pendant_length_bound,
        max_coordinate_passes=max_coordinate_passes,
        total_function_evaluation_count=total_function_evaluation_count,
        query_summaries=query_summaries,
        alternative_placements=alternative_rows,
    )


def place_queries_by_likelihood_from_alignment(
    reference_tree_path: Path,
    reference_alignment_path: Path,
    query_alignment_path: Path,
    *,
    model: str = "jc69",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    lower_pendant_length_bound: float = 1e-6,
    upper_pendant_length_bound: float = 5.0,
    max_coordinate_passes: int = 12,
) -> LikelihoodPlacementReport:
    """Place query sequences from FASTA files on one reference tree by likelihood."""
    return place_queries_by_likelihood(
        load_tree(reference_tree_path),
        load_fasta_alignment(reference_alignment_path),
        load_fasta_alignment(query_alignment_path),
        model=model,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        lower_pendant_length_bound=lower_pendant_length_bound,
        upper_pendant_length_bound=upper_pendant_length_bound,
        max_coordinate_passes=max_coordinate_passes,
    )


def _place_single_query(
    reference_tree: PhyloTree,
    reference_records: list[AlignmentRecord],
    query_record: AlignmentRecord,
    *,
    edge_rows: list[tuple[str, str | None, list[str], float]],
    model_name: str,
    root_prior_policy: str | None,
    root_prior: dict[str, float] | list[float] | tuple[float, ...] | None,
    fixed_root_state: str | None,
    lower_pendant_length_bound: float,
    upper_pendant_length_bound: float,
    max_coordinate_passes: int,
) -> tuple[list[LikelihoodPlacementAlternativeRow], LikelihoodPlacementQuerySummary]:
    combined_records = [*reference_records, query_record]
    compressed_patterns = compress_alignment_site_patterns_from_records(
        combined_records
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        combined_records,
        model_name=model_name,
        owner_name="JC69 likelihood placement",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )

    unresolved_rows: list[tuple[LikelihoodPlacementAlternativeRow, float]] = []
    for edge_id, child_name, descendant_taxa, original_branch_length in edge_rows:
        optimized = _optimize_edge_placement(
            reference_tree,
            compressed_patterns,
            query_id=query_record.identifier,
            edge_id=edge_id,
            original_branch_length=original_branch_length,
            model_name=specification.model_name,
            state_count=specification.state_count,
            observation_policy=specification.observation_policy,
            root_prior=specification.root_prior,
            transition_matrix_for_branch_length=(
                specification.transition_matrix_for_branch_length
            ),
            lower_pendant_length_bound=lower_pendant_length_bound,
            upper_pendant_length_bound=upper_pendant_length_bound,
            max_coordinate_passes=max_coordinate_passes,
        )
        unresolved_rows.append(
            (
                LikelihoodPlacementAlternativeRow(
                    query_id=query_record.identifier,
                    placement_rank=0,
                    edge_id=edge_id,
                    child_name=child_name,
                    descendant_taxa=descendant_taxa,
                    original_branch_length=original_branch_length,
                    optimized_proximal_length=optimized["proximal_length"],
                    optimized_distal_length=optimized["distal_length"],
                    optimized_pendant_length=optimized["pendant_length"],
                    log_likelihood=optimized["log_likelihood"],
                    likelihood_weight_ratio=0.0,
                    function_evaluation_count=optimized["function_evaluation_count"],
                    optimization_pass_count=optimized["optimization_pass_count"],
                    converged=optimized["converged"],
                    placed_tree_newick=optimized["placed_tree_newick"],
                ),
                optimized["log_likelihood"],
            )
        )

    ordered_rows = _rank_placement_rows(unresolved_rows)
    best_row = ordered_rows[0]
    equally_best_count = sum(
        1
        for row in ordered_rows
        if math.isclose(
            row.log_likelihood,
            best_row.log_likelihood,
            rel_tol=0.0,
            abs_tol=_PLACEMENT_LIKELIHOOD_EQUALITY_TOLERANCE,
        )
    )
    summary_row = LikelihoodPlacementQuerySummary(
        query_id=query_record.identifier,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        best_edge_id=best_row.edge_id,
        best_child_name=best_row.child_name,
        best_descendant_taxa=list(best_row.descendant_taxa),
        best_original_branch_length=best_row.original_branch_length,
        best_proximal_length=best_row.optimized_proximal_length,
        best_distal_length=best_row.optimized_distal_length,
        best_pendant_length=best_row.optimized_pendant_length,
        best_log_likelihood=best_row.log_likelihood,
        best_likelihood_weight_ratio=best_row.likelihood_weight_ratio,
        candidate_placement_count=len(ordered_rows),
        equally_best_placement_count=equally_best_count,
        best_tree_newick=best_row.placed_tree_newick,
    )
    return ordered_rows, summary_row


def _rank_placement_rows(
    unresolved_rows: list[tuple[LikelihoodPlacementAlternativeRow, float]],
) -> list[LikelihoodPlacementAlternativeRow]:
    max_log_likelihood = max(log_likelihood for _row, log_likelihood in unresolved_rows)
    scaled_weights = [
        math.exp(log_likelihood - max_log_likelihood)
        for _row, log_likelihood in unresolved_rows
    ]
    normalization = sum(scaled_weights)
    ordered = sorted(
        zip(unresolved_rows, scaled_weights, strict=True),
        key=lambda item: (-item[0][1], item[0][0].edge_id),
    )
    ranked_rows: list[LikelihoodPlacementAlternativeRow] = []
    for placement_rank, ((row, _log_likelihood), scaled_weight) in enumerate(
        ordered,
        start=1,
    ):
        ranked_rows.append(
            LikelihoodPlacementAlternativeRow(
                query_id=row.query_id,
                placement_rank=placement_rank,
                edge_id=row.edge_id,
                child_name=row.child_name,
                descendant_taxa=list(row.descendant_taxa),
                original_branch_length=row.original_branch_length,
                optimized_proximal_length=row.optimized_proximal_length,
                optimized_distal_length=row.optimized_distal_length,
                optimized_pendant_length=row.optimized_pendant_length,
                log_likelihood=row.log_likelihood,
                likelihood_weight_ratio=scaled_weight / normalization,
                function_evaluation_count=row.function_evaluation_count,
                optimization_pass_count=row.optimization_pass_count,
                converged=row.converged,
                placed_tree_newick=row.placed_tree_newick,
            )
        )
    return ranked_rows


def _optimize_edge_placement(
    reference_tree: PhyloTree,
    compressed_patterns,
    *,
    query_id: str,
    edge_id: str,
    original_branch_length: float,
    model_name: str,
    state_count: int,
    observation_policy: str,
    root_prior,
    transition_matrix_for_branch_length,
    lower_pendant_length_bound: float,
    upper_pendant_length_bound: float,
    max_coordinate_passes: int,
) -> dict[str, float | int | bool | str]:
    initial_pendant_length = min(
        max(
            lower_pendant_length_bound,
            max(original_branch_length, lower_pendant_length_bound),
        ),
        upper_pendant_length_bound,
    )
    total_function_evaluation_count = 0

    def evaluate_candidate(
        distal_length: float,
        pendant_length: float,
    ) -> tuple[str, float]:
        placed_tree = _build_placed_tree(
            reference_tree,
            query_id=query_id,
            edge_id=edge_id,
            distal_length=distal_length,
            pendant_length=pendant_length,
        )
        report = evaluate_selected_dna_site_log_likelihoods_from_patterns(
            placed_tree,
            compressed_patterns,
            model_name=model_name,
            state_count=state_count,
            observation_policy=observation_policy,
            root_prior=root_prior,
            parameter_values={},
            transition_matrix_for_child=lambda child: (
                transition_matrix_for_branch_length(
                    max(float(child.branch_length or 0.0), 0.0)
                )
            ),
        )
        return report.tree_newick, report.log_likelihood

    if original_branch_length <= 0.0:
        pendant_search = _optimize_fixed_distal_placement(
            evaluate_candidate,
            distal_length=0.0,
            lower_pendant_length_bound=lower_pendant_length_bound,
            upper_pendant_length_bound=upper_pendant_length_bound,
        )
        total_function_evaluation_count += pendant_search.function_evaluation_count
        return {
            "proximal_length": 0.0,
            "distal_length": 0.0,
            "pendant_length": float(pendant_search.parameter_value),
            "log_likelihood": float(pendant_search.objective_value),
            "function_evaluation_count": total_function_evaluation_count,
            "converged": pendant_search.converged,
            "placed_tree_newick": pendant_search.payload,
            _OPTIMIZATION_PASS_COUNT_FIELD: 1,
        }

    def evaluate_coordinate_values(
        candidate_values: dict[str, float],
    ) -> tuple[str, float]:
        return evaluate_candidate(
            float(candidate_values["distal_length"]),
            float(candidate_values["pendant_length"]),
        )

    interior_search = run_bounded_coordinate_likelihood_search(
        initial_values={
            "distal_length": original_branch_length / 2.0,
            "pendant_length": initial_pendant_length,
        },
        bounds_by_name={
            "distal_length": (0.0, original_branch_length),
            "pendant_length": (
                lower_pendant_length_bound,
                upper_pendant_length_bound,
            ),
        },
        evaluate=evaluate_coordinate_values,
        max_coordinate_passes=max_coordinate_passes,
    )
    total_function_evaluation_count += interior_search.function_evaluation_count
    best_candidate = {
        "distal_length": float(interior_search.parameter_values["distal_length"]),
        "pendant_length": float(interior_search.parameter_values["pendant_length"]),
        "log_likelihood": float(interior_search.objective_value),
        "optimization_pass_count": interior_search.optimization_pass_count,
        "converged": interior_search.converged,
        "placed_tree_newick": interior_search.payload,
    }

    for boundary_distal_length in (0.0, original_branch_length):
        boundary_search = _optimize_fixed_distal_placement(
            evaluate_candidate,
            distal_length=boundary_distal_length,
            lower_pendant_length_bound=lower_pendant_length_bound,
            upper_pendant_length_bound=upper_pendant_length_bound,
        )
        total_function_evaluation_count += boundary_search.function_evaluation_count
        if boundary_search.objective_value > best_candidate["log_likelihood"]:
            best_candidate = {
                "distal_length": boundary_distal_length,
                "pendant_length": float(boundary_search.parameter_value),
                "log_likelihood": float(boundary_search.objective_value),
                "converged": boundary_search.converged,
                "placed_tree_newick": boundary_search.payload,
            }
            best_candidate[_OPTIMIZATION_PASS_COUNT_FIELD] = 1

    distal_length = float(best_candidate["distal_length"])
    proximal_length = max(original_branch_length - distal_length, 0.0)
    return {
        "proximal_length": proximal_length,
        "distal_length": distal_length,
        "pendant_length": float(best_candidate["pendant_length"]),
        "log_likelihood": float(best_candidate["log_likelihood"]),
        "function_evaluation_count": total_function_evaluation_count,
        "optimization_pass_count": int(best_candidate[_OPTIMIZATION_PASS_COUNT_FIELD]),
        "converged": bool(best_candidate["converged"]),
        "placed_tree_newick": str(best_candidate["placed_tree_newick"]),
    }


def _optimize_fixed_distal_placement(
    evaluate_candidate,
    *,
    distal_length: float,
    lower_pendant_length_bound: float,
    upper_pendant_length_bound: float,
) -> BoundedLikelihoodSearchResult[str]:
    return run_bounded_likelihood_search(
        lower_bound=lower_pendant_length_bound,
        upper_bound=upper_pendant_length_bound,
        evaluate=lambda pendant_length: evaluate_candidate(
            distal_length,
            pendant_length,
        ),
    )


def _build_placed_tree(
    reference_tree: PhyloTree,
    *,
    query_id: str,
    edge_id: str,
    distal_length: float,
    pendant_length: float,
) -> PhyloTree:
    placed_tree = reference_tree.copy().refresh()
    child = placed_tree.node_by_id(edge_id)
    parent = child.parent
    if parent is None:
        raise _PlacementEvaluationError("placement edge cannot target the tree root")
    original_branch_length = float(child.branch_length or 0.0)
    if distal_length < 0.0 or distal_length > original_branch_length:
        raise _PlacementEvaluationError("distal length must lie on the selected edge")
    proximal_length = max(original_branch_length - distal_length, 0.0)
    child.branch_length = distal_length
    placement_node = TreeNode(
        branch_length=proximal_length,
        children=[
            child,
            TreeNode(name=query_id, branch_length=pendant_length),
        ],
    )
    replacement_children = list(parent.children)
    for child_index, observed_child in enumerate(replacement_children):
        if observed_child.node_id == edge_id:
            replacement_children[child_index] = placement_node
            break
    parent.replace_children(replacement_children)
    return placed_tree.refresh()


def _iter_reference_edges(
    tree: PhyloTree,
) -> list[tuple[str, str | None, list[str], float]]:
    rows: list[tuple[str, str | None, list[str], float]] = []
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise InvalidBranchLengthError(
                "likelihood placement requires stable edge node_id values"
            )
        rows.append(
            (
                child.node_id,
                child.name,
                child.descendant_taxa,
                float(child.branch_length or 0.0),
            )
        )
    return rows


def _validate_query_alignment_inputs(
    reference_tree: PhyloTree,
    reference_records: list[AlignmentRecord],
    query_records: list[AlignmentRecord],
) -> None:
    if not reference_records:
        raise InvalidAlignmentError(
            "likelihood placement requires at least one reference alignment record"
        )
    if not query_records:
        raise InvalidAlignmentError(
            "likelihood placement requires at least one query alignment record"
        )
    reference_length = len(reference_records[0].sequence)
    query_length = len(query_records[0].sequence)
    if query_length != reference_length:
        raise InvalidAlignmentError(
            "likelihood placement requires query alignments with the same aligned length as the reference alignment"
        )

    reference_ids = [record.identifier for record in reference_records]
    query_ids = [record.identifier for record in query_records]
    if len(set(query_ids)) != len(query_ids):
        raise InvalidAlignmentError(
            "likelihood placement requires uniquely named query alignment records"
        )
    overlapping_ids = sorted(set(reference_ids) & set(query_ids))
    if overlapping_ids:
        raise AlignmentTaxonMismatchError(
            "likelihood placement requires query identifiers that do not overlap the reference tree taxa"
            f" (overlap: {', '.join(overlapping_ids)})"
        )
    tree_taxa = [leaf.name for leaf in reference_tree.iter_leaves()]
    if any(name is None for name in tree_taxa):
        raise AlignmentTaxonMismatchError(
            "likelihood placement requires named reference-tree tips"
        )
    if set(reference_ids) != {name for name in tree_taxa if name is not None}:
        raise AlignmentTaxonMismatchError(
            "likelihood placement requires the reference alignment to match the reference tree taxa exactly"
        )


def write_likelihood_placement_summary_table(
    path: Path,
    report: LikelihoodPlacementReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "query_id",
                "site_count",
                "pattern_count",
                "best_edge_id",
                "best_child_name",
                "best_descendant_taxa",
                "best_original_branch_length",
                "best_proximal_length",
                "best_distal_length",
                "best_pendant_length",
                "best_log_likelihood",
                "best_likelihood_weight_ratio",
                "candidate_placement_count",
                "equally_best_placement_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.query_summaries:
            writer.writerow(
                {
                    "query_id": row.query_id,
                    "site_count": row.site_count,
                    "pattern_count": row.pattern_count,
                    "best_edge_id": row.best_edge_id,
                    "best_child_name": row.best_child_name or "",
                    "best_descendant_taxa": "|".join(row.best_descendant_taxa),
                    "best_original_branch_length": format(
                        row.best_original_branch_length,
                        ".15g",
                    ),
                    "best_proximal_length": format(row.best_proximal_length, ".15g"),
                    "best_distal_length": format(row.best_distal_length, ".15g"),
                    "best_pendant_length": format(row.best_pendant_length, ".15g"),
                    "best_log_likelihood": format(row.best_log_likelihood, ".15g"),
                    "best_likelihood_weight_ratio": format(
                        row.best_likelihood_weight_ratio,
                        ".15g",
                    ),
                    "candidate_placement_count": row.candidate_placement_count,
                    "equally_best_placement_count": row.equally_best_placement_count,
                }
            )
    return path


def write_likelihood_placement_alternative_table(
    path: Path,
    report: LikelihoodPlacementReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "query_id",
                "placement_rank",
                "edge_id",
                "child_name",
                "descendant_taxa",
                "original_branch_length",
                "optimized_proximal_length",
                "optimized_distal_length",
                "optimized_pendant_length",
                "log_likelihood",
                "likelihood_weight_ratio",
                "function_evaluation_count",
                "optimization_pass_count",
                "converged",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.alternative_placements:
            writer.writerow(
                {
                    "query_id": row.query_id,
                    "placement_rank": row.placement_rank,
                    "edge_id": row.edge_id,
                    "child_name": row.child_name or "",
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "original_branch_length": format(
                        row.original_branch_length,
                        ".15g",
                    ),
                    "optimized_proximal_length": format(
                        row.optimized_proximal_length,
                        ".15g",
                    ),
                    "optimized_distal_length": format(
                        row.optimized_distal_length,
                        ".15g",
                    ),
                    "optimized_pendant_length": format(
                        row.optimized_pendant_length,
                        ".15g",
                    ),
                    "log_likelihood": format(row.log_likelihood, ".15g"),
                    "likelihood_weight_ratio": format(
                        row.likelihood_weight_ratio,
                        ".15g",
                    ),
                    "function_evaluation_count": row.function_evaluation_count,
                    "optimization_pass_count": row.optimization_pass_count,
                    "converged": str(row.converged).lower(),
                }
            )
    return path


def write_likelihood_placement_tree_set(
    path: Path,
    report: LikelihoodPlacementReport,
) -> Path:
    return write_newick_tree_set(
        path,
        [loads_newick(row.best_tree_newick) for row in report.query_summaries],
    )


def write_likelihood_placement_run_json(
    path: Path,
    report: LikelihoodPlacementReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_likelihood_placement_artifacts(
    out_dir: Path,
    report: LikelihoodPlacementReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    best_tree_path = write_likelihood_placement_tree_set(
        out_dir / "best_placements.nwk",
        report,
    )
    summary_path = write_likelihood_placement_summary_table(
        out_dir / "summary.tsv",
        report,
    )
    alternative_path = write_likelihood_placement_alternative_table(
        out_dir / "alternative_placements.tsv",
        report,
    )
    run_json_path = write_likelihood_placement_run_json(out_dir / "run.json", report)
    return {
        "best_tree_path": best_tree_path,
        "summary_path": summary_path,
        "alternative_path": alternative_path,
        "run_json_path": run_json_path,
    }
