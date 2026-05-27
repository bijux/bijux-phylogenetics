from __future__ import annotations

import json
from pathlib import Path

import numpy

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodNniTraceRow,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    BranchReoptimizationResult,
    normalize_nucleotide_topology_search_records,
    optimize_selected_nucleotide_branch_lengths,
    prefer_higher_likelihood,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
    resolve_nucleotide_topology_search_tree,
    validate_branch_reoptimization_policy,
    validate_nucleotide_topology_search_tree,
)
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES = frozenset({"coordinate-branch-lengths"})


def search_nucleotide_likelihood_nni(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
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
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodNniSearchReport:
    """Search one rooted binary nucleotide tree by likelihood-improving rooted NNI moves."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = resolve_nucleotide_topology_search_records(
        records
    )
    validated_branch_reoptimization_policy = validate_branch_reoptimization_policy(
        branch_reoptimization_policy
    )
    validate_nucleotide_topology_search_tree(
        resolved_tree,
        workflow_name="nucleotide likelihood NNI search",
    )
    normalized_records, compressed_patterns = normalize_nucleotide_topology_search_records(
        resolved_records,
        owner_name="nucleotide likelihood NNI search",
    )
    resolved_surface = resolve_nucleotide_topology_search_surface(
        resolved_tree,
        normalized_records,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
    )
    input_tree_newick = dumps_newick(resolved_tree)
    start_result = _reoptimize_tree(
        resolved_tree,
        compressed_patterns=compressed_patterns,
        resolved_surface=resolved_surface,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    current_tree = start_result.optimized_tree.copy().refresh()
    current_log_likelihood = start_result.log_likelihood
    start_tree_newick = dumps_newick(current_tree)
    trace_rows = [
        NucleotideLikelihoodNniTraceRow(
            event_index=1,
            event_kind="start",
            iteration=0,
            log_likelihood_before=None,
            log_likelihood_after=current_log_likelihood,
            log_likelihood_delta=None,
            tree_before_newick=None,
            tree_after_newick=start_tree_newick,
            pivot_branch_id=None,
            sibling_clade_id=None,
            exchanged_clade_id=None,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            branch_optimization_pass_count=start_result.optimization_pass_count,
            branch_function_evaluation_count=start_result.function_evaluation_count,
            stopping_reason=None,
        )
    ]
    accepted_move_count = 0
    evaluated_neighbor_count = 0
    total_branch_optimization_pass_count = start_result.optimization_pass_count
    total_branch_function_evaluation_count = start_result.function_evaluation_count
    while True:
        improving_candidate = None
        improving_result: BranchReoptimizationResult | None = None
        improving_newick: str | None = None
        for candidate in iter_rooted_nni_move_candidates(current_tree):
            neighbor_tree = apply_rooted_nni_move(current_tree, candidate)
            neighbor_result = _reoptimize_tree(
                neighbor_tree,
                compressed_patterns=compressed_patterns,
                resolved_surface=resolved_surface,
                branch_reoptimization_policy=validated_branch_reoptimization_policy,
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
                improvement_tolerance=improvement_tolerance,
                max_coordinate_passes=max_coordinate_passes,
            )
            evaluated_neighbor_count += 1
            neighbor_newick = dumps_newick(neighbor_result.optimized_tree)
            total_branch_optimization_pass_count += (
                neighbor_result.optimization_pass_count
            )
            total_branch_function_evaluation_count += (
                neighbor_result.function_evaluation_count
            )
            if neighbor_result.log_likelihood <= current_log_likelihood:
                continue
            if (
                improving_result is None
                or improving_newick is None
                or prefer_higher_likelihood(
                    neighbor_result.log_likelihood,
                    neighbor_newick,
                    improving_result.log_likelihood,
                    improving_newick,
                )
            ):
                improving_candidate = candidate
                improving_result = neighbor_result
                improving_newick = neighbor_newick
        if improving_candidate is None or improving_result is None or improving_newick is None:
            stopping_reason = "no-improving-neighbor"
            break
        accepted_move_count += 1
        score_before = current_log_likelihood
        tree_before_newick = dumps_newick(current_tree)
        current_tree = improving_result.optimized_tree.copy().refresh()
        current_log_likelihood = improving_result.log_likelihood
        trace_rows.append(
            NucleotideLikelihoodNniTraceRow(
                event_index=len(trace_rows) + 1,
                event_kind="accepted-move",
                iteration=accepted_move_count,
                log_likelihood_before=score_before,
                log_likelihood_after=current_log_likelihood,
                log_likelihood_delta=current_log_likelihood - score_before,
                tree_before_newick=tree_before_newick,
                tree_after_newick=improving_newick,
                pivot_branch_id=improving_candidate.pivot_branch_id,
                sibling_clade_id=improving_candidate.sibling_clade_id,
                exchanged_clade_id=improving_candidate.exchanged_clade_id,
                branch_reoptimization_policy=validated_branch_reoptimization_policy,
                branch_optimization_pass_count=improving_result.optimization_pass_count,
                branch_function_evaluation_count=improving_result.function_evaluation_count,
                stopping_reason=None,
            )
        )
    trace_rows.append(
        NucleotideLikelihoodNniTraceRow(
            event_index=len(trace_rows) + 1,
            event_kind="final",
            iteration=accepted_move_count,
            log_likelihood_before=None,
            log_likelihood_after=current_log_likelihood,
            log_likelihood_delta=None,
            tree_before_newick=None,
            tree_after_newick=dumps_newick(current_tree),
            pivot_branch_id=None,
            sibling_clade_id=None,
            exchanged_clade_id=None,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            branch_optimization_pass_count=0,
            branch_function_evaluation_count=0,
            stopping_reason=stopping_reason,
        )
    )
    return NucleotideLikelihoodNniSearchReport(
        algorithm="nucleotide-likelihood-nni-search",
        model_name=resolved_surface.model_name,
        tree_path=None if resolved_tree_path is None else str(resolved_tree_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        input_tree_newick=input_tree_newick,
        start_tree_newick=start_tree_newick,
        start_log_likelihood=start_result.log_likelihood,
        final_tree_newick=dumps_newick(current_tree),
        final_log_likelihood=current_log_likelihood,
        accepted_move_count=accepted_move_count,
        evaluated_neighbor_count=evaluated_neighbor_count,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        substitution_parameter_policy=resolved_surface.substitution_parameter_policy,
        substitution_parameter_values=resolved_surface.substitution_parameter_values,
        substitution_parameter_warnings=resolved_surface.substitution_parameter_warnings,
        total_branch_optimization_pass_count=total_branch_optimization_pass_count,
        total_branch_function_evaluation_count=total_branch_function_evaluation_count,
        stopping_reason=stopping_reason,
        trace_rows=trace_rows,
    )


def search_nucleotide_likelihood_nni_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
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
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodNniSearchReport:
    """Search one rooted binary nucleotide tree path by likelihood-improving rooted NNI moves."""
    return search_nucleotide_likelihood_nni(
        tree_path,
        alignment_path,
        model_name=model_name,
        branch_reoptimization_policy=branch_reoptimization_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def write_nucleotide_likelihood_nni_trace_table(
    path: Path,
    report: NucleotideLikelihoodNniSearchReport,
) -> Path:
    """Write one deterministic rooted likelihood NNI trace table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "event_index",
        "event_kind",
        "iteration",
        "log_likelihood_before",
        "log_likelihood_after",
        "log_likelihood_delta",
        "tree_before_newick",
        "tree_after_newick",
        "pivot_branch_id",
        "sibling_clade_id",
        "exchanged_clade_id",
        "branch_reoptimization_policy",
        "branch_optimization_pass_count",
        "branch_function_evaluation_count",
        "stopping_reason",
    ]
    rows = ["\t".join(columns)]
    for row in report.trace_rows:
        payload = [
            row.event_index,
            row.event_kind,
            row.iteration,
            row.log_likelihood_before,
            row.log_likelihood_after,
            row.log_likelihood_delta,
            row.tree_before_newick,
            row.tree_after_newick,
            row.pivot_branch_id,
            row.sibling_clade_id,
            row.exchanged_clade_id,
            row.branch_reoptimization_policy,
            row.branch_optimization_pass_count,
            row.branch_function_evaluation_count,
            row.stopping_reason,
        ]
        rows.append(
            "\t".join(
                ""
                if value is None
                else format(value, ".15g")
                if isinstance(value, float)
                else str(value)
                for value in payload
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_nni_run_json(
    path: Path,
    report: NucleotideLikelihoodNniSearchReport,
) -> Path:
    """Write one machine-readable rooted likelihood NNI search payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "model_name": report.model_name,
        "tree_path": report.tree_path,
        "alignment_path": report.alignment_path,
        "taxon_count": report.taxon_count,
        "site_count": report.site_count,
        "pattern_count": report.pattern_count,
        "input_tree_newick": report.input_tree_newick,
        "start_tree_newick": report.start_tree_newick,
        "start_log_likelihood": report.start_log_likelihood,
        "final_tree_newick": report.final_tree_newick,
        "final_log_likelihood": report.final_log_likelihood,
        "accepted_move_count": report.accepted_move_count,
        "evaluated_neighbor_count": report.evaluated_neighbor_count,
        "branch_reoptimization_policy": report.branch_reoptimization_policy,
        "substitution_parameter_policy": report.substitution_parameter_policy,
        "substitution_parameter_values": report.substitution_parameter_values,
        "substitution_parameter_warnings": report.substitution_parameter_warnings,
        "total_branch_optimization_pass_count": report.total_branch_optimization_pass_count,
        "total_branch_function_evaluation_count": report.total_branch_function_evaluation_count,
        "stopping_reason": report.stopping_reason,
        "trace_rows": [
            {
                "event_index": row.event_index,
                "event_kind": row.event_kind,
                "iteration": row.iteration,
                "log_likelihood_before": row.log_likelihood_before,
                "log_likelihood_after": row.log_likelihood_after,
                "log_likelihood_delta": row.log_likelihood_delta,
                "tree_before_newick": row.tree_before_newick,
                "tree_after_newick": row.tree_after_newick,
                "pivot_branch_id": row.pivot_branch_id,
                "sibling_clade_id": row.sibling_clade_id,
                "exchanged_clade_id": row.exchanged_clade_id,
                "branch_reoptimization_policy": row.branch_reoptimization_policy,
                "branch_optimization_pass_count": row.branch_optimization_pass_count,
                "branch_function_evaluation_count": row.branch_function_evaluation_count,
                "stopping_reason": row.stopping_reason,
            }
            for row in report.trace_rows
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_nni_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodNniSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted likelihood NNI search run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    start_tree_path = write_newick(
        out_dir / "start_tree.nwk",
        loads_newick(report.start_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.final_tree_newick),
    )
    trace_path = write_nucleotide_likelihood_nni_trace_table(
        out_dir / "search_trace.tsv",
        report,
    )
    run_json_path = write_nucleotide_likelihood_nni_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "input_tree_path": input_tree_path,
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "trace_path": trace_path,
        "run_json_path": run_json_path,
    }


def _reoptimize_tree(
    tree: PhyloTree,
    *,
    compressed_patterns,
    resolved_surface,
    branch_reoptimization_policy: str,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> BranchReoptimizationResult:
    if branch_reoptimization_policy not in _SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES:
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
