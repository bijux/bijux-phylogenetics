from __future__ import annotations

import json
from pathlib import Path

import numpy

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.equal_best_topologies import (
    build_nucleotide_likelihood_equal_best_tree_report,
    initialize_nucleotide_likelihood_equal_best_topology_accumulator,
    record_nucleotide_likelihood_equal_best_topology,
    validate_nucleotide_likelihood_equal_best_likelihood_tolerance,
    validate_nucleotide_likelihood_equal_best_tree_cap,
)
from bijux_phylogenetics.phylo.likelihood.search_convergence import (
    resolve_nucleotide_likelihood_search_convergence_decision,
    validate_nucleotide_likelihood_search_improvement_tolerance,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodTbrSearchReport,
    NucleotideLikelihoodTbrTraceRow,
)
from bijux_phylogenetics.phylo.likelihood.search_artifacts import (
    write_nucleotide_likelihood_best_tree_set,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    BranchReoptimizationResult,
    initialize_generated_nucleotide_topology_search_tree,
    normalize_nucleotide_topology_search_records,
    prefer_higher_likelihood,
    reoptimize_nucleotide_topology_tree,
    resolve_reoptimized_branch_clade_ids,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
    resolve_nucleotide_topology_search_tree,
    validate_nucleotide_topology_search_tree,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.rooted_tbr import (
    RootedTbrMoveCandidate,
    apply_rooted_tbr_move,
    iter_rooted_tbr_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_TBR_BRANCH_REOPTIMIZATION_POLICIES = frozenset(
    {"coordinate-branch-lengths"}
)


def search_nucleotide_likelihood_tbr(
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
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    search_improvement_tolerance: float = 0.0,
    equal_best_likelihood_tolerance: float = 0.0,
    equal_best_tree_cap: int = 1,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodTbrSearchReport:
    """Search one rooted binary nucleotide tree by likelihood-improving rooted TBR moves."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = resolve_nucleotide_topology_search_records(
        records
    )
    validated_branch_reoptimization_policy = (
        validate_nucleotide_likelihood_tbr_branch_reoptimization_policy(
            branch_reoptimization_policy
        )
    )
    validated_search_improvement_tolerance = (
        validate_nucleotide_likelihood_search_improvement_tolerance(
            search_improvement_tolerance
        )
    )
    validated_equal_best_likelihood_tolerance = (
        validate_nucleotide_likelihood_equal_best_likelihood_tolerance(
            equal_best_likelihood_tolerance
        )
    )
    validated_equal_best_tree_cap = validate_nucleotide_likelihood_equal_best_tree_cap(
        equal_best_tree_cap
    )
    validate_nucleotide_topology_search_tree(
        resolved_tree,
        workflow_name="nucleotide likelihood TBR search",
    )
    normalized_records, compressed_patterns = normalize_nucleotide_topology_search_records(
        resolved_records,
        owner_name="nucleotide likelihood TBR search",
    )
    resolved_surface = resolve_nucleotide_topology_search_surface(
        resolved_tree,
        normalized_records,
        model_name=model_name,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    input_tree_newick = dumps_newick(resolved_tree)
    start_result = reoptimize_nucleotide_topology_tree(
        resolved_tree,
        compressed_patterns=compressed_patterns,
        resolved_surface=resolved_surface,
        branch_reoptimization_policy="coordinate-branch-lengths",
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    current_tree = start_result.optimized_tree.copy().refresh()
    current_log_likelihood = start_result.log_likelihood
    start_tree_newick = dumps_newick(current_tree)
    equal_best_accumulator = initialize_nucleotide_likelihood_equal_best_topology_accumulator()
    record_nucleotide_likelihood_equal_best_topology(
        equal_best_accumulator,
        tree_newick=start_tree_newick,
        log_likelihood=start_result.log_likelihood,
        likelihood_tolerance=validated_equal_best_likelihood_tolerance,
    )
    trace_rows = [
        NucleotideLikelihoodTbrTraceRow(
            event_index=1,
            event_kind="start",
            iteration=0,
            move_type="tbr",
            candidate_topology_fingerprint=rooted_topology_fingerprint(current_tree),
            log_likelihood_before=None,
            log_likelihood_after=current_log_likelihood,
            log_likelihood_delta=None,
            accepted_move=False,
            trace_reason="search-start",
            tree_before_newick=None,
            tree_after_newick=start_tree_newick,
            cut_edge_id=None,
            left_attachment_branch_id=None,
            right_attachment_branch_id=None,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            branch_reoptimization_scope="all-branches",
            optimized_branch_count=len(start_result.optimized_branch_ids),
            optimized_branch_clade_ids=resolve_reoptimized_branch_clade_ids(
                current_tree,
                start_result.optimized_branch_ids,
            ),
            branch_reoptimization_converged=start_result.converged,
            branch_optimization_pass_count=start_result.optimization_pass_count,
            branch_function_evaluation_count=start_result.function_evaluation_count,
            boundary_warning_messages=[
                *resolved_surface.substitution_parameter_warnings,
                *start_result.boundary_warning_messages,
            ],
            stopping_reason=None,
        )
    ]
    accepted_move_count = 0
    evaluated_neighbor_count = 0
    total_branch_optimization_pass_count = start_result.optimization_pass_count
    total_branch_function_evaluation_count = start_result.function_evaluation_count
    seen_topology_fingerprints = {rooted_topology_fingerprint(current_tree)}
    while True:
        improving_candidate: RootedTbrMoveCandidate | None = None
        improving_result: BranchReoptimizationResult | None = None
        improving_newick: str | None = None
        improving_topology_fingerprint: str | None = None
        best_positive_delta: float | None = None
        best_positive_newick: str | None = None
        seen_neighbor_fingerprints: set[str] = set()
        failure_detected = False
        for candidate in iter_rooted_tbr_move_candidates(current_tree):
            try:
                moved_tree = apply_rooted_tbr_move(current_tree, candidate)
            except Exception:
                failure_detected = True
                break
            moved_fingerprint = rooted_topology_fingerprint(moved_tree)
            if moved_fingerprint in seen_neighbor_fingerprints:
                continue
            seen_neighbor_fingerprints.add(moved_fingerprint)
            initialized_neighbor_tree = initialize_generated_nucleotide_topology_search_tree(
                moved_tree,
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
            )
            try:
                neighbor_result = reoptimize_nucleotide_topology_tree(
                    initialized_neighbor_tree,
                    compressed_patterns=compressed_patterns,
                    resolved_surface=resolved_surface,
                    branch_reoptimization_policy=validated_branch_reoptimization_policy,
                    lower_branch_length_bound=lower_branch_length_bound,
                    upper_branch_length_bound=upper_branch_length_bound,
                    improvement_tolerance=improvement_tolerance,
                    max_coordinate_passes=max_coordinate_passes,
                )
            except Exception:
                failure_detected = True
                break
            evaluated_neighbor_count += 1
            total_branch_optimization_pass_count += (
                neighbor_result.optimization_pass_count
            )
            total_branch_function_evaluation_count += (
                neighbor_result.function_evaluation_count
            )
            optimized_neighbor_newick = dumps_newick(neighbor_result.optimized_tree)
            record_nucleotide_likelihood_equal_best_topology(
                equal_best_accumulator,
                tree_newick=optimized_neighbor_newick,
                log_likelihood=neighbor_result.log_likelihood,
                likelihood_tolerance=validated_equal_best_likelihood_tolerance,
            )
            optimized_neighbor_topology_fingerprint = rooted_topology_fingerprint(
                neighbor_result.optimized_tree
            )
            log_likelihood_delta = neighbor_result.log_likelihood - current_log_likelihood
            if log_likelihood_delta > 0.0 and (
                best_positive_delta is None
                or best_positive_newick is None
                or prefer_higher_likelihood(
                    neighbor_result.log_likelihood,
                    optimized_neighbor_newick,
                    current_log_likelihood + best_positive_delta,
                    best_positive_newick,
                )
            ):
                best_positive_delta = log_likelihood_delta
                best_positive_newick = optimized_neighbor_newick
            if neighbor_result.log_likelihood <= current_log_likelihood:
                continue
            if log_likelihood_delta <= validated_search_improvement_tolerance:
                continue
            if (
                improving_result is None
                or improving_newick is None
                or prefer_higher_likelihood(
                    neighbor_result.log_likelihood,
                    optimized_neighbor_newick,
                    improving_result.log_likelihood,
                    improving_newick,
                )
            ):
                improving_candidate = candidate
                improving_result = neighbor_result
                improving_newick = optimized_neighbor_newick
                improving_topology_fingerprint = optimized_neighbor_topology_fingerprint
        convergence_decision = resolve_nucleotide_likelihood_search_convergence_decision(
            best_improving_delta=(
                best_positive_delta
                if improving_result is None
                else improving_result.log_likelihood - current_log_likelihood
            ),
            improvement_tolerance=validated_search_improvement_tolerance,
            candidate_topology_fingerprint=improving_topology_fingerprint,
            seen_topology_fingerprints=seen_topology_fingerprints,
            failure_detected=failure_detected,
        )
        if convergence_decision.should_stop:
            stopping_reason = convergence_decision.stopping_reason or "search-failure"
            break
        accepted_move_count += 1
        log_likelihood_before = current_log_likelihood
        tree_before_newick = dumps_newick(current_tree)
        current_tree = improving_result.optimized_tree.copy().refresh()
        current_log_likelihood = improving_result.log_likelihood
        seen_topology_fingerprints.add(rooted_topology_fingerprint(current_tree))
        trace_rows.append(
            NucleotideLikelihoodTbrTraceRow(
                event_index=len(trace_rows) + 1,
                event_kind="accepted-move",
                iteration=accepted_move_count,
                move_type="tbr",
                candidate_topology_fingerprint=rooted_topology_fingerprint(current_tree),
                log_likelihood_before=log_likelihood_before,
                log_likelihood_after=current_log_likelihood,
                log_likelihood_delta=current_log_likelihood - log_likelihood_before,
                accepted_move=True,
                trace_reason="accepted-improving-move",
                tree_before_newick=tree_before_newick,
                tree_after_newick=improving_newick,
                cut_edge_id=improving_candidate.cut_edge_id,
                left_attachment_branch_id=improving_candidate.left_attachment_branch_id,
                right_attachment_branch_id=improving_candidate.right_attachment_branch_id,
                branch_reoptimization_policy=validated_branch_reoptimization_policy,
                branch_reoptimization_scope="all-branches",
                optimized_branch_count=len(improving_result.optimized_branch_ids),
                optimized_branch_clade_ids=resolve_reoptimized_branch_clade_ids(
                    current_tree,
                    improving_result.optimized_branch_ids,
                ),
                branch_reoptimization_converged=improving_result.converged,
                branch_optimization_pass_count=improving_result.optimization_pass_count,
                branch_function_evaluation_count=improving_result.function_evaluation_count,
                boundary_warning_messages=list(improving_result.boundary_warning_messages),
                stopping_reason=None,
            )
        )
    trace_rows.append(
        NucleotideLikelihoodTbrTraceRow(
            event_index=len(trace_rows) + 1,
            event_kind="final",
            iteration=accepted_move_count,
            move_type="tbr",
            candidate_topology_fingerprint=rooted_topology_fingerprint(current_tree),
            log_likelihood_before=None,
            log_likelihood_after=current_log_likelihood,
            log_likelihood_delta=None,
            accepted_move=False,
            trace_reason=stopping_reason,
            tree_before_newick=None,
            tree_after_newick=dumps_newick(current_tree),
            cut_edge_id=None,
            left_attachment_branch_id=None,
            right_attachment_branch_id=None,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            branch_reoptimization_scope="none",
            optimized_branch_count=0,
            optimized_branch_clade_ids=[],
            branch_reoptimization_converged=None,
            branch_optimization_pass_count=0,
            branch_function_evaluation_count=0,
            boundary_warning_messages=[],
            stopping_reason=stopping_reason,
        )
    )
    return NucleotideLikelihoodTbrSearchReport(
        algorithm="nucleotide-likelihood-tbr-search",
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
        equal_best_tree_report=build_nucleotide_likelihood_equal_best_tree_report(
            equal_best_accumulator,
            likelihood_tolerance=validated_equal_best_likelihood_tolerance,
            retention_cap=validated_equal_best_tree_cap,
        ),
        trace_rows=trace_rows,
    )


def search_nucleotide_likelihood_tbr_from_alignment(
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
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    search_improvement_tolerance: float = 0.0,
    equal_best_likelihood_tolerance: float = 0.0,
    equal_best_tree_cap: int = 1,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodTbrSearchReport:
    """Search one rooted binary nucleotide tree path by likelihood-improving rooted TBR moves."""
    return search_nucleotide_likelihood_tbr(
        tree_path,
        alignment_path,
        model_name=model_name,
        branch_reoptimization_policy=branch_reoptimization_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        search_improvement_tolerance=search_improvement_tolerance,
        equal_best_likelihood_tolerance=equal_best_likelihood_tolerance,
        equal_best_tree_cap=equal_best_tree_cap,
        max_coordinate_passes=max_coordinate_passes,
    )


def validate_nucleotide_likelihood_tbr_branch_reoptimization_policy(policy: str) -> str:
    """Validate one rooted likelihood TBR branch-reoptimization policy."""
    normalized_policy = policy.strip().lower()
    if normalized_policy not in _SUPPORTED_TBR_BRANCH_REOPTIMIZATION_POLICIES:
        raise ValueError(
            "branch_reoptimization_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_TBR_BRANCH_REOPTIMIZATION_POLICIES))
        )
    return normalized_policy


def write_nucleotide_likelihood_tbr_trace_table(
    path: Path,
    report: NucleotideLikelihoodTbrSearchReport,
) -> Path:
    """Write one deterministic rooted likelihood TBR trace table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "event_index",
        "event_kind",
        "iteration",
        "move_type",
        "candidate_topology_fingerprint",
        "log_likelihood_before",
        "log_likelihood_after",
        "log_likelihood_delta",
        "accepted_move",
        "trace_reason",
        "tree_before_newick",
        "tree_after_newick",
        "cut_edge_id",
        "left_attachment_branch_id",
        "right_attachment_branch_id",
        "branch_reoptimization_policy",
        "branch_reoptimization_scope",
        "optimized_branch_count",
        "optimized_branch_clade_ids",
        "branch_reoptimization_converged",
        "branch_optimization_pass_count",
        "branch_function_evaluation_count",
        "boundary_warning_messages",
        "stopping_reason",
    ]
    rows = ["\t".join(columns)]
    for row in report.trace_rows:
        payload = [
            row.event_index,
            row.event_kind,
            row.iteration,
            row.move_type,
            row.candidate_topology_fingerprint,
            row.log_likelihood_before,
            row.log_likelihood_after,
            row.log_likelihood_delta,
            str(row.accepted_move).lower(),
            row.trace_reason,
            row.tree_before_newick,
            row.tree_after_newick,
            row.cut_edge_id,
            row.left_attachment_branch_id,
            row.right_attachment_branch_id,
            row.branch_reoptimization_policy,
            row.branch_reoptimization_scope,
            row.optimized_branch_count,
            ",".join(row.optimized_branch_clade_ids),
            row.branch_reoptimization_converged,
            row.branch_optimization_pass_count,
            row.branch_function_evaluation_count,
            ",".join(row.boundary_warning_messages),
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


def write_nucleotide_likelihood_tbr_run_json(
    path: Path,
    report: NucleotideLikelihoodTbrSearchReport,
) -> Path:
    """Write one machine-readable rooted likelihood TBR search payload."""
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
        "equal_best_tree_report": {
            "likelihood_tolerance": report.equal_best_tree_report.likelihood_tolerance,
            "retention_cap": report.equal_best_tree_report.retention_cap,
            "retained_tree_count": report.equal_best_tree_report.retained_tree_count,
            "omitted_tree_count": report.equal_best_tree_report.omitted_tree_count,
            "best_log_likelihood": report.equal_best_tree_report.best_log_likelihood,
            "consensus_method": report.equal_best_tree_report.consensus_method,
            "consensus_newick": report.equal_best_tree_report.consensus_newick,
            "rows": [
                {
                    "retained_rank": row.retained_rank,
                    "topology_fingerprint": row.topology_fingerprint,
                    "tree_newick": row.tree_newick,
                    "log_likelihood": row.log_likelihood,
                }
                for row in report.equal_best_tree_report.rows
            ],
        },
        "trace_rows": [
            {
                "event_index": row.event_index,
                "event_kind": row.event_kind,
                "iteration": row.iteration,
                "move_type": row.move_type,
                "candidate_topology_fingerprint": row.candidate_topology_fingerprint,
                "log_likelihood_before": row.log_likelihood_before,
                "log_likelihood_after": row.log_likelihood_after,
                "log_likelihood_delta": row.log_likelihood_delta,
                "accepted_move": row.accepted_move,
                "trace_reason": row.trace_reason,
                "tree_before_newick": row.tree_before_newick,
                "tree_after_newick": row.tree_after_newick,
                "cut_edge_id": row.cut_edge_id,
                "left_attachment_branch_id": row.left_attachment_branch_id,
                "right_attachment_branch_id": row.right_attachment_branch_id,
                "branch_reoptimization_policy": row.branch_reoptimization_policy,
                "branch_reoptimization_scope": row.branch_reoptimization_scope,
                "optimized_branch_count": row.optimized_branch_count,
                "optimized_branch_clade_ids": row.optimized_branch_clade_ids,
                "branch_reoptimization_converged": row.branch_reoptimization_converged,
                "branch_optimization_pass_count": row.branch_optimization_pass_count,
                "branch_function_evaluation_count": row.branch_function_evaluation_count,
                "boundary_warning_messages": row.boundary_warning_messages,
                "stopping_reason": row.stopping_reason,
            }
            for row in report.trace_rows
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_tbr_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodTbrSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted likelihood TBR search run."""
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
    best_tree_path = write_nucleotide_likelihood_best_tree_set(
        out_dir / "best_trees.nwk",
        report.equal_best_tree_report,
    )
    trace_path = write_nucleotide_likelihood_tbr_trace_table(
        out_dir / "search_trace.tsv",
        report,
    )
    run_json_path = write_nucleotide_likelihood_tbr_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "input_tree_path": input_tree_path,
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "best_tree_path": best_tree_path,
        "trace_path": trace_path,
        "run_json_path": run_json_path,
    }
