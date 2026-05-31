from __future__ import annotations

import json
from pathlib import Path

import numpy

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.topology.affected_subtrees import summarize_affected_subtrees
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodSprSearchReport,
    NucleotideLikelihoodSprTraceRow,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    BranchReoptimizationResult,
    normalize_nucleotide_topology_search_records,
    prefer_higher_likelihood,
    reoptimize_nucleotide_topology_tree_branch_subset,
    reoptimize_nucleotide_topology_tree,
    resolve_reoptimized_branch_clade_ids,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
    resolve_nucleotide_topology_search_tree,
    validate_nucleotide_topology_search_tree,
)
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveCandidate,
    apply_rooted_spr_move,
    iter_rooted_spr_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, descendant_taxa

_SUPPORTED_SPR_BRANCH_REOPTIMIZATION_POLICIES = frozenset(
    {"coordinate-branch-lengths", "spr-local-affected-branches"}
)


def search_nucleotide_likelihood_spr(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
    evaluation_budget: int | None = None,
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
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodSprSearchReport:
    """Search one rooted binary nucleotide tree by likelihood-improving rooted SPR moves."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = resolve_nucleotide_topology_search_records(
        records
    )
    validated_branch_reoptimization_policy = validate_nucleotide_likelihood_spr_branch_reoptimization_policy(
        branch_reoptimization_policy
    )
    validated_evaluation_budget = validate_likelihood_spr_evaluation_budget(
        evaluation_budget
    )
    validate_nucleotide_topology_search_tree(
        resolved_tree,
        workflow_name="nucleotide likelihood SPR search",
    )
    normalized_records, compressed_patterns = normalize_nucleotide_topology_search_records(
        resolved_records,
        owner_name="nucleotide likelihood SPR search",
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
    trace_rows = [
        NucleotideLikelihoodSprTraceRow(
            event_index=1,
            event_kind="start",
            iteration=0,
            log_likelihood_before=None,
            log_likelihood_after=current_log_likelihood,
            log_likelihood_delta=None,
            tree_before_newick=None,
            tree_after_newick=start_tree_newick,
            pruned_clade_id=None,
            regraft_target_branch_id=None,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            branch_reoptimization_scope="all-branches",
            affected_branch_clade_ids=[],
            optimized_branch_count=len(start_result.optimized_branch_ids),
            optimized_branch_clade_ids=resolve_reoptimized_branch_clade_ids(
                current_tree,
                start_result.optimized_branch_ids,
            ),
            branch_optimization_pass_count=start_result.optimization_pass_count,
            branch_function_evaluation_count=start_result.function_evaluation_count,
            stopping_reason=None,
        )
    ]
    accepted_move_count = 0
    evaluated_neighbor_count = 0
    total_branch_optimization_pass_count = start_result.optimization_pass_count
    total_branch_function_evaluation_count = start_result.function_evaluation_count
    stopping_reason = "no-improving-neighbor"
    while True:
        improving_candidate: RootedSprMoveCandidate | None = None
        improving_result: BranchReoptimizationResult | None = None
        improving_newick: str | None = None
        budget_exhausted = False
        seen_neighbor_newicks: set[str] = set()
        for candidate in iter_rooted_spr_move_candidates(current_tree):
            if (
                validated_evaluation_budget is not None
                and evaluated_neighbor_count >= validated_evaluation_budget
            ):
                budget_exhausted = True
                break
            neighbor_tree = apply_rooted_spr_move(current_tree, candidate)
            neighbor_newick = dumps_newick(neighbor_tree)
            if neighbor_newick in seen_neighbor_newicks:
                continue
            seen_neighbor_newicks.add(neighbor_newick)
            affected_branch_clade_ids = resolve_spr_affected_branch_clade_ids(
                current_tree,
                neighbor_tree,
                candidate,
            )
            if validated_branch_reoptimization_policy == "coordinate-branch-lengths":
                neighbor_result = reoptimize_nucleotide_topology_tree(
                    neighbor_tree,
                    compressed_patterns=compressed_patterns,
                    resolved_surface=resolved_surface,
                    branch_reoptimization_policy=validated_branch_reoptimization_policy,
                    lower_branch_length_bound=lower_branch_length_bound,
                    upper_branch_length_bound=upper_branch_length_bound,
                    improvement_tolerance=improvement_tolerance,
                    max_coordinate_passes=max_coordinate_passes,
                )
            else:
                neighbor_result = reoptimize_nucleotide_topology_tree_branch_subset(
                    neighbor_tree,
                    compressed_patterns=compressed_patterns,
                    resolved_surface=resolved_surface,
                    optimized_branch_ids=resolve_spr_local_branch_ids(
                        neighbor_tree,
                        affected_branch_clade_ids,
                    ),
                    lower_branch_length_bound=lower_branch_length_bound,
                    upper_branch_length_bound=upper_branch_length_bound,
                    improvement_tolerance=improvement_tolerance,
                    max_coordinate_passes=max_coordinate_passes,
                )
            evaluated_neighbor_count += 1
            total_branch_optimization_pass_count += (
                neighbor_result.optimization_pass_count
            )
            total_branch_function_evaluation_count += (
                neighbor_result.function_evaluation_count
            )
            optimized_neighbor_newick = dumps_newick(neighbor_result.optimized_tree)
            if neighbor_result.log_likelihood <= current_log_likelihood:
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
        if improving_candidate is None or improving_result is None or improving_newick is None:
            stopping_reason = (
                "evaluation-budget-exhausted"
                if budget_exhausted
                else "no-improving-neighbor"
            )
            break
        accepted_move_count += 1
        log_likelihood_before = current_log_likelihood
        tree_before = current_tree.copy().refresh()
        tree_before_newick = dumps_newick(tree_before)
        current_tree = improving_result.optimized_tree.copy().refresh()
        current_log_likelihood = improving_result.log_likelihood
        trace_rows.append(
            NucleotideLikelihoodSprTraceRow(
                event_index=len(trace_rows) + 1,
                event_kind="accepted-move",
                iteration=accepted_move_count,
                log_likelihood_before=log_likelihood_before,
                log_likelihood_after=current_log_likelihood,
                log_likelihood_delta=current_log_likelihood - log_likelihood_before,
                tree_before_newick=tree_before_newick,
                tree_after_newick=improving_newick,
                pruned_clade_id=improving_candidate.pruned_clade_id,
                regraft_target_branch_id=improving_candidate.regraft_target_branch_id,
                branch_reoptimization_policy=validated_branch_reoptimization_policy,
                branch_reoptimization_scope=resolve_spr_branch_reoptimization_scope(
                    validated_branch_reoptimization_policy
                ),
                affected_branch_clade_ids=resolve_spr_affected_branch_clade_ids(
                    tree_before,
                    current_tree,
                    improving_candidate,
                ),
                optimized_branch_count=len(improving_result.optimized_branch_ids),
                optimized_branch_clade_ids=resolve_reoptimized_branch_clade_ids(
                    current_tree,
                    improving_result.optimized_branch_ids,
                ),
                branch_optimization_pass_count=improving_result.optimization_pass_count,
                branch_function_evaluation_count=improving_result.function_evaluation_count,
                stopping_reason=None,
            )
        )
        if budget_exhausted:
            stopping_reason = "evaluation-budget-exhausted"
            break
    trace_rows.append(
        NucleotideLikelihoodSprTraceRow(
            event_index=len(trace_rows) + 1,
            event_kind="final",
            iteration=accepted_move_count,
            log_likelihood_before=None,
            log_likelihood_after=current_log_likelihood,
            log_likelihood_delta=None,
            tree_before_newick=None,
            tree_after_newick=dumps_newick(current_tree),
            pruned_clade_id=None,
            regraft_target_branch_id=None,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            branch_reoptimization_scope="none",
            affected_branch_clade_ids=[],
            optimized_branch_count=0,
            optimized_branch_clade_ids=[],
            branch_optimization_pass_count=0,
            branch_function_evaluation_count=0,
            stopping_reason=stopping_reason,
        )
    )
    return NucleotideLikelihoodSprSearchReport(
        algorithm="nucleotide-likelihood-spr-search",
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
        evaluation_budget=validated_evaluation_budget,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        substitution_parameter_policy=resolved_surface.substitution_parameter_policy,
        substitution_parameter_values=resolved_surface.substitution_parameter_values,
        substitution_parameter_warnings=resolved_surface.substitution_parameter_warnings,
        total_branch_optimization_pass_count=total_branch_optimization_pass_count,
        total_branch_function_evaluation_count=total_branch_function_evaluation_count,
        stopping_reason=stopping_reason,
        trace_rows=trace_rows,
    )


def search_nucleotide_likelihood_spr_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    branch_reoptimization_policy: str = "coordinate-branch-lengths",
    evaluation_budget: int | None = None,
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
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodSprSearchReport:
    """Search one rooted binary nucleotide tree path by likelihood-improving rooted SPR moves."""
    return search_nucleotide_likelihood_spr(
        tree_path,
        alignment_path,
        model_name=model_name,
        branch_reoptimization_policy=branch_reoptimization_policy,
        evaluation_budget=evaluation_budget,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def validate_nucleotide_likelihood_spr_branch_reoptimization_policy(policy: str) -> str:
    """Validate one rooted likelihood SPR branch-reoptimization policy."""
    normalized_policy = policy.strip().lower()
    if normalized_policy not in _SUPPORTED_SPR_BRANCH_REOPTIMIZATION_POLICIES:
        raise ValueError(
            "branch_reoptimization_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_SPR_BRANCH_REOPTIMIZATION_POLICIES))
        )
    return normalized_policy


def resolve_spr_branch_reoptimization_scope(policy: str) -> str:
    """Render one durable recomputation-scope label for rooted likelihood SPR moves."""
    if policy == "coordinate-branch-lengths":
        return "all-branches"
    if policy == "spr-local-affected-branches":
        return "local-affected-branches"
    raise ValueError(f"unsupported rooted likelihood SPR reoptimization policy '{policy}'")


def resolve_spr_affected_branch_clade_ids(
    original_tree: PhyloTree,
    moved_tree: PhyloTree,
    candidate: RootedSprMoveCandidate,
) -> list[str]:
    """Resolve the affected branch-clade ledger for one rooted SPR move."""
    affected_branch_clade_ids = list(
        summarize_affected_subtrees(original_tree, moved_tree).affected_branch_clade_ids
    )
    affected_branch_clade_ids.append(candidate.pruned_clade_id)
    if candidate.regraft_target_branch_id != "root":
        affected_branch_clade_ids.append(candidate.regraft_target_branch_id)
    return sorted(set(affected_branch_clade_ids))


def resolve_spr_local_branch_ids(
    moved_tree: PhyloTree,
    affected_branch_clade_ids: list[str],
) -> list[str]:
    """Resolve moved-tree branch identifiers for one local rooted SPR reoptimization set."""
    branch_ids_by_clade_id = {
        canonical_clade_id(frozenset(descendant_taxa(node))): node.node_id or ""
        for node in moved_tree.iter_nodes(order="preorder")
        if node is not moved_tree.root
    }
    return [
        branch_ids_by_clade_id[clade_id]
        for clade_id in affected_branch_clade_ids
        if clade_id in branch_ids_by_clade_id
    ]


def write_nucleotide_likelihood_spr_trace_table(
    path: Path,
    report: NucleotideLikelihoodSprSearchReport,
) -> Path:
    """Write one deterministic rooted likelihood SPR trace table."""
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
        "pruned_clade_id",
        "regraft_target_branch_id",
        "branch_reoptimization_policy",
        "branch_reoptimization_scope",
        "affected_branch_clade_ids",
        "optimized_branch_count",
        "optimized_branch_clade_ids",
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
            row.pruned_clade_id,
            row.regraft_target_branch_id,
            row.branch_reoptimization_policy,
            row.branch_reoptimization_scope,
            ",".join(row.affected_branch_clade_ids),
            row.optimized_branch_count,
            ",".join(row.optimized_branch_clade_ids),
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


def write_nucleotide_likelihood_spr_run_json(
    path: Path,
    report: NucleotideLikelihoodSprSearchReport,
) -> Path:
    """Write one machine-readable rooted likelihood SPR search payload."""
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
        "evaluation_budget": report.evaluation_budget,
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
                "pruned_clade_id": row.pruned_clade_id,
                "regraft_target_branch_id": row.regraft_target_branch_id,
                "branch_reoptimization_policy": row.branch_reoptimization_policy,
                "branch_reoptimization_scope": row.branch_reoptimization_scope,
                "affected_branch_clade_ids": row.affected_branch_clade_ids,
                "optimized_branch_count": row.optimized_branch_count,
                "optimized_branch_clade_ids": row.optimized_branch_clade_ids,
                "branch_optimization_pass_count": row.branch_optimization_pass_count,
                "branch_function_evaluation_count": row.branch_function_evaluation_count,
                "stopping_reason": row.stopping_reason,
            }
            for row in report.trace_rows
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_spr_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodSprSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted likelihood SPR search run."""
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
    trace_path = write_nucleotide_likelihood_spr_trace_table(
        out_dir / "search_trace.tsv",
        report,
    )
    run_json_path = write_nucleotide_likelihood_spr_run_json(
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


def validate_likelihood_spr_evaluation_budget(
    evaluation_budget: int | None,
) -> int | None:
    """Validate the optional rooted SPR evaluation budget."""
    if evaluation_budget is None:
        return None
    if evaluation_budget < 1:
        raise ValueError("evaluation_budget must be at least one when provided")
    return evaluation_budget
