from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key
import json
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.newick import (
    dumps_newick,
    loads_newick,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodMultiStartRunSummary,
    NucleotideLikelihoodMultiStartSearchReport,
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodSprSearchReport,
)
from bijux_phylogenetics.phylo.likelihood.nni_search import (
    search_nucleotide_likelihood_nni,
)
from bijux_phylogenetics.phylo.likelihood.spr_search import (
    search_nucleotide_likelihood_spr,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_generation import (
    build_random_likelihood_start_tree,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    build_nucleotide_likelihood_starting_tree_pool,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_selection import (
    select_nucleotide_likelihood_starting_tree_pool,
    validate_nucleotide_likelihood_starting_tree_selection_policy,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_validation import (
    validate_nucleotide_likelihood_starting_tree,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    normalize_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_tree,
    validate_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_LIKELIHOOD_MULTI_START_METHODS = frozenset({"nni", "spr"})
_SUPPORTED_START_TREE_SOURCE_POLICIES = frozenset({"input-tree-plus-random-tree"})
_MAX_RANDOM_START_ATTEMPTS = 10_000


@dataclass(frozen=True, slots=True)
class _StartTreeCandidate:
    """One prepared likelihood multi-start tree with an explicit provenance label."""

    source_kind: str
    source_label: str
    generation_seed: int | None
    tree: PhyloTree


def search_nucleotide_likelihood_multi_start(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    local_search_method: str = "nni",
    start_tree_count: int = 4,
    start_tree_source_policy: str = "input-tree-plus-random-tree",
    starting_tree_selection_policy: str | None = None,
    selected_start_tree_count: int | None = None,
    starting_tree_selection_seed: int = 1,
    starting_tree_strategy_priority: tuple[str, ...] | list[str] | None = None,
    start_tree_seed: int = 1,
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
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodMultiStartSearchReport:
    """Search a nucleotide likelihood surface from multiple independent rooted starts."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = (
        resolve_nucleotide_topology_search_records(records)
    )
    normalized_local_search_method = validate_likelihood_multi_start_method(
        local_search_method
    )
    validated_start_tree_count = validate_likelihood_multi_start_start_tree_count(
        start_tree_count
    )
    validated_start_tree_source_policy = validate_likelihood_multi_start_source_policy(
        start_tree_source_policy
    )
    validated_starting_tree_selection_policy = (
        None
        if starting_tree_selection_policy is None
        else validate_nucleotide_likelihood_starting_tree_selection_policy(
            starting_tree_selection_policy
        )
    )
    validated_branch_reoptimization_policy = validate_branch_reoptimization_policy(
        branch_reoptimization_policy
    )
    validated_evaluation_budget = validate_likelihood_multi_start_evaluation_budget(
        normalized_local_search_method,
        evaluation_budget,
    )
    _normalized_records, compressed_patterns = (
        normalize_nucleotide_topology_search_records(
            resolved_records,
            owner_name="nucleotide likelihood multi-start search",
        )
    )
    validate_nucleotide_likelihood_starting_tree(
        resolved_tree,
        compressed_patterns,
        model_name=model_name,
        workflow_name="nucleotide likelihood multi-start search",
    )
    if validated_starting_tree_selection_policy is None:
        start_tree_candidates = build_likelihood_multi_start_candidates(
            resolved_tree,
            start_tree_count=validated_start_tree_count,
            start_tree_source_policy=validated_start_tree_source_policy,
            start_tree_seed=start_tree_seed,
        )
        available_start_tree_count = len(start_tree_candidates)
        report_start_tree_source_policy = validated_start_tree_source_policy
    else:
        starting_tree_pool_report = (
            _build_scored_starting_tree_pool_for_multi_start_search(
                resolved_tree,
                resolved_records,
                model_name=model_name,
                random_start_tree_count=max(1, validated_start_tree_count - 2),
                random_start_tree_seed=start_tree_seed,
            )
        )
        selected_start_trees = select_nucleotide_likelihood_starting_tree_pool(
            starting_tree_pool_report,
            starting_tree_selection_policy=validated_starting_tree_selection_policy,
            selected_start_tree_count=selected_start_tree_count,
            selection_seed=starting_tree_selection_seed,
            strategy_priority=starting_tree_strategy_priority,
        )
        start_tree_candidates = [
            _StartTreeCandidate(
                source_kind=row.source_strategy,
                source_label=row.tree_id,
                generation_seed=row.generation_seed,
                tree=loads_newick(row.tree_newick),
            )
            for row in selected_start_trees
        ]
        available_start_tree_count = len(
            starting_tree_pool_report.starting_tree_summaries
        )
        report_start_tree_source_policy = "scored-starting-tree-pool"
    run_summaries: list[NucleotideLikelihoodMultiStartRunSummary] = []
    local_reports: list[
        NucleotideLikelihoodNniSearchReport | NucleotideLikelihoodSprSearchReport
    ] = []
    for candidate in start_tree_candidates:
        validate_nucleotide_likelihood_starting_tree(
            candidate.tree,
            compressed_patterns,
            model_name=model_name,
            workflow_name="nucleotide likelihood multi-start search",
        )
        local_report = _run_local_search(
            candidate.tree,
            resolved_records,
            model_name=model_name,
            local_search_method=normalized_local_search_method,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            evaluation_budget=validated_evaluation_budget,
            kappa=kappa,
            base_frequencies=base_frequencies,
            exchangeabilities=exchangeabilities,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
        local_reports.append(local_report)
        run_summaries.append(
            NucleotideLikelihoodMultiStartRunSummary(
                search_algorithm=local_report.algorithm,
                start_tree_source_kind=candidate.source_kind,
                start_tree_source_label=candidate.source_label,
                start_tree_generation_seed=candidate.generation_seed,
                start_tree_newick=local_report.start_tree_newick,
                start_log_likelihood=local_report.start_log_likelihood,
                final_tree_newick=local_report.final_tree_newick,
                final_log_likelihood=local_report.final_log_likelihood,
                final_topology_fingerprint=rooted_topology_fingerprint_from_newick(
                    local_report.final_tree_newick
                ),
                search_iteration_count=resolve_multi_start_local_search_iteration_count(
                    local_report
                ),
                accepted_move_count=local_report.accepted_move_count,
                evaluated_neighbor_count=local_report.evaluated_neighbor_count,
                final_likelihood_rank=0,
                branch_reoptimization_policy=local_report.branch_reoptimization_policy,
                substitution_parameter_policy=local_report.substitution_parameter_policy,
                substitution_parameter_values=dict(
                    local_report.substitution_parameter_values
                ),
                substitution_parameter_warnings=list(
                    local_report.substitution_parameter_warnings
                ),
                total_branch_optimization_pass_count=local_report.total_branch_optimization_pass_count,
                total_branch_function_evaluation_count=local_report.total_branch_function_evaluation_count,
                stopping_reason=local_report.stopping_reason,
                best_run=False,
            )
        )
    ranked_run_indices = rank_likelihood_multi_start_runs(run_summaries)
    for final_likelihood_rank, ranked_run_index in enumerate(
        ranked_run_indices,
        start=1,
    ):
        run_summaries[ranked_run_index].final_likelihood_rank = final_likelihood_rank
    best_run_index = ranked_run_indices[0]
    run_summaries[best_run_index].best_run = True
    best_run = run_summaries[best_run_index]
    first_local_report = local_reports[0]
    return NucleotideLikelihoodMultiStartSearchReport(
        algorithm="nucleotide-likelihood-multi-start-search",
        model_name=first_local_report.model_name,
        local_search_method=normalized_local_search_method,
        tree_path=None if resolved_tree_path is None else str(resolved_tree_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        input_tree_newick=dumps_newick(resolved_tree),
        start_tree_source_policy=report_start_tree_source_policy,
        starting_tree_selection_policy=validated_starting_tree_selection_policy,
        input_tree_included=any(
            row.start_tree_source_kind == "input-tree" for row in run_summaries
        ),
        available_start_tree_count=available_start_tree_count,
        generated_start_tree_count=sum(
            1 for row in run_summaries if row.start_tree_source_kind != "input-tree"
        ),
        start_tree_count=len(run_summaries),
        start_tree_seed=start_tree_seed,
        evaluation_budget=validated_evaluation_budget,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        best_run_source_label=best_run.start_tree_source_label,
        best_final_tree_newick=best_run.final_tree_newick,
        best_final_log_likelihood=best_run.final_log_likelihood,
        best_final_topology_fingerprint=best_run.final_topology_fingerprint,
        run_summaries=run_summaries,
    )


def search_nucleotide_likelihood_multi_start_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    local_search_method: str = "nni",
    start_tree_count: int = 4,
    start_tree_source_policy: str = "input-tree-plus-random-tree",
    starting_tree_selection_policy: str | None = None,
    selected_start_tree_count: int | None = None,
    starting_tree_selection_seed: int = 1,
    starting_tree_strategy_priority: tuple[str, ...] | list[str] | None = None,
    start_tree_seed: int = 1,
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
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodMultiStartSearchReport:
    """Search one tree path and alignment path by independent likelihood restarts."""
    return search_nucleotide_likelihood_multi_start(
        tree_path,
        alignment_path,
        model_name=model_name,
        local_search_method=local_search_method,
        start_tree_count=start_tree_count,
        start_tree_source_policy=start_tree_source_policy,
        starting_tree_selection_policy=starting_tree_selection_policy,
        selected_start_tree_count=selected_start_tree_count,
        starting_tree_selection_seed=starting_tree_selection_seed,
        starting_tree_strategy_priority=starting_tree_strategy_priority,
        start_tree_seed=start_tree_seed,
        branch_reoptimization_policy=branch_reoptimization_policy,
        evaluation_budget=evaluation_budget,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def validate_likelihood_multi_start_method(local_search_method: str) -> str:
    """Validate the bounded set of local search engines used inside restarts."""
    normalized_local_search_method = local_search_method.strip().lower()
    if normalized_local_search_method not in _SUPPORTED_LIKELIHOOD_MULTI_START_METHODS:
        raise ValueError(
            "local_search_method must be one of "
            + ", ".join(sorted(_SUPPORTED_LIKELIHOOD_MULTI_START_METHODS))
        )
    return normalized_local_search_method


def validate_likelihood_multi_start_start_tree_count(start_tree_count: int) -> int:
    """Require at least one input tree plus one independent generated restart."""
    if start_tree_count < 2:
        raise ValueError("start_tree_count must be at least two for multi-start search")
    return start_tree_count


def validate_likelihood_multi_start_source_policy(start_tree_source_policy: str) -> str:
    """Validate the supported start-tree generation policy."""
    normalized_policy = start_tree_source_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_START_TREE_SOURCE_POLICIES:
        raise ValueError(
            "start_tree_source_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_START_TREE_SOURCE_POLICIES))
        )
    return normalized_policy


def validate_likelihood_multi_start_evaluation_budget(
    local_search_method: str,
    evaluation_budget: int | None,
) -> int | None:
    """Validate restart-level forwarding of optional SPR neighbor budgets."""
    if local_search_method == "nni":
        if evaluation_budget is not None:
            raise ValueError(
                "evaluation_budget is supported only when local_search_method is 'spr'"
            )
        return None
    if evaluation_budget is None:
        return None
    if evaluation_budget < 1:
        raise ValueError("evaluation_budget must be at least one when provided")
    return evaluation_budget


def build_likelihood_multi_start_candidates(
    input_tree: PhyloTree,
    *,
    start_tree_count: int,
    start_tree_source_policy: str,
    start_tree_seed: int,
) -> list[_StartTreeCandidate]:
    """Build one deterministic distinct start tree family over one input taxon scope."""
    validated_policy = validate_likelihood_multi_start_source_policy(
        start_tree_source_policy
    )
    if validated_policy != "input-tree-plus-random-tree":
        raise ValueError(
            "start_tree_source_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_START_TREE_SOURCE_POLICIES))
        )
    ordered_taxa = sorted(input_tree.tip_names)
    candidates = [
        _StartTreeCandidate(
            source_kind="input-tree",
            source_label="input-tree",
            generation_seed=None,
            tree=input_tree.copy().refresh(),
        )
    ]
    seen_start_tree_newicks = {dumps_newick(candidates[0].tree)}
    next_seed = start_tree_seed
    attempt_count = 0
    while len(candidates) < start_tree_count:
        attempt_count += 1
        if attempt_count > _MAX_RANDOM_START_ATTEMPTS:
            raise ValueError(
                "could not generate enough distinct random start trees for multi-start search"
            )
        random_tree = build_random_likelihood_start_tree(ordered_taxa, seed=next_seed)
        next_seed += 1
        candidate_newick = dumps_newick(random_tree)
        if candidate_newick in seen_start_tree_newicks:
            continue
        seen_start_tree_newicks.add(candidate_newick)
        candidates.append(
            _StartTreeCandidate(
                source_kind="random-tree",
                source_label=f"random-tree-seed-{next_seed - 1}",
                generation_seed=next_seed - 1,
                tree=random_tree,
            )
        )
    return candidates


def select_best_likelihood_multi_start_run(
    run_summaries: list[NucleotideLikelihoodMultiStartRunSummary],
) -> int:
    """Choose the deterministic best run by likelihood, topology fingerprint, and tree text."""
    if not run_summaries:
        raise ValueError("run_summaries must not be empty")
    best_index = 0
    best_run = run_summaries[0]
    for index, candidate in enumerate(run_summaries[1:], start=1):
        if _prefer_multi_start_run(candidate, best_run):
            best_index = index
            best_run = candidate
    return best_index


def rank_likelihood_multi_start_runs(
    run_summaries: list[NucleotideLikelihoodMultiStartRunSummary],
) -> list[int]:
    """Rank multi-start runs by final likelihood with deterministic tie-breaking."""
    if not run_summaries:
        raise ValueError("run_summaries must not be empty")
    return sorted(
        range(len(run_summaries)),
        key=cmp_to_key(
            lambda left_index, right_index: _compare_multi_start_runs(
                run_summaries[left_index],
                run_summaries[right_index],
            )
        ),
    )


def _build_scored_starting_tree_pool_for_multi_start_search(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    random_start_tree_count: int,
    random_start_tree_seed: int,
):
    candidate_seed = random_start_tree_seed
    for _attempt in range(_MAX_RANDOM_START_ATTEMPTS):
        try:
            return build_nucleotide_likelihood_starting_tree_pool(
                tree,
                records,
                model_name=model_name,
                random_start_tree_count=random_start_tree_count,
                random_start_tree_seed=candidate_seed,
            )
        except ValueError as error:
            if "duplicate topology hash" not in str(error):
                raise
            candidate_seed += 1
    raise ValueError(
        "could not build enough distinct scored start trees for multi-start search"
    )


def rooted_topology_fingerprint_from_newick(tree_newick: str) -> str:
    """Fingerprint one rooted topology from a deterministic Newick record."""
    return rooted_topology_fingerprint(loads_newick(tree_newick))


def _prefer_multi_start_run(
    left: NucleotideLikelihoodMultiStartRunSummary,
    right: NucleotideLikelihoodMultiStartRunSummary,
) -> bool:
    if left.final_log_likelihood > right.final_log_likelihood and not math.isclose(
        left.final_log_likelihood, right.final_log_likelihood
    ):
        return True
    if right.final_log_likelihood > left.final_log_likelihood and not math.isclose(
        left.final_log_likelihood, right.final_log_likelihood
    ):
        return False
    if left.final_topology_fingerprint != right.final_topology_fingerprint:
        return left.final_topology_fingerprint < right.final_topology_fingerprint
    if left.final_tree_newick != right.final_tree_newick:
        return left.final_tree_newick < right.final_tree_newick
    return left.start_tree_source_label < right.start_tree_source_label


def _compare_multi_start_runs(
    left: NucleotideLikelihoodMultiStartRunSummary,
    right: NucleotideLikelihoodMultiStartRunSummary,
) -> int:
    if _prefer_multi_start_run(left, right):
        return -1
    if _prefer_multi_start_run(right, left):
        return 1
    return 0


def resolve_multi_start_local_search_iteration_count(
    local_report: NucleotideLikelihoodNniSearchReport
    | NucleotideLikelihoodSprSearchReport,
) -> int:
    """Resolve one comparable local-search iteration count for a multi-start run."""
    return local_report.iteration_count


def _run_local_search(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    local_search_method: str,
    branch_reoptimization_policy: str,
    evaluation_budget: int | None,
    kappa: float | None,
    base_frequencies: dict[str, float] | numpy.ndarray | None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ),
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> NucleotideLikelihoodNniSearchReport | NucleotideLikelihoodSprSearchReport:
    if local_search_method == "nni":
        return search_nucleotide_likelihood_nni(
            tree,
            records,
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
    return search_nucleotide_likelihood_spr(
        tree,
        records,
        model_name=model_name,
        branch_reoptimization_policy=branch_reoptimization_policy,
        evaluation_budget=evaluation_budget,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def write_nucleotide_likelihood_multi_start_summary_table(
    path: Path,
    report: NucleotideLikelihoodMultiStartSearchReport,
) -> Path:
    """Write one deterministic multi-start likelihood restart summary ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "start_tree_source_kind",
        "start_tree_source_label",
        "start_tree_generation_seed",
        "search_algorithm",
        "start_log_likelihood",
        "final_log_likelihood",
        "final_likelihood_rank",
        "final_topology_fingerprint",
        "search_iteration_count",
        "accepted_move_count",
        "evaluated_neighbor_count",
        "branch_reoptimization_policy",
        "substitution_parameter_policy",
        "stopping_reason",
        "best_run",
        "start_tree_newick",
        "final_tree_newick",
    ]
    rows = ["\t".join(columns)]
    for row in report.run_summaries:
        payload = [
            row.start_tree_source_kind,
            row.start_tree_source_label,
            row.start_tree_generation_seed,
            row.search_algorithm,
            row.start_log_likelihood,
            row.final_log_likelihood,
            row.final_likelihood_rank,
            row.final_topology_fingerprint,
            row.search_iteration_count,
            row.accepted_move_count,
            row.evaluated_neighbor_count,
            row.branch_reoptimization_policy,
            row.substitution_parameter_policy,
            row.stopping_reason,
            row.best_run,
            row.start_tree_newick,
            row.final_tree_newick,
        ]
        rows.append("\t".join("" if value is None else str(value) for value in payload))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_multi_start_run_json(
    path: Path,
    report: NucleotideLikelihoodMultiStartSearchReport,
) -> Path:
    """Write one machine-readable rooted likelihood multi-start payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "model_name": report.model_name,
        "local_search_method": report.local_search_method,
        "tree_path": report.tree_path,
        "alignment_path": report.alignment_path,
        "taxon_count": report.taxon_count,
        "site_count": report.site_count,
        "pattern_count": report.pattern_count,
        "input_tree_newick": report.input_tree_newick,
        "start_tree_source_policy": report.start_tree_source_policy,
        "starting_tree_selection_policy": report.starting_tree_selection_policy,
        "input_tree_included": report.input_tree_included,
        "available_start_tree_count": report.available_start_tree_count,
        "generated_start_tree_count": report.generated_start_tree_count,
        "start_tree_count": report.start_tree_count,
        "start_tree_seed": report.start_tree_seed,
        "evaluation_budget": report.evaluation_budget,
        "branch_reoptimization_policy": report.branch_reoptimization_policy,
        "best_run_source_label": report.best_run_source_label,
        "best_final_tree_newick": report.best_final_tree_newick,
        "best_final_log_likelihood": report.best_final_log_likelihood,
        "best_final_topology_fingerprint": report.best_final_topology_fingerprint,
        "run_summaries": [
            {
                "search_algorithm": row.search_algorithm,
                "start_tree_source_kind": row.start_tree_source_kind,
                "start_tree_source_label": row.start_tree_source_label,
                "start_tree_generation_seed": row.start_tree_generation_seed,
                "start_tree_newick": row.start_tree_newick,
                "start_log_likelihood": row.start_log_likelihood,
                "final_tree_newick": row.final_tree_newick,
                "final_log_likelihood": row.final_log_likelihood,
                "final_likelihood_rank": row.final_likelihood_rank,
                "final_topology_fingerprint": row.final_topology_fingerprint,
                "search_iteration_count": row.search_iteration_count,
                "accepted_move_count": row.accepted_move_count,
                "evaluated_neighbor_count": row.evaluated_neighbor_count,
                "branch_reoptimization_policy": row.branch_reoptimization_policy,
                "substitution_parameter_policy": row.substitution_parameter_policy,
                "substitution_parameter_values": row.substitution_parameter_values,
                "substitution_parameter_warnings": row.substitution_parameter_warnings,
                "total_branch_optimization_pass_count": row.total_branch_optimization_pass_count,
                "total_branch_function_evaluation_count": row.total_branch_function_evaluation_count,
                "stopping_reason": row.stopping_reason,
                "best_run": row.best_run,
            }
            for row in report.run_summaries
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_nucleotide_likelihood_multi_start_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodMultiStartSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one rooted likelihood multi-start run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    start_tree_path = write_newick_tree_set(
        out_dir / "start_trees.nwk",
        [loads_newick(row.start_tree_newick) for row in report.run_summaries],
    )
    best_tree_path = write_newick(
        out_dir / "best_tree.nwk",
        loads_newick(report.best_final_tree_newick),
    )
    summary_path = write_nucleotide_likelihood_multi_start_summary_table(
        out_dir / "restart_summary.tsv",
        report,
    )
    run_json_path = write_nucleotide_likelihood_multi_start_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "input_tree_path": input_tree_path,
        "start_tree_path": start_tree_path,
        "best_tree_path": best_tree_path,
        "summary_path": summary_path,
        "run_json_path": run_json_path,
    }
