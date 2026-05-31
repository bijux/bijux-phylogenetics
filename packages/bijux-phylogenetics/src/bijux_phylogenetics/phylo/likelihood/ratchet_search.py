from __future__ import annotations

import json
from pathlib import Path
from random import Random

import numpy

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodRatchetBestTreeHistory,
    NucleotideLikelihoodRatchetCycle,
    NucleotideLikelihoodRatchetReport,
    NucleotideLikelihoodSprSearchReport,
    NucleotideLikelihoodTbrSearchReport,
)
from bijux_phylogenetics.phylo.likelihood.nni_search import (
    search_nucleotide_likelihood_nni,
    validate_nucleotide_likelihood_nni_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.spr_search import (
    search_nucleotide_likelihood_spr,
    validate_likelihood_spr_evaluation_budget,
    validate_nucleotide_likelihood_spr_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.tbr_search import (
    search_nucleotide_likelihood_tbr,
    validate_nucleotide_likelihood_tbr_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    normalize_nucleotide_topology_search_records,
    prefer_higher_likelihood,
    reoptimize_nucleotide_topology_tree,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
    resolve_nucleotide_topology_search_tree,
    validate_nucleotide_topology_search_tree,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_RATCHET_LOCAL_SEARCH_METHODS = frozenset({"nni", "spr", "tbr"})


def search_nucleotide_likelihood_ratchet(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    local_search_method: str = "nni",
    cycle_count: int = 4,
    perturbation_seed: int = 1,
    perturbed_site_count: int = 1,
    perturbation_factor: int = 2,
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
) -> NucleotideLikelihoodRatchetReport:
    """Run a native temporary-site-reweighting likelihood ratchet over one tree."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = resolve_nucleotide_topology_search_records(
        records
    )
    validated_local_search_method = validate_nucleotide_likelihood_ratchet_local_search_method(
        local_search_method
    )
    validated_cycle_count = validate_nucleotide_likelihood_ratchet_cycle_count(
        cycle_count
    )
    validated_perturbation_factor = (
        validate_nucleotide_likelihood_ratchet_perturbation_factor(
            perturbation_factor
        )
    )
    validated_branch_reoptimization_policy = (
        validate_nucleotide_likelihood_ratchet_branch_reoptimization_policy(
            branch_reoptimization_policy,
            local_search_method=validated_local_search_method,
        )
    )
    validated_evaluation_budget = validate_nucleotide_likelihood_ratchet_evaluation_budget(
        evaluation_budget,
        local_search_method=validated_local_search_method,
    )
    validate_nucleotide_topology_search_tree(
        resolved_tree,
        workflow_name="nucleotide likelihood ratchet search",
    )
    normalized_records, compressed_patterns = normalize_nucleotide_topology_search_records(
        resolved_records,
        owner_name="nucleotide likelihood ratchet search",
    )
    validated_perturbed_site_count = validate_nucleotide_likelihood_ratchet_perturbed_site_count(
        perturbed_site_count,
        site_count=compressed_patterns.alignment_length,
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
    current_tree_newick = dumps_newick(start_result.optimized_tree)
    current_score = start_result.log_likelihood
    best_tree_newick = current_tree_newick
    best_score = current_score
    best_tree_history_rows = [
        NucleotideLikelihoodRatchetBestTreeHistory(
            history_index=1,
            cycle_index=0,
            best_log_likelihood=best_score,
            best_tree_newick=best_tree_newick,
            best_topology_fingerprint=rooted_topology_fingerprint(
                start_result.optimized_tree
            ),
        )
    ]
    cycle_rows: list[NucleotideLikelihoodRatchetCycle] = []
    rng = Random(perturbation_seed)
    latest_restored_report: (
        NucleotideLikelihoodNniSearchReport
        | NucleotideLikelihoodSprSearchReport
        | NucleotideLikelihoodTbrSearchReport
        | None
    ) = None
    for cycle_index in range(1, validated_cycle_count + 1):
        cycle_start_tree_newick = current_tree_newick
        cycle_start_score = current_score
        reweighted_site_positions = sorted(
            rng.sample(
                range(1, compressed_patterns.alignment_length + 1),
                validated_perturbed_site_count,
            )
        )
        perturbed_records = build_nucleotide_likelihood_ratchet_records(
            normalized_records,
            reweighted_site_positions=reweighted_site_positions,
            perturbation_factor=validated_perturbation_factor,
        )
        perturbed_patterns = compress_alignment_site_patterns_from_records(
            perturbed_records
        )
        perturbed_report = _run_nucleotide_likelihood_ratchet_local_search(
            loads_newick(current_tree_newick),
            perturbed_records,
            model_name=model_name,
            local_search_method=validated_local_search_method,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            evaluation_budget=validated_evaluation_budget,
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
        restored_report = _run_nucleotide_likelihood_ratchet_local_search(
            loads_newick(perturbed_report.final_tree_newick),
            normalized_records,
            model_name=model_name,
            local_search_method=validated_local_search_method,
            branch_reoptimization_policy=validated_branch_reoptimization_policy,
            evaluation_budget=validated_evaluation_budget,
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
        latest_restored_report = restored_report
        current_tree_newick = restored_report.final_tree_newick
        current_score = restored_report.final_log_likelihood
        best_tree_improved = prefer_higher_likelihood(
            current_score,
            current_tree_newick,
            best_score,
            best_tree_newick,
        )
        if best_tree_improved:
            best_score = current_score
            best_tree_newick = current_tree_newick
            best_tree_history_rows.append(
                NucleotideLikelihoodRatchetBestTreeHistory(
                    history_index=len(best_tree_history_rows) + 1,
                    cycle_index=cycle_index,
                    best_log_likelihood=best_score,
                    best_tree_newick=best_tree_newick,
                    best_topology_fingerprint=rooted_topology_fingerprint(
                        loads_newick(best_tree_newick)
                    ),
                )
            )
        cycle_rows.append(
            NucleotideLikelihoodRatchetCycle(
                cycle_index=cycle_index,
                start_log_likelihood=cycle_start_score,
                start_tree_newick=cycle_start_tree_newick,
                reweighted_site_positions=reweighted_site_positions,
                temporary_site_weights={
                    position: validated_perturbation_factor
                    for position in reweighted_site_positions
                },
                perturbation_factor=validated_perturbation_factor,
                perturbed_alignment_length=perturbed_patterns.alignment_length,
                perturbed_pattern_count=perturbed_patterns.pattern_count,
                perturbed_search_algorithm=perturbed_report.algorithm,
                perturbed_score=perturbed_report.final_log_likelihood,
                perturbed_tree_newick=perturbed_report.final_tree_newick,
                perturbed_accepted_move_count=perturbed_report.accepted_move_count,
                perturbed_evaluated_neighbor_count=perturbed_report.evaluated_neighbor_count,
                perturbed_stopping_reason=perturbed_report.stopping_reason,
                restored_search_algorithm=restored_report.algorithm,
                restored_score=restored_report.final_log_likelihood,
                restored_tree_newick=restored_report.final_tree_newick,
                restored_accepted_move_count=restored_report.accepted_move_count,
                restored_evaluated_neighbor_count=restored_report.evaluated_neighbor_count,
                restored_stopping_reason=restored_report.stopping_reason,
                best_score_after_cycle=best_score,
                best_tree_after_cycle=best_tree_newick,
                best_tree_improved=best_tree_improved,
            )
        )
    final_substitution_parameter_policy = (
        resolved_surface.substitution_parameter_policy
        if latest_restored_report is None
        else latest_restored_report.substitution_parameter_policy
    )
    final_substitution_parameter_values = (
        dict(resolved_surface.substitution_parameter_values)
        if latest_restored_report is None
        else dict(latest_restored_report.substitution_parameter_values)
    )
    final_substitution_parameter_warnings = (
        list(resolved_surface.substitution_parameter_warnings)
        if latest_restored_report is None
        else list(latest_restored_report.substitution_parameter_warnings)
    )
    return NucleotideLikelihoodRatchetReport(
        algorithm="nucleotide-likelihood-ratchet-search",
        model_name=resolved_surface.model_name,
        local_search_method=validated_local_search_method,
        tree_path=None if resolved_tree_path is None else str(resolved_tree_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        cycle_count=validated_cycle_count,
        perturbation_seed=perturbation_seed,
        perturbed_site_count=validated_perturbed_site_count,
        perturbation_factor=validated_perturbation_factor,
        input_tree_newick=dumps_newick(resolved_tree),
        start_tree_newick=dumps_newick(start_result.optimized_tree),
        start_log_likelihood=start_result.log_likelihood,
        final_tree_newick=current_tree_newick,
        final_log_likelihood=current_score,
        best_tree_newick=best_tree_newick,
        best_log_likelihood=best_score,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        evaluation_budget=validated_evaluation_budget,
        substitution_parameter_policy=final_substitution_parameter_policy,
        substitution_parameter_values=final_substitution_parameter_values,
        substitution_parameter_warnings=final_substitution_parameter_warnings,
        cycle_rows=cycle_rows,
        best_tree_history_rows=best_tree_history_rows,
    )


def search_nucleotide_likelihood_ratchet_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    local_search_method: str = "nni",
    cycle_count: int = 4,
    perturbation_seed: int = 1,
    perturbed_site_count: int = 1,
    perturbation_factor: int = 2,
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
) -> NucleotideLikelihoodRatchetReport:
    """Run a native likelihood ratchet from one tree path and one alignment path."""
    return search_nucleotide_likelihood_ratchet(
        tree_path,
        alignment_path,
        model_name=model_name,
        local_search_method=local_search_method,
        cycle_count=cycle_count,
        perturbation_seed=perturbation_seed,
        perturbed_site_count=perturbed_site_count,
        perturbation_factor=perturbation_factor,
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


def validate_nucleotide_likelihood_ratchet_local_search_method(
    local_search_method: str,
) -> str:
    """Validate one native likelihood ratchet local search method."""
    normalized_local_search_method = local_search_method.strip().lower()
    if normalized_local_search_method not in _SUPPORTED_RATCHET_LOCAL_SEARCH_METHODS:
        raise ValueError(
            "local_search_method must be one of "
            + ", ".join(sorted(_SUPPORTED_RATCHET_LOCAL_SEARCH_METHODS))
        )
    return normalized_local_search_method


def validate_nucleotide_likelihood_ratchet_cycle_count(cycle_count: int) -> int:
    """Require at least one likelihood ratchet cycle."""
    if cycle_count < 1:
        raise ValueError("cycle_count must be at least one")
    return cycle_count


def validate_nucleotide_likelihood_ratchet_perturbation_factor(
    perturbation_factor: int,
) -> int:
    """Require a real temporary reweighting factor greater than one."""
    if perturbation_factor < 2:
        raise ValueError("perturbation_factor must be at least two")
    return perturbation_factor


def validate_nucleotide_likelihood_ratchet_perturbed_site_count(
    perturbed_site_count: int,
    *,
    site_count: int,
) -> int:
    """Require one valid positive subset of aligned site positions."""
    if perturbed_site_count < 1:
        raise ValueError("perturbed_site_count must be at least one")
    if perturbed_site_count > site_count:
        raise ValueError(
            "perturbed_site_count must not exceed the alignment site count"
        )
    return perturbed_site_count


def validate_nucleotide_likelihood_ratchet_branch_reoptimization_policy(
    branch_reoptimization_policy: str,
    *,
    local_search_method: str,
) -> str:
    """Validate one ratchet local-search branch reoptimization policy."""
    if local_search_method == "nni":
        return validate_nucleotide_likelihood_nni_branch_reoptimization_policy(
            branch_reoptimization_policy
        )
    if local_search_method == "spr":
        return validate_nucleotide_likelihood_spr_branch_reoptimization_policy(
            branch_reoptimization_policy
        )
    return validate_nucleotide_likelihood_tbr_branch_reoptimization_policy(
        branch_reoptimization_policy
    )


def validate_nucleotide_likelihood_ratchet_evaluation_budget(
    evaluation_budget: int | None,
    *,
    local_search_method: str,
) -> int | None:
    """Validate one optional ratchet evaluation budget."""
    if local_search_method != "spr":
        if evaluation_budget is not None:
            raise ValueError(
                "evaluation_budget is supported only when local_search_method is 'spr'"
            )
        return None
    return validate_likelihood_spr_evaluation_budget(evaluation_budget)


def build_nucleotide_likelihood_ratchet_records(
    records: list[AlignmentRecord],
    *,
    reweighted_site_positions: list[int],
    perturbation_factor: int,
) -> list[AlignmentRecord]:
    """Duplicate selected aligned site columns to realize temporary integer weights."""
    selected_positions = set(reweighted_site_positions)
    reweighted_records: list[AlignmentRecord] = []
    for record in records:
        sequence_fragments: list[str] = []
        for site_position, state in enumerate(record.sequence, start=1):
            sequence_fragments.append(state)
            if site_position in selected_positions:
                sequence_fragments.extend(state for _ in range(perturbation_factor - 1))
        reweighted_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence="".join(sequence_fragments),
            )
        )
    return reweighted_records


def _run_nucleotide_likelihood_ratchet_local_search(
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
    root_prior_policy: str | None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None,
    fixed_root_state: str | None,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> (
    NucleotideLikelihoodNniSearchReport
    | NucleotideLikelihoodSprSearchReport
    | NucleotideLikelihoodTbrSearchReport
):
    """Run one native local tree search for one ratchet cycle phase."""
    if local_search_method == "nni":
        return search_nucleotide_likelihood_nni(
            tree,
            records,
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
            max_coordinate_passes=max_coordinate_passes,
        )
    if local_search_method == "spr":
        return search_nucleotide_likelihood_spr(
            tree,
            records,
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
    return search_nucleotide_likelihood_tbr(
        tree,
        records,
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
        max_coordinate_passes=max_coordinate_passes,
    )


def write_nucleotide_likelihood_ratchet_cycle_table(
    path: Path,
    report: NucleotideLikelihoodRatchetReport,
) -> Path:
    """Write one deterministic likelihood ratchet cycle ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "cycle_index",
        "start_log_likelihood",
        "start_tree_newick",
        "reweighted_site_positions",
        "temporary_site_weights",
        "perturbation_factor",
        "perturbed_alignment_length",
        "perturbed_pattern_count",
        "perturbed_search_algorithm",
        "perturbed_score",
        "perturbed_tree_newick",
        "perturbed_accepted_move_count",
        "perturbed_evaluated_neighbor_count",
        "perturbed_stopping_reason",
        "restored_search_algorithm",
        "restored_score",
        "restored_tree_newick",
        "restored_accepted_move_count",
        "restored_evaluated_neighbor_count",
        "restored_stopping_reason",
        "best_score_after_cycle",
        "best_tree_after_cycle",
        "best_tree_improved",
    ]
    rows = ["\t".join(columns)]
    for row in report.cycle_rows:
        payload = [
            row.cycle_index,
            row.start_log_likelihood,
            row.start_tree_newick,
            ",".join(str(position) for position in row.reweighted_site_positions),
            ",".join(
                f"{position}:{weight}"
                for position, weight in sorted(row.temporary_site_weights.items())
            ),
            row.perturbation_factor,
            row.perturbed_alignment_length,
            row.perturbed_pattern_count,
            row.perturbed_search_algorithm,
            row.perturbed_score,
            row.perturbed_tree_newick,
            row.perturbed_accepted_move_count,
            row.perturbed_evaluated_neighbor_count,
            row.perturbed_stopping_reason,
            row.restored_search_algorithm,
            row.restored_score,
            row.restored_tree_newick,
            row.restored_accepted_move_count,
            row.restored_evaluated_neighbor_count,
            row.restored_stopping_reason,
            row.best_score_after_cycle,
            row.best_tree_after_cycle,
            str(row.best_tree_improved).lower(),
        ]
        rows.append(
            "\t".join(
                format(value, ".15g") if isinstance(value, float) else str(value)
                for value in payload
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_ratchet_best_tree_history_table(
    path: Path,
    report: NucleotideLikelihoodRatchetReport,
) -> Path:
    """Write one deterministic likelihood ratchet best-tree history ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "history_index",
        "cycle_index",
        "best_log_likelihood",
        "best_topology_fingerprint",
        "best_tree_newick",
    ]
    rows = ["\t".join(columns)]
    for row in report.best_tree_history_rows:
        payload = [
            row.history_index,
            row.cycle_index,
            row.best_log_likelihood,
            row.best_topology_fingerprint,
            row.best_tree_newick,
        ]
        rows.append(
            "\t".join(
                format(value, ".15g") if isinstance(value, float) else str(value)
                for value in payload
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_ratchet_run_json(
    path: Path,
    report: NucleotideLikelihoodRatchetReport,
) -> Path:
    """Write one machine-readable likelihood ratchet payload."""
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
        "cycle_count": report.cycle_count,
        "perturbation_seed": report.perturbation_seed,
        "perturbed_site_count": report.perturbed_site_count,
        "perturbation_factor": report.perturbation_factor,
        "input_tree_newick": report.input_tree_newick,
        "start_tree_newick": report.start_tree_newick,
        "start_log_likelihood": report.start_log_likelihood,
        "final_tree_newick": report.final_tree_newick,
        "final_log_likelihood": report.final_log_likelihood,
        "best_tree_newick": report.best_tree_newick,
        "best_log_likelihood": report.best_log_likelihood,
        "branch_reoptimization_policy": report.branch_reoptimization_policy,
        "evaluation_budget": report.evaluation_budget,
        "substitution_parameter_policy": report.substitution_parameter_policy,
        "substitution_parameter_values": report.substitution_parameter_values,
        "substitution_parameter_warnings": report.substitution_parameter_warnings,
        "cycle_rows": [
            {
                "cycle_index": row.cycle_index,
                "start_log_likelihood": row.start_log_likelihood,
                "start_tree_newick": row.start_tree_newick,
                "reweighted_site_positions": row.reweighted_site_positions,
                "temporary_site_weights": row.temporary_site_weights,
                "perturbation_factor": row.perturbation_factor,
                "perturbed_alignment_length": row.perturbed_alignment_length,
                "perturbed_pattern_count": row.perturbed_pattern_count,
                "perturbed_search_algorithm": row.perturbed_search_algorithm,
                "perturbed_score": row.perturbed_score,
                "perturbed_tree_newick": row.perturbed_tree_newick,
                "perturbed_accepted_move_count": row.perturbed_accepted_move_count,
                "perturbed_evaluated_neighbor_count": row.perturbed_evaluated_neighbor_count,
                "perturbed_stopping_reason": row.perturbed_stopping_reason,
                "restored_search_algorithm": row.restored_search_algorithm,
                "restored_score": row.restored_score,
                "restored_tree_newick": row.restored_tree_newick,
                "restored_accepted_move_count": row.restored_accepted_move_count,
                "restored_evaluated_neighbor_count": row.restored_evaluated_neighbor_count,
                "restored_stopping_reason": row.restored_stopping_reason,
                "best_score_after_cycle": row.best_score_after_cycle,
                "best_tree_after_cycle": row.best_tree_after_cycle,
                "best_tree_improved": row.best_tree_improved,
            }
            for row in report.cycle_rows
        ],
        "best_tree_history_rows": [
            {
                "history_index": row.history_index,
                "cycle_index": row.cycle_index,
                "best_log_likelihood": row.best_log_likelihood,
                "best_tree_newick": row.best_tree_newick,
                "best_topology_fingerprint": row.best_topology_fingerprint,
            }
            for row in report.best_tree_history_rows
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_ratchet_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodRatchetReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one likelihood ratchet run."""
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
    cycle_table_path = write_nucleotide_likelihood_ratchet_cycle_table(
        out_dir / "cycle_table.tsv",
        report,
    )
    best_tree_history_path = write_nucleotide_likelihood_ratchet_best_tree_history_table(
        out_dir / "best_tree_history.tsv",
        report,
    )
    run_json_path = write_nucleotide_likelihood_ratchet_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "input_tree_path": input_tree_path,
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "cycle_table_path": cycle_table_path,
        "best_tree_history_path": best_tree_history_path,
        "run_json_path": run_json_path,
    }
