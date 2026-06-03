from __future__ import annotations

from json import dumps as json_dumps
from pathlib import Path
from random import Random

import numpy

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodNniSearchReport,
    NucleotideLikelihoodSprSearchReport,
    NucleotideLikelihoodStochasticTopologyPerturbationSearchReport,
    NucleotideLikelihoodTopologyPerturbationStep,
)
from bijux_phylogenetics.phylo.likelihood.nni_search import (
    search_nucleotide_likelihood_nni,
    validate_nucleotide_likelihood_nni_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.spr_search import (
    search_nucleotide_likelihood_spr,
    validate_nucleotide_likelihood_spr_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    initialize_generated_nucleotide_topology_search_tree,
    normalize_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_tree,
    validate_nucleotide_topology_search_tree,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
)
from bijux_phylogenetics.phylo.topology.rooted_spr import (
    apply_rooted_spr_move,
    iter_rooted_spr_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_PERTURBATION_MOVE_FAMILIES = frozenset({"nni", "spr"})
_SUPPORTED_LOCAL_SEARCH_METHODS = frozenset({"nni", "spr"})


def search_nucleotide_likelihood_stochastic_topology_perturbation(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    perturbation_move_family: str = "nni",
    local_search_method: str = "nni",
    perturbation_seed: int = 1,
    perturbation_move_count: int = 1,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodStochasticTopologyPerturbationSearchReport:
    """Perturb topology by random legal rooted moves before one native local search."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = (
        resolve_nucleotide_topology_search_records(records)
    )
    validated_move_family = validate_nucleotide_likelihood_perturbation_move_family(
        perturbation_move_family
    )
    validated_local_search_method = validate_nucleotide_likelihood_local_search_method(
        local_search_method
    )
    validated_perturbation_move_count = (
        validate_nucleotide_likelihood_perturbation_move_count(perturbation_move_count)
    )
    validated_branch_reoptimization_policy = (
        validate_nucleotide_likelihood_perturbation_branch_reoptimization_policy(
            branch_reoptimization_policy,
            local_search_method=validated_local_search_method,
        )
    )
    validate_nucleotide_topology_search_tree(
        resolved_tree,
        workflow_name="nucleotide likelihood stochastic topology perturbation search",
    )
    normalized_records, compressed_patterns = (
        normalize_nucleotide_topology_search_records(
            resolved_records,
            owner_name="nucleotide likelihood stochastic topology perturbation search",
        )
    )
    rng = Random(perturbation_seed)  # nosec B311
    perturbed_tree = resolved_tree.copy().refresh()
    perturbation_steps: list[NucleotideLikelihoodTopologyPerturbationStep] = []
    for step_index in range(1, validated_perturbation_move_count + 1):
        perturbed_tree, step = apply_random_nucleotide_likelihood_topology_perturbation(
            perturbed_tree,
            move_family=validated_move_family,
            step_index=step_index,
            rng=rng,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
        )
        perturbation_steps.append(step)
    local_search_report = run_nucleotide_likelihood_topology_perturbation_local_search(
        perturbed_tree,
        normalized_records,
        model_name=model_name,
        local_search_method=validated_local_search_method,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
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
    return NucleotideLikelihoodStochasticTopologyPerturbationSearchReport(
        algorithm="nucleotide-likelihood-stochastic-topology-perturbation-search",
        model_name=local_search_report.model_name,
        perturbation_move_family=validated_move_family,
        local_search_method=validated_local_search_method,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        tree_path=None if resolved_tree_path is None else str(resolved_tree_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        perturbation_seed=perturbation_seed,
        perturbation_move_count_requested=validated_perturbation_move_count,
        perturbation_move_count_applied=len(perturbation_steps),
        input_tree_newick=dumps_newick(resolved_tree),
        perturbed_tree_newick=dumps_newick(perturbed_tree),
        perturbed_topology_fingerprint=rooted_topology_fingerprint(perturbed_tree),
        local_search_algorithm=local_search_report.algorithm,
        local_search_start_tree_newick=local_search_report.start_tree_newick,
        local_search_start_log_likelihood=local_search_report.start_log_likelihood,
        final_tree_newick=local_search_report.final_tree_newick,
        final_log_likelihood=local_search_report.final_log_likelihood,
        final_topology_fingerprint=rooted_topology_fingerprint(
            loads_newick(local_search_report.final_tree_newick)
        ),
        local_search_accepted_move_count=local_search_report.accepted_move_count,
        local_search_evaluated_neighbor_count=local_search_report.evaluated_neighbor_count,
        local_search_stopping_reason=local_search_report.stopping_reason,
        substitution_parameter_policy=local_search_report.substitution_parameter_policy,
        substitution_parameter_values=dict(
            local_search_report.substitution_parameter_values
        ),
        substitution_parameter_warnings=list(
            local_search_report.substitution_parameter_warnings
        ),
        perturbation_steps=perturbation_steps,
    )


def search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    perturbation_move_family: str = "nni",
    local_search_method: str = "nni",
    perturbation_seed: int = 1,
    perturbation_move_count: int = 1,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> NucleotideLikelihoodStochasticTopologyPerturbationSearchReport:
    """Run one path-based stochastic topology perturbation search from file inputs."""
    return search_nucleotide_likelihood_stochastic_topology_perturbation(
        tree_path,
        alignment_path,
        model_name=model_name,
        perturbation_move_family=perturbation_move_family,
        local_search_method=local_search_method,
        perturbation_seed=perturbation_seed,
        perturbation_move_count=perturbation_move_count,
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


def validate_nucleotide_likelihood_perturbation_move_family(move_family: str) -> str:
    """Validate one random topology-perturbation move family."""
    normalized_move_family = move_family.strip().lower()
    if normalized_move_family not in _SUPPORTED_PERTURBATION_MOVE_FAMILIES:
        raise ValueError(
            "perturbation_move_family must be one of "
            + ", ".join(sorted(_SUPPORTED_PERTURBATION_MOVE_FAMILIES))
        )
    return normalized_move_family


def validate_nucleotide_likelihood_local_search_method(local_search_method: str) -> str:
    """Validate one local search method used after topology perturbation."""
    normalized_local_search_method = local_search_method.strip().lower()
    if normalized_local_search_method not in _SUPPORTED_LOCAL_SEARCH_METHODS:
        raise ValueError(
            "local_search_method must be one of "
            + ", ".join(sorted(_SUPPORTED_LOCAL_SEARCH_METHODS))
        )
    return normalized_local_search_method


def validate_nucleotide_likelihood_perturbation_move_count(
    perturbation_move_count: int,
) -> int:
    """Validate the number of random legal topology perturbations."""
    if perturbation_move_count < 1:
        raise ValueError("perturbation_move_count must be at least one")
    return perturbation_move_count


def validate_nucleotide_likelihood_perturbation_branch_reoptimization_policy(
    branch_reoptimization_policy: str,
    *,
    local_search_method: str,
) -> str:
    """Validate one local-search branch reoptimization policy by chosen method."""
    if local_search_method == "nni":
        return validate_nucleotide_likelihood_nni_branch_reoptimization_policy(
            branch_reoptimization_policy
        )
    return validate_nucleotide_likelihood_spr_branch_reoptimization_policy(
        branch_reoptimization_policy
    )


def apply_random_nucleotide_likelihood_topology_perturbation(
    tree: PhyloTree,
    *,
    move_family: str,
    step_index: int,
    rng: Random,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
) -> tuple[PhyloTree, NucleotideLikelihoodTopologyPerturbationStep]:
    """Apply one random legal rooted topology perturbation and record its path row."""
    if move_family == "nni":
        return _apply_random_nni_perturbation(
            tree,
            step_index=step_index,
            rng=rng,
        )
    if move_family == "spr":
        return _apply_random_spr_perturbation(
            tree,
            step_index=step_index,
            rng=rng,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
        )
    raise ValueError(f"unsupported perturbation move family '{move_family}'")


def run_nucleotide_likelihood_topology_perturbation_local_search(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    local_search_method: str,
    branch_reoptimization_policy: str,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None,
    fixed_root_state: str | None,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> NucleotideLikelihoodNniSearchReport | NucleotideLikelihoodSprSearchReport:
    """Run one native local search from a pre-perturbed topology."""
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
    return search_nucleotide_likelihood_spr(
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


def _apply_random_nni_perturbation(
    tree: PhyloTree,
    *,
    step_index: int,
    rng: Random,
) -> tuple[PhyloTree, NucleotideLikelihoodTopologyPerturbationStep]:
    candidates = [
        candidate
        for candidate in iter_rooted_nni_move_candidates(tree)
        if rooted_topology_fingerprint(apply_rooted_nni_move(tree, candidate))
        != rooted_topology_fingerprint(tree)
    ]
    if not candidates:
        raise ValueError(
            "nucleotide likelihood stochastic topology perturbation search "
            "requires at least one topology-changing rooted NNI perturbation move"
        )
    selected_candidate = candidates[rng.randrange(len(candidates))]
    tree_before_newick = dumps_newick(tree)
    moved_tree = apply_rooted_nni_move(tree, selected_candidate)
    tree_after_newick = dumps_newick(moved_tree)
    return moved_tree, NucleotideLikelihoodTopologyPerturbationStep(
        step_index=step_index,
        move_family="nni",
        tree_before_newick=tree_before_newick,
        tree_after_newick=tree_after_newick,
        topology_fingerprint_before=rooted_topology_fingerprint(tree),
        topology_fingerprint_after=rooted_topology_fingerprint(moved_tree),
        pivot_branch_id=selected_candidate.pivot_branch_id,
        sibling_clade_id=selected_candidate.sibling_clade_id,
        exchanged_clade_id=selected_candidate.exchanged_clade_id,
        pruned_clade_id=None,
        regraft_target_branch_id=None,
    )


def _apply_random_spr_perturbation(
    tree: PhyloTree,
    *,
    step_index: int,
    rng: Random,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
) -> tuple[PhyloTree, NucleotideLikelihoodTopologyPerturbationStep]:
    topology_fingerprint_before = rooted_topology_fingerprint(tree)
    candidate_rows = [
        (
            candidate,
            initialize_generated_nucleotide_topology_search_tree(
                apply_rooted_spr_move(tree, candidate),
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
            ),
        )
        for candidate in iter_rooted_spr_move_candidates(tree)
    ]
    topology_changing_rows = [
        (candidate, moved_tree)
        for candidate, moved_tree in candidate_rows
        if rooted_topology_fingerprint(moved_tree) != topology_fingerprint_before
    ]
    if not topology_changing_rows:
        raise ValueError(
            "nucleotide likelihood stochastic topology perturbation search "
            "requires at least one topology-changing rooted SPR perturbation move"
        )
    selected_candidate, initialized_moved_tree = topology_changing_rows[
        rng.randrange(len(topology_changing_rows))
    ]
    tree_before_newick = dumps_newick(tree)
    tree_after_newick = dumps_newick(initialized_moved_tree)
    return initialized_moved_tree, NucleotideLikelihoodTopologyPerturbationStep(
        step_index=step_index,
        move_family="spr",
        tree_before_newick=tree_before_newick,
        tree_after_newick=tree_after_newick,
        topology_fingerprint_before=topology_fingerprint_before,
        topology_fingerprint_after=rooted_topology_fingerprint(initialized_moved_tree),
        pivot_branch_id=None,
        sibling_clade_id=None,
        exchanged_clade_id=None,
        pruned_clade_id=selected_candidate.pruned_clade_id,
        regraft_target_branch_id=selected_candidate.regraft_target_branch_id,
    )


def write_nucleotide_likelihood_topology_perturbation_trace_table(
    path: Path,
    report: NucleotideLikelihoodStochasticTopologyPerturbationSearchReport,
) -> Path:
    """Write one deterministic topology-perturbation path table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "step_index",
        "move_family",
        "tree_before_newick",
        "tree_after_newick",
        "topology_fingerprint_before",
        "topology_fingerprint_after",
        "pivot_branch_id",
        "sibling_clade_id",
        "exchanged_clade_id",
        "pruned_clade_id",
        "regraft_target_branch_id",
    ]
    rows = ["\t".join(columns)]
    for row in report.perturbation_steps:
        payload = [
            row.step_index,
            row.move_family,
            row.tree_before_newick,
            row.tree_after_newick,
            row.topology_fingerprint_before,
            row.topology_fingerprint_after,
            row.pivot_branch_id,
            row.sibling_clade_id,
            row.exchanged_clade_id,
            row.pruned_clade_id,
            row.regraft_target_branch_id,
        ]
        rows.append("\t".join("" if value is None else str(value) for value in payload))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_stochastic_topology_perturbation_run_json(
    path: Path,
    report: NucleotideLikelihoodStochasticTopologyPerturbationSearchReport,
) -> Path:
    """Write one machine-readable stochastic topology perturbation payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "model_name": report.model_name,
        "perturbation_move_family": report.perturbation_move_family,
        "local_search_method": report.local_search_method,
        "branch_reoptimization_policy": report.branch_reoptimization_policy,
        "tree_path": report.tree_path,
        "alignment_path": report.alignment_path,
        "taxon_count": report.taxon_count,
        "site_count": report.site_count,
        "pattern_count": report.pattern_count,
        "perturbation_seed": report.perturbation_seed,
        "perturbation_move_count_requested": report.perturbation_move_count_requested,
        "perturbation_move_count_applied": report.perturbation_move_count_applied,
        "input_tree_newick": report.input_tree_newick,
        "perturbed_tree_newick": report.perturbed_tree_newick,
        "perturbed_topology_fingerprint": report.perturbed_topology_fingerprint,
        "local_search_algorithm": report.local_search_algorithm,
        "local_search_start_tree_newick": report.local_search_start_tree_newick,
        "local_search_start_log_likelihood": report.local_search_start_log_likelihood,
        "final_tree_newick": report.final_tree_newick,
        "final_log_likelihood": report.final_log_likelihood,
        "final_topology_fingerprint": report.final_topology_fingerprint,
        "local_search_accepted_move_count": report.local_search_accepted_move_count,
        "local_search_evaluated_neighbor_count": report.local_search_evaluated_neighbor_count,
        "local_search_stopping_reason": report.local_search_stopping_reason,
        "substitution_parameter_policy": report.substitution_parameter_policy,
        "substitution_parameter_values": report.substitution_parameter_values,
        "substitution_parameter_warnings": report.substitution_parameter_warnings,
        "perturbation_steps": [
            {
                "step_index": row.step_index,
                "move_family": row.move_family,
                "tree_before_newick": row.tree_before_newick,
                "tree_after_newick": row.tree_after_newick,
                "topology_fingerprint_before": row.topology_fingerprint_before,
                "topology_fingerprint_after": row.topology_fingerprint_after,
                "pivot_branch_id": row.pivot_branch_id,
                "sibling_clade_id": row.sibling_clade_id,
                "exchanged_clade_id": row.exchanged_clade_id,
                "pruned_clade_id": row.pruned_clade_id,
                "regraft_target_branch_id": row.regraft_target_branch_id,
            }
            for row in report.perturbation_steps
        ],
    }
    path.write_text(
        json_dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_nucleotide_likelihood_stochastic_topology_perturbation_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodStochasticTopologyPerturbationSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one stochastic topology perturbation run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    input_tree_path = write_newick(
        out_dir / "input_tree.nwk",
        loads_newick(report.input_tree_newick),
    )
    perturbed_tree_path = write_newick(
        out_dir / "perturbed_tree.nwk",
        loads_newick(report.perturbed_tree_newick),
    )
    final_tree_path = write_newick(
        out_dir / "final_tree.nwk",
        loads_newick(report.final_tree_newick),
    )
    perturbation_trace_path = (
        write_nucleotide_likelihood_topology_perturbation_trace_table(
            out_dir / "perturbation_trace.tsv",
            report,
        )
    )
    run_json_path = (
        write_nucleotide_likelihood_stochastic_topology_perturbation_run_json(
            out_dir / "run.json",
            report,
        )
    )
    return {
        "input_tree_path": input_tree_path,
        "perturbed_tree_path": perturbed_tree_path,
        "final_tree_path": final_tree_path,
        "perturbation_trace_path": perturbation_trace_path,
        "run_json_path": run_json_path,
    }
