from __future__ import annotations

import json
import math
from pathlib import Path
from random import Random

import numpy

from bijux_phylogenetics.io.newick import dumps_newick, loads_newick, write_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodSimulatedAnnealingSearchReport,
    NucleotideLikelihoodSimulatedAnnealingTraceRow,
)
from bijux_phylogenetics.phylo.likelihood.nni_search import (
    resolve_nni_branch_reoptimization_scope,
    resolve_nni_local_branch_ids,
    validate_nucleotide_likelihood_nni_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.spr_search import (
    resolve_spr_affected_branch_clade_ids,
    resolve_spr_branch_reoptimization_scope,
    resolve_spr_local_branch_ids,
    validate_nucleotide_likelihood_spr_branch_reoptimization_policy,
)
from bijux_phylogenetics.phylo.likelihood.topology_search import (
    initialize_generated_nucleotide_topology_search_tree,
    normalize_nucleotide_topology_search_records,
    prefer_higher_likelihood,
    reoptimize_nucleotide_topology_tree,
    reoptimize_nucleotide_topology_tree_branch_subset,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
    resolve_nucleotide_topology_search_tree,
    resolve_reoptimized_branch_clade_ids,
    validate_nucleotide_topology_search_tree,
)
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveCandidate,
    apply_rooted_spr_move,
    iter_rooted_spr_move_candidates,
    rooted_topology_fingerprint,
)
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    RootedNniMoveCandidate,
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_SUPPORTED_ANNEALING_MOVE_FAMILIES = frozenset({"nni", "spr"})


def search_nucleotide_likelihood_simulated_annealing(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    proposal_move_family: str = "nni",
    annealing_seed: int = 1,
    iteration_count: int = 8,
    initial_temperature: float = 10.0,
    cooling_rate: float = 0.85,
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
) -> NucleotideLikelihoodSimulatedAnnealingSearchReport:
    """Search one rooted nucleotide tree by simulated annealing over legal topology moves."""
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = (
        resolve_nucleotide_topology_search_records(records)
    )
    validated_move_family = validate_nucleotide_likelihood_annealing_move_family(
        proposal_move_family
    )
    validated_iteration_count = (
        validate_nucleotide_likelihood_annealing_iteration_count(iteration_count)
    )
    validated_initial_temperature = (
        validate_nucleotide_likelihood_annealing_initial_temperature(
            initial_temperature
        )
    )
    validated_cooling_rate = validate_nucleotide_likelihood_annealing_cooling_rate(
        cooling_rate
    )
    validated_branch_reoptimization_policy = (
        validate_nucleotide_likelihood_annealing_branch_reoptimization_policy(
            branch_reoptimization_policy,
            proposal_move_family=validated_move_family,
        )
    )
    validate_nucleotide_topology_search_tree(
        resolved_tree,
        workflow_name="nucleotide likelihood simulated annealing search",
    )
    normalized_records, compressed_patterns = (
        normalize_nucleotide_topology_search_records(
            resolved_records,
            owner_name="nucleotide likelihood simulated annealing search",
        )
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
    current_tree = start_result.optimized_tree.copy().refresh()
    current_log_likelihood = start_result.log_likelihood
    current_tree_newick = dumps_newick(current_tree)
    best_tree = current_tree.copy().refresh()
    best_tree_newick = current_tree_newick
    best_log_likelihood = current_log_likelihood
    accepted_move_count = 0
    rejected_move_count = 0
    accepted_worse_move_count = 0
    total_branch_optimization_pass_count = start_result.optimization_pass_count
    total_branch_function_evaluation_count = start_result.function_evaluation_count
    trace_rows: list[NucleotideLikelihoodSimulatedAnnealingTraceRow] = []
    rng = Random(annealing_seed)  # nosec B311
    for iteration in range(1, validated_iteration_count + 1):
        temperature = validated_initial_temperature * validated_cooling_rate ** (
            iteration - 1
        )
        current_log_likelihood_before = current_log_likelihood
        current_tree_before_newick = current_tree_newick
        proposal = select_random_nucleotide_likelihood_annealing_proposal(
            current_tree,
            move_family=validated_move_family,
            rng=rng,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
        )
        proposal_result, optimized_branch_clade_ids = (
            evaluate_nucleotide_likelihood_annealing_proposal(
                current_tree,
                proposal,
                compressed_patterns=compressed_patterns,
                resolved_surface=resolved_surface,
                branch_reoptimization_policy=validated_branch_reoptimization_policy,
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
                improvement_tolerance=improvement_tolerance,
                max_coordinate_passes=max_coordinate_passes,
            )
        )
        total_branch_optimization_pass_count += proposal_result.optimization_pass_count
        total_branch_function_evaluation_count += (
            proposal_result.function_evaluation_count
        )
        proposed_tree_newick = dumps_newick(proposal_result.optimized_tree)
        log_likelihood_delta = proposal_result.log_likelihood - current_log_likelihood
        acceptance_probability = (
            resolve_nucleotide_likelihood_annealing_acceptance_probability(
                log_likelihood_delta,
                temperature=temperature,
            )
        )
        acceptance_uniform_draw = rng.random()
        accepted_move = acceptance_uniform_draw <= acceptance_probability
        acceptance_decision = resolve_nucleotide_likelihood_annealing_decision(
            log_likelihood_delta,
            accepted_move=accepted_move,
        )
        if validated_move_family == "nni":
            branch_reoptimization_scope = resolve_nni_branch_reoptimization_scope(
                validated_branch_reoptimization_policy
            )
        else:
            branch_reoptimization_scope = resolve_spr_branch_reoptimization_scope(
                validated_branch_reoptimization_policy
            )
        best_tree_improved = False
        if accepted_move:
            accepted_move_count += 1
            if log_likelihood_delta < 0.0:
                accepted_worse_move_count += 1
            current_tree = proposal_result.optimized_tree.copy().refresh()
            current_log_likelihood = proposal_result.log_likelihood
            current_tree_newick = proposed_tree_newick
            best_tree_improved = prefer_higher_likelihood(
                current_log_likelihood,
                current_tree_newick,
                best_log_likelihood,
                best_tree_newick,
            )
            if best_tree_improved:
                best_tree = current_tree.copy().refresh()
                best_tree_newick = current_tree_newick
                best_log_likelihood = current_log_likelihood
        else:
            rejected_move_count += 1
        trace_rows.append(
            NucleotideLikelihoodSimulatedAnnealingTraceRow(
                iteration=iteration,
                temperature=temperature,
                move_family=validated_move_family,
                current_log_likelihood_before=current_log_likelihood_before,
                proposed_log_likelihood=proposal_result.log_likelihood,
                log_likelihood_delta=log_likelihood_delta,
                acceptance_probability=acceptance_probability,
                acceptance_uniform_draw=acceptance_uniform_draw,
                acceptance_decision=acceptance_decision,
                accepted_move=accepted_move,
                best_tree_improved=best_tree_improved,
                current_tree_before_newick=current_tree_before_newick,
                proposed_tree_newick=proposed_tree_newick,
                current_tree_after_newick=current_tree_newick,
                pivot_branch_id=proposal.pivot_branch_id,
                sibling_clade_id=proposal.sibling_clade_id,
                exchanged_clade_id=proposal.exchanged_clade_id,
                pruned_clade_id=proposal.pruned_clade_id,
                regraft_target_branch_id=proposal.regraft_target_branch_id,
                branch_reoptimization_policy=validated_branch_reoptimization_policy,
                branch_reoptimization_scope=branch_reoptimization_scope,
                optimized_branch_count=len(proposal_result.optimized_branch_ids),
                optimized_branch_clade_ids=optimized_branch_clade_ids,
                branch_reoptimization_converged=bool(proposal_result.converged),
                branch_optimization_pass_count=proposal_result.optimization_pass_count,
                branch_function_evaluation_count=proposal_result.function_evaluation_count,
                boundary_warning_messages=list(
                    proposal_result.boundary_warning_messages
                ),
            )
        )
    return NucleotideLikelihoodSimulatedAnnealingSearchReport(
        algorithm="nucleotide-likelihood-simulated-annealing-search",
        model_name=resolved_surface.model_name,
        proposal_move_family=validated_move_family,
        branch_reoptimization_policy=validated_branch_reoptimization_policy,
        tree_path=None if resolved_tree_path is None else str(resolved_tree_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        annealing_seed=annealing_seed,
        iteration_count_requested=validated_iteration_count,
        iteration_count_completed=len(trace_rows),
        initial_temperature=validated_initial_temperature,
        cooling_rate=validated_cooling_rate,
        input_tree_newick=dumps_newick(resolved_tree),
        start_tree_newick=dumps_newick(start_result.optimized_tree),
        start_log_likelihood=start_result.log_likelihood,
        final_tree_newick=current_tree_newick,
        final_log_likelihood=current_log_likelihood,
        best_tree_newick=best_tree_newick,
        best_log_likelihood=best_log_likelihood,
        best_topology_fingerprint=rooted_topology_fingerprint(best_tree),
        accepted_move_count=accepted_move_count,
        rejected_move_count=rejected_move_count,
        accepted_worse_move_count=accepted_worse_move_count,
        substitution_parameter_policy=resolved_surface.substitution_parameter_policy,
        substitution_parameter_values=dict(
            resolved_surface.substitution_parameter_values
        ),
        substitution_parameter_warnings=list(
            resolved_surface.substitution_parameter_warnings
        ),
        total_branch_optimization_pass_count=total_branch_optimization_pass_count,
        total_branch_function_evaluation_count=total_branch_function_evaluation_count,
        stopping_reason="temperature-schedule-exhausted",
        trace_rows=trace_rows,
    )


def search_nucleotide_likelihood_simulated_annealing_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    proposal_move_family: str = "nni",
    annealing_seed: int = 1,
    iteration_count: int = 8,
    initial_temperature: float = 10.0,
    cooling_rate: float = 0.85,
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
) -> NucleotideLikelihoodSimulatedAnnealingSearchReport:
    """Run one path-based nucleotide likelihood simulated-annealing search."""
    return search_nucleotide_likelihood_simulated_annealing(
        tree_path,
        alignment_path,
        model_name=model_name,
        proposal_move_family=proposal_move_family,
        annealing_seed=annealing_seed,
        iteration_count=iteration_count,
        initial_temperature=initial_temperature,
        cooling_rate=cooling_rate,
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


def validate_nucleotide_likelihood_annealing_move_family(move_family: str) -> str:
    """Validate one topology-move family for native likelihood simulated annealing."""
    normalized_move_family = move_family.strip().lower()
    if normalized_move_family not in _SUPPORTED_ANNEALING_MOVE_FAMILIES:
        raise ValueError(
            "proposal_move_family must be one of "
            + ", ".join(sorted(_SUPPORTED_ANNEALING_MOVE_FAMILIES))
        )
    return normalized_move_family


def validate_nucleotide_likelihood_annealing_iteration_count(
    iteration_count: int,
) -> int:
    """Require at least one simulated-annealing proposal iteration."""
    if iteration_count < 1:
        raise ValueError("iteration_count must be at least one")
    return iteration_count


def validate_nucleotide_likelihood_annealing_initial_temperature(
    initial_temperature: float,
) -> float:
    """Require one strictly positive initial temperature."""
    if initial_temperature <= 0.0:
        raise ValueError("initial_temperature must be greater than zero")
    return initial_temperature


def validate_nucleotide_likelihood_annealing_cooling_rate(
    cooling_rate: float,
) -> float:
    """Require one geometric cooling rate in the open interval (0, 1]."""
    if cooling_rate <= 0.0 or cooling_rate > 1.0:
        raise ValueError("cooling_rate must be greater than zero and at most one")
    return cooling_rate


def validate_nucleotide_likelihood_annealing_branch_reoptimization_policy(
    branch_reoptimization_policy: str,
    *,
    proposal_move_family: str,
) -> str:
    """Validate one annealing branch-reoptimization policy against the move family."""
    if proposal_move_family == "nni":
        return validate_nucleotide_likelihood_nni_branch_reoptimization_policy(
            branch_reoptimization_policy
        )
    return validate_nucleotide_likelihood_spr_branch_reoptimization_policy(
        branch_reoptimization_policy
    )


class _AnnealingProposal:
    def __init__(
        self,
        *,
        move_family: str,
        tree_before_newick: str,
        proposed_tree: PhyloTree,
        branch_reoptimization_scope: str,
        pivot_branch_id: str | None = None,
        sibling_clade_id: str | None = None,
        exchanged_clade_id: str | None = None,
        pruned_clade_id: str | None = None,
        regraft_target_branch_id: str | None = None,
        nni_candidate: RootedNniMoveCandidate | None = None,
        spr_candidate: RootedSprMoveCandidate | None = None,
    ) -> None:
        self.move_family = move_family
        self.tree_before_newick = tree_before_newick
        self.proposed_tree = proposed_tree
        self.branch_reoptimization_scope = branch_reoptimization_scope
        self.pivot_branch_id = pivot_branch_id
        self.sibling_clade_id = sibling_clade_id
        self.exchanged_clade_id = exchanged_clade_id
        self.pruned_clade_id = pruned_clade_id
        self.regraft_target_branch_id = regraft_target_branch_id
        self.nni_candidate = nni_candidate
        self.spr_candidate = spr_candidate


def select_random_nucleotide_likelihood_annealing_proposal(
    tree: PhyloTree,
    *,
    move_family: str,
    rng: Random,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
) -> _AnnealingProposal:
    """Select one random topology-changing annealing proposal from the current tree."""
    if move_family == "nni":
        return _select_random_nni_annealing_proposal(tree, rng=rng)
    return _select_random_spr_annealing_proposal(
        tree,
        rng=rng,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
    )


def _select_random_nni_annealing_proposal(
    tree: PhyloTree,
    *,
    rng: Random,
) -> _AnnealingProposal:
    tree_before_newick = dumps_newick(tree)
    candidate_rows = []
    seen_topologies: set[str] = set()
    current_topology = rooted_topology_fingerprint(tree)
    for candidate in iter_rooted_nni_move_candidates(tree):
        proposed_tree = apply_rooted_nni_move(tree, candidate)
        proposal_topology = rooted_topology_fingerprint(proposed_tree)
        if (
            proposal_topology == current_topology
            or proposal_topology in seen_topologies
        ):
            continue
        seen_topologies.add(proposal_topology)
        candidate_rows.append(
            _AnnealingProposal(
                move_family="nni",
                tree_before_newick=tree_before_newick,
                proposed_tree=proposed_tree,
                branch_reoptimization_scope="candidate-dependent",
                pivot_branch_id=candidate.pivot_branch_id,
                sibling_clade_id=candidate.sibling_clade_id,
                exchanged_clade_id=candidate.exchanged_clade_id,
                nni_candidate=candidate,
            )
        )
    if not candidate_rows:
        raise ValueError(
            "nucleotide likelihood simulated annealing search requires at least one "
            "topology-changing rooted NNI proposal"
        )
    return candidate_rows[rng.randrange(len(candidate_rows))]


def _select_random_spr_annealing_proposal(
    tree: PhyloTree,
    *,
    rng: Random,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
) -> _AnnealingProposal:
    tree_before_newick = dumps_newick(tree)
    current_topology = rooted_topology_fingerprint(tree)
    candidate_rows = []
    seen_topologies: set[str] = set()
    for candidate in iter_rooted_spr_move_candidates(tree):
        proposed_tree = initialize_generated_nucleotide_topology_search_tree(
            apply_rooted_spr_move(tree, candidate),
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
        )
        proposal_topology = rooted_topology_fingerprint(proposed_tree)
        if (
            proposal_topology == current_topology
            or proposal_topology in seen_topologies
        ):
            continue
        seen_topologies.add(proposal_topology)
        candidate_rows.append(
            _AnnealingProposal(
                move_family="spr",
                tree_before_newick=tree_before_newick,
                proposed_tree=proposed_tree,
                branch_reoptimization_scope="candidate-dependent",
                pruned_clade_id=candidate.pruned_clade_id,
                regraft_target_branch_id=candidate.regraft_target_branch_id,
                spr_candidate=candidate,
            )
        )
    if not candidate_rows:
        raise ValueError(
            "nucleotide likelihood simulated annealing search requires at least one "
            "topology-changing rooted SPR proposal"
        )
    return candidate_rows[rng.randrange(len(candidate_rows))]


def evaluate_nucleotide_likelihood_annealing_proposal(
    current_tree: PhyloTree,
    proposal: _AnnealingProposal,
    *,
    compressed_patterns,
    resolved_surface,
    branch_reoptimization_policy: str,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
):
    """Evaluate one simulated-annealing topology proposal under the chosen policy."""
    if proposal.move_family == "nni":
        if branch_reoptimization_policy == "coordinate-branch-lengths":
            result = reoptimize_nucleotide_topology_tree(
                proposal.proposed_tree,
                compressed_patterns=compressed_patterns,
                resolved_surface=resolved_surface,
                branch_reoptimization_policy=branch_reoptimization_policy,
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
                improvement_tolerance=improvement_tolerance,
                max_coordinate_passes=max_coordinate_passes,
            )
        else:
            result = reoptimize_nucleotide_topology_tree_branch_subset(
                proposal.proposed_tree,
                compressed_patterns=compressed_patterns,
                resolved_surface=resolved_surface,
                optimized_branch_ids=resolve_nni_local_branch_ids(
                    current_tree,
                    proposal.proposed_tree,
                    proposal.nni_candidate,
                ),
                lower_branch_length_bound=lower_branch_length_bound,
                upper_branch_length_bound=upper_branch_length_bound,
                improvement_tolerance=improvement_tolerance,
                max_coordinate_passes=max_coordinate_passes,
            )
        return (
            result,
            resolve_reoptimized_branch_clade_ids(
                result.optimized_tree,
                result.optimized_branch_ids,
            ),
        )
    affected_branch_clade_ids = resolve_spr_affected_branch_clade_ids(
        current_tree,
        proposal.proposed_tree,
        proposal.spr_candidate,
    )
    if branch_reoptimization_policy == "coordinate-branch-lengths":
        result = reoptimize_nucleotide_topology_tree(
            proposal.proposed_tree,
            compressed_patterns=compressed_patterns,
            resolved_surface=resolved_surface,
            branch_reoptimization_policy=branch_reoptimization_policy,
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
    else:
        result = reoptimize_nucleotide_topology_tree_branch_subset(
            proposal.proposed_tree,
            compressed_patterns=compressed_patterns,
            resolved_surface=resolved_surface,
            optimized_branch_ids=resolve_spr_local_branch_ids(
                proposal.proposed_tree,
                affected_branch_clade_ids,
            ),
            lower_branch_length_bound=lower_branch_length_bound,
            upper_branch_length_bound=upper_branch_length_bound,
            improvement_tolerance=improvement_tolerance,
            max_coordinate_passes=max_coordinate_passes,
        )
    return (
        result,
        resolve_reoptimized_branch_clade_ids(
            result.optimized_tree,
            result.optimized_branch_ids,
        ),
    )


def resolve_nucleotide_likelihood_annealing_acceptance_probability(
    log_likelihood_delta: float,
    *,
    temperature: float,
) -> float:
    """Resolve one Metropolis acceptance probability for a likelihood delta."""
    if log_likelihood_delta >= 0.0:
        return 1.0
    return min(1.0, math.exp(log_likelihood_delta / temperature))


def resolve_nucleotide_likelihood_annealing_decision(
    log_likelihood_delta: float,
    *,
    accepted_move: bool,
) -> str:
    """Render one durable annealing acceptance decision label."""
    if log_likelihood_delta > 0.0:
        return "accepted-improving-move"
    if log_likelihood_delta == 0.0:
        return "accepted-equal-move"
    if accepted_move:
        return "accepted-worse-move"
    return "rejected-worse-move"


def write_nucleotide_likelihood_simulated_annealing_trace_table(
    path: Path,
    report: NucleotideLikelihoodSimulatedAnnealingSearchReport,
) -> Path:
    """Write one deterministic simulated-annealing proposal trace table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "iteration",
        "temperature",
        "move_family",
        "current_log_likelihood_before",
        "proposed_log_likelihood",
        "log_likelihood_delta",
        "acceptance_probability",
        "acceptance_uniform_draw",
        "acceptance_decision",
        "accepted_move",
        "best_tree_improved",
        "current_tree_before_newick",
        "proposed_tree_newick",
        "current_tree_after_newick",
        "pivot_branch_id",
        "sibling_clade_id",
        "exchanged_clade_id",
        "pruned_clade_id",
        "regraft_target_branch_id",
        "branch_reoptimization_policy",
        "branch_reoptimization_scope",
        "optimized_branch_count",
        "optimized_branch_clade_ids",
        "branch_reoptimization_converged",
        "branch_optimization_pass_count",
        "branch_function_evaluation_count",
        "boundary_warning_messages",
    ]
    rows = ["\t".join(columns)]
    for row in report.trace_rows:
        payload = [
            row.iteration,
            row.temperature,
            row.move_family,
            row.current_log_likelihood_before,
            row.proposed_log_likelihood,
            row.log_likelihood_delta,
            row.acceptance_probability,
            row.acceptance_uniform_draw,
            row.acceptance_decision,
            str(row.accepted_move).lower(),
            str(row.best_tree_improved).lower(),
            row.current_tree_before_newick,
            row.proposed_tree_newick,
            row.current_tree_after_newick,
            row.pivot_branch_id,
            row.sibling_clade_id,
            row.exchanged_clade_id,
            row.pruned_clade_id,
            row.regraft_target_branch_id,
            row.branch_reoptimization_policy,
            row.branch_reoptimization_scope,
            row.optimized_branch_count,
            ",".join(row.optimized_branch_clade_ids),
            str(row.branch_reoptimization_converged).lower(),
            row.branch_optimization_pass_count,
            row.branch_function_evaluation_count,
            "|".join(row.boundary_warning_messages),
        ]
        rows.append(
            "\t".join(
                ""
                if value is None
                else (format(value, ".15g") if isinstance(value, float) else str(value))
                for value in payload
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def write_nucleotide_likelihood_simulated_annealing_run_json(
    path: Path,
    report: NucleotideLikelihoodSimulatedAnnealingSearchReport,
) -> Path:
    """Write one machine-readable simulated-annealing payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "model_name": report.model_name,
        "proposal_move_family": report.proposal_move_family,
        "branch_reoptimization_policy": report.branch_reoptimization_policy,
        "tree_path": report.tree_path,
        "alignment_path": report.alignment_path,
        "taxon_count": report.taxon_count,
        "site_count": report.site_count,
        "pattern_count": report.pattern_count,
        "annealing_seed": report.annealing_seed,
        "iteration_count_requested": report.iteration_count_requested,
        "iteration_count_completed": report.iteration_count_completed,
        "initial_temperature": report.initial_temperature,
        "cooling_rate": report.cooling_rate,
        "input_tree_newick": report.input_tree_newick,
        "start_tree_newick": report.start_tree_newick,
        "start_log_likelihood": report.start_log_likelihood,
        "final_tree_newick": report.final_tree_newick,
        "final_log_likelihood": report.final_log_likelihood,
        "best_tree_newick": report.best_tree_newick,
        "best_log_likelihood": report.best_log_likelihood,
        "best_topology_fingerprint": report.best_topology_fingerprint,
        "accepted_move_count": report.accepted_move_count,
        "rejected_move_count": report.rejected_move_count,
        "accepted_worse_move_count": report.accepted_worse_move_count,
        "substitution_parameter_policy": report.substitution_parameter_policy,
        "substitution_parameter_values": report.substitution_parameter_values,
        "substitution_parameter_warnings": report.substitution_parameter_warnings,
        "total_branch_optimization_pass_count": report.total_branch_optimization_pass_count,
        "total_branch_function_evaluation_count": report.total_branch_function_evaluation_count,
        "stopping_reason": report.stopping_reason,
        "trace_rows": [
            {
                "iteration": row.iteration,
                "temperature": row.temperature,
                "move_family": row.move_family,
                "current_log_likelihood_before": row.current_log_likelihood_before,
                "proposed_log_likelihood": row.proposed_log_likelihood,
                "log_likelihood_delta": row.log_likelihood_delta,
                "acceptance_probability": row.acceptance_probability,
                "acceptance_uniform_draw": row.acceptance_uniform_draw,
                "acceptance_decision": row.acceptance_decision,
                "accepted_move": row.accepted_move,
                "best_tree_improved": row.best_tree_improved,
                "current_tree_before_newick": row.current_tree_before_newick,
                "proposed_tree_newick": row.proposed_tree_newick,
                "current_tree_after_newick": row.current_tree_after_newick,
                "pivot_branch_id": row.pivot_branch_id,
                "sibling_clade_id": row.sibling_clade_id,
                "exchanged_clade_id": row.exchanged_clade_id,
                "pruned_clade_id": row.pruned_clade_id,
                "regraft_target_branch_id": row.regraft_target_branch_id,
                "branch_reoptimization_policy": row.branch_reoptimization_policy,
                "branch_reoptimization_scope": row.branch_reoptimization_scope,
                "optimized_branch_count": row.optimized_branch_count,
                "optimized_branch_clade_ids": row.optimized_branch_clade_ids,
                "branch_reoptimization_converged": row.branch_reoptimization_converged,
                "branch_optimization_pass_count": row.branch_optimization_pass_count,
                "branch_function_evaluation_count": row.branch_function_evaluation_count,
                "boundary_warning_messages": row.boundary_warning_messages,
            }
            for row in report.trace_rows
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_nucleotide_likelihood_simulated_annealing_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodSimulatedAnnealingSearchReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one simulated-annealing search."""
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
    best_tree_path = write_newick(
        out_dir / "best_tree.nwk",
        loads_newick(report.best_tree_newick),
    )
    trace_table_path = write_nucleotide_likelihood_simulated_annealing_trace_table(
        out_dir / "search_trace.tsv",
        report,
    )
    run_json_path = write_nucleotide_likelihood_simulated_annealing_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "input_tree_path": input_tree_path,
        "start_tree_path": start_tree_path,
        "final_tree_path": final_tree_path,
        "best_tree_path": best_tree_path,
        "trace_table_path": trace_table_path,
        "run_json_path": run_json_path,
    }
