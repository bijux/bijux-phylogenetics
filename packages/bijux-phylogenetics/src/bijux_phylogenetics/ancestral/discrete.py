from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import json

from bijux_phylogenetics.ancestral.common import (
    dump_pruned_tree,
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
    write_ancestral_rows,
)
import numpy


@dataclass(slots=True)
class DiscreteAncestralEstimate:
    """One discrete ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    state_set: list[str]
    most_likely_state: str
    state_probabilities: dict[str, float]
    ambiguous: bool
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class DiscreteAncestralReport:
    """Discrete ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    minimal_change_count: int | None
    parsimonious_root_state_count: int | None
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    estimates: list[DiscreteAncestralEstimate]


@dataclass(slots=True)
class DiscreteAncestralSummary:
    """Reviewer-facing summary for one discrete ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    unstable_node_count: int
    weak_support_node_count: int
    observed_state_count: int
    sparse_state_count: int
    minimal_change_count: int | None
    parsimonious_root_state_count: int | None
    root_node: str
    root_most_likely_state: str
    root_confidence: float
    warning_count: int


@dataclass(slots=True)
class DiscreteAncestralExclusion:
    """One excluded tip from a discrete ancestral reconstruction."""

    taxon: str
    reason: str


def reconstruct_discrete_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
) -> DiscreteAncestralReport:
    """Reconstruct discrete ancestral states under Fitch or Mk likelihood models."""
    resolved_model = _resolve_discrete_model_name(model)
    if resolved_model == "fitch" and state_ordering != "unordered":
        raise ValueError(
            "ordered discrete ancestral reconstruction requires a likelihood model"
        )
    if resolved_model == "fitch" and (
        root_prior_mode != "equal" or fixed_root_state is not None
    ):
        raise ValueError(
            "fitch discrete ancestral reconstruction does not support root-prior assumptions"
        )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if resolved_model != "fitch":
        estimates, resolved_ordered_states = _reconstruct_likelihood_estimates(
            dataset,
            model=resolved_model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
        )
        unstable_nodes = [
            estimate.node
            for estimate in estimates
            if estimate.unstable and not estimate.is_tip
        ]
        weak_support_nodes = [
            estimate.node
            for estimate in estimates
            if not estimate.is_tip and estimate.confidence < 0.75
        ]
        warnings = list(dataset.warnings)
        if unstable_nodes:
            warnings.append(
                "one or more discrete ancestral nodes remain unstable across candidate states"
            )
        if weak_support_nodes:
            warnings.append(
                "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions"
            )
        return DiscreteAncestralReport(
            tree_path=tree_path,
            traits_path=traits_path,
            taxon_column=dataset.taxon_column,
            trait=trait,
            model=resolved_model,
            state_ordering=state_ordering,
            ordered_states=resolved_ordered_states,
            taxon_count=len(dataset.taxa),
            observed_states=dataset.observed_states,
            state_counts=dataset.state_counts,
            sparse_states=dataset.sparse_states,
            analysis_tree_newick=dump_pruned_tree(dataset.tree),
            dropped_missing_taxa=dataset.dropped_missing_taxa,
            minimal_change_count=None,
            parsimonious_root_state_count=None,
            warnings=warnings,
            unstable_nodes=unstable_nodes,
            weak_support_nodes=weak_support_nodes,
            estimates=estimates,
        )
    estimates: list[DiscreteAncestralEstimate] = []

    candidate_sets: dict[str, set[str]] = {}
    minimal_change_count = 0

    def record_candidate_sets(node) -> tuple[set[str], int]:
        if node.is_leaf():
            state = dataset.states_by_taxon[node.name]
            candidate_set = {state}
            candidate_sets[node_signature(node)] = candidate_set
            return candidate_set, 0
        child_results = [record_candidate_sets(child) for child in node.children]
        candidate = set(child_results[0][0])
        minimal_changes = sum(result[1] for result in child_results)
        for child_set, _ in child_results[1:]:
            intersection = candidate & child_set
            if intersection:
                candidate = intersection
            else:
                candidate |= child_set
                minimal_changes += 1
        candidate_sets[node_signature(node)] = candidate
        return candidate, minimal_changes

    _, minimal_change_count = record_candidate_sets(dataset.tree.root)

    for node in dataset.tree.iter_nodes():
        signature = node_signature(node)
        if node.is_leaf():
            resolved_state = dataset.states_by_taxon[node.name]
            probabilities = {resolved_state: 1.0}
            state_set = [resolved_state]
        else:
            state_set = sorted(candidate_sets[signature])
            probabilities = {state: 1.0 / len(state_set) for state in state_set}
            resolved_state = state_set[0]
        estimates.append(
            _build_discrete_estimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                state_set=state_set,
                most_likely_state=resolved_state,
                state_probabilities=probabilities,
            )
        )
    unstable_nodes = [
        estimate.node
        for estimate in estimates
        if estimate.unstable and not estimate.is_tip
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more discrete ancestral nodes remain unstable across candidate states"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions"
        )

    return DiscreteAncestralReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=resolved_model,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        taxon_count=len(dataset.taxa),
        observed_states=dataset.observed_states,
        state_counts=dataset.state_counts,
        sparse_states=dataset.sparse_states,
        analysis_tree_newick=dump_pruned_tree(dataset.tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        minimal_change_count=minimal_change_count,
        parsimonious_root_state_count=len(
            candidate_sets[node_signature(dataset.tree.root)]
        ),
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        estimates=estimates,
    )


def summarize_discrete_ancestral_report(
    report: DiscreteAncestralReport,
) -> DiscreteAncestralSummary:
    """Summarize the main review facts for one discrete ancestral report."""
    internal_estimates = [estimate for estimate in report.estimates if not estimate.is_tip]
    if not internal_estimates:
        raise ValueError(
            "discrete ancestral summary requires at least one internal-node estimate"
        )
    ambiguous_internal_node_count = sum(
        1 for estimate in internal_estimates if estimate.ambiguous
    )
    root_estimate = max(
        internal_estimates,
        key=lambda estimate: (
            len(estimate.descendant_taxa),
            estimate.node,
        ),
    )
    return DiscreteAncestralSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.dropped_missing_taxa),
        internal_node_count=len(internal_estimates),
        ambiguous_internal_node_count=ambiguous_internal_node_count,
        unstable_node_count=len(report.unstable_nodes),
        weak_support_node_count=len(report.weak_support_nodes),
        observed_state_count=len(report.observed_states),
        sparse_state_count=len(report.sparse_states),
        minimal_change_count=report.minimal_change_count,
        parsimonious_root_state_count=report.parsimonious_root_state_count,
        root_node=root_estimate.node,
        root_most_likely_state=root_estimate.most_likely_state,
        root_confidence=root_estimate.confidence,
        warning_count=len(report.warnings),
    )


def discrete_ancestral_exclusions(
    report: DiscreteAncestralReport,
) -> list[DiscreteAncestralExclusion]:
    """Return one explicit exclusion row per dropped tip taxon."""
    return [
        DiscreteAncestralExclusion(
            taxon=taxon,
            reason="missing_discrete_trait_state",
        )
        for taxon in report.dropped_missing_taxa
    ]


def write_discrete_ancestral_summary_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one summary ledger for a discrete ancestral reconstruction."""
    summary = summarize_discrete_ancestral_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "unstable_node_count",
            "weak_support_node_count",
            "observed_state_count",
            "sparse_state_count",
            "minimal_change_count",
            "parsimonious_root_state_count",
            "root_node",
            "root_most_likely_state",
            "root_confidence",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "unstable_node_count": str(summary.unstable_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "observed_state_count": str(summary.observed_state_count),
                "sparse_state_count": str(summary.sparse_state_count),
                "minimal_change_count": _format_optional_int(
                    summary.minimal_change_count
                ),
                "parsimonious_root_state_count": _format_optional_int(
                    summary.parsimonious_root_state_count
                ),
                "root_node": summary.root_node,
                "root_most_likely_state": summary.root_most_likely_state,
                "root_confidence": str(summary.root_confidence),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_discrete_ancestral_probability_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one internal-node marginal-probability ledger for a discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_state",
            "state_set",
            "state_probabilities",
            "confidence",
            "ambiguous",
            "unstable",
            "interpretation",
        ],
        rows=[
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "most_likely_state": estimate.most_likely_state,
                "state_set": ",".join(estimate.state_set),
                "state_probabilities": json.dumps(
                    estimate.state_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(estimate.confidence),
                "ambiguous": str(estimate.ambiguous).lower(),
                "unstable": str(estimate.unstable).lower(),
                "interpretation": estimate.interpretation,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def write_discrete_ancestral_exclusion_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one explicit excluded-tip ledger for a discrete reconstruction."""
    exclusions = discrete_ancestral_exclusions(report)
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in exclusions
        ],
    )


def _resolve_discrete_model_name(model: str) -> str:
    aliases = {
        "fitch": "fitch",
        "equal-rates": "equal-rates",
        "er": "equal-rates",
        "symmetric": "symmetric",
        "sym": "symmetric",
        "all-rates-different": "all-rates-different",
        "ard": "all-rates-different",
    }
    resolved = aliases.get(model)
    if resolved is None:
        raise ValueError(f"unsupported discrete ancestral model: {model}")
    return resolved


def _reconstruct_likelihood_estimates(
    dataset,
    *,
    model: str,
    state_ordering: str,
    ordered_states: list[str] | None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
) -> tuple[list[DiscreteAncestralEstimate], list[str]]:
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    rate_matrix, default_root_prior = _fit_discrete_mk_model(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
    )
    root_prior = _resolve_root_prior(
        state_order,
        state_counts=dataset.state_counts,
        mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
        default_root_prior=default_root_prior,
    )
    posterior_by_node = _estimate_marginal_state_probabilities(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    estimates: list[DiscreteAncestralEstimate] = []
    for node in dataset.tree.iter_nodes():
        signature = node_signature(node)
        probabilities = posterior_by_node[signature]
        material_states = _material_state_set(probabilities)
        estimates.append(
            _build_discrete_estimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                state_set=material_states,
                most_likely_state=max(
                    sorted(probabilities),
                    key=lambda state: probabilities[state],
                ),
                state_probabilities=probabilities,
            )
        )
    return estimates, (state_order if state_ordering == "ordered" else [])


def _resolve_state_order(
    observed_states: list[str],
    *,
    state_ordering: str,
    ordered_states: list[str] | None,
) -> list[str]:
    if state_ordering == "unordered":
        return list(observed_states)
    if ordered_states is None:
        return list(observed_states)
    missing_states = sorted(set(observed_states) - set(ordered_states))
    if missing_states:
        raise ValueError(
            "ordered discrete ancestral reconstruction is missing observed states: "
            + ", ".join(missing_states)
        )
    return list(ordered_states)


def _fit_discrete_mk_model(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    parameter_count = _parameter_count(
        len(state_order),
        model=model,
        state_ordering=state_ordering,
    )
    initial_candidates = [
        numpy.full(parameter_count, math.log(scale), dtype=float)
        for scale in (0.1, 1.0, 3.0)
    ]
    best_log_parameters: numpy.ndarray | None = None
    best_log_likelihood = float("-inf")
    for initial in initial_candidates:
        candidate, candidate_score = _optimize_log_parameters(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            initial_log_parameters=initial,
        )
        if candidate_score > best_log_likelihood:
            best_log_parameters = candidate
            best_log_likelihood = candidate_score
    assert best_log_parameters is not None
    rate_matrix = _rate_matrix_from_log_parameters(
        best_log_parameters,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
    )
    root_prior = _uniform_root_prior(len(state_order))
    return rate_matrix, root_prior


def _optimize_log_parameters(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    initial_log_parameters: numpy.ndarray,
) -> tuple[numpy.ndarray, float]:
    simplex = [numpy.clip(initial_log_parameters.copy(), -10.0, 5.0)]
    for index in range(initial_log_parameters.size):
        vertex = simplex[0].copy()
        vertex[index] += 0.75
        simplex.append(numpy.clip(vertex, -10.0, 5.0))
    scores = [
        _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            log_parameters=vertex,
        )
        for vertex in simplex
    ]
    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5
    for _ in range(240):
        ordering = sorted(
            range(len(simplex)),
            key=lambda index: scores[index],
            reverse=True,
        )
        simplex = [simplex[index] for index in ordering]
        scores = [scores[index] for index in ordering]
        if (
            max(numpy.linalg.norm(vertex - simplex[0]) for vertex in simplex[1:]) < 1e-3
            and max(abs(score - scores[0]) for score in scores[1:]) < 1e-5
        ):
            break
        centroid = numpy.mean(simplex[:-1], axis=0)
        reflected = numpy.clip(centroid + alpha * (centroid - simplex[-1]), -10.0, 5.0)
        reflected_score = _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            log_parameters=reflected,
        )
        if scores[0] >= reflected_score > scores[-2]:
            simplex[-1] = reflected
            scores[-1] = reflected_score
            continue
        if reflected_score > scores[0]:
            expanded = numpy.clip(
                centroid + gamma * (reflected - centroid),
                -10.0,
                5.0,
            )
            expanded_score = _evaluate_log_likelihood(
                tree,
                states_by_taxon,
                state_order=state_order,
                model=model,
                state_ordering=state_ordering,
                log_parameters=expanded,
            )
            if expanded_score > reflected_score:
                simplex[-1] = expanded
                scores[-1] = expanded_score
            else:
                simplex[-1] = reflected
                scores[-1] = reflected_score
            continue
        if reflected_score > scores[-1]:
            contracted = numpy.clip(
                centroid + rho * (reflected - centroid),
                -10.0,
                5.0,
            )
        else:
            contracted = numpy.clip(
                centroid + rho * (simplex[-1] - centroid),
                -10.0,
                5.0,
            )
        contracted_score = _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            log_parameters=contracted,
        )
        if contracted_score > scores[-1]:
            simplex[-1] = contracted
            scores[-1] = contracted_score
            continue
        best_vertex = simplex[0]
        new_simplex = [best_vertex]
        new_scores = [scores[0]]
        for vertex in simplex[1:]:
            shrunk = numpy.clip(
                best_vertex + sigma * (vertex - best_vertex),
                -10.0,
                5.0,
            )
            new_simplex.append(shrunk)
            new_scores.append(
                _evaluate_log_likelihood(
                    tree,
                    states_by_taxon,
                    state_order=state_order,
                    model=model,
                    state_ordering=state_ordering,
                    log_parameters=shrunk,
                )
            )
        simplex = new_simplex
        scores = new_scores
    best_index = max(range(len(scores)), key=lambda index: scores[index])
    return simplex[best_index], scores[best_index]


def _evaluate_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    log_parameters: numpy.ndarray,
) -> float:
    rate_matrix = _rate_matrix_from_log_parameters(
        log_parameters,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
    )
    return _tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=_uniform_root_prior(len(state_order)),
    )


def _parameter_count(
    state_count: int,
    *,
    model: str,
    state_ordering: str,
) -> int:
    if model == "equal-rates":
        return 1
    if state_ordering == "ordered":
        edge_count = max(state_count - 1, 0)
        if model == "symmetric":
            return edge_count
        return edge_count * 2
    if model == "symmetric":
        return state_count * max(state_count - 1, 0) // 2
    return state_count * max(state_count - 1, 0)


def _rate_matrix_from_log_parameters(
    log_parameters: numpy.ndarray,
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
) -> numpy.ndarray:
    state_count = len(state_order)
    rate_matrix = numpy.zeros((state_count, state_count), dtype=float)
    parameter_index = 0
    if model == "equal-rates":
        rate = math.exp(float(log_parameters[0]))
        for left_index in range(state_count):
            for right_index in range(state_count):
                if _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                ):
                    rate_matrix[left_index, right_index] = rate
    elif model == "symmetric":
        for left_index in range(state_count):
            for right_index in range(left_index + 1, state_count):
                if not _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                ):
                    continue
                rate = math.exp(float(log_parameters[parameter_index]))
                parameter_index += 1
                rate_matrix[left_index, right_index] = rate
                rate_matrix[right_index, left_index] = rate
    else:
        for left_index in range(state_count):
            for right_index in range(state_count):
                if not _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                ):
                    continue
                rate_matrix[left_index, right_index] = math.exp(
                    float(log_parameters[parameter_index])
                )
                parameter_index += 1
    for state_index in range(state_count):
        rate_matrix[state_index, state_index] = -float(
            numpy.sum(rate_matrix[state_index, :])
        )
    return rate_matrix


def _transition_allowed(
    left_index: int,
    right_index: int,
    *,
    state_count: int,
    state_ordering: str,
) -> bool:
    if left_index == right_index:
        return False
    if state_ordering == "unordered":
        return True
    return abs(left_index - right_index) == 1 and max(left_index, right_index) < state_count


def _tree_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> float:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = _transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def visit(node) -> tuple[numpy.ndarray, float]:
        if node.is_leaf():
            likelihood = numpy.zeros(len(state_order), dtype=float)
            likelihood[state_index[states_by_taxon[node.name]]] = 1.0
            return likelihood, 0.0
        partial = numpy.ones(len(state_order), dtype=float)
        log_scale = 0.0
        for child in node.children:
            child_partial, child_scale = visit(child)
            branch_transition = transition(_branch_length(child))
            partial *= branch_transition @ child_partial
            log_scale += child_scale
        scale = float(partial.sum())
        if scale <= 0.0:
            return partial, float("-inf")
        partial /= scale
        return partial, log_scale + math.log(scale)

    root_partial, subtree_log_scale = visit(tree.root)
    root_weight = root_prior * root_partial
    root_scale = float(root_weight.sum())
    if root_scale <= 0.0:
        return float("-inf")
    return subtree_log_scale + math.log(root_scale)


def _estimate_marginal_state_probabilities(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    state_index = {state: index for index, state in enumerate(state_order)}
    partial_by_node: dict[str, numpy.ndarray] = {}
    child_contributions: dict[str, dict[str, numpy.ndarray]] = {}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = _transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def postorder(node) -> numpy.ndarray:
        signature = node_signature(node)
        if node.is_leaf():
            partial = numpy.zeros(len(state_order), dtype=float)
            partial[state_index[states_by_taxon[node.name]]] = 1.0
            partial_by_node[signature] = partial
            child_contributions[signature] = {}
            return partial
        partial = numpy.ones(len(state_order), dtype=float)
        contribution_by_child: dict[str, numpy.ndarray] = {}
        for child in node.children:
            child_partial = postorder(child)
            contribution = transition(_branch_length(child)) @ child_partial
            contribution_by_child[node_signature(child)] = contribution
            partial *= contribution
        scale = float(partial.sum())
        if scale > 0.0:
            partial /= scale
        partial_by_node[signature] = partial
        child_contributions[signature] = contribution_by_child
        return partial

    postorder(tree.root)
    posterior_by_node: dict[str, numpy.ndarray] = {}
    root_signature = node_signature(tree.root)
    posterior_by_node[root_signature] = _normalize_array(
        root_prior * partial_by_node[root_signature]
    )

    def preorder(node, down_message: numpy.ndarray) -> None:
        parent_signature = node_signature(node)
        if not node.is_leaf():
            for child in node.children:
                sibling_support = down_message.copy()
                child_signature = node_signature(child)
                for sibling in node.children:
                    if sibling is child:
                        continue
                    sibling_support *= child_contributions[parent_signature][
                        node_signature(sibling)
                    ]
                branch_transition = transition(_branch_length(child))
                child_down = _normalize_array(sibling_support @ branch_transition)
                posterior_by_node[child_signature] = _normalize_array(
                    child_down * partial_by_node[child_signature]
                )
                preorder(child, child_down)

    preorder(tree.root, _normalize_array(root_prior))
    return {
        node: {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                _normalize_array(probabilities),
                strict=True,
            )
        }
        for node, probabilities in posterior_by_node.items()
    }


def _transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    if branch_length <= 0.0:
        return numpy.eye(rate_matrix.shape[0], dtype=float)
    eigenvalues, eigenvectors = numpy.linalg.eig(rate_matrix)
    inverse_vectors = numpy.linalg.inv(eigenvectors)
    diagonal = numpy.diag(numpy.exp(eigenvalues * branch_length))
    transition = eigenvectors @ diagonal @ inverse_vectors
    transition = numpy.real_if_close(transition, tol=1000).astype(float)
    transition[transition < 0.0] = 0.0
    row_sums = transition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return transition / row_sums


def _normalize_array(values: numpy.ndarray) -> numpy.ndarray:
    total = float(values.sum())
    if total <= 0.0:
        return numpy.full(values.shape[0], 1.0 / values.shape[0], dtype=float)
    return values / total


def _uniform_root_prior(state_count: int) -> numpy.ndarray:
    return numpy.full(state_count, 1.0 / state_count, dtype=float)


def _empirical_root_prior(
    state_order: list[str], state_counts: dict[str, int]
) -> numpy.ndarray:
    return _normalize_array(
        numpy.array(
            [float(state_counts.get(state, 0)) for state in state_order],
            dtype=float,
        )
    )


def _fixed_root_prior(state_order: list[str], fixed_root_state: str) -> numpy.ndarray:
    if fixed_root_state not in state_order:
        raise ValueError(
            "fixed root state is not available in the analyzed state vocabulary: "
            f"{fixed_root_state}"
        )
    prior = numpy.zeros(len(state_order), dtype=float)
    prior[state_order.index(fixed_root_state)] = 1.0
    return prior


def _resolve_root_prior(
    state_order: list[str],
    *,
    state_counts: dict[str, int],
    mode: str,
    fixed_root_state: str | None,
    default_root_prior: numpy.ndarray | None = None,
) -> numpy.ndarray:
    if mode == "equal":
        if fixed_root_state is not None:
            raise ValueError(
                "fixed_root_state requires root_prior_mode 'fixed'"
            )
        if default_root_prior is not None:
            return default_root_prior
        return _uniform_root_prior(len(state_order))
    if mode == "empirical":
        if fixed_root_state is not None:
            raise ValueError(
                "fixed_root_state requires root_prior_mode 'fixed'"
            )
        return _empirical_root_prior(state_order, state_counts)
    if mode == "fixed":
        if fixed_root_state is None:
            raise ValueError(
                "root_prior_mode 'fixed' requires a fixed_root_state"
            )
        return _fixed_root_prior(state_order, fixed_root_state)
    raise ValueError(f"unsupported discrete ancestral root prior mode: {mode}")


def _branch_length(node) -> float:
    if node.branch_length is None:
        return 1.0
    return max(float(node.branch_length), 0.0)


def _material_state_set(state_probabilities: dict[str, float]) -> list[str]:
    return sorted(
        state
        for state, probability in state_probabilities.items()
        if probability >= 0.1
    ) or [
        max(sorted(state_probabilities), key=lambda state: state_probabilities[state])
    ]


def _build_discrete_estimate(
    *,
    node: str,
    node_name: str | None,
    is_tip: bool,
    descendant_taxa: list[str],
    state_set: list[str] | None = None,
    most_likely_state: str,
    state_probabilities: dict[str, float],
) -> DiscreteAncestralEstimate:
    resolved_state_set = sorted(state_set or state_probabilities)
    ordered_probabilities = sorted(state_probabilities.values(), reverse=True)
    confidence = ordered_probabilities[0] if ordered_probabilities else 0.0
    runner_up = ordered_probabilities[1] if len(ordered_probabilities) > 1 else 0.0
    unstable = not is_tip and ((confidence - runner_up) < 0.15 or confidence < 0.7)
    if is_tip:
        interpretation = "observed tip state"
    elif unstable:
        interpretation = "unstable node state"
    elif confidence >= 0.9:
        interpretation = "strongly supported node state"
    else:
        interpretation = "moderately supported node state"
    return DiscreteAncestralEstimate(
        node=node,
        node_name=node_name,
        is_tip=is_tip,
        descendant_taxa=descendant_taxa,
        state_set=resolved_state_set,
        most_likely_state=most_likely_state,
        state_probabilities=state_probabilities,
        ambiguous=len(resolved_state_set) > 1,
        confidence=confidence,
        interpretation=interpretation,
        unstable=unstable,
        downstream_risks=_discrete_downstream_risks(unstable),
    )


def _discrete_downstream_risks(unstable: bool) -> list[str]:
    if not unstable:
        return []
    return [
        "transition counts and inferred ancestral geography may change under alternative state models",
        "biological narratives about ancestral states should be treated as provisional for this node",
    ]


def _format_optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)
