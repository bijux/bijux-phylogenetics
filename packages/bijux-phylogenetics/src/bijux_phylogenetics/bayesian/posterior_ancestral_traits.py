from __future__ import annotations

from dataclasses import dataclass
import math
from statistics import NormalDist, mean, median

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    matrix_multiply,
    matrix_vector_multiply,
    stable_covariance,
    transpose,
)
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .brownian_continuous_trait import BrownianContinuousTraitRunReport
from .discrete_trait_mk import (
    DiscreteTraitMkNodeStateSummary,
    DiscreteTraitMkRunReport,
)
from .ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitRunReport,
)
from .posterior_sets.diagnostics import highest_posterior_density_interval

_CONTINUOUS_POSTERIOR_MIXTURE_QUANTILE_COUNT = 32


@dataclass(frozen=True, slots=True)
class _PosteriorCladeMetadata:
    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _NormalizedDiscreteNodeStateSummary:
    clade_id: str
    node_id: str
    node_name: str | None
    descendant_taxa: tuple[str, ...]
    most_likely_state: str
    state_probabilities: dict[str, float]


@dataclass(frozen=True, slots=True)
class _PosteriorContinuousNodeDistribution:
    clade_id: str
    node_id: str
    node_name: str | None
    descendant_taxa: tuple[str, ...]
    conditional_mean: float
    conditional_standard_deviation: float


@dataclass(frozen=True, slots=True)
class _TreeDepthIndex:
    depth_by_node_id: dict[str, float]
    ancestor_depths_by_node_id: dict[str, dict[str, float]]


@dataclass(frozen=True, slots=True)
class PosteriorDiscreteTraitStateProbabilityRow:
    """One clade-state posterior probability aggregated across sampled states."""

    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    state: str
    clade_posterior_probability: float
    conditional_posterior_probability: float
    marginal_posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorDiscreteTraitNodeSummaryRow:
    """One clade-level discrete-trait posterior summary across sampled states."""

    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_state: str
    clade_posterior_probability: float
    max_conditional_posterior_probability: float
    max_marginal_posterior_probability: float
    conditional_posterior_entropy: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorDiscreteTraitReport:
    """Posterior ancestral discrete-trait summary across one Bayesian Mk chain."""

    sample_count: int
    distinct_topology_count: int
    sampled_transition_models: list[str]
    state_order: list[str]
    tree_uncertainty_policy: str
    state_probability_rows: list[PosteriorDiscreteTraitStateProbabilityRow]
    node_summary_rows: list[PosteriorDiscreteTraitNodeSummaryRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorContinuousTraitNodeSummaryRow:
    """One clade-level continuous-trait posterior summary across sampled states."""

    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    clade_posterior_probability: float
    conditional_posterior_mean: float
    conditional_posterior_median: float
    conditional_hpd_95_lower: float
    conditional_hpd_95_upper: float
    mean_conditional_standard_deviation: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorContinuousTraitReport:
    """Posterior ancestral continuous-trait summary across Bayesian samples."""

    sample_count: int
    distinct_topology_count: int
    sampled_trait_models: list[str]
    tree_uncertainty_policy: str
    node_summary_rows: list[PosteriorContinuousTraitNodeSummaryRow]
    warnings: list[str]


def summarize_discrete_trait_mk_posterior_ancestral_states(
    run_report: DiscreteTraitMkRunReport,
) -> PosteriorDiscreteTraitReport:
    """Summarize posterior ancestral discrete-trait states from one Mk chain."""
    if not isinstance(run_report, DiscreteTraitMkRunReport):
        raise PhylogeneticsError(
            "posterior ancestral discrete-trait summary requires one DiscreteTraitMkRunReport",
            code="posterior_ancestral_discrete_trait_run_report_type_invalid",
        )
    if not run_report.posterior_rows:
        raise PhylogeneticsError(
            "posterior ancestral discrete-trait summary requires at least one posterior row",
            code="posterior_ancestral_discrete_trait_posterior_rows_empty",
        )
    sample_count = len(run_report.posterior_rows)
    topology_ids = {row.topology_id for row in run_report.posterior_rows}
    clade_metadata_by_id: dict[str, _PosteriorCladeMetadata] = {}
    clade_presence_count_by_id: dict[str, int] = {}
    state_probability_sum_by_key: dict[tuple[str, str], float] = {}
    sampled_transition_models = {
        row.transition_model_name for row in run_report.posterior_rows
    }
    for row in run_report.posterior_rows:
        present_clade_ids: set[str] = set()
        for node_summary in row.node_state_summaries:
            normalized_summary = _normalize_discrete_node_state_summary(node_summary)
            clade_metadata_by_id.setdefault(
                normalized_summary.clade_id,
                _PosteriorCladeMetadata(
                    clade_id=normalized_summary.clade_id,
                    representative_node_id=normalized_summary.node_id,
                    node_name=normalized_summary.node_name,
                    descendant_taxa=tuple(normalized_summary.descendant_taxa),
                ),
            )
            present_clade_ids.add(normalized_summary.clade_id)
            for state in run_report.state_order:
                if state not in normalized_summary.state_probabilities:
                    raise PhylogeneticsError(
                        "posterior ancestral discrete-trait summary requires every node-state summary to report every modeled discrete state",
                        code="posterior_ancestral_discrete_trait_state_probability_missing",
                        details={
                            "clade_id": normalized_summary.clade_id,
                            "missing_state": state,
                            "reported_states": sorted(
                                normalized_summary.state_probabilities
                            ),
                        },
                    )
                accumulation_key = (normalized_summary.clade_id, state)
                state_probability_sum_by_key[accumulation_key] = float(
                    format(
                        state_probability_sum_by_key.get(accumulation_key, 0.0)
                        + normalized_summary.state_probabilities[state],
                        ".15g",
                    )
                )
        for clade_id in present_clade_ids:
            clade_presence_count_by_id[clade_id] = (
                clade_presence_count_by_id.get(clade_id, 0) + 1
            )
    state_probability_rows = _build_discrete_state_probability_rows(
        clade_metadata_by_id=clade_metadata_by_id,
        clade_presence_count_by_id=clade_presence_count_by_id,
        state_probability_sum_by_key=state_probability_sum_by_key,
        sample_count=sample_count,
        state_order=run_report.state_order,
    )
    return PosteriorDiscreteTraitReport(
        sample_count=sample_count,
        distinct_topology_count=len(topology_ids),
        sampled_transition_models=sorted(sampled_transition_models),
        state_order=list(run_report.state_order),
        tree_uncertainty_policy=_resolve_tree_uncertainty_policy(len(topology_ids)),
        state_probability_rows=state_probability_rows,
        node_summary_rows=_build_discrete_node_summary_rows(
            state_probability_rows=state_probability_rows,
            state_order=run_report.state_order,
        ),
        warnings=_build_summary_warnings(
            distinct_topology_count=len(topology_ids),
            sample_count=sample_count,
            clade_presence_count_by_id=clade_presence_count_by_id,
        ),
    )


def _normalize_discrete_node_state_summary(
    node_summary: DiscreteTraitMkNodeStateSummary,
) -> _NormalizedDiscreteNodeStateSummary:
    canonical_descendant_taxa = tuple(sorted(node_summary.descendant_taxa))
    return _NormalizedDiscreteNodeStateSummary(
        clade_id=canonical_clade_id(frozenset(canonical_descendant_taxa)),
        node_id=node_summary.node_id,
        node_name=node_summary.node_name,
        descendant_taxa=canonical_descendant_taxa,
        most_likely_state=node_summary.most_likely_state,
        state_probabilities={
            state: float(format(probability, ".15g"))
            for state, probability in sorted(node_summary.state_probabilities.items())
        },
    )


def _build_discrete_state_probability_rows(
    *,
    clade_metadata_by_id: dict[str, _PosteriorCladeMetadata],
    clade_presence_count_by_id: dict[str, int],
    state_probability_sum_by_key: dict[tuple[str, str], float],
    sample_count: int,
    state_order: list[str],
) -> list[PosteriorDiscreteTraitStateProbabilityRow]:
    rows: list[PosteriorDiscreteTraitStateProbabilityRow] = []
    for clade_id, metadata in sorted(
        clade_metadata_by_id.items(),
        key=lambda item: (len(item[1].descendant_taxa), item[1].descendant_taxa),
    ):
        supporting_sample_count = clade_presence_count_by_id[clade_id]
        clade_probability = float(
            format(supporting_sample_count / sample_count, ".15g")
        )
        for state in state_order:
            posterior_sum = state_probability_sum_by_key.get((clade_id, state), 0.0)
            conditional_probability = (
                float(format(posterior_sum / supporting_sample_count, ".15g"))
                if supporting_sample_count > 0
                else 0.0
            )
            marginal_probability = float(format(posterior_sum / sample_count, ".15g"))
            rows.append(
                PosteriorDiscreteTraitStateProbabilityRow(
                    clade_id=clade_id,
                    representative_node_id=metadata.representative_node_id,
                    node_name=metadata.node_name,
                    descendant_taxa=list(metadata.descendant_taxa),
                    state=state,
                    clade_posterior_probability=clade_probability,
                    conditional_posterior_probability=conditional_probability,
                    marginal_posterior_probability=marginal_probability,
                    supporting_sample_count=supporting_sample_count,
                    total_sample_count=sample_count,
                )
            )
    return rows


def _build_discrete_node_summary_rows(
    *,
    state_probability_rows: list[PosteriorDiscreteTraitStateProbabilityRow],
    state_order: list[str],
) -> list[PosteriorDiscreteTraitNodeSummaryRow]:
    grouped_rows: dict[str, dict[str, PosteriorDiscreteTraitStateProbabilityRow]] = {}
    for row in state_probability_rows:
        grouped_rows.setdefault(row.clade_id, {})[row.state] = row
    summary_rows: list[PosteriorDiscreteTraitNodeSummaryRow] = []
    for clade_id in sorted(grouped_rows):
        state_rows = grouped_rows[clade_id]
        ordered_rows = [state_rows[state] for state in state_order]
        best_row = max(
            ordered_rows,
            key=lambda row: row.marginal_posterior_probability,
        )
        conditional_entropy = -math.fsum(
            row.conditional_posterior_probability
            * math.log(row.conditional_posterior_probability)
            for row in ordered_rows
            if row.conditional_posterior_probability > 0.0
        )
        summary_rows.append(
            PosteriorDiscreteTraitNodeSummaryRow(
                clade_id=clade_id,
                representative_node_id=best_row.representative_node_id,
                node_name=best_row.node_name,
                descendant_taxa=list(best_row.descendant_taxa),
                most_likely_state=best_row.state,
                clade_posterior_probability=best_row.clade_posterior_probability,
                max_conditional_posterior_probability=max(
                    row.conditional_posterior_probability for row in ordered_rows
                ),
                max_marginal_posterior_probability=best_row.marginal_posterior_probability,
                conditional_posterior_entropy=float(
                    format(conditional_entropy, ".15g")
                ),
                supporting_sample_count=best_row.supporting_sample_count,
                total_sample_count=best_row.total_sample_count,
            )
        )
    return summary_rows


def _build_summary_warnings(
    *,
    distinct_topology_count: int,
    sample_count: int,
    clade_presence_count_by_id: dict[str, int],
) -> list[str]:
    warnings: list[str] = []
    if distinct_topology_count > 1:
        warnings.append(
            "posterior ancestral trait summary aggregated comparable clades across multiple sampled topologies"
        )
    if any(count < sample_count for count in clade_presence_count_by_id.values()):
        warnings.append(
            "one or more clades were absent from part of the posterior sample and are summarized conditionally on clade presence"
        )
    return warnings


def _resolve_tree_uncertainty_policy(distinct_topology_count: int) -> str:
    if distinct_topology_count <= 1:
        return "fixed-topology-posterior-aggregation"
    return "clade-marginal-posterior-aggregation-across-sampled-trees"


def summarize_brownian_continuous_trait_posterior_ancestral_states(
    run_report: BrownianContinuousTraitRunReport,
) -> PosteriorContinuousTraitReport:
    """Summarize posterior ancestral continuous traits from one Brownian chain."""
    if not isinstance(run_report, BrownianContinuousTraitRunReport):
        raise PhylogeneticsError(
            "posterior ancestral Brownian trait summary requires one BrownianContinuousTraitRunReport",
            code="posterior_ancestral_brownian_trait_run_report_type_invalid",
        )
    return _summarize_continuous_trait_posterior_states(
        sampled_states=run_report.chain_report.sampled_states,
        tip_values=run_report.tip_values,
        sampled_trait_models=[row.model_name for row in run_report.posterior_rows],
        distribution_builder=_evaluate_brownian_node_distributions,
    )


def summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states(
    run_report: OrnsteinUhlenbeckContinuousTraitRunReport,
) -> PosteriorContinuousTraitReport:
    """Summarize posterior ancestral continuous traits from one OU chain."""
    if not isinstance(run_report, OrnsteinUhlenbeckContinuousTraitRunReport):
        raise PhylogeneticsError(
            "posterior ancestral OU trait summary requires one OrnsteinUhlenbeckContinuousTraitRunReport",
            code="posterior_ancestral_ou_trait_run_report_type_invalid",
        )
    return _summarize_continuous_trait_posterior_states(
        sampled_states=run_report.chain_report.sampled_states,
        tip_values=run_report.tip_values,
        sampled_trait_models=[row.model_name for row in run_report.posterior_rows],
        distribution_builder=_evaluate_ou_node_distributions,
    )


def summarize_continuous_trait_posterior_ancestral_states(
    run_report: BrownianContinuousTraitRunReport
    | OrnsteinUhlenbeckContinuousTraitRunReport,
) -> PosteriorContinuousTraitReport:
    """Summarize posterior ancestral continuous traits from one supported chain."""
    if isinstance(run_report, BrownianContinuousTraitRunReport):
        return summarize_brownian_continuous_trait_posterior_ancestral_states(
            run_report
        )
    if isinstance(run_report, OrnsteinUhlenbeckContinuousTraitRunReport):
        return summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states(
            run_report
        )
    raise PhylogeneticsError(
        "posterior ancestral continuous-trait summary requires one BrownianContinuousTraitRunReport or OrnsteinUhlenbeckContinuousTraitRunReport",
        code="posterior_ancestral_continuous_trait_run_report_type_invalid",
    )


def _summarize_continuous_trait_posterior_states(
    *,
    sampled_states,
    tip_values: dict[str, float],
    sampled_trait_models: list[str],
    distribution_builder,
) -> PosteriorContinuousTraitReport:
    if not sampled_states:
        raise PhylogeneticsError(
            "posterior ancestral continuous-trait summary requires at least one sampled Bayesian phylogenetic state",
            code="posterior_ancestral_continuous_trait_sampled_states_empty",
        )
    sample_count = len(sampled_states)
    topology_ids: set[str] = set()
    clade_metadata_by_id: dict[str, _PosteriorCladeMetadata] = {}
    clade_presence_count_by_id: dict[str, int] = {}
    distributions_by_clade_id: dict[
        str, list[_PosteriorContinuousNodeDistribution]
    ] = {}
    for sampled_state in sampled_states:
        topology_ids.add(sampled_state.tree.topology_id)
        sampled_tree = sampled_state.tree.to_tree()
        present_clade_ids: set[str] = set()
        for distribution in distribution_builder(
            sampled_state=sampled_state,
            tree=sampled_tree,
            tip_values=tip_values,
        ):
            clade_metadata_by_id.setdefault(
                distribution.clade_id,
                _PosteriorCladeMetadata(
                    clade_id=distribution.clade_id,
                    representative_node_id=distribution.node_id,
                    node_name=distribution.node_name,
                    descendant_taxa=distribution.descendant_taxa,
                ),
            )
            distributions_by_clade_id.setdefault(distribution.clade_id, []).append(
                distribution
            )
            present_clade_ids.add(distribution.clade_id)
        for clade_id in present_clade_ids:
            clade_presence_count_by_id[clade_id] = (
                clade_presence_count_by_id.get(clade_id, 0) + 1
            )
    return PosteriorContinuousTraitReport(
        sample_count=sample_count,
        distinct_topology_count=len(topology_ids),
        sampled_trait_models=sorted(set(sampled_trait_models)),
        tree_uncertainty_policy=_resolve_tree_uncertainty_policy(len(topology_ids)),
        node_summary_rows=_build_continuous_node_summary_rows(
            clade_metadata_by_id=clade_metadata_by_id,
            clade_presence_count_by_id=clade_presence_count_by_id,
            distributions_by_clade_id=distributions_by_clade_id,
            sample_count=sample_count,
        ),
        warnings=_build_summary_warnings(
            distinct_topology_count=len(topology_ids),
            sample_count=sample_count,
            clade_presence_count_by_id=clade_presence_count_by_id,
        ),
    )


def _build_continuous_node_summary_rows(
    *,
    clade_metadata_by_id: dict[str, _PosteriorCladeMetadata],
    clade_presence_count_by_id: dict[str, int],
    distributions_by_clade_id: dict[str, list[_PosteriorContinuousNodeDistribution]],
    sample_count: int,
) -> list[PosteriorContinuousTraitNodeSummaryRow]:
    summary_rows: list[PosteriorContinuousTraitNodeSummaryRow] = []
    for clade_id, metadata in sorted(
        clade_metadata_by_id.items(),
        key=lambda item: (len(item[1].descendant_taxa), item[1].descendant_taxa),
    ):
        distributions = distributions_by_clade_id[clade_id]
        supporting_sample_count = clade_presence_count_by_id[clade_id]
        mixture_draws = _build_continuous_mixture_draws(distributions)
        hpd_95_lower, hpd_95_upper = highest_posterior_density_interval(mixture_draws)
        summary_rows.append(
            PosteriorContinuousTraitNodeSummaryRow(
                clade_id=clade_id,
                representative_node_id=metadata.representative_node_id,
                node_name=metadata.node_name,
                descendant_taxa=list(metadata.descendant_taxa),
                clade_posterior_probability=float(
                    format(supporting_sample_count / sample_count, ".15g")
                ),
                conditional_posterior_mean=float(
                    format(
                        mean(
                            distribution.conditional_mean
                            for distribution in distributions
                        ),
                        ".15g",
                    )
                ),
                conditional_posterior_median=float(
                    format(median(mixture_draws), ".15g")
                ),
                conditional_hpd_95_lower=float(format(hpd_95_lower, ".15g")),
                conditional_hpd_95_upper=float(format(hpd_95_upper, ".15g")),
                mean_conditional_standard_deviation=float(
                    format(
                        mean(
                            distribution.conditional_standard_deviation
                            for distribution in distributions
                        ),
                        ".15g",
                    )
                ),
                supporting_sample_count=supporting_sample_count,
                total_sample_count=sample_count,
            )
        )
    return summary_rows


def _build_continuous_mixture_draws(
    distributions: list[_PosteriorContinuousNodeDistribution],
) -> list[float]:
    quantile_grid = [
        NormalDist().inv_cdf(
            (index + 0.5) / _CONTINUOUS_POSTERIOR_MIXTURE_QUANTILE_COUNT
        )
        for index in range(_CONTINUOUS_POSTERIOR_MIXTURE_QUANTILE_COUNT)
    ]
    draws: list[float] = []
    for distribution in distributions:
        if distribution.conditional_standard_deviation == 0.0:
            draws.extend(
                [distribution.conditional_mean]
                * _CONTINUOUS_POSTERIOR_MIXTURE_QUANTILE_COUNT
            )
            continue
        draws.extend(
            [
                float(
                    format(
                        distribution.conditional_mean
                        + (distribution.conditional_standard_deviation * quantile),
                        ".15g",
                    )
                )
                for quantile in quantile_grid
            ]
        )
    return draws


def _evaluate_brownian_node_distributions(
    *,
    sampled_state,
    tree: PhyloTree,
    tip_values: dict[str, float],
) -> list[_PosteriorContinuousNodeDistribution]:
    root_state = float(
        sampled_state.model_parameters.scalar_parameters.get("root-state")
    )
    sigma_squared = float(
        sampled_state.model_parameters.scalar_parameters.get("sigma-squared")
    )
    return _evaluate_continuous_node_distributions(
        tree=tree,
        tip_values=tip_values,
        location=root_state,
        scale=sigma_squared,
        covariance_evaluator=_evaluate_brownian_covariance,
    )


def _evaluate_ou_node_distributions(
    *,
    sampled_state,
    tree: PhyloTree,
    tip_values: dict[str, float],
) -> list[_PosteriorContinuousNodeDistribution]:
    optimum = float(sampled_state.model_parameters.scalar_parameters.get("optimum"))
    alpha = float(sampled_state.model_parameters.scalar_parameters.get("alpha"))
    sigma_squared = float(
        sampled_state.model_parameters.scalar_parameters.get("sigma-squared")
    )
    return _evaluate_continuous_node_distributions(
        tree=tree,
        tip_values=tip_values,
        location=optimum,
        scale=sigma_squared,
        covariance_evaluator=lambda left_node, right_node, depth_index: (
            _evaluate_ou_covariance(
                left_node=left_node,
                right_node=right_node,
                depth_index=depth_index,
                alpha=alpha,
            )
        ),
    )


def _evaluate_continuous_node_distributions(
    *,
    tree: PhyloTree,
    tip_values: dict[str, float],
    location: float,
    scale: float,
    covariance_evaluator,
) -> list[_PosteriorContinuousNodeDistribution]:
    if tree.rooted is not True:
        raise PhylogeneticsError(
            "posterior ancestral continuous-trait summary requires rooted sampled trees",
            code="posterior_ancestral_continuous_trait_tree_rooting_invalid",
        )
    ordered_taxa = sorted(tip_values)
    if set(tree.tip_names) != set(ordered_taxa):
        raise PhylogeneticsError(
            "posterior ancestral continuous-trait summary requires each sampled tree to match the observed tip-value taxon set",
            code="posterior_ancestral_continuous_trait_tip_value_taxa_mismatch",
            details={
                "tree_taxa": sorted(tree.tip_names),
                "tip_value_taxa": ordered_taxa,
            },
        )
    depth_index = _build_tree_depth_index(tree)
    tip_lookup = {
        node.name: node for node in tree.iter_leaves() if node.name is not None
    }
    ordered_tip_nodes = [tip_lookup[taxon] for taxon in ordered_taxa]
    ordered_tip_values = [tip_values[taxon] for taxon in ordered_taxa]
    root_distribution = _build_posterior_continuous_node_distribution(
        node=tree.root,
        conditional_mean=location,
        conditional_standard_deviation=0.0,
    )
    internal_nodes = [
        node
        for node in tree.iter_internal_nodes(order="preorder")
        if node is not tree.root
    ]
    if not internal_nodes:
        return [root_distribution]
    tip_covariance = stable_covariance(
        [
            [
                covariance_evaluator(left_node, right_node, depth_index)
                for right_node in ordered_tip_nodes
            ]
            for left_node in ordered_tip_nodes
        ]
    )
    inverse_tip_covariance = invert_matrix(tip_covariance)
    residual_vector = [tip_value - location for tip_value in ordered_tip_values]
    projected_residuals = matrix_vector_multiply(
        inverse_tip_covariance,
        residual_vector,
    )
    cross_covariance = [
        [
            covariance_evaluator(internal_node, tip_node, depth_index)
            for tip_node in ordered_tip_nodes
        ]
        for internal_node in internal_nodes
    ]
    conditional_mean_offsets = matrix_vector_multiply(
        cross_covariance,
        projected_residuals,
    )
    internal_covariance = [
        [
            covariance_evaluator(left_node, right_node, depth_index)
            for right_node in internal_nodes
        ]
        for left_node in internal_nodes
    ]
    projected_cross_covariance = matrix_multiply(
        cross_covariance,
        inverse_tip_covariance,
    )
    conditional_covariance = _subtract_matrices(
        internal_covariance,
        matrix_multiply(projected_cross_covariance, transpose(cross_covariance)),
    )
    return [
        root_distribution,
        *[
            _build_posterior_continuous_node_distribution(
                node=internal_node,
                conditional_mean=location + conditional_mean_offsets[index],
                conditional_standard_deviation=math.sqrt(
                    max(conditional_covariance[index][index], 0.0) * scale
                ),
            )
            for index, internal_node in enumerate(internal_nodes)
        ],
    ]


def _build_tree_depth_index(tree: PhyloTree) -> _TreeDepthIndex:
    depth_by_node_id: dict[str, float] = {}
    ancestor_depths_by_node_id: dict[str, dict[str, float]] = {}

    def visit(node: TreeNode, ancestors: dict[str, float], depth: float) -> None:
        if node.node_id is None:
            raise PhylogeneticsError(
                "posterior ancestral trait summary requires stable sampled tree node identifiers",
                code="posterior_ancestral_trait_node_id_missing",
            )
        if node is not tree.root:
            if node.branch_length is None:
                raise PhylogeneticsError(
                    "posterior ancestral trait summary requires complete branch lengths on every sampled tree",
                    code="posterior_ancestral_trait_branch_length_missing",
                )
            depth += node.branch_length
        current_ancestors = dict(ancestors)
        current_ancestors[node.node_id] = depth
        depth_by_node_id[node.node_id] = depth
        ancestor_depths_by_node_id[node.node_id] = current_ancestors
        for child in node.children:
            visit(child, current_ancestors, depth)

    visit(tree.root, {}, 0.0)
    return _TreeDepthIndex(
        depth_by_node_id=depth_by_node_id,
        ancestor_depths_by_node_id=ancestor_depths_by_node_id,
    )


def _evaluate_brownian_covariance(
    left_node: TreeNode,
    right_node: TreeNode,
    depth_index: _TreeDepthIndex,
) -> float:
    return _shared_ancestor_depth(left_node, right_node, depth_index)


def _evaluate_ou_covariance(
    *,
    left_node: TreeNode,
    right_node: TreeNode,
    depth_index: _TreeDepthIndex,
    alpha: float,
) -> float:
    left_depth = depth_index.depth_by_node_id[left_node.node_id or ""]
    right_depth = depth_index.depth_by_node_id[right_node.node_id or ""]
    shared_depth = _shared_ancestor_depth(left_node, right_node, depth_index)
    if left_node.node_id == right_node.node_id:
        return (1.0 - math.exp(-2.0 * alpha * left_depth)) / (2.0 * alpha)
    return (
        math.exp(-alpha * ((left_depth - shared_depth) + (right_depth - shared_depth)))
        * (1.0 - math.exp(-2.0 * alpha * shared_depth))
        / (2.0 * alpha)
    )


def _shared_ancestor_depth(
    left_node: TreeNode,
    right_node: TreeNode,
    depth_index: _TreeDepthIndex,
) -> float:
    left_node_id = left_node.node_id or ""
    right_node_id = right_node.node_id or ""
    left_ancestors = depth_index.ancestor_depths_by_node_id[left_node_id]
    right_ancestors = depth_index.ancestor_depths_by_node_id[right_node_id]
    shared_ancestor_ids = set(left_ancestors) & set(right_ancestors)
    return max(left_ancestors[node_id] for node_id in shared_ancestor_ids)


def _build_posterior_continuous_node_distribution(
    *,
    node: TreeNode,
    conditional_mean: float,
    conditional_standard_deviation: float,
) -> _PosteriorContinuousNodeDistribution:
    canonical_descendant_taxa = tuple(sorted(node.descendant_taxa))
    if node.node_id is None:
        raise PhylogeneticsError(
            "posterior ancestral continuous-trait summary requires stable node identifiers",
            code="posterior_ancestral_continuous_trait_node_id_missing",
        )
    return _PosteriorContinuousNodeDistribution(
        clade_id=canonical_clade_id(frozenset(canonical_descendant_taxa)),
        node_id=node.node_id,
        node_name=node.name,
        descendant_taxa=canonical_descendant_taxa,
        conditional_mean=float(format(conditional_mean, ".15g")),
        conditional_standard_deviation=float(
            format(conditional_standard_deviation, ".15g")
        ),
    )


def _subtract_matrices(
    left: list[list[float]],
    right: list[list[float]],
) -> list[list[float]]:
    return [
        [
            float(format(left_value - right_value, ".15g"))
            for left_value, right_value in zip(left_row, right_row, strict=True)
        ]
        for left_row, right_row in zip(left, right, strict=True)
    ]
