from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .discrete_trait_mk import (
    DiscreteTraitMkNodeStateSummary,
    DiscreteTraitMkRunReport,
)


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
