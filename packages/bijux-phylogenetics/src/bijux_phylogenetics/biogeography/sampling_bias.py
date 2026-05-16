from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.biogeography.geographic_states import (
    GeographicExcludedTaxonRow,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.comparative.discrete_evolution import (
    NodeStateEstimate,
    TransitionEvent,
    TransitionSupportRow,
    _estimate_node_states,
    _estimate_transition_support_rows,
    _fit_transition_matrix,
    _fitch_candidate_sets,
    _resolve_state_order,
    _root_prior,
    _stationary_frequencies,
    _transition_events,
    audit_discrete_state_coding,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

_MODEL_ALIAS_TO_INTERNAL = {
    "er": "equal-rates",
    "sym": "symmetric",
    "ard": "all-rates-different",
    "equal-rates": "equal-rates",
    "symmetric": "symmetric",
    "all-rates-different": "all-rates-different",
}

_INTERNAL_MODEL_TO_ALIAS = {
    "equal-rates": "er",
    "symmetric": "sym",
    "all-rates-different": "ard",
}

_DOMINANT_REGION_THRESHOLD = 0.8


@dataclass(frozen=True, slots=True)
class GeographicSamplingBiasSummary:
    """One summary row for region-sampling bias review on one tree."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    weighting_mode: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_region_count: int
    region_dominated: bool
    dominant_region: str
    dominant_region_fraction: float
    weighted_region_dominated: bool
    weighted_dominant_region: str
    weighted_dominant_region_fraction: float
    root_region_unweighted: str
    root_region_weighted: str
    root_region_changed: bool
    compared_internal_node_count: int
    changed_internal_node_count: int
    compared_transition_count: int
    changed_transition_count: int
    warning_count: int


@dataclass(frozen=True, slots=True)
class GeographicSamplingCountRow:
    """One observed-region sample-count and weight row."""

    region: str
    sample_count: int
    sample_fraction: float
    applied_weight: float
    weighted_sample_count: float
    weighted_sample_fraction: float
    dominant_unweighted: bool
    dominant_weighted: bool


@dataclass(frozen=True, slots=True)
class GeographicSamplingBiasNodeRow:
    """One weighted-versus-unweighted node region comparison row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    is_root: bool
    unweighted_region: str
    weighted_region: str
    unweighted_confidence: float
    weighted_confidence: float
    confidence_delta: float
    changed: bool
    unweighted_region_probabilities: dict[str, float]
    weighted_region_probabilities: dict[str, float]


@dataclass(frozen=True, slots=True)
class GeographicSamplingBiasTransitionRow:
    """One weighted-versus-unweighted branch transition comparison row."""

    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    unweighted_source_region: str
    unweighted_target_region: str
    weighted_source_region: str
    weighted_target_region: str
    unweighted_transition: str
    weighted_transition: str
    unweighted_changed: bool
    weighted_changed: bool
    changed_by_weighting: bool
    unweighted_support: float
    weighted_support: float


@dataclass(slots=True)
class GeographicSamplingBiasReport:
    """Owned review surface for region-sampling bias correction on one tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    weighting_mode: str
    summary: GeographicSamplingBiasSummary
    count_rows: list[GeographicSamplingCountRow]
    node_rows: list[GeographicSamplingBiasNodeRow]
    transition_rows: list[GeographicSamplingBiasTransitionRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class _SamplingBiasModelSurface:
    taxon_column: str
    estimates: list[NodeStateEstimate]
    events: list[TransitionEvent]
    support_rows: list[TransitionSupportRow]


def summarize_geographic_sampling_bias(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
    weights_path: Path | None = None,
    region_column: str = "region",
    weight_column: str = "weight",
) -> GeographicSamplingBiasReport:
    """Review how explicit region weights change biogeographic state inference."""
    internal_model = _resolve_internal_model(model)
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_regions,
    )
    exclusion_rows = _build_exclusion_rows(audit)
    included_counts = _included_region_counts(audit)
    weights, weighting_mode = _resolve_region_weights(
        included_counts,
        weights_path=weights_path,
        region_column=region_column,
        weight_column=weight_column,
    )
    baseline = _run_sampling_bias_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_regions=allowed_regions,
        region_weights=None,
    )
    weighted = _run_sampling_bias_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_regions=allowed_regions,
        region_weights=weights,
    )
    count_rows = _build_count_rows(included_counts, weights)
    node_rows = _build_node_rows(baseline, weighted)
    transition_rows = _build_transition_rows(baseline, weighted)
    warnings = _build_warnings(count_rows, weighting_mode, node_rows)
    summary = _build_summary(
        trait=trait,
        taxon_column=baseline.taxon_column,
        model=_INTERNAL_MODEL_TO_ALIAS[internal_model],
        internal_model=internal_model,
        analyzed_taxon_count=sum(included_counts.values()),
        excluded_taxon_count=len(exclusion_rows),
        weighting_mode=weighting_mode,
        count_rows=count_rows,
        node_rows=node_rows,
        transition_rows=transition_rows,
        warning_count=len(warnings),
    )
    return GeographicSamplingBiasReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=baseline.taxon_column,
        model=summary.model,
        internal_model=internal_model,
        weighting_mode=weighting_mode,
        summary=summary,
        count_rows=count_rows,
        node_rows=node_rows,
        transition_rows=transition_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )


def write_geographic_sampling_bias_summary_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one summary ledger for geographic sampling-bias review."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "weighting_mode",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_region_count",
            "region_dominated",
            "dominant_region",
            "dominant_region_fraction",
            "weighted_region_dominated",
            "weighted_dominant_region",
            "weighted_dominant_region_fraction",
            "root_region_unweighted",
            "root_region_weighted",
            "root_region_changed",
            "compared_internal_node_count",
            "changed_internal_node_count",
            "compared_transition_count",
            "changed_transition_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "weighting_mode": summary.weighting_mode,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_region_count": str(summary.observed_region_count),
                "region_dominated": str(summary.region_dominated).lower(),
                "dominant_region": summary.dominant_region,
                "dominant_region_fraction": str(summary.dominant_region_fraction),
                "weighted_region_dominated": str(
                    summary.weighted_region_dominated
                ).lower(),
                "weighted_dominant_region": summary.weighted_dominant_region,
                "weighted_dominant_region_fraction": str(
                    summary.weighted_dominant_region_fraction
                ),
                "root_region_unweighted": summary.root_region_unweighted,
                "root_region_weighted": summary.root_region_weighted,
                "root_region_changed": str(summary.root_region_changed).lower(),
                "compared_internal_node_count": str(
                    summary.compared_internal_node_count
                ),
                "changed_internal_node_count": str(summary.changed_internal_node_count),
                "compared_transition_count": str(summary.compared_transition_count),
                "changed_transition_count": str(summary.changed_transition_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_geographic_sampling_count_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one region sample-count and weighting ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "region",
            "sample_count",
            "sample_fraction",
            "applied_weight",
            "weighted_sample_count",
            "weighted_sample_fraction",
            "dominant_unweighted",
            "dominant_weighted",
        ],
        rows=[
            {
                "region": row.region,
                "sample_count": str(row.sample_count),
                "sample_fraction": str(row.sample_fraction),
                "applied_weight": str(row.applied_weight),
                "weighted_sample_count": str(row.weighted_sample_count),
                "weighted_sample_fraction": str(row.weighted_sample_fraction),
                "dominant_unweighted": str(row.dominant_unweighted).lower(),
                "dominant_weighted": str(row.dominant_weighted).lower(),
            }
            for row in report.count_rows
        ],
    )


def write_geographic_sampling_bias_node_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one weighted-versus-unweighted node comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "is_root",
            "unweighted_region",
            "weighted_region",
            "unweighted_confidence",
            "weighted_confidence",
            "confidence_delta",
            "changed",
            "unweighted_region_probabilities",
            "weighted_region_probabilities",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "is_root": str(row.is_root).lower(),
                "unweighted_region": row.unweighted_region,
                "weighted_region": row.weighted_region,
                "unweighted_confidence": str(row.unweighted_confidence),
                "weighted_confidence": str(row.weighted_confidence),
                "confidence_delta": str(row.confidence_delta),
                "changed": str(row.changed).lower(),
                "unweighted_region_probabilities": json.dumps(
                    row.unweighted_region_probabilities,
                    sort_keys=True,
                ),
                "weighted_region_probabilities": json.dumps(
                    row.weighted_region_probabilities,
                    sort_keys=True,
                ),
            }
            for row in report.node_rows
        ],
    )


def write_geographic_sampling_bias_transition_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one weighted-versus-unweighted transition comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "unweighted_source_region",
            "unweighted_target_region",
            "weighted_source_region",
            "weighted_target_region",
            "unweighted_transition",
            "weighted_transition",
            "unweighted_changed",
            "weighted_changed",
            "changed_by_weighting",
            "unweighted_support",
            "weighted_support",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "unweighted_source_region": row.unweighted_source_region,
                "unweighted_target_region": row.unweighted_target_region,
                "weighted_source_region": row.weighted_source_region,
                "weighted_target_region": row.weighted_target_region,
                "unweighted_transition": row.unweighted_transition,
                "weighted_transition": row.weighted_transition,
                "unweighted_changed": str(row.unweighted_changed).lower(),
                "weighted_changed": str(row.weighted_changed).lower(),
                "changed_by_weighting": str(row.changed_by_weighting).lower(),
                "unweighted_support": str(row.unweighted_support),
                "weighted_support": str(row.weighted_support),
            }
            for row in report.transition_rows
        ],
    )


def write_geographic_sampling_bias_exclusion_table(
    path: Path,
    report: GeographicSamplingBiasReport,
) -> Path:
    """Write one excluded-taxa ledger for geographic sampling-bias review."""
    return write_taxon_rows(
        path,
        columns=["taxon", "raw_region", "normalized_region", "reason", "note"],
        rows=[
            {
                "taxon": row.taxon,
                "raw_region": row.raw_region,
                "normalized_region": row.normalized_region or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )


def _included_region_counts(audit) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in audit.rows:
        if not row.included or row.normalized_state is None:
            continue
        counts[row.normalized_state] = counts.get(row.normalized_state, 0) + 1
    if not counts:
        raise AncestralReconstructionError(
            "geographic sampling-bias review requires at least one usable observed region"
        )
    return dict(sorted(counts.items()))


def _resolve_region_weights(
    counts: dict[str, int],
    *,
    weights_path: Path | None,
    region_column: str,
    weight_column: str,
) -> tuple[dict[str, float], str]:
    if weights_path is None:
        total = sum(counts.values())
        region_count = len(counts)
        return (
            {
                region: stable_value(total / max(region_count * count, 1))
                for region, count in counts.items()
            },
            "inverse-frequency",
        )
    table = load_taxon_table(weights_path, taxon_column=region_column)
    if weight_column not in table.columns:
        raise AncestralReconstructionError(
            f"region weight table does not contain column '{weight_column}'"
        )
    weights: dict[str, float] = {}
    for row in table.rows:
        region = row[table.taxon_column]
        raw_weight = row[weight_column]
        try:
            weight = float(raw_weight)
        except ValueError as error:
            raise AncestralReconstructionError(
                f"region weight for '{region}' is not numeric: {raw_weight!r}"
            ) from error
        if not math.isfinite(weight) or weight <= 0.0:
            raise AncestralReconstructionError(
                f"region weight for '{region}' must be finite and positive"
            )
        weights[region] = stable_value(weight)
    missing = sorted(region for region in counts if region not in weights)
    if missing:
        raise AncestralReconstructionError(
            "region weight table is missing observed regions: " + ", ".join(missing)
        )
    return ({region: weights[region] for region in sorted(counts)}, "explicit")


def _build_count_rows(
    counts: dict[str, int],
    weights: dict[str, float],
) -> list[GeographicSamplingCountRow]:
    total = sum(counts.values())
    weighted_counts = {
        region: stable_value(count * weights[region])
        for region, count in counts.items()
    }
    weighted_total = sum(weighted_counts.values())
    dominant_fraction = max(counts.values()) / max(total, 1)
    weighted_dominant_fraction = max(weighted_counts.values()) / max(
        weighted_total, 1.0
    )
    return [
        GeographicSamplingCountRow(
            region=region,
            sample_count=count,
            sample_fraction=stable_value(count / max(total, 1)),
            applied_weight=weights[region],
            weighted_sample_count=weighted_counts[region],
            weighted_sample_fraction=stable_value(
                weighted_counts[region] / max(weighted_total, 1.0)
            ),
            dominant_unweighted=stable_value(count / max(total, 1))
            == stable_value(dominant_fraction),
            dominant_weighted=stable_value(
                weighted_counts[region] / max(weighted_total, 1.0)
            )
            == stable_value(weighted_dominant_fraction),
        )
        for region, count in sorted(counts.items())
    ]


def _build_node_rows(
    baseline: _SamplingBiasModelSurface,
    weighted: _SamplingBiasModelSurface,
) -> list[GeographicSamplingBiasNodeRow]:
    weighted_by_node = {estimate.node: estimate for estimate in weighted.estimates}
    internal_estimates = [
        estimate for estimate in baseline.estimates if not estimate.is_tip
    ]
    root_node = max(
        internal_estimates,
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    ).node
    rows: list[GeographicSamplingBiasNodeRow] = []
    for estimate in internal_estimates:
        weighted_estimate = weighted_by_node[estimate.node]
        weighted_region = weighted_estimate.most_likely_state
        weighted_confidence = stable_value(
            weighted_estimate.state_probabilities[weighted_region]
        )
        unweighted_confidence = stable_value(
            estimate.state_probabilities[estimate.most_likely_state]
        )
        rows.append(
            GeographicSamplingBiasNodeRow(
                node=estimate.node,
                node_name=estimate.node_name,
                descendant_taxa=list(estimate.descendant_taxa),
                is_root=estimate.node == root_node,
                unweighted_region=estimate.most_likely_state,
                weighted_region=weighted_region,
                unweighted_confidence=unweighted_confidence,
                weighted_confidence=weighted_confidence,
                confidence_delta=stable_value(
                    weighted_confidence - unweighted_confidence
                ),
                changed=weighted_region != estimate.most_likely_state,
                unweighted_region_probabilities=dict(estimate.state_probabilities),
                weighted_region_probabilities=dict(
                    weighted_estimate.state_probabilities
                ),
            )
        )
    return rows


def _build_transition_rows(
    baseline: _SamplingBiasModelSurface,
    weighted: _SamplingBiasModelSurface,
) -> list[GeographicSamplingBiasTransitionRow]:
    baseline_support = {
        (row.parent_node, row.child_node): row for row in baseline.support_rows
    }
    weighted_support = {
        (row.parent_node, row.child_node): row for row in weighted.support_rows
    }
    weighted_event_by_branch = {
        (event.parent_node, event.child_node): event for event in weighted.events
    }
    rows: list[GeographicSamplingBiasTransitionRow] = []
    for event in baseline.events:
        weighted_event = weighted_event_by_branch[(event.parent_node, event.child_node)]
        rows.append(
            GeographicSamplingBiasTransitionRow(
                parent_node=event.parent_node,
                child_node=event.child_node,
                child_descendant_taxa=event.child_node.split("|"),
                unweighted_source_region=event.source_state,
                unweighted_target_region=event.target_state,
                weighted_source_region=weighted_event.source_state,
                weighted_target_region=weighted_event.target_state,
                unweighted_transition=f"{event.source_state}->{event.target_state}",
                weighted_transition=(
                    f"{weighted_event.source_state}->{weighted_event.target_state}"
                ),
                unweighted_changed=event.changed,
                weighted_changed=weighted_event.changed,
                changed_by_weighting=(
                    event.source_state != weighted_event.source_state
                    or event.target_state != weighted_event.target_state
                ),
                unweighted_support=baseline_support[
                    (event.parent_node, event.child_node)
                ].support,
                weighted_support=weighted_support[
                    (event.parent_node, event.child_node)
                ].support,
            )
        )
    return rows


def _build_warnings(
    count_rows: list[GeographicSamplingCountRow],
    weighting_mode: str,
    node_rows: list[GeographicSamplingBiasNodeRow],
) -> list[str]:
    warnings: list[str] = []
    dominant_row = max(count_rows, key=lambda row: (row.sample_fraction, row.region))
    if dominant_row.sample_fraction >= _DOMINANT_REGION_THRESHOLD:
        warnings.append(
            "observed regions are dominated by one sampled region and the baseline reconstruction may reflect sample imbalance"
        )
    weighted_dominant_row = max(
        count_rows,
        key=lambda row: (row.weighted_sample_fraction, row.region),
    )
    if weighting_mode == "inverse-frequency":
        warnings.append(
            "inverse-frequency region weights rebalance observed regions to equal weighted mass before comparing ancestral conclusions"
        )
    else:
        warnings.append(
            "explicit region weights reweight ancestral region probabilities and branchwise transitions for sampling-bias review"
        )
    if weighted_dominant_row.weighted_sample_fraction >= _DOMINANT_REGION_THRESHOLD:
        warnings.append(
            "weighted region counts remain dominated by one region after correction"
        )
    if any(row.changed and row.is_root for row in node_rows):
        warnings.append(
            "the weighted correction changes the most likely root region relative to the baseline reconstruction"
        )
    return warnings


def _build_summary(
    *,
    trait: str,
    taxon_column: str,
    model: str,
    internal_model: str,
    analyzed_taxon_count: int,
    excluded_taxon_count: int,
    weighting_mode: str,
    count_rows: list[GeographicSamplingCountRow],
    node_rows: list[GeographicSamplingBiasNodeRow],
    transition_rows: list[GeographicSamplingBiasTransitionRow],
    warning_count: int,
) -> GeographicSamplingBiasSummary:
    root_row = next(row for row in node_rows if row.is_root)
    dominant_row = max(count_rows, key=lambda row: (row.sample_fraction, row.region))
    weighted_dominant_row = max(
        count_rows,
        key=lambda row: (row.weighted_sample_fraction, row.region),
    )
    changed_node_count = sum(1 for row in node_rows if row.changed)
    changed_transition_count = sum(
        1
        for row in transition_rows
        if row.unweighted_transition != row.weighted_transition
    )
    return GeographicSamplingBiasSummary(
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        internal_model=internal_model,
        weighting_mode=weighting_mode,
        analyzed_taxon_count=analyzed_taxon_count,
        excluded_taxon_count=excluded_taxon_count,
        observed_region_count=len(count_rows),
        region_dominated=dominant_row.sample_fraction >= _DOMINANT_REGION_THRESHOLD,
        dominant_region=dominant_row.region,
        dominant_region_fraction=dominant_row.sample_fraction,
        weighted_region_dominated=(
            weighted_dominant_row.weighted_sample_fraction >= _DOMINANT_REGION_THRESHOLD
        ),
        weighted_dominant_region=weighted_dominant_row.region,
        weighted_dominant_region_fraction=weighted_dominant_row.weighted_sample_fraction,
        root_region_unweighted=root_row.unweighted_region,
        root_region_weighted=root_row.weighted_region,
        root_region_changed=root_row.changed,
        compared_internal_node_count=len(node_rows),
        changed_internal_node_count=changed_node_count,
        compared_transition_count=len(transition_rows),
        changed_transition_count=changed_transition_count,
        warning_count=warning_count,
    )


def _normalize_probabilities(probabilities: dict[str, float]) -> dict[str, float]:
    total = sum(probabilities.values())
    if total <= 0.0:
        uniform = 1.0 / max(len(probabilities), 1)
        return {region: stable_value(uniform) for region in sorted(probabilities)}
    return {
        region: stable_value(probability / total)
        for region, probability in sorted(probabilities.items())
    }


def _most_likely_region(probabilities: dict[str, float]) -> str:
    return max(sorted(probabilities), key=lambda region: probabilities[region])


def _resolve_internal_model(model: str) -> str:
    if model not in _MODEL_ALIAS_TO_INTERNAL:
        raise ValueError(f"unsupported biogeography model alias: {model}")
    return _MODEL_ALIAS_TO_INTERNAL[model]


def _build_exclusion_rows(audit) -> list[GeographicExcludedTaxonRow]:
    return [
        GeographicExcludedTaxonRow(
            taxon=row.taxon,
            raw_region=row.raw_state,
            normalized_region=row.normalized_state,
            reason=row.issue_code or "excluded",
            note=row.note,
        )
        for row in audit.rows
        if not row.included
    ]


def _run_sampling_bias_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
    model: str,
    allowed_regions: list[str] | None,
    region_weights: dict[str, float] | None,
) -> _SamplingBiasModelSurface:
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if len(dataset.observed_states) < 2:
        raise AncestralReconstructionError(
            "geographic sampling-bias review requires at least two observed regions"
        )
    state_order = _resolve_state_order(
        dataset.observed_states,
        allowed_states=allowed_regions,
        ordered_states=None,
        state_ordering="unordered",
    )
    candidate_sets = _fitch_candidate_sets(dataset.tree, dataset.states_by_taxon)
    stationary = (
        _weighted_stationary_frequencies(
            dataset.states_by_taxon, state_order, region_weights
        )
        if region_weights is not None
        else _stationary_frequencies(dataset.states_by_taxon, state_order)
    )
    priority_weights = (
        {state: region_weights.get(state, 1.0) for state in state_order}
        if region_weights is not None
        else dict.fromkeys(state_order, 1.0)
    )
    er_resolved = _resolve_biased_er_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        priority_weights,
    )
    er_events = _transition_events(dataset.tree, er_resolved)
    matrix = _fit_transition_matrix(
        model,
        state_order,
        stationary,
        er_events,
        state_ordering="unordered",
    )
    root_prior = _root_prior(
        model,
        stationary,
        candidate_sets[node_signature(dataset.tree.root)],
    )
    estimates = _estimate_node_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        matrix,
        root_prior,
        state_ordering="unordered",
    )
    resolved_states = {
        estimate.node: estimate.most_likely_state for estimate in estimates
    }
    events = _transition_events(dataset.tree, resolved_states)
    support_rows = _estimate_transition_support_rows(
        estimates=estimates,
        events=events,
        transition_matrix=matrix,
    )
    return _SamplingBiasModelSurface(
        taxon_column=dataset.taxon_column,
        estimates=estimates,
        events=events,
        support_rows=support_rows,
    )


def _weighted_stationary_frequencies(
    states_by_taxon: dict[str, str],
    state_order: list[str],
    region_weights: dict[str, float],
) -> dict[str, float]:
    weighted_counts = {
        state: stable_value(
            sum(
                region_weights.get(observed_state, 1.0)
                for observed_state in states_by_taxon.values()
                if observed_state == state
            )
        )
        for state in state_order
    }
    return _normalize_probabilities(weighted_counts)


def _resolve_biased_er_states(
    tree,
    candidate_sets: dict[str, list[str]],
    states_by_taxon: dict[str, str],
    state_order: list[str],
    priority_weights: dict[str, float],
) -> dict[str, str]:
    support_by_node = {
        node_signature(node): _state_support(node, states_by_taxon, state_order)
        for node in tree.iter_nodes()
    }
    resolved: dict[str, str] = {}
    root_signature = node_signature(tree.root)
    root_candidates = candidate_sets[root_signature]
    resolved[root_signature] = _best_supported_state(
        root_candidates,
        support_by_node[root_signature],
        priority_weights,
    )

    def visit(node, parent_state: str) -> None:
        for child in node.children:
            signature = node_signature(child)
            candidates = candidate_sets[signature]
            if parent_state in candidates:
                chosen = parent_state
            else:
                chosen = _best_supported_state(
                    candidates,
                    support_by_node[signature],
                    priority_weights,
                )
            resolved[signature] = chosen
            if not child.is_leaf():
                visit(child, chosen)

    visit(tree.root, resolved[root_signature])
    return resolved


def _state_support(
    node,
    states_by_taxon: dict[str, str],
    state_order: list[str],
) -> dict[str, int]:
    counts = dict.fromkeys(state_order, 0)
    for taxon in node_descendant_taxa(node):
        state = states_by_taxon.get(taxon)
        if state is not None:
            counts[state] += 1
    return counts


def _best_supported_state(
    candidate_states: list[str],
    support_counts: dict[str, int],
    priority_weights: dict[str, float],
) -> str:
    return max(
        sorted(candidate_states),
        key=lambda state: (
            support_counts.get(state, 0),
            priority_weights.get(state, 0.0),
            state,
        ),
    )
