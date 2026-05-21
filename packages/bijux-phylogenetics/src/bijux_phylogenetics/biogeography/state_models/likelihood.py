from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
    estimate_ancestral_geographic_states,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

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


@dataclass(frozen=True, slots=True)
class GeographicStateSummary:
    """One summary row for geographic state modeling on one tree."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_region_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    changed_branch_count: int
    strongly_supported_transition_count: int
    transition_rate_row_count: int
    root_region: str
    root_region_probability: float
    warning_count: int


@dataclass(frozen=True, slots=True)
class GeographicRegionProbabilityRow:
    """One internal-node ancestral region probability row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_region: str
    region_probabilities: dict[str, float]
    confidence: float
    ambiguous: bool
    is_root: bool


@dataclass(frozen=True, slots=True)
class GeographicTransitionRateRow:
    """One pairwise geographic transition-rate row."""

    source_region: str
    target_region: str
    rate: float
    lower_95_interval: float
    upper_95_interval: float
    effective_transition_count: float


@dataclass(frozen=True, slots=True)
class GeographicTransitionEventRow:
    """One branchwise inferred geographic transition row."""

    parent_node: str
    child_node: str
    source_region: str
    target_region: str
    changed: bool
    support: float
    strongly_supported: bool


@dataclass(frozen=True, slots=True)
class GeographicExcludedTaxonRow:
    """One excluded metadata row from the geographic workflow."""

    taxon: str
    raw_region: str
    normalized_region: str | None
    reason: str
    note: str


@dataclass(slots=True)
class GeographicStateModelReport:
    """Owned biogeography review surface for one ancestral geographic model."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    summary: GeographicStateSummary
    node_rows: list[GeographicRegionProbabilityRow]
    transition_rate_rows: list[GeographicTransitionRateRow]
    transition_event_rows: list[GeographicTransitionEventRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


def summarize_geographic_state_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
) -> GeographicStateModelReport:
    """Build a governed biogeography review surface for one tree and region table."""
    internal_model = _resolve_internal_model(model)
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_regions,
    )
    reconstruction = estimate_ancestral_geographic_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_states=allowed_regions,
    )
    node_rows = _build_node_rows(reconstruction)
    transition_rate_rows = _build_transition_rate_rows(reconstruction)
    transition_event_rows = _build_transition_event_rows(reconstruction)
    exclusion_rows = _build_exclusion_rows(audit)
    summary = _build_summary(
        reconstruction,
        model=_INTERNAL_MODEL_TO_ALIAS[internal_model],
        exclusion_rows=exclusion_rows,
        node_rows=node_rows,
        transition_rate_rows=transition_rate_rows,
    )
    return GeographicStateModelReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=reconstruction.taxon_column,
        model=summary.model,
        internal_model=internal_model,
        likelihood_method=reconstruction.likelihood_method,
        summary=summary,
        node_rows=node_rows,
        transition_rate_rows=transition_rate_rows,
        transition_event_rows=transition_event_rows,
        exclusion_rows=exclusion_rows,
        warnings=list(reconstruction.warnings),
    )


def write_geographic_state_summary_table(
    path: Path,
    report: GeographicStateModelReport,
) -> Path:
    """Write one summary ledger for a geographic state model."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "likelihood_method",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_region_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "changed_branch_count",
            "strongly_supported_transition_count",
            "transition_rate_row_count",
            "root_region",
            "root_region_probability",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "likelihood_method": summary.likelihood_method,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_region_count": str(summary.observed_region_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "changed_branch_count": str(summary.changed_branch_count),
                "strongly_supported_transition_count": str(
                    summary.strongly_supported_transition_count
                ),
                "transition_rate_row_count": str(summary.transition_rate_row_count),
                "root_region": summary.root_region,
                "root_region_probability": str(summary.root_region_probability),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_geographic_region_probability_table(
    path: Path,
    report: GeographicStateModelReport,
) -> Path:
    """Write one internal-node ancestral region probability ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_region",
            "region_probabilities",
            "confidence",
            "ambiguous",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "most_likely_region": row.most_likely_region,
                "region_probabilities": json.dumps(
                    row.region_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_geographic_transition_rate_table(
    path: Path,
    report: GeographicStateModelReport,
) -> Path:
    """Write one pairwise geographic transition-rate ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "source_region",
            "target_region",
            "rate",
            "lower_95_interval",
            "upper_95_interval",
            "effective_transition_count",
        ],
        rows=[
            {
                "source_region": row.source_region,
                "target_region": row.target_region,
                "rate": str(row.rate),
                "lower_95_interval": str(row.lower_95_interval),
                "upper_95_interval": str(row.upper_95_interval),
                "effective_transition_count": str(row.effective_transition_count),
            }
            for row in report.transition_rate_rows
        ],
    )


def write_geographic_transition_event_table(
    path: Path,
    report: GeographicStateModelReport,
) -> Path:
    """Write one branchwise geographic transition-event ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "source_region",
            "target_region",
            "changed",
            "support",
            "strongly_supported",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "source_region": row.source_region,
                "target_region": row.target_region,
                "changed": str(row.changed).lower(),
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
            }
            for row in report.transition_event_rows
        ],
    )


def write_geographic_exclusion_table(
    path: Path,
    report: GeographicStateModelReport,
) -> Path:
    """Write one excluded-taxa ledger for geographic state modeling."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_region",
            "normalized_region",
            "reason",
            "note",
        ],
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


def _resolve_internal_model(model: str) -> str:
    try:
        return _MODEL_ALIAS_TO_INTERNAL[model]
    except KeyError as error:
        raise ValueError(
            "biogeography model must be one of er, sym, ard, equal-rates, symmetric, or all-rates-different"
        ) from error


def _build_node_rows(reconstruction) -> list[GeographicRegionProbabilityRow]:
    internal_estimates = [
        estimate for estimate in reconstruction.estimates if not estimate.is_tip
    ]
    root_node = max(
        internal_estimates,
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    ).node
    rows: list[GeographicRegionProbabilityRow] = []
    for estimate in internal_estimates:
        confidence = max(estimate.state_probabilities.values(), default=0.0)
        rows.append(
            GeographicRegionProbabilityRow(
                node=estimate.node,
                node_name=estimate.node_name,
                descendant_taxa=estimate.descendant_taxa,
                most_likely_region=estimate.most_likely_state,
                region_probabilities=estimate.state_probabilities,
                confidence=confidence,
                ambiguous=estimate.ambiguous,
                is_root=estimate.node == root_node,
            )
        )
    return rows


def _build_transition_rate_rows(reconstruction) -> list[GeographicTransitionRateRow]:
    uncertainty_by_pair = {
        (row.source_state, row.target_state): row
        for row in reconstruction.transition_model.uncertainty.rows
    }
    rows: list[GeographicTransitionRateRow] = []
    for source_row in reconstruction.transition_model.transition_matrix:
        for target_region, rate in sorted(source_row.target_rates.items()):
            uncertainty = uncertainty_by_pair[(source_row.source_state, target_region)]
            rows.append(
                GeographicTransitionRateRow(
                    source_region=source_row.source_state,
                    target_region=target_region,
                    rate=rate,
                    lower_95_interval=uncertainty.lower_95_interval,
                    upper_95_interval=uncertainty.upper_95_interval,
                    effective_transition_count=uncertainty.effective_transition_count,
                )
            )
    return rows


def _build_transition_event_rows(reconstruction) -> list[GeographicTransitionEventRow]:
    support_by_branch = {
        (row.parent_node, row.child_node): row
        for row in reconstruction.transition_summary.support_rows
    }
    return [
        GeographicTransitionEventRow(
            parent_node=event.parent_node,
            child_node=event.child_node,
            source_region=event.source_state,
            target_region=event.target_state,
            changed=event.changed,
            support=support_by_branch[(event.parent_node, event.child_node)].support,
            strongly_supported=support_by_branch[
                (event.parent_node, event.child_node)
            ].strongly_supported,
        )
        for event in reconstruction.transition_summary.events
    ]


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


def _build_summary(
    reconstruction,
    *,
    model: str,
    exclusion_rows: list[GeographicExcludedTaxonRow],
    node_rows: list[GeographicRegionProbabilityRow],
    transition_rate_rows: list[GeographicTransitionRateRow],
) -> GeographicStateSummary:
    root_probabilities = reconstruction.transition_model.root_state_probabilities
    root_region = max(
        root_probabilities,
        key=lambda state: (root_probabilities[state], state),
    )
    return GeographicStateSummary(
        trait=reconstruction.trait,
        taxon_column=reconstruction.taxon_column,
        model=model,
        internal_model=reconstruction.model,
        likelihood_method=reconstruction.likelihood_method,
        analyzed_taxon_count=reconstruction.taxon_count,
        excluded_taxon_count=len(exclusion_rows),
        observed_region_count=len(reconstruction.observed_states),
        internal_node_count=len(node_rows),
        ambiguous_internal_node_count=sum(row.ambiguous for row in node_rows),
        changed_branch_count=reconstruction.transition_summary.transition_count,
        strongly_supported_transition_count=(
            reconstruction.transition_summary.strongly_supported_transition_count
        ),
        transition_rate_row_count=len(transition_rate_rows),
        root_region=root_region,
        root_region_probability=root_probabilities[root_region],
        warning_count=len(reconstruction.warnings),
    )
