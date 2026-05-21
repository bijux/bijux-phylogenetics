from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .likelihood import GeographicExcludedTaxonRow

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

_TRUE_MATRIX_VALUES = {"1", "true", "yes", "allowed", "x"}
_FALSE_MATRIX_VALUES = {"", "0", "false", "no", "forbidden"}
_UNSUPPORTED_CLAIM_COMPETITOR_RATIO = 0.9


@dataclass(frozen=True, slots=True)
class ConstrainedGeographicFitRow:
    """One geographic fit under one transition-constraint regime."""

    constraint_mode: str
    model: str
    analyzed_taxon_count: int
    log_likelihood: float
    parameter_count: int
    aic: float
    root_region: str
    root_confidence: float


@dataclass(frozen=True, slots=True)
class ConstrainedGeographicTransitionRow:
    """One directed geographic transition compared across two regimes."""

    source_region: str
    target_region: str
    transition_allowed: bool
    unconstrained_rate: float
    constrained_rate: float
    rate_delta: float


@dataclass(frozen=True, slots=True)
class UnsupportedGeographicTransitionClaimRow:
    """One unconstrained branchwise transition claim forbidden by adjacency."""

    parent_node: str
    child_node: str
    descendant_taxa: list[str]
    unconstrained_source_region: str
    unconstrained_target_region: str
    unconstrained_support: float
    constrained_source_region: str
    constrained_target_region: str
    constrained_support: float
    claim_resolved: bool


@dataclass(frozen=True, slots=True)
class ConstrainedGeographicSummary:
    """Reviewer-facing summary for constrained geography modeling."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_region_count: int
    allowed_transition_count: int
    forbidden_transition_count: int
    constrained_log_likelihood: float
    unconstrained_log_likelihood: float
    likelihood_difference: float
    constrained_parameter_count: int
    unconstrained_parameter_count: int
    constrained_aic: float
    unconstrained_aic: float
    delta_aic: float
    preferred_constraint: str
    unsupported_transition_claim_count: int
    warning_count: int


@dataclass(slots=True)
class ConstrainedGeographicReport:
    """Comparison between constrained and unconstrained geographic fits."""

    tree_path: Path
    traits_path: Path
    adjacency_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    allowed_transition_pairs: list[tuple[str, str]]
    fit_rows: list[ConstrainedGeographicFitRow]
    transition_rows: list[ConstrainedGeographicTransitionRow]
    unsupported_claim_rows: list[UnsupportedGeographicTransitionClaimRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class _GeographicAdjacencyMatrix:
    path: Path
    region_order: list[str]
    allowed_transition_pairs: list[tuple[str, str]]


def summarize_constrained_geographic_model(
    tree_path: Path,
    traits_path: Path,
    adjacency_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "ard",
) -> ConstrainedGeographicReport:
    """Compare constrained and unconstrained geographic fits under one adjacency matrix."""
    adjacency = _load_geographic_adjacency_matrix(adjacency_path)
    internal_model = _resolve_internal_model(model)
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=adjacency.region_order,
    )
    exclusion_rows = _build_exclusion_rows(audit)
    _validate_in_tree_region_labels(audit)
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    observed_states = set(dataset.observed_states)
    filtered_allowed_transition_pairs = [
        (source_region, target_region)
        for source_region, target_region in adjacency.allowed_transition_pairs
        if source_region in observed_states and target_region in observed_states
    ]
    constrained_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_transition_pairs=filtered_allowed_transition_pairs,
    )
    unconstrained_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
    )
    fit_rows = [
        _build_fit_row(constrained_report, constraint_mode="constrained"),
        _build_fit_row(unconstrained_report, constraint_mode="unconstrained"),
    ]
    transition_rows = _build_transition_rows(
        dataset.observed_states,
        constrained_report=constrained_report,
        unconstrained_report=unconstrained_report,
        allowed_transition_pairs=filtered_allowed_transition_pairs,
    )
    unsupported_claim_rows = _build_unsupported_claim_rows(
        dataset.tree.root,
        constrained_report=constrained_report,
        unconstrained_report=unconstrained_report,
        allowed_transition_pairs=set(filtered_allowed_transition_pairs),
    )
    warnings = list(constrained_report.warnings)
    if unsupported_claim_rows:
        warnings.append(
            "one or more unconstrained geographic transition claims violate the supplied region adjacency matrix"
        )
    return ConstrainedGeographicReport(
        tree_path=tree_path,
        traits_path=traits_path,
        adjacency_path=adjacency_path,
        trait=trait,
        taxon_column=dataset.taxon_column,
        model=_INTERNAL_MODEL_TO_ALIAS[internal_model],
        internal_model=internal_model,
        allowed_transition_pairs=filtered_allowed_transition_pairs,
        fit_rows=fit_rows,
        transition_rows=transition_rows,
        unsupported_claim_rows=unsupported_claim_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )


def summarize_constrained_geographic_report(
    report: ConstrainedGeographicReport,
) -> ConstrainedGeographicSummary:
    """Summarize the main review facts for constrained geography modeling."""
    constrained_fit = next(
        row for row in report.fit_rows if row.constraint_mode == "constrained"
    )
    unconstrained_fit = next(
        row for row in report.fit_rows if row.constraint_mode == "unconstrained"
    )
    delta_aic = stable_value(constrained_fit.aic - unconstrained_fit.aic)
    if abs(delta_aic) < 1e-9:
        preferred_constraint = "tie"
    elif delta_aic < 0.0:
        preferred_constraint = "constrained"
    else:
        preferred_constraint = "unconstrained"
    observed_regions = sorted(
        {row.source_region for row in report.transition_rows}
        | {row.target_region for row in report.transition_rows}
    )
    total_directed_pairs = len(observed_regions) * max(len(observed_regions) - 1, 0)
    return ConstrainedGeographicSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        internal_model=report.internal_model,
        analyzed_taxon_count=constrained_fit.analyzed_taxon_count,
        excluded_taxon_count=len(report.exclusion_rows),
        observed_region_count=len(observed_regions),
        allowed_transition_count=len(report.allowed_transition_pairs),
        forbidden_transition_count=total_directed_pairs
        - len(report.allowed_transition_pairs),
        constrained_log_likelihood=constrained_fit.log_likelihood,
        unconstrained_log_likelihood=unconstrained_fit.log_likelihood,
        likelihood_difference=stable_value(
            constrained_fit.log_likelihood - unconstrained_fit.log_likelihood
        ),
        constrained_parameter_count=constrained_fit.parameter_count,
        unconstrained_parameter_count=unconstrained_fit.parameter_count,
        constrained_aic=constrained_fit.aic,
        unconstrained_aic=unconstrained_fit.aic,
        delta_aic=delta_aic,
        preferred_constraint=preferred_constraint,
        unsupported_transition_claim_count=len(report.unsupported_claim_rows),
        warning_count=len(report.warnings),
    )


def write_constrained_geographic_summary_table(
    path: Path,
    report: ConstrainedGeographicReport,
) -> Path:
    """Write one overall summary ledger for constrained geography modeling."""
    summary = summarize_constrained_geographic_report(report)
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_region_count",
            "allowed_transition_count",
            "forbidden_transition_count",
            "constrained_log_likelihood",
            "unconstrained_log_likelihood",
            "likelihood_difference",
            "constrained_parameter_count",
            "unconstrained_parameter_count",
            "constrained_aic",
            "unconstrained_aic",
            "delta_aic",
            "preferred_constraint",
            "unsupported_transition_claim_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_region_count": str(summary.observed_region_count),
                "allowed_transition_count": str(summary.allowed_transition_count),
                "forbidden_transition_count": str(summary.forbidden_transition_count),
                "constrained_log_likelihood": str(summary.constrained_log_likelihood),
                "unconstrained_log_likelihood": str(
                    summary.unconstrained_log_likelihood
                ),
                "likelihood_difference": str(summary.likelihood_difference),
                "constrained_parameter_count": str(summary.constrained_parameter_count),
                "unconstrained_parameter_count": str(
                    summary.unconstrained_parameter_count
                ),
                "constrained_aic": str(summary.constrained_aic),
                "unconstrained_aic": str(summary.unconstrained_aic),
                "delta_aic": str(summary.delta_aic),
                "preferred_constraint": summary.preferred_constraint,
                "unsupported_transition_claim_count": str(
                    summary.unsupported_transition_claim_count
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_constrained_geographic_fit_table(
    path: Path,
    report: ConstrainedGeographicReport,
) -> Path:
    """Write one constrained-versus-unconstrained fit ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "constraint_mode",
            "model",
            "analyzed_taxon_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "root_region",
            "root_confidence",
        ],
        rows=[
            {
                "constraint_mode": row.constraint_mode,
                "model": row.model,
                "analyzed_taxon_count": str(row.analyzed_taxon_count),
                "log_likelihood": str(row.log_likelihood),
                "parameter_count": str(row.parameter_count),
                "aic": str(row.aic),
                "root_region": row.root_region,
                "root_confidence": str(row.root_confidence),
            }
            for row in report.fit_rows
        ],
    )


def write_constrained_geographic_transition_table(
    path: Path,
    report: ConstrainedGeographicReport,
) -> Path:
    """Write one transition-by-transition constrained comparison ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "source_region",
            "target_region",
            "transition_allowed",
            "unconstrained_rate",
            "constrained_rate",
            "rate_delta",
        ],
        rows=[
            {
                "source_region": row.source_region,
                "target_region": row.target_region,
                "transition_allowed": str(row.transition_allowed).lower(),
                "unconstrained_rate": str(row.unconstrained_rate),
                "constrained_rate": str(row.constrained_rate),
                "rate_delta": str(row.rate_delta),
            }
            for row in report.transition_rows
        ],
    )


def write_unsupported_geographic_transition_claim_table(
    path: Path,
    report: ConstrainedGeographicReport,
) -> Path:
    """Write one branchwise ledger of forbidden unconstrained transition claims."""
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "descendant_taxa",
            "unconstrained_source_region",
            "unconstrained_target_region",
            "unconstrained_support",
            "constrained_source_region",
            "constrained_target_region",
            "constrained_support",
            "claim_resolved",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "unconstrained_source_region": row.unconstrained_source_region,
                "unconstrained_target_region": row.unconstrained_target_region,
                "unconstrained_support": str(row.unconstrained_support),
                "constrained_source_region": row.constrained_source_region,
                "constrained_target_region": row.constrained_target_region,
                "constrained_support": str(row.constrained_support),
                "claim_resolved": str(row.claim_resolved).lower(),
            }
            for row in report.unsupported_claim_rows
        ],
    )


def write_constrained_geographic_exclusion_table(
    path: Path,
    report: ConstrainedGeographicReport,
) -> Path:
    """Write one excluded-taxa ledger for constrained geography modeling."""
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


def _build_fit_row(
    report: DiscreteAncestralReport,
    *,
    constraint_mode: str,
) -> ConstrainedGeographicFitRow:
    root_estimate = max(
        (estimate for estimate in report.estimates if not estimate.is_tip),
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    )
    return ConstrainedGeographicFitRow(
        constraint_mode=constraint_mode,
        model=report.model,
        analyzed_taxon_count=report.taxon_count,
        log_likelihood=stable_value(report.log_likelihood or 0.0),
        parameter_count=int(report.parameter_count or 0),
        aic=stable_value(report.aic or 0.0),
        root_region=root_estimate.most_likely_state,
        root_confidence=stable_value(root_estimate.confidence),
    )


def _build_transition_rows(
    observed_states: list[str],
    *,
    constrained_report: DiscreteAncestralReport,
    unconstrained_report: DiscreteAncestralReport,
    allowed_transition_pairs: list[tuple[str, str]],
) -> list[ConstrainedGeographicTransitionRow]:
    constrained_lookup = {
        (row.source_state, row.target_state): stable_value(row.rate)
        for row in constrained_report.transition_rate_rows
    }
    unconstrained_lookup = {
        (row.source_state, row.target_state): stable_value(row.rate)
        for row in unconstrained_report.transition_rate_rows
    }
    allowed_set = set(allowed_transition_pairs)
    rows: list[ConstrainedGeographicTransitionRow] = []
    for source_region in sorted(observed_states):
        for target_region in sorted(observed_states):
            if source_region == target_region:
                continue
            unconstrained_rate = unconstrained_lookup.get(
                (source_region, target_region),
                0.0,
            )
            constrained_rate = constrained_lookup.get(
                (source_region, target_region),
                0.0,
            )
            rows.append(
                ConstrainedGeographicTransitionRow(
                    source_region=source_region,
                    target_region=target_region,
                    transition_allowed=(source_region, target_region) in allowed_set,
                    unconstrained_rate=stable_value(unconstrained_rate),
                    constrained_rate=stable_value(constrained_rate),
                    rate_delta=stable_value(constrained_rate - unconstrained_rate),
                )
            )
    return rows


def _build_unsupported_claim_rows(
    root,
    *,
    constrained_report: DiscreteAncestralReport,
    unconstrained_report: DiscreteAncestralReport,
    allowed_transition_pairs: set[tuple[str, str]],
) -> list[UnsupportedGeographicTransitionClaimRow]:
    constrained_by_node = {
        estimate.node: estimate for estimate in constrained_report.estimates
    }
    unconstrained_by_node = {
        estimate.node: estimate for estimate in unconstrained_report.estimates
    }
    rows: list[UnsupportedGeographicTransitionClaimRow] = []

    def visit(parent, child) -> None:
        parent_signature = node_signature(parent)
        child_signature = node_signature(child)
        unconstrained_parent = unconstrained_by_node[parent_signature]
        unconstrained_child = unconstrained_by_node[child_signature]
        transition_pair = _select_unsupported_transition_pair(
            unconstrained_parent,
            unconstrained_child,
            allowed_transition_pairs=allowed_transition_pairs,
        )
        if transition_pair is not None:
            (
                unconstrained_source_region,
                unconstrained_target_region,
                unconstrained_support,
            ) = transition_pair
            constrained_parent = constrained_by_node[parent_signature]
            constrained_child = constrained_by_node[child_signature]
            constrained_transition_pair = (
                constrained_parent.most_likely_state,
                constrained_child.most_likely_state,
            )
            rows.append(
                UnsupportedGeographicTransitionClaimRow(
                    parent_node=parent_signature,
                    child_node=child_signature,
                    descendant_taxa=node_descendant_taxa(child),
                    unconstrained_source_region=unconstrained_source_region,
                    unconstrained_target_region=unconstrained_target_region,
                    unconstrained_support=unconstrained_support,
                    constrained_source_region=constrained_parent.most_likely_state,
                    constrained_target_region=constrained_child.most_likely_state,
                    constrained_support=stable_value(
                        min(
                            constrained_parent.confidence,
                            constrained_child.confidence,
                        )
                    ),
                    claim_resolved=constrained_transition_pair
                    != (
                        unconstrained_source_region,
                        unconstrained_target_region,
                    ),
                )
            )
        for grandchild in child.children:
            visit(child, grandchild)

    for child in root.children:
        visit(root, child)
    return rows


def _select_unsupported_transition_pair(
    parent_estimate,
    child_estimate,
    *,
    allowed_transition_pairs: set[tuple[str, str]],
) -> tuple[str, str, float] | None:
    primary_transition_pair = (
        parent_estimate.most_likely_state,
        child_estimate.most_likely_state,
    )
    if primary_transition_pair[0] == primary_transition_pair[1]:
        return None
    primary_support = stable_value(
        min(
            parent_estimate.state_probabilities.get(primary_transition_pair[0], 0.0),
            child_estimate.state_probabilities.get(primary_transition_pair[1], 0.0),
        )
    )
    if primary_transition_pair not in allowed_transition_pairs:
        return (
            primary_transition_pair[0],
            primary_transition_pair[1],
            primary_support,
        )
    strongest_forbidden_pair: tuple[str, str] | None = None
    strongest_forbidden_support = 0.0
    for (
        source_region,
        source_probability,
    ) in parent_estimate.state_probabilities.items():
        if source_probability <= 0.0:
            continue
        for (
            target_region,
            target_probability,
        ) in child_estimate.state_probabilities.items():
            if target_probability <= 0.0 or source_region == target_region:
                continue
            candidate_pair = (source_region, target_region)
            if candidate_pair in allowed_transition_pairs:
                continue
            candidate_support = stable_value(
                min(source_probability, target_probability)
            )
            if candidate_support > strongest_forbidden_support or (
                candidate_support == strongest_forbidden_support
                and strongest_forbidden_pair is not None
                and candidate_pair < strongest_forbidden_pair
            ):
                strongest_forbidden_pair = candidate_pair
                strongest_forbidden_support = candidate_support
    if strongest_forbidden_pair is None:
        return None
    if strongest_forbidden_support < (
        primary_support * _UNSUPPORTED_CLAIM_COMPETITOR_RATIO
    ):
        return None
    return (
        strongest_forbidden_pair[0],
        strongest_forbidden_pair[1],
        strongest_forbidden_support,
    )


def _load_geographic_adjacency_matrix(path: Path) -> _GeographicAdjacencyMatrix:
    if not path.exists():
        raise FileNotFoundError(f"geographic adjacency matrix file not found: {path}")
    delimiter = _detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = [row for row in reader if any(cell.strip() for cell in row)]
    if not rows:
        raise ValueError("geographic adjacency matrix is empty")
    header = [cell.strip() for cell in rows[0]]
    if len(header) < 2:
        raise ValueError(
            "geographic adjacency matrix must include one row-label column and at least one region column"
        )
    region_order = header[1:]
    if len(set(region_order)) != len(region_order):
        raise ValueError(
            "geographic adjacency matrix contains duplicate region columns"
        )
    row_lookup: dict[str, list[str]] = {}
    for row in rows[1:]:
        normalized = [cell.strip() for cell in row]
        if len(normalized) != len(region_order) + 1:
            raise ValueError(
                "geographic adjacency matrix rows must match the declared region columns"
            )
        source_region = normalized[0]
        if not source_region:
            raise ValueError("geographic adjacency matrix contains an empty row label")
        if source_region not in region_order:
            raise ValueError(
                "geographic adjacency matrix row label is not present in the header: "
                f"{source_region}"
            )
        if source_region in row_lookup:
            raise ValueError(
                "geographic adjacency matrix contains a duplicate row label: "
                f"{source_region}"
            )
        row_lookup[source_region] = normalized[1:]
    missing_rows = [region for region in region_order if region not in row_lookup]
    if missing_rows:
        raise ValueError(
            "geographic adjacency matrix is missing rows for: "
            + ", ".join(sorted(missing_rows))
        )
    allowed_transition_pairs: list[tuple[str, str]] = []
    for source_region in region_order:
        values = row_lookup[source_region]
        for target_region, raw_value in zip(region_order, values, strict=True):
            allowed = _parse_adjacency_value(
                raw_value,
                source_region=source_region,
                target_region=target_region,
            )
            if allowed:
                allowed_transition_pairs.append((source_region, target_region))
    if not allowed_transition_pairs:
        raise ValueError(
            "geographic adjacency matrix does not allow any directed transitions"
        )
    return _GeographicAdjacencyMatrix(
        path=path,
        region_order=region_order,
        allowed_transition_pairs=allowed_transition_pairs,
    )


def _parse_adjacency_value(
    raw_value: str,
    *,
    source_region: str,
    target_region: str,
) -> bool:
    normalized = raw_value.strip().lower()
    if source_region == target_region:
        if normalized in _FALSE_MATRIX_VALUES:
            return False
        raise ValueError(
            "geographic adjacency matrix diagonal must be empty or false-like for "
            f"{source_region}"
        )
    if normalized in _TRUE_MATRIX_VALUES:
        return True
    if normalized in _FALSE_MATRIX_VALUES:
        return False
    raise ValueError(
        "geographic adjacency matrix values must be 0/1 or false/true-like, got "
        f"'{raw_value}' for {source_region}->{target_region}"
    )


def _detect_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    header = path.read_text(encoding="utf-8").splitlines()[0] if path.exists() else ""
    return "\t" if "\t" in header else ","


def _resolve_internal_model(model: str) -> str:
    try:
        return _MODEL_ALIAS_TO_INTERNAL[model]
    except KeyError as error:
        raise ValueError(
            "constrained geography model must be one of er, sym, ard, equal-rates, symmetric, or all-rates-different"
        ) from error


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


def _validate_in_tree_region_labels(
    audit,
) -> None:
    invalid_rows = [
        row
        for row in audit.rows
        if row.in_tree and not row.included and row.issue_code not in {"missing-state"}
    ]
    if not invalid_rows:
        return
    invalid_labels = sorted(
        {row.normalized_state or row.raw_state for row in invalid_rows}
    )
    raise AncestralReconstructionError(
        "geographic adjacency matrix does not define one or more analyzed region labels: "
        + ", ".join(invalid_labels)
    )
