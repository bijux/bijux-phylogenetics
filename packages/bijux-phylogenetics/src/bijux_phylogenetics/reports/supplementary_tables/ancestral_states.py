from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.ancestral import (
    ContinuousAncestralReport,
    DiscreteAncestralReport,
    continuous_ancestral_exclusions,
    discrete_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
)

from .columns import ancestral_state_table_columns
from .models import (
    SupplementaryAncestralStateRow,
    SupplementaryAncestralStateTableResult,
)
from .shared import stringify_list, write_dict_rows


def _serialize_ancestral_state_row(
    row: SupplementaryAncestralStateRow,
) -> dict[str, str]:
    return {
        "tree_source": row.tree_source,
        "traits_source": row.traits_source,
        "trait": row.trait,
        "reconstruction_kind": row.reconstruction_kind,
        "model": row.model,
        "estimator": "" if row.estimator is None else row.estimator,
        "state_ordering": "" if row.state_ordering is None else row.state_ordering,
        "root_prior_mode": ("" if row.root_prior_mode is None else row.root_prior_mode),
        "fixed_root_state": (
            "" if row.fixed_root_state is None else row.fixed_root_state
        ),
        "alpha": "" if row.alpha is None else str(row.alpha),
        "analysis_taxon_count": str(row.analysis_taxon_count),
        "excluded_taxon_count": str(row.excluded_taxon_count),
        "excluded_taxa": stringify_list(row.excluded_taxa),
        "warning_count": str(row.warning_count),
        "warnings": stringify_list(row.warnings),
        "node": row.node,
        "node_name": "" if row.node_name is None else row.node_name,
        "descendant_taxa": stringify_list(row.descendant_taxa),
        "descendant_taxon_count": str(row.descendant_taxon_count),
        "estimate_value": "" if row.estimate_value is None else str(row.estimate_value),
        "most_likely_state": (
            "" if row.most_likely_state is None else row.most_likely_state
        ),
        "state_set": stringify_list(row.state_set),
        "state_probabilities": json.dumps(row.state_probabilities, sort_keys=True),
        "standard_error": "" if row.standard_error is None else str(row.standard_error),
        "lower_95_interval": (
            "" if row.lower_95_interval is None else str(row.lower_95_interval)
        ),
        "upper_95_interval": (
            "" if row.upper_95_interval is None else str(row.upper_95_interval)
        ),
        "confidence": str(row.confidence),
        "ambiguous": "" if row.ambiguous is None else str(row.ambiguous).lower(),
        "unstable": str(row.unstable).lower(),
        "interpretation": row.interpretation,
        "downstream_risks": stringify_list(row.downstream_risks),
    }


def _write_ancestral_state_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryAncestralStateRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_ancestral_state_row(row) for row in rows],
    )


def _build_continuous_ancestral_state_rows(
    report: ContinuousAncestralReport,
) -> list[SupplementaryAncestralStateRow]:
    excluded_taxa = [
        f"{row.taxon}:{row.reason}" for row in continuous_ancestral_exclusions(report)
    ]
    return [
        SupplementaryAncestralStateRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            trait=report.trait,
            reconstruction_kind="continuous",
            model=report.model,
            estimator=report.estimator,
            state_ordering=None,
            root_prior_mode=None,
            fixed_root_state=None,
            alpha=report.alpha,
            analysis_taxon_count=report.taxon_count,
            excluded_taxon_count=len(excluded_taxa),
            excluded_taxa=excluded_taxa,
            warning_count=len(report.warnings),
            warnings=list(report.warnings),
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            descendant_taxon_count=len(estimate.descendant_taxa),
            estimate_value=estimate.estimate,
            most_likely_state=None,
            state_set=[],
            state_probabilities={},
            standard_error=estimate.standard_error,
            lower_95_interval=estimate.lower_95_interval,
            upper_95_interval=estimate.upper_95_interval,
            confidence=estimate.confidence,
            ambiguous=None,
            unstable=estimate.unstable,
            interpretation=estimate.interpretation,
            downstream_risks=list(estimate.downstream_risks),
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _build_discrete_ancestral_state_rows(
    report: DiscreteAncestralReport,
) -> list[SupplementaryAncestralStateRow]:
    excluded_taxa = [
        f"{row.taxon}:{row.reason}" for row in discrete_ancestral_exclusions(report)
    ]
    return [
        SupplementaryAncestralStateRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            trait=report.trait,
            reconstruction_kind="discrete",
            model=report.model,
            estimator=None,
            state_ordering=report.state_ordering,
            root_prior_mode=report.root_prior_mode,
            fixed_root_state=report.fixed_root_state,
            alpha=None,
            analysis_taxon_count=report.taxon_count,
            excluded_taxon_count=len(excluded_taxa),
            excluded_taxa=excluded_taxa,
            warning_count=len(report.warnings),
            warnings=list(report.warnings),
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            descendant_taxon_count=len(estimate.descendant_taxa),
            estimate_value=None,
            most_likely_state=estimate.most_likely_state,
            state_set=list(estimate.state_set),
            state_probabilities=dict(estimate.state_probabilities),
            standard_error=None,
            lower_95_interval=None,
            upper_95_interval=None,
            confidence=estimate.confidence,
            ambiguous=estimate.ambiguous,
            unstable=estimate.unstable,
            interpretation=estimate.interpretation,
            downstream_risks=list(estimate.downstream_risks),
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def write_supplementary_ancestral_state_table(
    path: Path,
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    reconstruction_kind: str,
    taxon_column: str | None = None,
    model: str | None = None,
    estimator: str | None = None,
    alpha: float = 1.0,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
) -> SupplementaryAncestralStateTableResult:
    """Write one supplementary internal-node ancestral-state table."""
    if reconstruction_kind == "continuous":
        report = reconstruct_continuous_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model="brownian" if model is None else model,
            estimator=estimator,
            alpha=alpha,
        )
        rows = _build_continuous_ancestral_state_rows(report)
        resolved_model = report.model
        unstable_node_count = len(report.unstable_nodes)
    elif reconstruction_kind == "discrete":
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model="fitch" if model is None else model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
        )
        rows = _build_discrete_ancestral_state_rows(report)
        resolved_model = report.model
        unstable_node_count = len(report.unstable_nodes)
    else:
        raise ValueError("reconstruction_kind must be 'continuous' or 'discrete'")
    columns = ancestral_state_table_columns()
    _write_ancestral_state_rows(path, columns=columns, rows=rows)
    return SupplementaryAncestralStateTableResult(
        output_path=path,
        row_count=len(rows),
        reconstruction_kind=reconstruction_kind,
        model=resolved_model,
        analysis_taxon_count=0 if not rows else rows[0].analysis_taxon_count,
        excluded_taxon_count=0 if not rows else rows[0].excluded_taxon_count,
        unstable_node_count=unstable_node_count,
        columns=columns,
        rows=rows,
    )
