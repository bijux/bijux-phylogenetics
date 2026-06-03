from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.reports.publication.tree import (
    TreeSupportRow,
    summarize_tree_support,
)
from bijux_phylogenetics.trees import (
    TreeSetCladeSupportReport,
    TreeSetCladeSupportRow,
    compute_reference_tree_clade_support,
    extract_tree_clades,
)

from .columns import clade_support_table_columns
from .models import SupplementaryCladeSupportRow, SupplementaryCladeSupportTableResult
from .shared import stringify_list, write_dict_rows


def _serialize_clade_support_row(
    row: SupplementaryCladeSupportRow,
) -> dict[str, object]:
    return {
        "tree_source": row.tree_source,
        "comparison_tree_set_source": ""
        if row.comparison_tree_set_source is None
        else row.comparison_tree_set_source,
        "clade_id": row.clade_id,
        "node_kind": row.node_kind,
        "node_label": "" if row.node_label is None else row.node_label,
        "descendant_taxa": stringify_list(row.descendant_taxa),
        "support": "" if row.support is None else row.support,
        "support_fraction": ""
        if row.support_fraction is None
        else row.support_fraction,
        "support_class": row.support_class,
        "support_method": row.support_method,
        "branch_length": "" if row.branch_length is None else row.branch_length,
        "root_depth": "" if row.root_depth is None else row.root_depth,
        "supporting_tree_count": ""
        if row.supporting_tree_count is None
        else row.supporting_tree_count,
        "clade_frequency": "" if row.clade_frequency is None else row.clade_frequency,
        "support_percent": "" if row.support_percent is None else row.support_percent,
        "frequency_method": ""
        if row.frequency_method is None
        else row.frequency_method,
        "frequency_status": ""
        if row.frequency_status is None
        else row.frequency_status,
        "frequency_explanation": ""
        if row.frequency_explanation is None
        else row.frequency_explanation,
    }


def _clade_support_row_lookup(
    rows: list[TreeSetCladeSupportRow],
) -> dict[tuple[str, ...], TreeSetCladeSupportRow]:
    return {tuple(row.descendant_taxa): row for row in rows}


def _build_clade_support_row(
    *,
    tree_path: Path,
    comparison_tree_set_path: Path | None,
    support_row: TreeSupportRow,
    frequency_row: TreeSetCladeSupportRow | None,
) -> SupplementaryCladeSupportRow:
    return SupplementaryCladeSupportRow(
        tree_source=str(tree_path),
        comparison_tree_set_source=(
            None if comparison_tree_set_path is None else str(comparison_tree_set_path)
        ),
        clade_id=support_row.node,
        node_kind=support_row.node_kind,
        node_label=support_row.node_label,
        descendant_taxa=list(support_row.descendant_taxa),
        support=support_row.support,
        support_fraction=support_row.support_fraction,
        support_class=support_row.support_class,
        support_method="tree-label",
        branch_length=support_row.branch_length,
        root_depth=support_row.root_depth,
        supporting_tree_count=(
            None if frequency_row is None else frequency_row.supporting_tree_count
        ),
        clade_frequency=None
        if frequency_row is None
        else frequency_row.clade_frequency,
        support_percent=None
        if frequency_row is None
        else frequency_row.support_percent,
        frequency_method=(
            None if frequency_row is None else "reference-tree-clade-frequency"
        ),
        frequency_status=None
        if frequency_row is None
        else frequency_row.support_status,
        frequency_explanation=(
            None if frequency_row is None else frequency_row.explanation
        ),
    )


def _write_clade_support_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryCladeSupportRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_clade_support_row(row) for row in rows],
    )


def write_supplementary_clade_support_table(
    path: Path,
    *,
    tree_path: Path,
    comparison_tree_set_path: Path | None = None,
) -> SupplementaryCladeSupportTableResult:
    """Write one supplementary clade-support table from a reference tree and optional tree set."""
    clades = extract_tree_clades(tree_path)
    support_rows = summarize_tree_support(clades)
    frequency_report: TreeSetCladeSupportReport | None = None
    frequency_rows: dict[tuple[str, ...], TreeSetCladeSupportRow] = {}
    if comparison_tree_set_path is not None:
        frequency_report = compute_reference_tree_clade_support(
            tree_path,
            comparison_tree_set_path,
        )
        frequency_rows = _clade_support_row_lookup(frequency_report.rows)
    rows = [
        _build_clade_support_row(
            tree_path=tree_path,
            comparison_tree_set_path=comparison_tree_set_path,
            support_row=support_row,
            frequency_row=frequency_rows.get(tuple(support_row.descendant_taxa)),
        )
        for support_row in support_rows
    ]
    columns = clade_support_table_columns()
    _write_clade_support_rows(path, columns=columns, rows=rows)
    return SupplementaryCladeSupportTableResult(
        output_path=path,
        row_count=len(rows),
        supported_clade_count=sum(1 for row in rows if row.support is not None),
        frequency_scored_clade_count=sum(
            1 for row in rows if row.clade_frequency is not None
        ),
        frequency_partial_support_count=sum(
            1 for row in rows if row.frequency_status == "partial-support"
        ),
        frequency_absent_clade_count=sum(
            1 for row in rows if row.frequency_status == "absent"
        ),
        frequency_unscored_clade_count=sum(
            1 for row in rows if row.frequency_status == "not-counted"
        ),
        columns=columns,
        rows=rows,
    )
