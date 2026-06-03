from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import ComparativeMethodReport

from .contracts import (
    ComparativeAnalysisSummaryRow,
    ComparativeAuditTableRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
)


def _table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _write_rows(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter=_table_delimiter(path)
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_comparative_summary_table(
    path: Path, row: ComparativeAnalysisSummaryRow
) -> Path:
    return _write_rows(path, list(asdict(row).keys()), [asdict(row)])


def write_comparative_coefficient_table(
    path: Path, rows: list[ComparativeCoefficientTableRow]
) -> Path:
    return _write_rows(
        path, list(asdict(rows[0]).keys()), [asdict(row) for row in rows]
    )


def write_comparative_residual_table(
    path: Path, rows: list[ComparativeResidualTableRow]
) -> Path:
    rendered = []
    for row in rows:
        payload = asdict(row)
        payload["outlier_taxa"] = "|".join(row.outlier_taxa)
        payload["high_leverage_taxa"] = "|".join(row.high_leverage_taxa)
        payload["warnings"] = "|".join(row.warnings)
        rendered.append(payload)
    return _write_rows(path, list(rendered[0].keys()), rendered)


def write_comparative_signal_table(path: Path, row: ComparativeSignalTableRow) -> Path:
    return _write_rows(path, list(asdict(row).keys()), [asdict(row)])


def write_comparative_model_comparison_table(
    path: Path, report: ComparativeMethodReport
) -> Path:
    rows = [asdict(row) for row in report.snapshot.model_comparison.rows]
    return _write_rows(path, list(rows[0].keys()), rows)


def write_comparative_interpretation_table(
    path: Path, rows: list[ComparativeInterpretationRow]
) -> Path:
    return _write_rows(
        path, list(asdict(rows[0]).keys()), [asdict(row) for row in rows]
    )


def write_comparative_audit_table(
    path: Path, rows: list[ComparativeAuditTableRow]
) -> Path:
    rendered = []
    for row in rows:
        payload = asdict(row)
        payload["taxa_used"] = "|".join(row.taxa_used)
        payload["traits_used"] = "|".join(row.traits_used)
        payload["excluded_taxa"] = "|".join(row.excluded_taxa)
        payload["assumptions"] = "|".join(row.assumptions)
        payload["warnings"] = "|".join(row.warnings)
        rendered.append(payload)
    return _write_rows(path, list(rendered[0].keys()), rendered)


def write_comparative_contrast_table(
    path: Path, report: ComparativeMethodReport
) -> Path:
    rows = [
        {
            "node": row.node,
            "left_taxa": "|".join(row.left_taxa),
            "right_taxa": "|".join(row.right_taxa),
            "contrast": row.contrast,
            "expected_variance": row.expected_variance,
            "ancestral_value": row.ancestral_value,
        }
        for row in report.snapshot.contrasts.contrasts
    ]
    return _write_rows(path, list(rows[0].keys()), rows)
