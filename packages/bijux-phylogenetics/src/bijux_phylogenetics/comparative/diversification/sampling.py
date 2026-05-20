from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree

from .models import SamplingFractionIssue, SamplingFractionReport

_SAMPLING_COLUMNS = (
    "sampling_fraction",
    "sampling_proportion",
    "sampling_probability",
    "sampling_prob",
)


def _sampling_fraction_from_rows(rows: list[float]) -> float:
    return float(format(sum(rows) / max(len(rows), 1), ".15g"))


def _resolve_sampling_column(columns: list[str], requested: str | None) -> str | None:
    if requested is not None:
        return requested if requested in columns else None
    return next((column for column in _SAMPLING_COLUMNS if column in columns), None)


def detect_incomplete_taxon_sampling_metadata(
    tree_path: Path,
    metadata_path: Path,
    *,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
) -> SamplingFractionReport:
    """Inspect taxon sampling fractions keyed to the tree tips."""
    tree = load_tree(tree_path)
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    resolved_sampling_column = _resolve_sampling_column(table.columns, sampling_column)
    if resolved_sampling_column is None:
        return SamplingFractionReport(
            tree_path=tree_path,
            metadata_path=metadata_path,
            taxon_column=table.taxon_column,
            sampling_column=None,
            complete=False,
            matched_taxa=[],
            missing_taxa=sorted(tree.tip_names),
            invalid_rows=[],
            sampling_fraction=None,
            heterogeneous_values=False,
            warnings=["metadata does not declare a sampling-fraction column"],
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    matched_taxa: list[str] = []
    missing_taxa: list[str] = []
    invalid_rows: list[SamplingFractionIssue] = []
    fractions: list[float] = []
    for taxon in sorted(tree.tip_names):
        row = rows_by_taxon.get(taxon)
        if row is None:
            missing_taxa.append(taxon)
            continue
        matched_taxa.append(taxon)
        raw_value = row[resolved_sampling_column].strip()
        if not raw_value:
            invalid_rows.append(
                SamplingFractionIssue(
                    taxon=taxon,
                    code="missing-sampling-fraction",
                    raw_value=raw_value,
                    message="sampling fraction is missing for this taxon",
                )
            )
            continue
        try:
            value = float(raw_value)
        except ValueError:
            invalid_rows.append(
                SamplingFractionIssue(
                    taxon=taxon,
                    code="invalid-sampling-fraction",
                    raw_value=raw_value,
                    message="sampling fraction must be numeric",
                )
            )
            continue
        if value <= 0.0 or value > 1.0:
            invalid_rows.append(
                SamplingFractionIssue(
                    taxon=taxon,
                    code="out-of-range-sampling-fraction",
                    raw_value=raw_value,
                    message="sampling fraction must be greater than 0 and at most 1",
                )
            )
            continue
        fractions.append(value)

    heterogeneous_values = len({format(value, ".9g") for value in fractions}) > 1
    warnings: list[str] = []
    if missing_taxa:
        warnings.append("sampling metadata does not cover every tree tip")
    if invalid_rows:
        warnings.append(
            "one or more sampling fractions are missing, invalid, or out of range"
        )
    if heterogeneous_values:
        warnings.append(
            "sampling fractions vary across taxa; mean sampling fraction will be used when correction is applied"
        )
    return SamplingFractionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=table.taxon_column,
        sampling_column=resolved_sampling_column,
        complete=not missing_taxa and not invalid_rows,
        matched_taxa=matched_taxa,
        missing_taxa=missing_taxa,
        invalid_rows=invalid_rows,
        sampling_fraction=_sampling_fraction_from_rows(fractions)
        if fractions
        else None,
        heterogeneous_values=heterogeneous_values,
        warnings=warnings,
    )
