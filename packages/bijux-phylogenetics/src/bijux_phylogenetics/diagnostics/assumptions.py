from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path


@dataclass(frozen=True, slots=True)
class StandardizedSupportLabel:
    """One internal support-like label normalized into explicit support fields."""

    node: str
    raw_label: str
    raw_value: float
    scale: str
    support_fraction: float | None
    support_percent: float | None
    normalized_probability: float | None
    confidence_of_inference: str


@dataclass(slots=True)
class BranchLengthUnitReport:
    """Declared branch-length unit summary extracted from metadata when present."""

    metadata_path: Path
    column_name: str | None
    declared_unit: str | None
    normalized_unit: str | None
    compatible_with_time_tree: bool | None
    compatible_with_substitution_tree: bool | None
    conflicting_values: list[str]


@dataclass(slots=True)
class TreeAssumptionReport:
    """Assumption-oriented tree diagnostics with optional metadata-aware unit checks."""

    tree_path: Path
    metadata_path: Path | None
    standardized_support_labels: list[StandardizedSupportLabel]
    branch_length_units: BranchLengthUnitReport | None
    time_tree_compatible: bool
    substitution_tree_compatible: bool
    blockers: list[str]
    warnings: list[str]


_UNIT_COLUMNS = (
    "branch_length_units",
    "tree_branch_length_units",
    "branch_units",
    "units",
)
_TIME_UNITS = {"years", "year", "mya", "ma", "million-years", "million_years"}
_SUBSTITUTION_UNITS = {
    "substitutions-per-site",
    "substitutions_per_site",
    "subs/site",
    "expected-substitutions",
    "expected_substitutions",
}


def standardize_support_labels(tree_path: Path) -> list[StandardizedSupportLabel]:
    """Convert support-like internal labels into normalized support fields."""
    inspection = inspect_tree_path(tree_path)
    standardized: list[StandardizedSupportLabel] = []
    support_values = [
        row.numeric_value
        for row in inspection.likely_support_labels
        if row.numeric_value is not None
    ]
    has_fraction_scale = any(0.0 <= value <= 1.0 for value in support_values)
    has_percent_scale = any(1.0 < value <= 100.0 for value in support_values)
    for row in inspection.likely_support_labels:
        if row.numeric_value is None:
            continue
        normalized_probability: float | None = None
        confidence = "low"
        if 0.0 <= row.numeric_value <= 1.0:
            scale = "fraction"
            support_fraction = row.numeric_value
            support_percent = round(row.numeric_value * 100.0, 15)
            normalized_probability = support_fraction
            confidence = "medium" if has_percent_scale else "high"
            if row.numeric_value in {0.0, 1.0}:
                confidence = "medium" if not has_percent_scale else "low"
        elif 1.0 < row.numeric_value <= 100.0:
            scale = "percentage"
            support_fraction = round(row.numeric_value / 100.0, 15)
            support_percent = row.numeric_value
            normalized_probability = support_fraction
            confidence = "medium" if has_fraction_scale else "high"
        else:
            scale = "out-of-range"
            support_fraction = None
            support_percent = None
        standardized.append(
            StandardizedSupportLabel(
                node=row.node,
                raw_label=row.label,
                raw_value=row.numeric_value,
                scale=scale,
                support_fraction=support_fraction,
                support_percent=support_percent,
                normalized_probability=normalized_probability,
                confidence_of_inference=confidence,
            )
        )
    return standardized


def inspect_branch_length_units(
    metadata_path: Path,
    *,
    taxon_column: str | None = None,
) -> BranchLengthUnitReport:
    """Inspect metadata for declared branch-length units."""
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    column_name = next(
        (column for column in _UNIT_COLUMNS if column in table.columns), None
    )
    if column_name is None:
        return BranchLengthUnitReport(
            metadata_path=metadata_path,
            column_name=None,
            declared_unit=None,
            normalized_unit=None,
            compatible_with_time_tree=None,
            compatible_with_substitution_tree=None,
            conflicting_values=[],
        )

    observed_values = sorted(
        {row[column_name].strip() for row in table.rows if row[column_name].strip()}
    )
    if not observed_values:
        return BranchLengthUnitReport(
            metadata_path=metadata_path,
            column_name=column_name,
            declared_unit=None,
            normalized_unit=None,
            compatible_with_time_tree=None,
            compatible_with_substitution_tree=None,
            conflicting_values=[],
        )
    normalized_values = [value.lower().replace(" ", "-") for value in observed_values]
    if len(set(normalized_values)) > 1:
        return BranchLengthUnitReport(
            metadata_path=metadata_path,
            column_name=column_name,
            declared_unit=None,
            normalized_unit=None,
            compatible_with_time_tree=None,
            compatible_with_substitution_tree=None,
            conflicting_values=observed_values,
        )

    normalized_unit = normalized_values[0]
    return BranchLengthUnitReport(
        metadata_path=metadata_path,
        column_name=column_name,
        declared_unit=observed_values[0],
        normalized_unit=normalized_unit,
        compatible_with_time_tree=normalized_unit in _TIME_UNITS,
        compatible_with_substitution_tree=normalized_unit in _SUBSTITUTION_UNITS,
        conflicting_values=[],
    )


def assess_tree_assumptions(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
) -> TreeAssumptionReport:
    """Check whether a tree is compatible with common time-tree or substitution-tree assumptions."""
    inspection = inspect_tree_path(tree_path)
    standardized_support = standardize_support_labels(tree_path)
    unit_report = (
        inspect_branch_length_units(metadata_path, taxon_column=taxon_column)
        if metadata_path is not None
        else None
    )

    blockers: list[str] = []
    warnings: list[str] = []

    if inspection.branch_length_status != "complete":
        blockers.append("tree requires complete branch lengths")
    if inspection.zero_length_branch_count:
        warnings.append("tree contains zero-length branches")
    if any(
        warning.code == "negative_branch_lengths"
        for warning in inspection.tree_quality_warnings
    ):
        blockers.append("tree contains negative branch lengths")
    if unit_report is not None and unit_report.conflicting_values:
        blockers.append("metadata declares conflicting branch-length units")

    time_tree_compatible = inspection.is_ultrametric is True and not blockers
    substitution_tree_compatible = (
        inspection.branch_length_status == "complete"
        and not any(
            warning.code == "negative_branch_lengths"
            for warning in inspection.tree_quality_warnings
        )
    )

    if unit_report is not None and unit_report.normalized_unit is not None:
        if unit_report.compatible_with_time_tree is False:
            time_tree_compatible = False
            warnings.append("metadata declares non-time branch-length units")
        if unit_report.compatible_with_substitution_tree is False:
            substitution_tree_compatible = False
            warnings.append("metadata declares non-substitution branch-length units")
    elif metadata_path is not None:
        warnings.append("metadata does not declare branch-length units")

    if inspection.is_ultrametric is not True:
        time_tree_compatible = False
        warnings.append("tree is not ultrametric")

    return TreeAssumptionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        standardized_support_labels=standardized_support,
        branch_length_units=unit_report,
        time_tree_compatible=time_tree_compatible,
        substitution_tree_compatible=substitution_tree_compatible,
        blockers=blockers,
        warnings=warnings,
    )
