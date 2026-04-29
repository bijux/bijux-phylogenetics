from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.common import descendant_taxa
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class CalibrationValidationIssue:
    calibration_id: str
    code: str
    message: str


@dataclass(slots=True)
class ValidatedCalibration:
    calibration_id: str
    target_kind: str
    target_label: str
    taxa: list[str]
    minimum_age: float | None
    maximum_age: float | None
    distribution: str
    valid: bool


@dataclass(slots=True)
class FossilCalibrationValidationReport:
    tree_path: Path
    calibration_path: Path
    tree_taxa: list[str]
    calibration_count: int
    valid_calibration_count: int
    invalid_calibration_count: int
    calibrations: list[ValidatedCalibration]
    issues: list[CalibrationValidationIssue]


@dataclass(slots=True)
class ImpossibleCalibrationConstraintReport:
    tree_path: Path
    calibration_path: Path
    impossible_calibration_ids: list[str]
    issues: list[CalibrationValidationIssue]


def _read_delimited_rows(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"calibration table contains no header: {path}")
        return [{key: value or "" for key, value in row.items() if key is not None} for row in reader]


def _named_clades(tree: PhyloTree) -> dict[str, list[str]]:
    named: dict[str, list[str]] = {}
    for node in tree.iter_nodes():
        if node.name:
            named[node.name] = descendant_taxa(node)
    return named


def _clade_taxon_sets(tree: PhyloTree) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()
    for node in tree.iter_nodes():
        taxa = descendant_taxa(node)
        if taxa:
            clades.add(frozenset(taxa))
    return clades


def _parse_target_taxa(raw: str) -> list[str]:
    if not raw.strip():
        return []
    normalized = raw.replace(",", "|").replace(";", "|")
    return sorted({token.strip() for token in normalized.split("|") if token.strip()})


def _parse_age(raw: str, *, calibration_id: str, field_name: str, issues: list[CalibrationValidationIssue]) -> float | None:
    if not raw.strip():
        return None
    try:
        value = float(raw)
    except ValueError:
        issues.append(
            CalibrationValidationIssue(
                calibration_id=calibration_id,
                code="invalid-age",
                message=f"{field_name} must be numeric when provided",
            )
        )
        return None
    return value


def validate_fossil_calibration_table(tree_path: Path, calibration_path: Path) -> FossilCalibrationValidationReport:
    """Validate fossil calibration targets and age constraints against one tree."""
    tree = load_tree(tree_path)
    rows = _read_delimited_rows(calibration_path)
    named_clades = _named_clades(tree)
    clade_sets = _clade_taxon_sets(tree)
    tree_taxa = sorted(tree.tip_names)
    tree_taxa_set = set(tree_taxa)
    calibrations: list[ValidatedCalibration] = []
    issues: list[CalibrationValidationIssue] = []

    for index, row in enumerate(rows, start=1):
        calibration_id = row.get("calibration_id", "").strip() or f"calibration-{index}"
        clade_name = row.get("clade_name", "").strip()
        taxa = _parse_target_taxa(row.get("taxa", ""))
        distribution = row.get("distribution", "").strip() or "uniform"
        minimum_age = _parse_age(row.get("minimum_age", ""), calibration_id=calibration_id, field_name="minimum_age", issues=issues)
        maximum_age = _parse_age(row.get("maximum_age", ""), calibration_id=calibration_id, field_name="maximum_age", issues=issues)

        target_kind = "taxa"
        target_label = clade_name or "|".join(taxa)
        resolved_taxa: list[str] = []
        valid = True

        if clade_name:
            target_kind = "named-clade"
            resolved_taxa = named_clades.get(clade_name, [])
            if not resolved_taxa:
                valid = False
                issues.append(
                    CalibrationValidationIssue(
                        calibration_id=calibration_id,
                        code="unknown-clade-name",
                        message=f"named clade '{clade_name}' is not present in the tree",
                    )
                )
        elif taxa:
            missing_taxa = sorted(set(taxa) - tree_taxa_set)
            if missing_taxa:
                valid = False
                issues.append(
                    CalibrationValidationIssue(
                        calibration_id=calibration_id,
                        code="unknown-taxa",
                        message=f"calibration taxa are absent from the tree: {', '.join(missing_taxa)}",
                    )
                )
            elif frozenset(taxa) not in clade_sets:
                valid = False
                issues.append(
                    CalibrationValidationIssue(
                        calibration_id=calibration_id,
                        code="non-monophyletic-target",
                        message="calibration taxa do not map to a named clade, single taxon, or monophyletic descendant set",
                    )
                )
            resolved_taxa = taxa
        else:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="missing-target",
                    message="calibration must provide either a clade_name or taxa column",
                )
            )

        if minimum_age is not None and minimum_age < 0.0:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="negative-minimum-age",
                    message="minimum calibration age cannot be negative",
                )
            )
        if maximum_age is not None and maximum_age < 0.0:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="negative-maximum-age",
                    message="maximum calibration age cannot be negative",
                )
            )
        if minimum_age is not None and maximum_age is not None and minimum_age > maximum_age:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="minimum-exceeds-maximum",
                    message="minimum calibration age cannot exceed maximum calibration age",
                )
            )
        if minimum_age is None and maximum_age is None:
            valid = False
            issues.append(
                CalibrationValidationIssue(
                    calibration_id=calibration_id,
                    code="missing-age-bounds",
                    message="calibration must provide at least one age bound",
                )
            )

        calibrations.append(
            ValidatedCalibration(
                calibration_id=calibration_id,
                target_kind=target_kind,
                target_label=target_label,
                taxa=resolved_taxa,
                minimum_age=minimum_age,
                maximum_age=maximum_age,
                distribution=distribution,
                valid=valid,
            )
        )

    valid_count = sum(1 for calibration in calibrations if calibration.valid)
    return FossilCalibrationValidationReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        tree_taxa=tree_taxa,
        calibration_count=len(calibrations),
        valid_calibration_count=valid_count,
        invalid_calibration_count=len(calibrations) - valid_count,
        calibrations=calibrations,
        issues=issues,
    )


def detect_impossible_calibration_constraints(tree_path: Path, calibration_path: Path) -> ImpossibleCalibrationConstraintReport:
    """Return calibrations with impossible target or age constraints."""
    report = validate_fossil_calibration_table(tree_path, calibration_path)
    impossible_codes = {
        "unknown-clade-name",
        "unknown-taxa",
        "non-monophyletic-target",
        "missing-target",
        "negative-minimum-age",
        "negative-maximum-age",
        "minimum-exceeds-maximum",
        "missing-age-bounds",
    }
    impossible_ids = sorted(
        {
            issue.calibration_id
            for issue in report.issues
            if issue.code in impossible_codes
        }
    )
    return ImpossibleCalibrationConstraintReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        impossible_calibration_ids=impossible_ids,
        issues=[issue for issue in report.issues if issue.calibration_id in impossible_ids],
    )
