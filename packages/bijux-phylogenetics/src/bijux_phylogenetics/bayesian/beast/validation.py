from __future__ import annotations

from ._shared import (
    Path,
    _clade_taxon_sets,
    _ensure_inference_ready_alignment,
    _named_clades,
    _parse_age,
    _parse_target_taxa,
    _read_delimited_rows,
    _tree_root_age,
    load_fasta_alignment,
    load_taxon_table,
    load_tree,
    validate_tree_path,
)
from .models import (
    CalibrationDominanceObservation,
    CalibrationDominanceReport,
    CalibrationValidationIssue,
    FossilCalibrationValidationReport,
    ImpossibleCalibrationConstraintReport,
    TimeTreeReadinessReport,
    TipDatingValidationIssue,
    TipDatingValidationReport,
    ValidatedCalibration,
    ValidatedTipDate,
)


def validate_fossil_calibration_table(
    tree_path: Path, calibration_path: Path
) -> FossilCalibrationValidationReport:
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
        minimum_age = _parse_age(
            row.get("minimum_age", ""),
            calibration_id=calibration_id,
            field_name="minimum_age",
            issues=issues,
        )
        maximum_age = _parse_age(
            row.get("maximum_age", ""),
            calibration_id=calibration_id,
            field_name="maximum_age",
            issues=issues,
        )

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
        if (
            minimum_age is not None
            and maximum_age is not None
            and minimum_age > maximum_age
        ):
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


def detect_impossible_calibration_constraints(
    tree_path: Path, calibration_path: Path
) -> ImpossibleCalibrationConstraintReport:
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
        issues=[
            issue for issue in report.issues if issue.calibration_id in impossible_ids
        ],
    )


def assess_calibration_dominance(
    tree_path: Path, calibration_path: Path
) -> CalibrationDominanceReport:
    """Flag cases where one calibration contributes a disproportionate share of the dated-tree age range."""
    report = validate_fossil_calibration_table(tree_path, calibration_path)
    root_age = _tree_root_age(tree_path)
    observations: list[CalibrationDominanceObservation] = []
    dominant_calibration_ids: list[str] = []
    warnings: list[str] = []
    for calibration in report.calibrations:
        if not calibration.valid:
            continue
        span_fraction: float | None = None
        dominates_root_age = False
        warning: str | None = None
        if (
            calibration.minimum_age is not None
            and calibration.maximum_age is not None
            and root_age > 0.0
        ):
            span_fraction = round(
                (calibration.maximum_age - calibration.minimum_age) / root_age, 15
            )
            dominates_root_age = span_fraction >= 0.5
            if dominates_root_age:
                warning = "calibration age span covers at least half of the current root-age scale"
        elif root_age > 0.0 and (
            (
                calibration.minimum_age is not None
                and calibration.minimum_age / root_age >= 0.8
            )
            or (
                calibration.maximum_age is not None
                and calibration.maximum_age / root_age >= 0.8
            )
        ):
            dominates_root_age = True
            warning = "single-sided calibration bound lies close to the current root-age scale"
        if dominates_root_age:
            dominant_calibration_ids.append(calibration.calibration_id)
            if warning is not None:
                warnings.append(f"{calibration.calibration_id}: {warning}")
        observations.append(
            CalibrationDominanceObservation(
                calibration_id=calibration.calibration_id,
                target_label=calibration.target_label,
                bounded_span_fraction=span_fraction,
                dominates_root_age=dominates_root_age,
                warning=warning,
            )
        )
    if report.valid_calibration_count == 1 and report.calibrations:
        only_calibration = next(
            (calibration for calibration in report.calibrations if calibration.valid),
            None,
        )
        if (
            only_calibration is not None
            and only_calibration.calibration_id not in dominant_calibration_ids
        ):
            dominant_calibration_ids.append(only_calibration.calibration_id)
        warnings.append(
            "only one valid calibration remains, so time estimates are effectively driven by a single calibration target"
        )
    return CalibrationDominanceReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        root_age=root_age,
        valid_calibration_count=report.valid_calibration_count,
        dominant_calibration_ids=sorted(dominant_calibration_ids),
        observations=sorted(observations, key=lambda row: row.calibration_id),
        warnings=sorted(dict.fromkeys(warnings)),
    )


def validate_tip_dating_metadata(
    tree_path: Path,
    tip_dates_path: Path,
    *,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> TipDatingValidationReport:
    """Validate that dated tips resolve cleanly against the tree and optional alignment."""
    tree = load_tree(tree_path)
    tree_taxa = set(tree.tip_names)
    alignment_taxa: set[str] | None = None
    if alignment_path is not None:
        _ensure_inference_ready_alignment(alignment_path)
        alignment_taxa = {
            record.identifier for record in load_fasta_alignment(alignment_path)
        }
    table = load_taxon_table(tip_dates_path, taxon_column=taxon_column)
    if date_column not in table.columns:
        raise ValueError(f"tip-dating table does not contain column '{date_column}'")

    issues: list[TipDatingValidationIssue] = []
    tip_dates: list[ValidatedTipDate] = []
    extra_tip_taxa = sorted(set(table.taxa) - tree_taxa)
    missing_tree_taxa = sorted(tree_taxa - set(table.taxa))
    extra_alignment_taxa = (
        sorted(set(table.taxa) - alignment_taxa) if alignment_taxa is not None else []
    )
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_date = row.get(date_column, "")
        valid = True
        parsed_date: float | None = None
        if taxon not in tree_taxa:
            valid = False
            issues.append(
                TipDatingValidationIssue(
                    taxon=taxon,
                    code="taxon-missing-from-tree",
                    message="dated tip is absent from the tree",
                )
            )
        if alignment_taxa is not None and taxon not in alignment_taxa:
            valid = False
            issues.append(
                TipDatingValidationIssue(
                    taxon=taxon,
                    code="taxon-missing-from-alignment",
                    message="dated tip is absent from the alignment",
                )
            )
        if not raw_date.strip():
            valid = False
            issues.append(
                TipDatingValidationIssue(
                    taxon=taxon,
                    code="missing-date",
                    message="dated tip requires a numeric date value",
                )
            )
        else:
            try:
                parsed_date = float(raw_date)
            except ValueError:
                valid = False
                issues.append(
                    TipDatingValidationIssue(
                        taxon=taxon,
                        code="invalid-date",
                        message="dated tip value must be numeric",
                    )
                )
        tip_dates.append(ValidatedTipDate(taxon=taxon, date=parsed_date, valid=valid))

    valid_tip_count = sum(1 for tip in tip_dates if tip.valid)
    return TipDatingValidationReport(
        tree_path=tree_path,
        tip_dates_path=tip_dates_path,
        alignment_path=alignment_path,
        taxon_column=table.taxon_column,
        date_column=date_column,
        valid_tip_count=valid_tip_count,
        invalid_tip_count=len(tip_dates) - valid_tip_count,
        missing_tree_taxa=missing_tree_taxa,
        extra_tip_taxa=extra_tip_taxa,
        extra_alignment_taxa=extra_alignment_taxa,
        tip_dates=tip_dates,
        issues=issues,
    )


def assess_time_tree_readiness(
    tree_path: Path,
    *,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    alignment_path: Path | None = None,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> TimeTreeReadinessReport:
    """Decide whether a dataset is suitable for dated phylogenetics."""
    validation = validate_tree_path(tree_path)
    blockers: list[str] = []
    warnings: list[str] = []
    if not validation.rooted:
        blockers.append("time-tree analysis requires a rooted tree")
    if validation.ultrametric is not True:
        blockers.append("time-tree analysis requires ultrametric branch lengths")
    if validation.branch_length_status != "complete":
        blockers.append("time-tree analysis requires complete branch lengths")
    calibration_report = None
    calibration_dominance = None
    if calibration_path is not None:
        calibration_report = validate_fossil_calibration_table(
            tree_path, calibration_path
        )
        if calibration_report.invalid_calibration_count:
            blockers.append(
                "calibration table contains invalid fossil calibration targets or ages"
            )
        calibration_dominance = assess_calibration_dominance(
            tree_path, calibration_path
        )
        warnings.extend(calibration_dominance.warnings)
    tip_date_report = None
    if tip_dates_path is not None:
        tip_date_report = validate_tip_dating_metadata(
            tree_path,
            tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
            date_column=date_column,
        )
        if tip_date_report.invalid_tip_count:
            blockers.append(
                "tip-date table contains missing, invalid, or mismatched dated taxa"
            )
    if calibration_path is None and tip_dates_path is None:
        warnings.append(
            "no calibrations or tip dates were supplied, so the tree cannot be dated by the current workflow"
        )
    decision = "ready"
    if blockers:
        decision = "blocked"
    elif warnings:
        decision = "risky"
    return TimeTreeReadinessReport(
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        decision=decision,
        rooted=validation.rooted,
        ultrametric=validation.ultrametric is True,
        branch_length_status=validation.branch_length_status,
        blockers=sorted(dict.fromkeys(blockers)),
        warnings=sorted(dict.fromkeys([*warnings, *validation.warnings])),
        calibration_report=calibration_report,
        tip_date_report=tip_date_report,
        calibration_dominance=calibration_dominance,
    )
