from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.diagnostics import TraceConvergenceReport, summarize_trace_convergence
from bijux_phylogenetics.bayesian.posterior import summarize_maximum_clade_credibility_tree
from bijux_phylogenetics.comparative.common import descendant_taxa
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.engines.workflows import _ensure_inference_ready_alignment


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


@dataclass(slots=True)
class ValidatedTipDate:
    taxon: str
    date: float | None
    valid: bool


@dataclass(slots=True)
class TipDatingValidationIssue:
    taxon: str
    code: str
    message: str


@dataclass(slots=True)
class TipDatingValidationReport:
    tree_path: Path
    tip_dates_path: Path
    alignment_path: Path | None
    taxon_column: str
    date_column: str
    valid_tip_count: int
    invalid_tip_count: int
    missing_tree_taxa: list[str]
    extra_tip_taxa: list[str]
    extra_alignment_taxa: list[str]
    tip_dates: list[ValidatedTipDate]
    issues: list[TipDatingValidationIssue]


@dataclass(slots=True)
class CalibrationDominanceObservation:
    calibration_id: str
    target_label: str
    bounded_span_fraction: float | None
    dominates_root_age: bool
    warning: str | None


@dataclass(slots=True)
class CalibrationDominanceReport:
    tree_path: Path
    calibration_path: Path
    root_age: float
    valid_calibration_count: int
    dominant_calibration_ids: list[str]
    observations: list[CalibrationDominanceObservation]
    warnings: list[str]


@dataclass(slots=True)
class TimeTreeReadinessReport:
    tree_path: Path
    calibration_path: Path | None
    tip_dates_path: Path | None
    decision: str
    rooted: bool
    ultrametric: bool
    branch_length_status: str
    blockers: list[str]
    warnings: list[str]
    calibration_report: FossilCalibrationValidationReport | None
    tip_date_report: TipDatingValidationReport | None
    calibration_dominance: CalibrationDominanceReport | None


@dataclass(slots=True)
class BeastPreparationReport:
    alignment_path: Path
    output_path: Path
    tree_path: Path | None
    calibration_path: Path | None
    tip_dates_path: Path | None
    taxon_count: int
    character_count: int
    inferred_alphabet: str
    clock_model: str
    tree_prior: str
    chain_length: int
    log_every: int
    calibration_count: int
    tip_date_count: int


@dataclass(slots=True)
class BeastLogRow:
    state: int
    values: dict[str, float]


@dataclass(slots=True)
class BeastLogReport:
    path: Path
    row_count: int
    columns: list[str]
    rows: list[BeastLogRow]


@dataclass(slots=True)
class BeastLogValidationIssue:
    code: str
    message: str
    row: int | None = None
    column: str | None = None


@dataclass(slots=True)
class BeastPosteriorLogValidationReport:
    path: Path
    row_count: int
    state_count: int
    required_columns: list[str]
    observed_columns: list[str]
    missing_columns: list[str]
    issues: list[BeastLogValidationIssue]
    valid: bool


@dataclass(slots=True)
class BeastBurninSensitivitySlice:
    burnin_fraction: float
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    selected_tree_index: int
    clade_credibility_score: float
    posterior_mean: float | None
    likelihood_mean: float | None
    tree_height_mean: float | None


@dataclass(slots=True)
class BeastBurninSensitivityReport:
    posterior_tree_path: Path
    log_path: Path | None
    slices: list[BeastBurninSensitivitySlice]
    changed_mcc_count: int
    warnings: list[str]


@dataclass(slots=True)
class BeastChainMixingIssue:
    path: Path | None
    parameter: str
    code: str
    message: str
    observed_value: float
    threshold: float


@dataclass(slots=True)
class BeastChainMixingReport:
    log_paths: list[Path]
    chain_count: int
    converged: bool
    issues: list[BeastChainMixingIssue]
    chain_summaries: list[BeastConvergenceReport]


@dataclass(slots=True)
class BeastConvergenceReport:
    path: Path
    sample_count: int
    converged: bool
    ess_threshold: float
    mean_shift_threshold: float
    warnings: list[dict[str, object]]
    parameter_summaries: list[dict[str, object]]


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


def _tree_root_age(tree_path: Path) -> float:
    tree = load_tree(tree_path)
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    if not lengths:
        return 0.0
    return round(max(float(length) for length in lengths), 15)


def assess_calibration_dominance(tree_path: Path, calibration_path: Path) -> CalibrationDominanceReport:
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
        if calibration.minimum_age is not None and calibration.maximum_age is not None and root_age > 0.0:
            span_fraction = round((calibration.maximum_age - calibration.minimum_age) / root_age, 15)
            dominates_root_age = span_fraction >= 0.5
            if dominates_root_age:
                warning = "calibration age span covers at least half of the current root-age scale"
        elif root_age > 0.0 and (
            (calibration.minimum_age is not None and calibration.minimum_age / root_age >= 0.8)
            or (calibration.maximum_age is not None and calibration.maximum_age / root_age >= 0.8)
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
        only_calibration = next((calibration for calibration in report.calibrations if calibration.valid), None)
        if only_calibration is not None and only_calibration.calibration_id not in dominant_calibration_ids:
            dominant_calibration_ids.append(only_calibration.calibration_id)
        warnings.append("only one valid calibration remains, so time estimates are effectively driven by a single calibration target")
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
        alignment_taxa = {record.identifier for record in load_fasta_alignment(alignment_path)}
    table = load_taxon_table(tip_dates_path, taxon_column=taxon_column)
    if date_column not in table.columns:
        raise ValueError(f"tip-dating table does not contain column '{date_column}'")

    issues: list[TipDatingValidationIssue] = []
    tip_dates: list[ValidatedTipDate] = []
    extra_tip_taxa = sorted(set(table.taxa) - tree_taxa)
    missing_tree_taxa = sorted(tree_taxa - set(table.taxa))
    extra_alignment_taxa = sorted(set(table.taxa) - alignment_taxa) if alignment_taxa is not None else []
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_date = row.get(date_column, "")
        valid = True
        parsed_date: float | None = None
        if taxon not in tree_taxa:
            valid = False
            issues.append(TipDatingValidationIssue(taxon=taxon, code="taxon-missing-from-tree", message="dated tip is absent from the tree"))
        if alignment_taxa is not None and taxon not in alignment_taxa:
            valid = False
            issues.append(TipDatingValidationIssue(taxon=taxon, code="taxon-missing-from-alignment", message="dated tip is absent from the alignment"))
        if not raw_date.strip():
            valid = False
            issues.append(TipDatingValidationIssue(taxon=taxon, code="missing-date", message="dated tip requires a numeric date value"))
        else:
            try:
                parsed_date = float(raw_date)
            except ValueError:
                valid = False
                issues.append(TipDatingValidationIssue(taxon=taxon, code="invalid-date", message="dated tip value must be numeric"))
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
        calibration_report = validate_fossil_calibration_table(tree_path, calibration_path)
        if calibration_report.invalid_calibration_count:
            blockers.append("calibration table contains invalid fossil calibration targets or ages")
        calibration_dominance = assess_calibration_dominance(tree_path, calibration_path)
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
            blockers.append("tip-date table contains missing, invalid, or mismatched dated taxa")
    if calibration_path is None and tip_dates_path is None:
        warnings.append("no calibrations or tip dates were supplied, so the tree cannot be dated by the current workflow")
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


def prepare_beast_time_tree_analysis(
    alignment_path: Path,
    output_path: Path,
    *,
    tree_path: Path | None = None,
    calibration_path: Path | None = None,
    tip_dates_path: Path | None = None,
    clock_model: str = "strict",
    tree_prior: str = "yule",
    chain_length: int = 1000000,
    log_every: int = 1000,
    taxon_column: str | None = None,
    date_column: str = "date",
) -> BeastPreparationReport:
    """Prepare a deterministic BEAST-style XML configuration from alignment and dating inputs."""
    _ensure_inference_ready_alignment(alignment_path)
    records = load_fasta_alignment(alignment_path)
    inferred_alphabet = infer_alignment_alphabet(records)
    calibration_report = (
        validate_fossil_calibration_table(tree_path, calibration_path)
        if tree_path is not None and calibration_path is not None
        else None
    )
    if calibration_report is not None and calibration_report.invalid_calibration_count:
        raise ValueError("BEAST preparation requires all fossil calibrations to validate successfully")
    tip_date_report = (
        validate_tip_dating_metadata(
            tree_path or alignment_path,
            tip_dates_path,
            alignment_path=alignment_path,
            taxon_column=taxon_column,
            date_column=date_column,
        )
        if tip_dates_path is not None and tree_path is not None
        else None
    )
    if tip_dates_path is not None and tree_path is None:
        raise ValueError("BEAST preparation requires tree_path when tip_dates_path is provided")
    if tip_date_report is not None and tip_date_report.invalid_tip_count:
        raise ValueError("BEAST preparation requires all tip dates to validate successfully")

    sequence_block = "\n".join(
        f'    <sequence taxon="{record.identifier}" value="{record.sequence}" />'
        for record in records
    )
    tip_date_block = ""
    if tip_date_report is not None:
        tip_date_lines = [
            f'    <date taxon="{tip.taxon}" value="{tip.date}" />'
            for tip in tip_date_report.tip_dates
            if tip.valid and tip.date is not None
        ]
        tip_date_block = "\n".join(["  <tipDates>", *tip_date_lines, "  </tipDates>"])
    calibration_block = ""
    if calibration_report is not None:
        calibration_lines = [
            (
                f'    <calibration id="{calibration.calibration_id}" kind="{calibration.target_kind}" '
                f'target="{calibration.target_label}" minimum="{calibration.minimum_age}" '
                f'maximum="{calibration.maximum_age}" distribution="{calibration.distribution}" />'
            )
            for calibration in calibration_report.calibrations
        ]
        calibration_block = "\n".join(["  <calibrations>", *calibration_lines, "  </calibrations>"])
    xml = "\n".join(
        [
            '<?xml version="1.0" encoding="UTF-8"?>',
            '<beast version="bijux-beast-prep-1">',
            f'  <alignment alphabet="{inferred_alphabet}">',
            sequence_block,
            "  </alignment>",
            f'  <clockModel name="{clock_model}" />',
            f'  <treePrior name="{tree_prior}" />',
            *( [tip_date_block] if tip_date_block else [] ),
            *( [calibration_block] if calibration_block else [] ),
            f'  <run chainLength="{chain_length}" logEvery="{log_every}" />',
            "</beast>",
            "",
        ]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(xml, encoding="utf-8")
    return BeastPreparationReport(
        alignment_path=alignment_path,
        output_path=output_path,
        tree_path=tree_path,
        calibration_path=calibration_path,
        tip_dates_path=tip_dates_path,
        taxon_count=len(records),
        character_count=len(records[0].sequence),
        inferred_alphabet=inferred_alphabet,
        clock_model=clock_model,
        tree_prior=tree_prior,
        chain_length=chain_length,
        log_every=log_every,
        calibration_count=0 if calibration_report is None else calibration_report.calibration_count,
        tip_date_count=0 if tip_date_report is None else tip_date_report.valid_tip_count,
    )


def parse_beast_log(path: Path) -> BeastLogReport:
    """Parse a BEAST-style log table into deterministic numeric rows."""
    with path.open(encoding="utf-8", newline="") as handle:
        filtered_lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"BEAST log contains no header: {path}")
    state_field = "state" if "state" in reader.fieldnames else "State"
    if state_field not in reader.fieldnames:
        raise ValueError(f"BEAST log lacks a state column: {path}")
    columns = [field for field in reader.fieldnames if field and field != state_field]
    rows: list[BeastLogRow] = []
    for row in reader:
        values = {
            column: float(row[column])
            for column in columns
            if row.get(column) not in {None, ""}
        }
        rows.append(BeastLogRow(state=int(float(row[state_field])), values=values))
    if not rows:
        raise ValueError(f"BEAST log contains no sampled rows: {path}")
    return BeastLogReport(path=path, row_count=len(rows), columns=columns, rows=rows)


def validate_beast_posterior_log(
    path: Path,
    *,
    required_columns: tuple[str, ...] = ("posterior", "likelihood"),
) -> BeastPosteriorLogValidationReport:
    """Validate that a BEAST posterior log contains the required fields and monotonic sampled states."""
    issues: list[BeastLogValidationIssue] = []
    with path.open(encoding="utf-8", newline="") as handle:
        filtered_lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(filtered_lines, delimiter="\t")
    if reader.fieldnames is None:
        issues.append(BeastLogValidationIssue(code="missing-header", message="BEAST log contains no header row"))
        return BeastPosteriorLogValidationReport(
            path=path,
            row_count=0,
            state_count=0,
            required_columns=list(required_columns),
            observed_columns=[],
            missing_columns=list(required_columns),
            issues=issues,
            valid=False,
        )
    state_field = "state" if "state" in reader.fieldnames else "State" if "State" in reader.fieldnames else None
    observed_columns = [field for field in reader.fieldnames if field]
    if state_field is None:
        issues.append(BeastLogValidationIssue(code="missing-state-column", message="BEAST log lacks a state column"))
    missing_columns = [column for column in required_columns if column not in observed_columns]
    for column in missing_columns:
        issues.append(BeastLogValidationIssue(code="missing-required-column", message=f"missing required BEAST log column '{column}'", column=column))
    previous_state: int | None = None
    row_count = 0
    state_count = 0
    for row_index, row in enumerate(reader, start=2):
        row_count += 1
        if state_field is not None:
            raw_state = row.get(state_field, "")
            if not raw_state:
                issues.append(BeastLogValidationIssue(code="missing-state-value", message="row is missing a sampled state", row=row_index, column=state_field))
            else:
                try:
                    state_value = int(float(raw_state))
                    state_count += 1
                    if previous_state is not None and state_value <= previous_state:
                        issues.append(BeastLogValidationIssue(code="nonmonotonic-state", message="sampled states must increase strictly through the log", row=row_index, column=state_field))
                    previous_state = state_value
                except ValueError:
                    issues.append(BeastLogValidationIssue(code="invalid-state-value", message="sampled state must be numeric", row=row_index, column=state_field))
        for column in observed_columns:
            if column == state_field:
                continue
            raw_value = row.get(column, "")
            if raw_value in {None, ""}:
                issues.append(BeastLogValidationIssue(code="missing-parameter-value", message=f"missing sampled value for '{column}'", row=row_index, column=column))
                continue
            try:
                float(raw_value)
            except ValueError:
                issues.append(BeastLogValidationIssue(code="invalid-parameter-value", message=f"sampled value for '{column}' must be numeric", row=row_index, column=column))
    if row_count == 0:
        issues.append(BeastLogValidationIssue(code="missing-rows", message="BEAST log contains no sampled rows"))
    return BeastPosteriorLogValidationReport(
        path=path,
        row_count=row_count,
        state_count=state_count,
        required_columns=list(required_columns),
        observed_columns=observed_columns,
        missing_columns=missing_columns,
        issues=issues,
        valid=not issues,
    )


def assess_beast_burnin_sensitivity(
    posterior_tree_path: Path,
    *,
    log_path: Path | None = None,
    burnin_fractions: tuple[float, ...] = (0.1, 0.25, 0.5),
) -> BeastBurninSensitivityReport:
    """Compare posterior summaries across multiple BEAST burn-in fractions."""
    if not burnin_fractions:
        raise ValueError("burnin_fractions must contain at least one value")
    log_report = parse_beast_log(log_path) if log_path is not None else None
    ordered_fractions = tuple(sorted(dict.fromkeys(burnin_fractions)))
    slices: list[BeastBurninSensitivitySlice] = []
    previous_newick: str | None = None
    changed_mcc_count = 0
    for fraction in ordered_fractions:
        _, mcc_report = summarize_maximum_clade_credibility_tree(
            posterior_tree_path,
            burnin_fraction=fraction,
        )
        posterior_mean = None
        likelihood_mean = None
        tree_height_mean = None
        if log_report is not None:
            burnin_row_count = int(log_report.row_count * fraction)
            kept_rows = log_report.rows[burnin_row_count:]
            if kept_rows:
                posterior_values = [row.values["posterior"] for row in kept_rows if "posterior" in row.values]
                likelihood_values = [row.values["likelihood"] for row in kept_rows if "likelihood" in row.values]
                tree_height_values = [row.values["treeHeight"] for row in kept_rows if "treeHeight" in row.values]
                posterior_mean = _mean_or_none(posterior_values)
                likelihood_mean = _mean_or_none(likelihood_values)
                tree_height_mean = _mean_or_none(tree_height_values)
        slices.append(
            BeastBurninSensitivitySlice(
                burnin_fraction=fraction,
                burnin_tree_count=mcc_report.burnin_tree_count,
                kept_tree_count=mcc_report.kept_tree_count,
                rooted_topology_count=mcc_report.rooted_topology_count,
                selected_tree_index=mcc_report.selected_tree_index,
                clade_credibility_score=mcc_report.clade_credibility_score,
                posterior_mean=posterior_mean,
                likelihood_mean=likelihood_mean,
                tree_height_mean=tree_height_mean,
            )
        )
        if previous_newick is not None and previous_newick != mcc_report.mcc_newick:
            changed_mcc_count += 1
        previous_newick = mcc_report.mcc_newick
    warnings: list[str] = []
    if changed_mcc_count:
        warnings.append("maximum clade credibility topology changes across tested burn-in fractions")
    if len({row.rooted_topology_count for row in slices}) > 1:
        warnings.append("rooted topology diversity changes across tested burn-in fractions")
    return BeastBurninSensitivityReport(
        posterior_tree_path=posterior_tree_path,
        log_path=log_path,
        slices=slices,
        changed_mcc_count=changed_mcc_count,
        warnings=warnings,
    )


def assess_beast_chain_mixing(
    log_paths: list[Path],
    *,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
    cross_chain_mean_shift_threshold: float = 0.75,
    stuck_parameter_span_threshold: float = 1e-9,
) -> BeastChainMixingReport:
    """Flag low ESS, mean drift, stuck parameters, and inconsistent chain means across BEAST logs."""
    if not log_paths:
        raise ValueError("assess_beast_chain_mixing requires at least one log path")
    chain_summaries = [
        assess_beast_convergence(
            path,
            ess_threshold=ess_threshold,
            mean_shift_threshold=mean_shift_threshold,
        )
        for path in log_paths
    ]
    issues: list[BeastChainMixingIssue] = []
    for summary in chain_summaries:
        for warning in summary.warnings:
            issues.append(
                BeastChainMixingIssue(
                    path=summary.path,
                    parameter=str(warning["parameter"]),
                    code=str(warning["code"]),
                    message=str(warning["message"]),
                    observed_value=float(warning["observed_value"]),
                    threshold=float(warning["threshold"]),
                )
            )
        for parameter_summary in summary.parameter_summaries:
            span = float(parameter_summary["maximum"]) - float(parameter_summary["minimum"])
            if span <= stuck_parameter_span_threshold:
                issues.append(
                    BeastChainMixingIssue(
                        path=summary.path,
                        parameter=str(parameter_summary["parameter"]),
                        code="stuck-parameter",
                        message="parameter shows effectively no movement across the sampled chain",
                        observed_value=span,
                        threshold=stuck_parameter_span_threshold,
                    )
                )
    parameter_to_means: dict[str, list[tuple[Path, float]]] = {}
    for summary in chain_summaries:
        for parameter_summary in summary.parameter_summaries:
            parameter_to_means.setdefault(str(parameter_summary["parameter"]), []).append(
                (summary.path, float(parameter_summary["mean"]))
            )
    for parameter, chain_means in parameter_to_means.items():
        if len(chain_means) < 2:
            continue
        mean_values = [value for _, value in chain_means]
        span = max(mean_values) - min(mean_values)
        if span > cross_chain_mean_shift_threshold:
            issues.append(
                BeastChainMixingIssue(
                    path=None,
                    parameter=parameter,
                    code="inconsistent-chains",
                    message="independent chains disagree more than the allowed mean-shift threshold",
                    observed_value=span,
                    threshold=cross_chain_mean_shift_threshold,
                )
            )
    return BeastChainMixingReport(
        log_paths=log_paths,
        chain_count=len(log_paths),
        converged=not issues,
        issues=issues,
        chain_summaries=chain_summaries,
    )


def assess_beast_convergence(
    path: Path,
    *,
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> BeastConvergenceReport:
    """Flag low-ESS or unstable BEAST trace parameters."""
    report = parse_beast_log(path)
    convergence = summarize_trace_convergence(
        path=path,
        rows=[row.values for row in report.rows],
        columns=report.columns,
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
    )
    return _build_beast_convergence_report(convergence)


def _build_beast_convergence_report(convergence: TraceConvergenceReport) -> BeastConvergenceReport:
    return BeastConvergenceReport(
        path=convergence.path,
        sample_count=convergence.sample_count,
        converged=convergence.converged,
        ess_threshold=convergence.ess_threshold,
        mean_shift_threshold=convergence.mean_shift_threshold,
        warnings=[
            {
                "parameter": warning.parameter,
                "code": warning.code,
                "message": warning.message,
                "observed_value": warning.observed_value,
                "threshold": warning.threshold,
            }
            for warning in convergence.warnings
        ],
        parameter_summaries=[
            {
                "parameter": summary.parameter,
                "sample_count": summary.sample_count,
                "effective_sample_size": summary.effective_sample_size,
                "mean": summary.mean,
                "minimum": summary.minimum,
                "maximum": summary.maximum,
                "first_half_mean": summary.first_half_mean,
                "second_half_mean": summary.second_half_mean,
                "standardized_mean_shift": summary.standardized_mean_shift,
            }
            for summary in convergence.series
        ],
    )


def _mean_or_none(values: list[float]) -> float | None:
    if not values:
        return None
    return round(sum(values) / len(values), 6)
