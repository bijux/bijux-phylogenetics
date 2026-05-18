from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from importlib import metadata
import json
import math
import os
from pathlib import Path
import shutil
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousModeSearchControls,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness

from .registry import GeigerParityCase, list_geiger_parity_cases


@dataclass(frozen=True, slots=True)
class GeigerParityObservation:
    """One live parity comparison between Bijux and `geiger`."""

    case_id: str
    fixture_id: str
    function_name: str
    python_function_name: str
    input_fixtures: tuple[Path, ...]
    model_name: str
    optimizer_settings: dict[str, object] | None
    tolerance: float
    r_version: str | None
    geiger_version: str | None
    bijux_version: str
    bijux_commit: str | None
    status: str
    passed: bool
    mismatch_reason: str | None
    reproducible_artifact_root: Path | None
    reference_summary: dict[str, object] | None
    bijux_summary: dict[str, object] | None
    reference_rows: list[dict[str, object]] | None
    bijux_rows: list[dict[str, object]] | None
    reference_error: dict[str, object] | None
    bijux_error: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class GeigerParitySummaryRow:
    """One function-level summary across governed `geiger` parity cases."""

    function_name: str
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int


@dataclass(slots=True)
class GeigerParityReport:
    """Aggregate report for governed live `geiger` parity cases."""

    observations: list[GeigerParityObservation]
    summary_rows: list[GeigerParitySummaryRow]
    case_count: int
    passed_case_count: int
    failed_case_count: int
    skipped_case_count: int
    all_passed: bool
    limitations: list[str]


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[6]


def _geiger_runner_path() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "resources"
        / "reference"
        / "geiger_parity_runner.R"
    )


def _failure_root() -> Path:
    return _repository_root() / "artifacts" / "geiger-parity-failures"


def _reference_environment() -> dict[str, str]:
    environment = dict(os.environ)
    r_library = _repository_root() / "artifacts" / "r-lib"
    if "R_LIBS_USER" not in environment and r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    return environment


def _bijux_version() -> str:
    try:
        return metadata.version("bijux-phylogenetics")
    except metadata.PackageNotFoundError:
        return "0.1.0"


def _bijux_commit() -> str | None:
    result = subprocess.run(  # nosec
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        check=False,
        cwd=_repository_root(),
        text=True,
    )
    if result.returncode != 0:
        return None
    commit = result.stdout.strip()
    return commit or None


def _selected_cases(case_ids: list[str] | None) -> list[GeigerParityCase]:
    registry = {case.case_id: case for case in list_geiger_parity_cases()}
    if case_ids is None:
        return list(registry.values())
    missing = [case_id for case_id in case_ids if case_id not in registry]
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"unknown geiger parity case id(s): {missing_text}")
    return [registry[case_id] for case_id in case_ids]


def _write_case_file(path: Path, case: GeigerParityCase) -> Path:
    payload = {
        "case_id": case.case_id,
        "fixture_id": case.fixture_id,
        "function_name": case.function_name,
        "python_function_name": case.python_function_name,
        "operation": case.operation,
        "model_name": case.model_name,
        "python_mode": case.python_mode,
        "input_fixtures": [str(item) for item in case.input_fixtures],
        "tolerance": case.tolerance,
        "trait_name": case.trait_name,
        "taxon_column": case.taxon_column,
        "optimizer_settings": case.optimizer_settings,
        "candidate_model_names": None
        if case.candidate_model_names is None
        else list(case.candidate_model_names),
        "reference_control": case.reference_control,
        "coarse_grid_point_count": case.coarse_grid_point_count,
        "fine_grid_point_count": case.fine_grid_point_count,
        "initial_parameter_value": case.initial_parameter_value,
        "comparison_fields": list(case.comparison_fields),
        "lambda_bounds": None
        if case.lambda_bounds is None
        else list(case.lambda_bounds),
        "kappa_bounds": None
        if case.kappa_bounds is None
        else list(case.kappa_bounds),
        "delta_bounds": None
        if case.delta_bounds is None
        else list(case.delta_bounds),
        "ou_bounds": None if case.ou_bounds is None else list(case.ou_bounds),
        "early_burst_bounds": None
        if case.early_burst_bounds is None
        else list(case.early_burst_bounds),
        "field_tolerances": case.field_tolerances,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_rows_table(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, object]] = []
        for row in reader:
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if value is None or value == "":
                    parsed[key] = ""
                    continue
                if key in {"parameter", "model", "comparability_note"}:
                    parsed[key] = value
                    continue
                if value.lower() in {"true", "false"}:
                    parsed[key] = value.lower() == "true"
                    continue
                try:
                    parsed[key] = int(value)
                    continue
                except ValueError:
                    pass
                try:
                    parsed[key] = float(value)
                    continue
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
        return rows


def _optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None


def _field_tolerance(case: GeigerParityCase, key: str) -> float:
    if case.field_tolerances and key in case.field_tolerances:
        return case.field_tolerances[key]
    return case.tolerance


def _isclose(left: object, right: object, *, tolerance: float) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return left == right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(
            float(left),
            float(right),
            rel_tol=tolerance,
            abs_tol=tolerance,
        )
    return left == right


def _mismatch_reason(
    case: GeigerParityCase,
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> str | None:
    if reference_summary is None or bijux_summary is None:
        return "summary_missing"
    for key in case.comparison_fields:
        if key not in reference_summary:
            return f"reference_summary_missing:{key}"
        if key not in bijux_summary:
            return f"bijux_summary_missing:{key}"
        if not _isclose(
            reference_summary[key],
            bijux_summary[key],
            tolerance=_field_tolerance(case, key),
        ):
            return f"summary_field_mismatch:{key}"
    return None


def _row_mismatch_reason(
    case: GeigerParityCase,
    *,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
) -> str | None:
    if reference_rows is None or bijux_rows is None:
        return "rows_missing"
    normalized_reference_rows = _normalized_rows(case, reference_rows)
    normalized_bijux_rows = _normalized_rows(case, bijux_rows)
    if len(normalized_reference_rows) != len(normalized_bijux_rows):
        return "row_count_mismatch"
    for reference_row, bijux_row in zip(
        normalized_reference_rows,
        normalized_bijux_rows,
        strict=True,
    ):
        reference_id = reference_row.get("parameter", reference_row.get("model"))
        bijux_id = bijux_row.get("parameter", bijux_row.get("model"))
        if reference_id != bijux_id:
            return "row_identifier_mismatch"
        if set(reference_row) != set(bijux_row):
            return "row_field_set_mismatch"
        for key in reference_row:
            if not _isclose(
                reference_row.get(key),
                bijux_row.get(key),
                tolerance=case.tolerance,
            ):
                return f"row_field_mismatch:{reference_id}:{key}"
    return None


def _normalized_rows(
    case: GeigerParityCase,
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    normalized = [dict(row) for row in rows]
    if case.operation != "compare-fitcontinuous-models":
        return normalized
    for row in normalized:
        row.pop("rank", None)
    normalized.sort(key=lambda row: str(row.get("model", "")))
    return normalized


def _parameter_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for parameter in ("root_state", "rate", "log_likelihood", "aic", "aicc"):
        value = summary.get(parameter)
        if value in {None, ""}:
            continue
        rows.append({"parameter": parameter, "value": value})
    parameter_name = summary.get("parameter_name")
    parameter_value = summary.get("parameter_value")
    if isinstance(parameter_name, str) and parameter_name and parameter_value not in {
        None,
        "",
    }:
        rows.append({"parameter": parameter_name, "value": parameter_value})
    return rows


def _comparison_rows(report) -> list[dict[str, object]]:
    return [
        {
            "model": row.model,
            "rank": "" if row.rank is None else row.rank,
            "parameter_count": row.parameter_count,
            "log_likelihood": row.log_likelihood,
            "aic": row.aic,
            "aicc": row.aicc,
            "delta_aic": row.delta_aic,
            "delta_aicc": row.delta_aicc,
            "selected": row.selected,
            "comparable": row.comparable,
            "comparability_note": "" if row.comparability_note is None else row.comparability_note,
        }
        for row in report.rows
    ]


def _comparison_modes(
    candidate_model_names: tuple[str, ...] | None,
) -> tuple[str, ...] | None:
    if candidate_model_names is None:
        return None
    mode_names = {
        "BM": "brownian",
        "white": "white-noise",
        "lambda": "pagel-lambda",
        "kappa": "pagel-kappa",
        "delta": "pagel-delta",
        "OU": "ornstein-uhlenbeck",
        "EB": "early-burst",
    }
    return tuple(mode_names.get(model_name, model_name) for model_name in candidate_model_names)


def _standard_error_policy() -> str:
    return "fitcontinuous-standard-error-explicitly-excluded-this-round"


def _missing_value_policy() -> str:
    return "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values"


def _parameter_bound_policy(case: GeigerParityCase) -> str:
    if case.python_mode in {"brownian", "white-noise"}:
        return "closed-form-without-parameter-bounds"
    return "governed-bounded-grid-search"


def _bijux_optimizer_result(
    case: GeigerParityCase,
    report,
) -> dict[str, object]:
    if report.optimizer_diagnostics is not None:
        diagnostics = report.optimizer_diagnostics
        return {
            "optimizer_name": diagnostics.optimizer_name,
            "parameter_search_strategy": diagnostics.parameter_search_strategy,
            "converged": diagnostics.converged,
            "lower_bound": diagnostics.lower_bound,
            "upper_bound": diagnostics.upper_bound,
            "starting_parameter_policy": diagnostics.starting_parameter_policy,
            "starting_parameter_value": diagnostics.starting_parameter_value,
            "starting_parameter_log_likelihood": (
                diagnostics.starting_parameter_log_likelihood
            ),
            "coarse_grid_point_count": diagnostics.coarse_grid_point_count,
            "fine_grid_point_count": diagnostics.fine_grid_point_count,
            "function_evaluation_count": diagnostics.function_evaluation_count,
            "coarse_best_parameter": diagnostics.coarse_best_parameter,
            "coarse_best_log_likelihood": diagnostics.coarse_best_log_likelihood,
            "fine_search_start": diagnostics.fine_search_start,
            "fine_search_stop": diagnostics.fine_search_stop,
            "hit_lower_boundary": diagnostics.hit_lower_boundary,
            "hit_upper_boundary": diagnostics.hit_upper_boundary,
        }
    if case.python_mode in {"brownian", "white-noise"}:
        return {
            "optimizer_name": "closed-form-profile-solution",
            "parameter_search": "none",
            "converged": True,
            "parameter_count": 2,
        }
    if case.python_mode == "ornstein-uhlenbeck":
        return {
            "optimizer_name": "governed-two-stage-grid-search",
            "parameter_search": "bounded-grid-search",
            "converged": True,
            "parameter_count": 3,
            "coarse_grid_point_count": 81,
            "fine_grid_point_count": 81,
        }
    return {
        "optimizer_name": "governed-two-stage-grid-search",
        "parameter_search": "bounded-grid-search",
        "converged": True,
        "parameter_count": 3,
        "coarse_grid_point_count": 81,
        "fine_grid_point_count": 81,
    }


def _build_bijux_case_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    if case.operation == "compare-fitcontinuous-models":
        return _build_bijux_model_comparison_payload(case)
    tree_path, traits_path = case.input_fixtures
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=case.trait_name,
        taxon_column=case.taxon_column,
    )
    search_controls = None
    if (
        case.coarse_grid_point_count is not None
        or case.fine_grid_point_count is not None
        or case.initial_parameter_value is not None
    ):
        search_controls = ContinuousModeSearchControls(
            coarse_grid_point_count=(
                81
                if case.coarse_grid_point_count is None
                else case.coarse_grid_point_count
            ),
            fine_grid_point_count=(
                81
                if case.fine_grid_point_count is None
                else case.fine_grid_point_count
            ),
            initial_parameter_value=case.initial_parameter_value,
        )
    report = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait=case.trait_name,
        mode=case.python_mode,
        taxon_column=case.taxon_column,
        search_controls=search_controls,
        lambda_bounds=(0.0, 1.0)
        if case.lambda_bounds is None
        else case.lambda_bounds,
        kappa_bounds=(0.0, 3.0)
        if case.kappa_bounds is None
        else case.kappa_bounds,
        delta_bounds=(0.0, 3.0)
        if case.delta_bounds is None
        else case.delta_bounds,
        ou_bounds=(0.0, 10.0) if case.ou_bounds is None else case.ou_bounds,
        early_burst_bounds=(0.0, 50.0)
        if case.early_burst_bounds is None
        else case.early_burst_bounds,
    )
    excluded_taxa = sorted(
        {
            *readiness.missing_from_traits,
            *readiness.pruned_missing_value_taxa,
            *readiness.pruned_non_numeric_taxa,
        }
    )
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case.model_name,
        "excluded_taxon_count": len(excluded_taxa),
        "excluded_taxa": excluded_taxa,
        "missing_value_taxa": list(readiness.pruned_missing_value_taxa),
        "non_numeric_taxa": list(readiness.pruned_non_numeric_taxa),
        "missing_from_traits": list(readiness.missing_from_traits),
        "extra_trait_taxa": list(readiness.extra_trait_taxa),
        "missing_value_policy": _missing_value_policy(),
        "standard_error_policy": _standard_error_policy(),
        "parameter_bound_policy": _parameter_bound_policy(case),
        "hit_lower_parameter_boundary": (
            False
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.hit_lower_boundary
        ),
        "hit_upper_parameter_boundary": (
            False
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.hit_upper_boundary
        ),
        "identifiability_warning_kinds": [
            warning.kind for warning in report.identifiability_warnings
        ],
        "identifiability_warning_count": len(report.identifiability_warnings),
        "root_state": report.root_state,
        "rate": report.rate,
        "log_likelihood": report.log_likelihood,
        "aic": report.aic,
        "aicc": report.aicc,
        "parameter_name": report.parameter_name,
        "parameter_value": report.parameter_value,
        "optimizer_settings": case.optimizer_settings,
        "optimizer_result": _bijux_optimizer_result(case, report),
    }
    return summary, _parameter_rows(summary)


def _build_bijux_model_comparison_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_path, traits_path = case.input_fixtures
    report = compare_fitcontinuous_model_ranking(
        tree_path,
        traits_path,
        trait=case.trait_name,
        taxon_column=case.taxon_column,
        modes=_comparison_modes(case.candidate_model_names),
        lambda_bounds=(0.0, 1.0)
        if case.lambda_bounds is None
        else case.lambda_bounds,
        kappa_bounds=(0.0, 3.0)
        if case.kappa_bounds is None
        else case.kappa_bounds,
        delta_bounds=(0.0, 3.0)
        if case.delta_bounds is None
        else case.delta_bounds,
        ou_bounds=(0.0, 10.0) if case.ou_bounds is None else case.ou_bounds,
        early_burst_bounds=(0.0, 50.0)
        if case.early_burst_bounds is None
        else case.early_burst_bounds,
    )
    runner_up_rows = [
        row
        for row in report.rows
        if row.model != report.better_model and row.comparable and row.rank is not None
    ]
    runner_up_row = runner_up_rows[0] if runner_up_rows else None
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case.model_name,
        "selected_model": report.better_model,
        "model_ranking": [row.model for row in report.rows],
        "comparable_model_count": sum(1 for row in report.rows if row.comparable),
        "noncomparable_model_count": sum(1 for row in report.rows if not row.comparable),
        "runner_up_model": None if runner_up_row is None else runner_up_row.model,
        "runner_up_aicc_delta": (
            math.nan if runner_up_row is None else runner_up_row.delta_aicc
        ),
        "warning_count": len(report.warnings),
        "optimizer_settings": case.optimizer_settings,
    }
    return summary, _comparison_rows(report)


def _write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys()) if rows else ["parameter", "value"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: GeigerParityCase,
    case_file: Path,
    execution_root: Path,
    execution_payload: dict[str, object] | None,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
    reference_error: dict[str, object] | None,
    bijux_error: dict[str, object] | None,
    mismatch_reason: str,
) -> Path:
    artifact_root = failure_root / case.case_id
    if artifact_root.exists():
        shutil.rmtree(artifact_root)
    artifact_root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case_file, artifact_root / "case.json")
    if execution_root.exists():
        shutil.copytree(
            execution_root,
            artifact_root / "reference",
            dirs_exist_ok=True,
        )
    if execution_payload is not None:
        (artifact_root / "reference-execution.json").write_text(
            json.dumps(execution_payload, indent=2),
            encoding="utf-8",
        )
    if reference_summary is not None:
        (artifact_root / "reference-summary.json").write_text(
            json.dumps(reference_summary, indent=2),
            encoding="utf-8",
        )
    if bijux_summary is not None:
        (artifact_root / "bijux-summary.json").write_text(
            json.dumps(bijux_summary, indent=2),
            encoding="utf-8",
        )
    if reference_rows is not None:
        _write_rows_table(artifact_root / "reference-parameters.tsv", reference_rows)
    if bijux_rows is not None:
        _write_rows_table(artifact_root / "bijux-parameters.tsv", bijux_rows)
    if reference_error is not None:
        (artifact_root / "reference-error.json").write_text(
            json.dumps(reference_error, indent=2),
            encoding="utf-8",
        )
    if bijux_error is not None:
        (artifact_root / "bijux-error.json").write_text(
            json.dumps(bijux_error, indent=2),
            encoding="utf-8",
        )
    (artifact_root / "mismatch-reason.txt").write_text(
        mismatch_reason,
        encoding="utf-8",
    )
    return artifact_root


def _summary_rows(
    observations: list[GeigerParityObservation],
) -> list[GeigerParitySummaryRow]:
    by_function: dict[str, list[GeigerParityObservation]] = {}
    for observation in observations:
        by_function.setdefault(observation.function_name, []).append(observation)
    return [
        GeigerParitySummaryRow(
            function_name=function_name,
            case_count=len(items),
            passed_case_count=sum(1 for item in items if item.status == "passed"),
            failed_case_count=sum(1 for item in items if item.status == "failed"),
            skipped_case_count=sum(1 for item in items if item.status == "skipped"),
        )
        for function_name, items in sorted(by_function.items())
    ]


def run_geiger_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
) -> GeigerParityReport:
    """Run governed live `geiger` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids)
    observations: list[GeigerParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-geiger-parity-{case.case_id}-"
        ) as tmpdir:
            working_root = Path(tmpdir)
            case_file = _write_case_file(working_root / "case.json", case)
            execution_root = working_root / "reference"
            execution_root.mkdir(parents=True, exist_ok=True)
            bijux_summary: dict[str, object] | None = None
            bijux_rows: list[dict[str, object]] | None = None
            bijux_error: dict[str, object] | None = None
            try:
                bijux_summary, bijux_rows = _build_bijux_case_payload(case)
            except Exception as error:
                bijux_error = {
                    "error_type": type(error).__name__,
                    "message": str(error),
                }
            execution_payload: dict[str, object] | None = None
            reference_summary: dict[str, object] | None = None
            reference_rows: list[dict[str, object]] | None = None
            reference_error: dict[str, object] | None = None
            status = "failed"
            mismatch_reason: str | None = None
            artifact_root: Path | None = None
            r_version: str | None = None
            geiger_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                process = subprocess.run(  # nosec
                    [
                        rscript_executable,
                        str(_geiger_runner_path()),
                        str(case_file),
                        str(execution_root),
                    ],
                    capture_output=True,
                    check=False,
                    cwd=_repository_root(),
                    env=_reference_environment(),
                    text=True,
                )
                process_stdout = process.stdout
                process_stderr = process.stderr
            except FileNotFoundError:
                process = None
                status = "skipped"
                mismatch_reason = "rscript_unavailable"
            if process is not None and process.returncode == 0:
                execution_path = execution_root / "reference-execution.json"
                summary_path = execution_root / "reference-summary.json"
                rows_path = execution_root / "reference-parameters.tsv"
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    geiger_version = _optional_payload_string(
                        execution_payload,
                        "geiger_version",
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "geiger_package_unavailable",
                            )
                        )
                    elif execution_status != "ok":
                        reference_error = {
                            "error_type": str(
                                execution_payload.get(
                                    "error_type",
                                    execution_payload.get(
                                        "mismatch_reason",
                                        "reference_execution_failed",
                                    ),
                                )
                            ),
                            "message": str(execution_payload.get("message", "")),
                        }
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "reference_execution_failed",
                            )
                        )
                    elif not summary_path.exists():
                        mismatch_reason = "reference_summary_missing"
                    else:
                        reference_summary = _load_json(summary_path)
                        if rows_path.exists():
                            reference_rows = _load_rows_table(rows_path)
                        mismatch_reason = _mismatch_reason(
                            case,
                            reference_summary=reference_summary,
                            bijux_summary=bijux_summary,
                        )
                        if mismatch_reason is None:
                            if reference_rows is None:
                                mismatch_reason = "reference_rows_missing"
                            else:
                                mismatch_reason = _row_mismatch_reason(
                                    case,
                                    reference_rows=reference_rows,
                                    bijux_rows=bijux_rows,
                                )
                        if mismatch_reason is None:
                            status = "passed"
            elif process is not None and process.returncode != 0:
                mismatch_reason = "reference_execution_failed"
            if status != "passed":
                artifact_root = _persist_failure_bundle(
                    failure_root=active_failure_root,
                    case=case,
                    case_file=case_file,
                    execution_root=execution_root,
                    execution_payload=execution_payload,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                    mismatch_reason=mismatch_reason or "reference_execution_failed",
                )
                if process_stdout:
                    (artifact_root / "reference-stdout.txt").write_text(
                        process_stdout,
                        encoding="utf-8",
                    )
                if process_stderr:
                    (artifact_root / "reference-stderr.txt").write_text(
                        process_stderr,
                        encoding="utf-8",
                    )
            observations.append(
                GeigerParityObservation(
                    case_id=case.case_id,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixtures=case.input_fixtures,
                    model_name=case.model_name,
                    optimizer_settings=case.optimizer_settings,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    geiger_version=geiger_version,
                    bijux_version=bijux_version,
                    bijux_commit=bijux_commit,
                    status=status,
                    passed=status == "passed",
                    mismatch_reason=mismatch_reason,
                    reproducible_artifact_root=artifact_root,
                    reference_summary=reference_summary,
                    bijux_summary=bijux_summary,
                    reference_rows=reference_rows,
                    bijux_rows=bijux_rows,
                    reference_error=reference_error,
                    bijux_error=bijux_error,
                )
            )
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return GeigerParityReport(
        observations=observations,
        summary_rows=_summary_rows(observations),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        skipped_case_count=skipped_case_count,
        all_passed=case_count > 0
        and passed_case_count == case_count
        and failed_case_count == 0
        and skipped_case_count == 0,
        limitations=[
            "The governed live `geiger` parity registry stays narrow in this harness tranche and only covers the currently owned continuous-mode fit surfaces.",
            "This harness requires Rscript plus the `geiger`, `ape`, and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_geiger_parity_summary_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write one row per governed `geiger` function summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "function_name",
                "case_count",
                "passed_case_count",
                "failed_case_count",
                "skipped_case_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.summary_rows:
            writer.writerow(asdict(row))
    return path


def write_geiger_parity_observation_table(
    path: Path,
    report: GeigerParityReport,
) -> Path:
    """Write one row per governed `geiger` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "fixture_id",
                "function_name",
                "python_function_name",
                "input_fixtures",
                "model_name",
                "optimizer_settings",
                "tolerance",
                "r_version",
                "geiger_version",
                "bijux_version",
                "bijux_commit",
                "status",
                "passed",
                "mismatch_reason",
                "reproducible_artifact_root",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for observation in report.observations:
            writer.writerow(
                {
                    "case_id": observation.case_id,
                    "fixture_id": observation.fixture_id,
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixtures": json.dumps(
                        [str(path) for path in observation.input_fixtures]
                    ),
                    "model_name": observation.model_name,
                    "optimizer_settings": json.dumps(
                        observation.optimizer_settings,
                        sort_keys=True,
                    ),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "geiger_version": observation.geiger_version or "",
                    "bijux_version": observation.bijux_version,
                    "bijux_commit": observation.bijux_commit or "",
                    "status": observation.status,
                    "passed": str(observation.passed).lower(),
                    "mismatch_reason": observation.mismatch_reason or "",
                    "reproducible_artifact_root": ""
                    if observation.reproducible_artifact_root is None
                    else str(observation.reproducible_artifact_root),
                }
            )
    return path
