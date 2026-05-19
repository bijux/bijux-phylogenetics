from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from bijux_phylogenetics.comparative.evolutionary_modes import (
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    ContinuousModeSearchControls,
    FITCONTINUOUS_MODEL_RANKING_POLICY,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness

from ..optimizer_triage import build_geiger_optimizer_triage_rows
from ..boundary_warning_registry import build_geiger_boundary_warning_rows
from ..likelihood_policy import build_geiger_likelihood_policy_rows
from ..model_confidence import build_geiger_model_confidence_rows
from ..parameterization_registry import build_geiger_parameterization_registry_rows
from ..registry import (
    GeigerParityCase,
    list_geiger_parity_cases as list_geiger_parity_cases,
)
from .models import (
    GeigerParityObservation,
    GeigerParityReport,
    GeigerParitySummaryRow,
)
from .comparison import (
    comparison_rows as _comparison_rows,
    load_json as _load_json,
    load_rows_table as _load_rows_table,
    mismatch_reason as _mismatch_reason,
    optional_payload_string as _optional_payload_string,
    parameter_rows as _parameter_rows,
    row_mismatch_reason as _row_mismatch_reason,
)
from .runtime import (
    bijux_commit as _bijux_commit,
    bijux_version as _bijux_version,
    failure_root as _failure_root,
    geiger_runner_path as _geiger_runner_path,
    reference_environment as _reference_environment,
    repository_root as _repository_root,
    selected_cases as _selected_cases,
    write_case_file as _write_case_file,
)


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


def _discrete_missing_value_policy() -> str:
    return "prune-overlapping-missing-values"


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
    if case.operation == "fit-discrete-mk":
        diagnostics = report.optimizer_diagnostics
        return {
            "optimizer_name": diagnostics.optimizer_name,
            "parameter_count": diagnostics.parameter_count,
            "initial_candidate_count": diagnostics.initial_candidate_count,
            "best_initial_scale": diagnostics.best_initial_scale,
            "converged": diagnostics.converged,
            "iteration_count": diagnostics.iteration_count,
            "function_evaluation_count": diagnostics.function_evaluation_count,
            "simplex_shrink_count": diagnostics.simplex_shrink_count,
            "hit_lower_parameter_bound": diagnostics.hit_lower_parameter_bound,
            "hit_upper_parameter_bound": diagnostics.hit_upper_parameter_bound,
        }
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
            "profile_rows": (
                None
                if report.optimizer_profile_rows is None
                else [
                    {
                        "parameter_name": report.parameter_name,
                        "parameter_value": row.parameter_value,
                        "log_likelihood": row.log_likelihood,
                    }
                    for row in report.optimizer_profile_rows
                ]
            ),
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
    if case.operation == "fit-discrete-mk":
        return _build_bijux_discrete_case_payload(case)
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
        kappa_bounds=(0.0, 1.0)
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
        "likelihood_constant_policy": report.likelihood_constant_policy,
        "likelihood_comparison_policy": report.likelihood_comparison_policy,
        "parameter_name": report.parameter_name,
        "parameter_value": report.parameter_value,
        "optimizer_settings": case.optimizer_settings,
        "optimizer_result": _bijux_optimizer_result(case, report),
    }
    return summary, _parameter_rows(summary)


def _discrete_rate_rows(report) -> list[dict[str, object]]:
    return [
        {
            "source_state": row.source_state,
            "target_state": row.target_state,
            "transition_allowed": row.transition_allowed,
            "step_distance": row.step_distance,
            "rate": row.rate,
        }
        for row in report.transition_rate_rows
    ]


def _build_bijux_discrete_case_payload(
    case: GeigerParityCase,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    tree_path, traits_path = case.input_fixtures
    report = fit_discrete_mk_model(
        tree_path,
        traits_path,
        trait=case.trait_name,
        taxon_column=case.taxon_column,
        model=case.python_mode,
        transform=case.discrete_transform_name,
        lambda_bounds=(0.0, 1.0)
        if case.lambda_bounds is None
        else case.lambda_bounds,
        kappa_bounds=(0.0, 1.0)
        if case.kappa_bounds is None
        else case.kappa_bounds,
        delta_bounds=(math.exp(-5.0), 3.0)
        if case.delta_bounds is None
        else case.delta_bounds,
        early_burst_bounds=(-10.0, 10.0)
        if case.early_burst_bounds is None
        else case.early_burst_bounds,
    )
    input_audit = report.input_audit
    missing_value_taxa = list(input_audit.pruned_missing_value_taxa)
    excluded_taxa = sorted(
        set(input_audit.missing_from_traits) | set(missing_value_taxa)
    )
    diagnostics = report.optimizer_diagnostics
    transform_fit = report.transform_fit
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case.model_name,
        "transform_name": (
            None
            if transform_fit is None
            else {
                "lambda": "pagel-lambda",
                "kappa": "pagel-kappa",
                "delta": "pagel-delta",
                "early-burst": "early-burst",
            }.get(transform_fit.transform_name, transform_fit.transform_name)
        ),
        "observed_state_count": len(input_audit.observed_states),
        "state_order": list(report.state_order),
        "excluded_taxon_count": len(excluded_taxa),
        "excluded_taxa": excluded_taxa,
        "missing_value_taxa": missing_value_taxa,
        "missing_from_traits": list(input_audit.missing_from_traits),
        "extra_trait_taxa": list(input_audit.extra_trait_taxa),
        "missing_value_policy": _discrete_missing_value_policy(),
        "log_likelihood": report.log_likelihood,
        "parameter_count": report.parameter_count,
        "aic": report.aic,
        "aicc": report.aicc,
        "parameter_name": None if transform_fit is None else transform_fit.parameter_name,
        "parameter_value": (
            None if transform_fit is None else transform_fit.parameter_value
        ),
        "hit_lower_parameter_boundary": (
            diagnostics.hit_lower_parameter_bound
            if transform_fit is None
            else transform_fit.hit_lower_parameter_boundary
        ),
        "hit_upper_parameter_boundary": (
            diagnostics.hit_upper_parameter_bound
            if transform_fit is None
            else transform_fit.hit_upper_parameter_boundary
        ),
        "optimizer_settings": case.optimizer_settings,
        "optimizer_result": {
            **_bijux_optimizer_result(case, report),
            "profile_rows": (
                None
                if transform_fit is None
                else [
                    {
                        "parameter_name": transform_fit.parameter_name,
                        "parameter_value": row.transform_parameter_value,
                        "log_likelihood": row.log_likelihood,
                    }
                    for row in transform_fit.profile_rows
                ]
            ),
        },
    }
    return summary, _discrete_rate_rows(report)


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
        "likelihood_constant_policy": (
            report.likelihood_constant_policy
            if report.likelihood_constant_policy is not None
            else CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY
        ),
        "likelihood_comparison_policy": FITCONTINUOUS_MODEL_RANKING_POLICY,
        "noncomparable_likelihood_models": list(report.noncomparable_likelihood_models),
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
        optimizer_triage_rows=build_geiger_optimizer_triage_rows(observations),
        boundary_warning_rows=build_geiger_boundary_warning_rows(observations),
        likelihood_policy_rows=build_geiger_likelihood_policy_rows(observations),
        model_confidence_rows=build_geiger_model_confidence_rows(observations),
        parameterization_registry_rows=build_geiger_parameterization_registry_rows(
            observations
        ),
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
