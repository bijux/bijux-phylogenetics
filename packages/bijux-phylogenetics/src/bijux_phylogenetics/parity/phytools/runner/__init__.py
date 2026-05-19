from __future__ import annotations

import csv
from dataclasses import asdict
import json
import math
from pathlib import Path
import subprocess  # nosec B404 - parity helpers invoke repository-owned reference commands
import tempfile

from ..registry import PhytoolsParityCase, list_phytools_parity_cases
from .comparative_payloads import build_comparative_case_payload
from .continuous_payloads import build_continuous_case_payload
from .discrete_payloads import build_discrete_case_payload
from .models import (
    PhytoolsParityObservation,
    PhytoolsParityReport,
    PhytoolsParitySummaryRow,
)
from .runtime import (
    bijux_commit as _bijux_commit,
    bijux_version as _bijux_version,
    failure_root as _failure_root,
    phytools_runner_path as _phytools_runner_path,
    reference_environment as _reference_environment,
    repository_root as _repository_root,
)
from .signal_payloads import build_signal_case_payload


def _optional_payload_string(payload: dict[str, object], key: str) -> str | None:
    value = payload.get(key)
    return value if isinstance(value, str) else None


def _selected_cases(case_ids: list[str] | None) -> list[PhytoolsParityCase]:
    registry = {case.case_id: case for case in list_phytools_parity_cases()}
    if case_ids is None:
        return list(registry.values())
    missing = [case_id for case_id in case_ids if case_id not in registry]
    if missing:
        missing_text = ", ".join(sorted(missing))
        raise ValueError(f"unknown phytools parity case id(s): {missing_text}")
    return [registry[case_id] for case_id in case_ids]


def _write_case_file(path: Path, case: PhytoolsParityCase) -> Path:
    payload = {
        "case_id": case.case_id,
        "fixture_id": case.fixture_id,
        "function_name": case.function_name,
        "operation": case.operation,
        "input_fixtures": [str(path) for path in case.input_fixtures],
        "trait_name": case.trait_name,
        "taxon_column": case.taxon_column,
        "discrete_model": case.discrete_model,
        "root_prior_mode": case.root_prior_mode,
        "tolerance": case.tolerance,
        "permutation_count": case.permutation_count,
        "permutation_seed": case.permutation_seed,
        "stochastic_map_replicate_count": case.stochastic_map_replicate_count,
        "stochastic_map_seed": case.stochastic_map_seed,
        "density_resolution": case.density_resolution,
        "focal_state": case.focal_state,
        "simulation_states": case.simulation_states,
        "simulation_rate_rows": (
            None
            if case.simulation_rate_rows is None
            else [asdict(row) for row in case.simulation_rate_rows]
        ),
        "simulation_root_state": case.simulation_root_state,
        "simulation_root_state_probabilities": case.simulation_root_state_probabilities,
        "simulation_replicate_count": case.simulation_replicate_count,
        "simulation_seed": case.simulation_seed,
        "continuous_root_state": case.continuous_root_state,
        "continuous_sigma_squared": case.continuous_sigma_squared,
        "continuous_replicate_count": case.continuous_replicate_count,
        "continuous_seed": case.continuous_seed,
        "continuous_trait_names": case.continuous_trait_names,
        "continuous_root_states": case.continuous_root_states,
        "continuous_covariance_matrix": case.continuous_covariance_matrix,
        "comparative_formula": case.comparative_formula,
        "comparative_predictors": case.comparative_predictors,
        "comparative_lambda_value": case.comparative_lambda_value,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def _build_bijux_case_payload(
    case: PhytoolsParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None]:
    tree_path = case.input_fixtures[0]
    traits_path = case.input_fixtures[1] if len(case.input_fixtures) > 1 else None
    comparative_payload = build_comparative_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if comparative_payload is not None:
        return comparative_payload
    discrete_payload = build_discrete_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if discrete_payload is not None:
        return discrete_payload
    continuous_payload = build_continuous_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if continuous_payload is not None:
        return continuous_payload
    signal_payload = build_signal_case_payload(
        case,
        tree_path=tree_path,
        traits_path=traits_path,
    )
    if signal_payload is not None:
        return signal_payload
    raise ValueError(f"unsupported phytools parity operation: {case.operation}")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_rows_table(path: Path) -> list[dict[str, object]]:
    string_identity_fields = {
        "row_kind",
        "label",
        "source_state",
        "target_state",
        "node",
        "state",
    }
    boolean_fields = {
        "transition_allowed",
    }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, object]] = []
        for row in reader:
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if value is None or value == "":
                    parsed[key] = ""
                    continue
                if key in string_identity_fields:
                    parsed[key] = value
                    continue
                if key in boolean_fields:
                    lowered = value.lower()
                    if lowered == "true":
                        parsed[key] = True
                        continue
                    if lowered == "false":
                        parsed[key] = False
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


def _field_tolerance(case: PhytoolsParityCase, key: str) -> float:
    if case.field_tolerances and key in case.field_tolerances:
        return case.field_tolerances[key]
    return case.tolerance


def _row_field_tolerance(case: PhytoolsParityCase, key: str) -> float:
    if case.row_field_tolerances and key in case.row_field_tolerances:
        return case.row_field_tolerances[key]
    return case.tolerance


def _mismatch_reason(
    case: PhytoolsParityCase,
    *,
    reference_summary: dict[str, object] | None,
    bijux_summary: dict[str, object] | None,
) -> str | None:
    if reference_summary is None or bijux_summary is None:
        return "summary_missing"
    if case.operation == "phylogenetic-signal-lambda":
        compare_keys = ("taxon_count", "trait_name", "lambda_value", "log_likelihood")
    elif case.operation == "phylogenetic-signal-k":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "k",
            "p_value",
            "permutation_count",
            "permutation_seed",
            "simulated_k_minimum",
            "simulated_k_mean",
        )
    elif case.operation == "discrete-fit-mk":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
        )
    elif case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
    }:
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
            "requested_replicate_count",
            "successful_replicate_count",
            "simulation_failure_count",
            "seed",
            "mean_total_transition_count",
            "lower_95_total_transition_count",
            "upper_95_total_transition_count",
        )
        if case.operation == "discrete-stochastic-map":
            compare_keys = (
                compare_keys[:16]
                + ("conditioned_on_node_estimates",)
                + compare_keys[16:]
            )
        elif case.operation == "discrete-stochastic-map-description":
            compare_keys = compare_keys + ("branch_count",)
    elif case.operation == "discrete-stochastic-map-density":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "overparameterized",
            "baseline_model",
            "preferred_model_by_aic",
            "requested_replicate_count",
            "successful_replicate_count",
            "simulation_failure_count",
            "seed",
            "branch_count",
            "focal_state",
            "baseline_state",
            "resolution",
            "total_tree_depth",
        )
    elif case.operation == "simulate-discrete-history":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "branch_count",
            "state_count",
            "requested_replicate_count",
            "successful_replicate_count",
            "fixed_root_state",
            "root_prior_probabilities",
            "seed",
            "mean_total_transition_count",
            "lower_95_total_transition_count",
            "upper_95_total_transition_count",
        )
    elif case.operation == "simulate-continuous-brownian":
        compare_keys = (
            "taxon_count",
            "branch_count",
            "requested_replicate_count",
            "successful_replicate_count",
            "seed",
            "root_state",
            "sigma_squared",
        )
    elif case.operation == "simulate-continuous-correlated-brownian":
        compare_keys = (
            "taxon_count",
            "branch_count",
            "trait_count",
            "requested_replicate_count",
            "successful_replicate_count",
            "seed",
        )
    elif case.operation == "comparative-pgls-brownian":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "formula",
            "analysis_taxon_count",
            "coefficient_count",
            "model_matrix_row_count",
            "model_matrix_column_count",
            "categorical_predictor_count",
            "interaction_term_count",
            "lambda_value",
            "lambda_estimation_mode",
            "log_likelihood",
            "aic",
        )
    elif case.operation == "phylogenetic-residuals":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "predictor_name",
            "method",
            "excluded_taxon_count",
            "excluded_taxa",
        )
        if (
            case.comparative_lambda_value is None
            or not math.isclose(case.comparative_lambda_value, 1.0, abs_tol=1e-12)
        ):
            compare_keys = compare_keys + ("lambda_value", "log_likelihood")
    elif case.operation == "phylogenetic-anova":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "group_column",
            "excluded_taxon_count",
            "excluded_taxa",
            "group_count",
            "simulation_count",
            "seed",
            "pairwise_adjustment_method",
            "brownian_sigma_squared",
            "sum_of_squares_between",
            "sum_of_squares_within",
            "mean_square_between",
            "mean_square_within",
            "f_statistic",
            "p_value",
            "low_sample_group_count",
        )
    elif case.operation == "discrete-ancestral-rerooting":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "excluded_taxon_count",
            "excluded_taxa",
            "model",
            "state_count",
            "internal_node_count",
            "root_prior_mode",
            "phytools_rerooting_method_comparable",
        )
    elif case.operation == "continuous-ancestral-fast-anc":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "internal_node_count",
            "excluded_taxon_count",
            "excluded_taxa",
            "tree_is_ultrametric",
        )
    elif case.operation == "continuous-ancestral-anc-ml":
        compare_keys = (
            "taxon_count",
            "trait_name",
            "internal_node_count",
            "excluded_taxon_count",
            "excluded_taxa",
            "tree_is_ultrametric",
            "sigma_squared",
            "log_likelihood",
        )
    else:
        return "unsupported_operation"
    for key in compare_keys:
        if key not in reference_summary or key not in bijux_summary:
            return f"summary_field_missing:{key}"
        if not _isclose(
            reference_summary[key],
            bijux_summary[key],
            tolerance=_field_tolerance(case, key),
        ):
            return f"summary_mismatch:{key}"
    return None


def _row_mismatch_reason(
    case: PhytoolsParityCase,
    *,
    reference_rows: list[dict[str, object]] | None,
    bijux_rows: list[dict[str, object]] | None,
) -> str | None:
    if not case.compare_rows:
        return None
    if case.operation not in {
        "discrete-fit-mk",
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "simulate-discrete-history",
        "simulate-continuous-brownian",
        "simulate-continuous-correlated-brownian",
        "comparative-pgls-brownian",
        "phylogenetic-residuals",
        "phylogenetic-anova",
        "discrete-ancestral-rerooting",
        "continuous-ancestral-fast-anc",
        "continuous-ancestral-anc-ml",
    }:
        return None
    if reference_rows is None or bijux_rows is None:
        return "rows_missing"
    if case.operation == "discrete-fit-mk":
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("source_state", "")),
                str(row.get("target_state", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("source_state", "")),
                str(row.get("target_state", "")),
            ),
        )
    elif case.operation == "discrete-ancestral-rerooting":
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("node", "")),
                str(row.get("state", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("node", "")),
                str(row.get("state", "")),
            ),
        )
    elif case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "simulate-discrete-history",
        "simulate-continuous-brownian",
        "phylogenetic-residuals",
        "phylogenetic-anova",
    }:
        reference_rows = sorted(
            reference_rows,
            key=lambda row: (
                str(row.get("row_kind", "")),
                str(row.get("label", "")),
            ),
        )
        bijux_rows = sorted(
            bijux_rows,
            key=lambda row: (
                str(row.get("row_kind", "")),
                str(row.get("label", "")),
            ),
        )
    else:
        reference_rows = sorted(
            reference_rows, key=lambda row: str(row.get("node", ""))
        )
        bijux_rows = sorted(bijux_rows, key=lambda row: str(row.get("node", "")))
    if len(reference_rows) != len(bijux_rows):
        return "row_count_mismatch"
    if case.operation == "discrete-fit-mk":
        compare_keys = (
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        )
    elif case.operation == "discrete-ancestral-rerooting":
        compare_keys = (
            "node",
            "state",
            "probability",
        )
    elif case.operation in {
        "discrete-stochastic-map",
        "discrete-stochastic-map-count",
        "discrete-stochastic-map-description",
        "simulate-discrete-history",
    }:
        compare_keys = (
            "row_kind",
            "label",
            "mean_value",
            "lower_95_interval",
            "upper_95_interval",
            "presence_fraction",
        )
    elif (
        case.operation == "simulate-continuous-brownian"
        or case.operation == "simulate-continuous-correlated-brownian"
    ):
        compare_keys = (
            "row_kind",
            "label",
            "mean_value",
            "standard_deviation",
            "minimum",
            "median",
            "maximum",
            "covariance",
            "correlation",
        )
    elif case.operation == "comparative-pgls-brownian":
        compare_keys = (
            "row_kind",
            "label",
            "value",
        )
    else:
        compare_keys = (
            ("node", "estimate", "standard_error")
            if case.operation == "continuous-ancestral-fast-anc"
            else (
                "node",
                "estimate",
                "standard_error",
                "lower_95_interval",
                "upper_95_interval",
            )
        )
    for reference_row, bijux_row in zip(reference_rows, bijux_rows, strict=True):
        if case.operation == "phylogenetic-residuals":
            row_kind = str(reference_row.get("row_kind", ""))
            if row_kind != str(bijux_row.get("row_kind", "")):
                return "row_mismatch:row_kind"
            compare_keys = (
                ("row_kind", "label", "value")
                if row_kind == "coefficient_estimate"
                else (
                    "row_kind",
                    "label",
                    "observed_value",
                    "fitted_value",
                    "residual",
                )
            )
        elif case.operation == "phylogenetic-anova":
            row_kind = str(reference_row.get("row_kind", ""))
            if row_kind != str(bijux_row.get("row_kind", "")):
                return "row_mismatch:row_kind"
            compare_keys = (
                (
                    "row_kind",
                    "label",
                    "taxon_count",
                    "taxa",
                    "mean",
                    "variance",
                    "minimum",
                    "maximum",
                )
                if row_kind == "group_summary"
                else (
                    "row_kind",
                    "label",
                    "left_taxon_count",
                    "right_taxon_count",
                    "observed_t_statistic",
                    "uncorrected_p_value",
                    "adjusted_p_value",
                )
            )
        for key in compare_keys:
            if key not in reference_row or key not in bijux_row:
                return f"row_field_missing:{key}"
            if not _isclose(
                reference_row[key],
                bijux_row[key],
                tolerance=_row_field_tolerance(case, key),
            ):
                return f"row_mismatch:{key}"
    return None


def _persist_failure_bundle(
    *,
    failure_root: Path,
    case: PhytoolsParityCase,
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
        for child in artifact_root.iterdir():
            if child.is_file():
                child.unlink()
    artifact_root.mkdir(parents=True, exist_ok=True)
    (artifact_root / "case.json").write_text(
        case_file.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    if execution_root.exists():
        for source in execution_root.iterdir():
            if source.is_file():
                (artifact_root / source.name).write_text(
                    source.read_text(encoding="utf-8"),
                    encoding="utf-8",
                )
    if execution_payload is not None:
        (artifact_root / "reference-execution-copy.json").write_text(
            json.dumps(execution_payload, indent=2),
            encoding="utf-8",
        )
    if reference_summary is not None:
        (artifact_root / "reference-summary-copy.json").write_text(
            json.dumps(reference_summary, indent=2),
            encoding="utf-8",
        )
    if bijux_summary is not None:
        (artifact_root / "bijux-summary.json").write_text(
            json.dumps(bijux_summary, indent=2),
            encoding="utf-8",
        )
    if reference_rows is not None:
        (artifact_root / "reference-rows.json").write_text(
            json.dumps(reference_rows, indent=2),
            encoding="utf-8",
        )
    if bijux_rows is not None:
        (artifact_root / "bijux-rows.json").write_text(
            json.dumps(bijux_rows, indent=2),
            encoding="utf-8",
        )
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
    observations: list[PhytoolsParityObservation],
) -> list[PhytoolsParitySummaryRow]:
    function_names = sorted({item.function_name for item in observations})
    rows: list[PhytoolsParitySummaryRow] = []
    for function_name in function_names:
        matching = [
            observation
            for observation in observations
            if observation.function_name == function_name
        ]
        rows.append(
            PhytoolsParitySummaryRow(
                function_name=function_name,
                case_count=len(matching),
                passed_case_count=sum(
                    1 for observation in matching if observation.status == "passed"
                ),
                failed_case_count=sum(
                    1 for observation in matching if observation.status == "failed"
                ),
                skipped_case_count=sum(
                    1 for observation in matching if observation.status == "skipped"
                ),
            )
        )
    return rows


def run_phytools_parity_cases(
    *,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
) -> PhytoolsParityReport:
    """Run governed live `phytools` parity cases through the checked-in R runner."""
    selected = _selected_cases(case_ids)
    observations: list[PhytoolsParityObservation] = []
    active_failure_root = _failure_root() if failure_root is None else failure_root
    bijux_version = _bijux_version()
    bijux_commit = _bijux_commit()
    for case in selected:
        with tempfile.TemporaryDirectory(
            prefix=f"bijux-phytools-parity-{case.case_id}-"
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
            phytools_version: str | None = None
            process_stdout = ""
            process_stderr = ""
            try:
                # Repository-owned R parity runner.
                process = subprocess.run(  # nosec
                    [
                        rscript_executable,
                        str(_phytools_runner_path()),
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
                if not execution_path.exists():
                    mismatch_reason = "reference_execution_failed"
                else:
                    execution_payload = _load_json(execution_path)
                    r_version = _optional_payload_string(execution_payload, "r_version")
                    phytools_version = _optional_payload_string(
                        execution_payload, "phytools_version"
                    )
                    execution_status = execution_payload.get("status")
                    if execution_status == "unavailable":
                        status = "skipped"
                        mismatch_reason = str(
                            execution_payload.get(
                                "mismatch_reason",
                                "phytools_package_unavailable",
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
                        mismatch_reason = _mismatch_reason(
                            case,
                            reference_summary=reference_summary,
                            bijux_summary=bijux_summary,
                        )
                        if (
                            mismatch_reason is None
                            and case.compare_rows
                            and case.operation
                            in {
                                "discrete-fit-mk",
                                "discrete-stochastic-map",
                                "discrete-stochastic-map-count",
                                "discrete-stochastic-map-description",
                                "discrete-stochastic-map-density",
                                "simulate-discrete-history",
                                "simulate-continuous-brownian",
                                "simulate-continuous-correlated-brownian",
                                "comparative-pgls-brownian",
                                "phylogenetic-residuals",
                                "phylogenetic-anova",
                                "discrete-ancestral-rerooting",
                                "continuous-ancestral-fast-anc",
                                "continuous-ancestral-anc-ml",
                            }
                        ):
                            rows_path = execution_root / (
                                "fitmk-rate-matrix.tsv"
                                if case.operation == "discrete-fit-mk"
                                else (
                                    "stochastic-map-summary-rows.tsv"
                                    if case.operation
                                    in {
                                        "discrete-stochastic-map",
                                        "discrete-stochastic-map-count",
                                        "discrete-stochastic-map-description",
                                        "discrete-stochastic-map-density",
                                        "simulate-discrete-history",
                                    }
                                    else (
                                        "fastbm-summary-rows.tsv"
                                        if case.operation
                                        == "simulate-continuous-brownian"
                                        else (
                                            "simcorrs-summary-rows.tsv"
                                            if case.operation
                                            == "simulate-continuous-correlated-brownian"
                                            else (
                                                "pgls-summary-rows.tsv"
                                                if case.operation
                                                == "comparative-pgls-brownian"
                                                else (
                                                    "phyl-resid-summary-rows.tsv"
                                                    if case.operation
                                                    == "phylogenetic-residuals"
                                                    else (
                                                        "phyl-anova-summary-rows.tsv"
                                                        if case.operation
                                                        == "phylogenetic-anova"
                                                        else (
                                                            "rerooting-method-node-probabilities.tsv"
                                                            if case.operation
                                                            == "discrete-ancestral-rerooting"
                                                            else (
                                                                "fast-anc-node-estimates.tsv"
                                                                if case.operation
                                                                == "continuous-ancestral-fast-anc"
                                                                else "anc-ml-node-estimates.tsv"
                                                            )
                                                        )
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                            if not rows_path.exists():
                                mismatch_reason = "reference_rows_missing"
                            else:
                                reference_rows = _load_rows_table(rows_path)
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
                PhytoolsParityObservation(
                    case_id=case.case_id,
                    fixture_id=case.fixture_id,
                    function_name=case.function_name,
                    python_function_name=case.python_function_name,
                    input_fixtures=case.input_fixtures,
                    tolerance=case.tolerance,
                    r_version=r_version,
                    phytools_version=phytools_version,
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
    return PhytoolsParityReport(
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
            "The governed live `phytools` parity registry is intentionally narrow until later rounds expand the comparative fixture surface.",
            "This harness requires Rscript plus the `phytools` and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_phytools_parity_summary_table(
    path: Path,
    report: PhytoolsParityReport,
) -> Path:
    """Write one row per governed `phytools` function summary."""
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


def write_phytools_parity_observation_table(
    path: Path,
    report: PhytoolsParityReport,
) -> Path:
    """Write one row per governed `phytools` parity observation."""
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
                "tolerance",
                "r_version",
                "phytools_version",
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
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "phytools_version": observation.phytools_version or "",
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
