from __future__ import annotations

from pathlib import Path

from tests.support.fake_external_engines import write_executable


def fake_geiger_rscript(
    path: Path,
    *,
    geiger_available: bool = True,
    summary_overrides: dict[str, dict[str, object]] | None = None,
    reference_payloads: dict[str, dict[str, object]] | None = None,
) -> Path:
    summary_payload = repr(summary_overrides or {})
    reference_payload_repr = repr(reference_payloads or {})
    script = (
        """#!/usr/bin/env python3
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SUMMARY_OVERRIDES = __SUMMARY_OVERRIDES__
REFERENCE_PAYLOADS = __REFERENCE_PAYLOADS__


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys()) if rows else ["parameter", "value"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\\t")
        writer.writeheader()
        writer.writerows(rows)


def find_package_root() -> Path:
    for candidate in (Path.cwd(), *Path.cwd().parents):
        package_root = candidate / "packages" / "bijux-phylogenetics"
        if package_root.is_dir():
            return package_root
    raise ValueError("could not locate bijux-phylogenetics package root")


def find_repo_python(package_root: Path) -> Path:
    candidates = (
        package_root / ".venv" / "bin" / "python",
        package_root.parent / ".venv" / "bin" / "python",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError("could not locate repo python for fake geiger parity payload")


def parameter_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    for parameter in ("root_state", "rate", "log_likelihood", "aic", "aicc"):
        value = summary.get(parameter)
        if value is None:
            continue
        rows.append({"parameter": parameter, "value": value})
    parameter_name = summary.get("parameter_name")
    parameter_value = summary.get("parameter_value")
    if parameter_name and parameter_value is not None:
        rows.append({"parameter": parameter_name, "value": parameter_value})
    return rows


def comparison_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    return [dict(row) for row in summary["rows"]]


def rate_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    return [dict(row) for row in summary["rows"]]


def standard_error_policy() -> str:
    return "fitcontinuous-standard-error-explicitly-excluded-this-round"


def discrete_missing_value_policy() -> str:
    return "prune-overlapping-missing-values"


def missing_value_policy() -> str:
    return "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values"


def build_reference_payload(case_payload: dict[str, object]) -> tuple[dict[str, object], list[dict[str, object]]]:
    if case_payload["case_id"] in REFERENCE_PAYLOADS:
        summary = dict(REFERENCE_PAYLOADS[case_payload["case_id"]]["summary"])
        rows = [dict(row) for row in REFERENCE_PAYLOADS[case_payload["case_id"]]["rows"]]
        return summary, rows
    package_root = find_package_root()
    repo_python = find_repo_python(package_root)
    payload_root = package_root / "artifacts" / "fake-geiger-cases"
    payload_root.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        suffix=".json",
        prefix="case-",
        dir=payload_root,
        delete=False,
    ) as handle:
        handle.write(json.dumps(case_payload))
        payload_path = Path(handle.name)
    inline_script = '''
import json
import math
from pathlib import Path
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.comparative.common import summarize_numeric_trait_readiness
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousModeSearchControls,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
)

case_payload = json.loads(Path(__PAYLOAD_PATH__).read_text(encoding="utf-8"))
mode_lookup = {
    "BM": "brownian",
    "white": "white-noise",
    "lambda": "pagel-lambda",
    "kappa": "pagel-kappa",
    "delta": "pagel-delta",
    "OU": "ornstein-uhlenbeck",
    "EB": "early-burst",
}
tree_path, traits_path = case_payload["input_fixtures"]
if case_payload["operation"] == "fit-discrete-mk":
    tree = load_tree(Path(tree_path))
    table = load_taxon_table(Path(traits_path), taxon_column=case_payload["taxon_column"])
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_only_taxa = sorted(set(tree.tip_names) - set(rows_by_taxon))
    extra_trait_taxa = sorted(set(rows_by_taxon) - set(tree.tip_names))
    report = fit_discrete_mk_model(
        Path(tree_path),
        Path(traits_path),
        trait=case_payload["trait_name"],
        taxon_column=case_payload["taxon_column"],
        model=case_payload["python_mode"],
        transform=case_payload.get("discrete_transform_name"),
        lambda_bounds=tuple(case_payload.get("lambda_bounds") or (0.0, 1.0)),
        kappa_bounds=tuple(case_payload.get("kappa_bounds") or (0.0, 1.0)),
        delta_bounds=tuple(case_payload.get("delta_bounds") or (math.exp(-5.0), 3.0)),
        early_burst_bounds=tuple(case_payload.get("early_burst_bounds") or (-10.0, 10.0)),
    )
    missing_value_taxa = sorted(
        set(report.input_audit.pruned_missing_value_taxa) - set(tree_only_taxa)
    )
    excluded_taxa = sorted(set(tree_only_taxa) | set(missing_value_taxa))
    transform_fit = report.transform_fit
    summary = {
        "taxon_count": report.taxon_count,
        "trait_name": report.trait,
        "model_name": case_payload["model_name"],
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
        "observed_state_count": len(report.state_order),
        "state_order": list(report.state_order),
        "excluded_taxon_count": len(excluded_taxa),
        "excluded_taxa": excluded_taxa,
        "missing_value_taxa": missing_value_taxa,
        "missing_from_traits": tree_only_taxa,
        "extra_trait_taxa": extra_trait_taxa,
        "missing_value_policy": "__DISCRETE_MISSING_VALUE_POLICY__",
        "log_likelihood": report.log_likelihood,
        "parameter_count": report.parameter_count,
        "aic": report.aic,
        "aicc": report.aicc,
        "parameter_name": None if transform_fit is None else transform_fit.parameter_name,
        "parameter_value": None if transform_fit is None else transform_fit.parameter_value,
        "optimizer_settings": case_payload["optimizer_settings"],
        "optimizer_result": {
            "convergence_code": 0,
            "message": "fake geiger parity runner",
        },
        "rows": [
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "transition_allowed": row.transition_allowed,
                "step_distance": row.step_distance,
                "rate": row.rate,
            }
            for row in report.transition_rate_rows
        ],
    }
elif case_payload["operation"] == "compare-fitcontinuous-models":
    comparison = compare_fitcontinuous_model_ranking(
        Path(tree_path),
        Path(traits_path),
        trait=case_payload["trait_name"],
        taxon_column=case_payload["taxon_column"],
        modes=tuple(
            mode_lookup[model_name]
            for model_name in case_payload["candidate_model_names"]
        ),
        lambda_bounds=tuple(case_payload.get("lambda_bounds") or (0.0, 1.0)),
        kappa_bounds=tuple(case_payload.get("kappa_bounds") or (0.0, 3.0)),
        delta_bounds=tuple(case_payload.get("delta_bounds") or (0.0, 3.0)),
        ou_bounds=tuple(case_payload.get("ou_bounds") or (0.0, 10.0)),
        early_burst_bounds=tuple(case_payload.get("early_burst_bounds") or (0.0, 50.0)),
    )
    summary = {
        "taxon_count": comparison.taxon_count,
        "trait_name": comparison.trait,
        "model_name": case_payload["model_name"],
        "selected_model": comparison.better_model,
        "model_ranking": [row.model for row in comparison.rows],
        "comparable_model_count": sum(1 for row in comparison.rows if row.comparable),
        "noncomparable_model_count": sum(1 for row in comparison.rows if not row.comparable),
        "runner_up_model": next((row.model for row in comparison.rows if row.rank == 2), None),
        "runner_up_aicc_delta": next((row.delta_aicc for row in comparison.rows if row.rank == 2), None),
        "warning_count": len(comparison.warnings),
        "optimizer_settings": case_payload["optimizer_settings"],
        "rows": [
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
                "likelihood_constant_policy": (
                    ""
                    if row.likelihood_constant_policy is None
                    else row.likelihood_constant_policy
                ),
                "comparability_note": "" if row.comparability_note is None else row.comparability_note,
            }
            for row in comparison.rows
        ],
    }
else:
    readiness = summarize_numeric_trait_readiness(
        Path(tree_path),
        Path(traits_path),
        trait=case_payload["trait_name"],
        taxon_column=case_payload["taxon_column"],
    )
    report = fit_continuous_evolutionary_mode(
        Path(tree_path),
        Path(traits_path),
        trait=case_payload["trait_name"],
        mode=mode_lookup[case_payload["model_name"]],
        taxon_column=case_payload["taxon_column"],
        search_controls=(
            ContinuousModeSearchControls(
                coarse_grid_point_count=case_payload.get("coarse_grid_point_count") or 81,
                fine_grid_point_count=case_payload.get("fine_grid_point_count") or 81,
                initial_parameter_value=case_payload.get("initial_parameter_value"),
            )
            if (
                case_payload.get("coarse_grid_point_count") is not None
                or case_payload.get("fine_grid_point_count") is not None
                or case_payload.get("initial_parameter_value") is not None
            )
            else None
        ),
        lambda_bounds=tuple(case_payload.get("lambda_bounds") or (0.0, 1.0)),
        kappa_bounds=tuple(case_payload.get("kappa_bounds") or (0.0, 3.0)),
        delta_bounds=tuple(case_payload.get("delta_bounds") or (0.0, 3.0)),
        ou_bounds=tuple(case_payload.get("ou_bounds") or (0.0, 10.0)),
        early_burst_bounds=tuple(case_payload.get("early_burst_bounds") or (0.0, 50.0)),
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
        "model_name": case_payload["model_name"],
        "excluded_taxon_count": len(excluded_taxa),
        "excluded_taxa": excluded_taxa,
        "missing_value_taxa": list(readiness.pruned_missing_value_taxa),
        "non_numeric_taxa": list(readiness.pruned_non_numeric_taxa),
        "missing_from_traits": list(readiness.missing_from_traits),
        "extra_trait_taxa": list(readiness.extra_trait_taxa),
        "missing_value_policy": "__MISSING_VALUE_POLICY__",
        "standard_error_policy": "__STANDARD_ERROR_POLICY__",
        "parameter_bound_policy": (
            "closed-form-without-parameter-bounds"
            if case_payload["model_name"] in {"BM", "white"}
            else "governed-bounded-grid-search"
        ),
        "hit_lower_parameter_boundary": (
            False if report.optimizer_diagnostics is None else report.optimizer_diagnostics.hit_lower_boundary
        ),
        "hit_upper_parameter_boundary": (
            False if report.optimizer_diagnostics is None else report.optimizer_diagnostics.hit_upper_boundary
        ),
        "identifiability_warning_kinds": [warning.kind for warning in report.identifiability_warnings],
        "identifiability_warning_count": len(report.identifiability_warnings),
        "root_state": report.root_state,
        "rate": report.rate,
        "log_likelihood": report.log_likelihood,
        "aic": report.aic,
        "aicc": report.aicc,
        "parameter_name": report.parameter_name,
        "parameter_value": report.parameter_value,
        "optimizer_settings": case_payload["optimizer_settings"],
        "optimizer_result": {
            "convergence_code": 0,
            "message": "fake geiger parity runner",
        },
    }
print(json.dumps(summary))
'''
    command = [
        str(repo_python),
        "-c",
        inline_script
        .replace("__PAYLOAD_PATH__", repr(str(payload_path)))
        .replace("__DISCRETE_MISSING_VALUE_POLICY__", discrete_missing_value_policy())
        .replace("__MISSING_VALUE_POLICY__", missing_value_policy())
        .replace("__STANDARD_ERROR_POLICY__", standard_error_policy()),
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_root / "src")
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(package_root),
            env=env,
        )
    finally:
        payload_path.unlink(missing_ok=True)
    summary = json.loads(result.stdout)
    if case_payload["case_id"] in SUMMARY_OVERRIDES:
        summary.update(SUMMARY_OVERRIDES[case_payload["case_id"]])
    if case_payload["operation"] == "compare-fitcontinuous-models":
        return summary, comparison_rows(summary)
    if case_payload["operation"] == "fit-discrete-mk":
        return summary, rate_rows(summary)
    return summary, parameter_rows(summary)


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit(2)
    case_path = Path(sys.argv[2])
    output_root = Path(sys.argv[3])
    output_root.mkdir(parents=True, exist_ok=True)
    execution_path = output_root / "reference-execution.json"
    summary_path = output_root / "reference-summary.json"
    rows_path = output_root / "reference-parameters.tsv"
    case_payload = json.loads(case_path.read_text(encoding="utf-8"))
    if not __GEIGER_AVAILABLE__:
        write_json(
            execution_path,
            {
                "status": "unavailable",
                "mismatch_reason": "geiger_package_unavailable",
                "r_version": "4.4.0",
                "geiger_version": None,
            },
        )
        return 0
    summary, rows = build_reference_payload(case_payload)
    write_json(summary_path, summary)
    write_rows_table(rows_path, rows)
    write_json(
        execution_path,
        {
            "status": "ok",
            "r_version": "4.4.0",
            "geiger_version": "2.0.11",
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
""".replace("__GEIGER_AVAILABLE__", "True" if geiger_available else "False")
        .replace("__SUMMARY_OVERRIDES__", summary_payload)
        .replace("__REFERENCE_PAYLOADS__", reference_payload_repr)
    )
    return write_executable(
        path,
        script,
    )
