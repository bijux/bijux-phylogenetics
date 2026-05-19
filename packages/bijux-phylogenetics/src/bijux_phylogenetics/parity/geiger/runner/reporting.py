from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path
import shutil

from bijux_phylogenetics.parity.geiger.registry import GeigerParityCase

from .models import (
    GeigerParityObservation,
    GeigerParityReport,
    GeigerParitySummaryRow,
)


def write_rows_table(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys()) if rows else ["parameter", "value"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def persist_failure_bundle(
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
        write_rows_table(artifact_root / "reference-parameters.tsv", reference_rows)
    if bijux_rows is not None:
        write_rows_table(artifact_root / "bijux-parameters.tsv", bijux_rows)
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


def summary_rows(
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
