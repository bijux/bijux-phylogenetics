from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path

from ..registry import PhytoolsParityCase
from .models import (
    PhytoolsParityObservation,
    PhytoolsParityReport,
    PhytoolsParitySummaryRow,
)


def persist_failure_bundle(
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


def summary_rows(
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


def build_phytools_parity_report(
    observations: list[PhytoolsParityObservation],
) -> PhytoolsParityReport:
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return PhytoolsParityReport(
        observations=observations,
        summary_rows=summary_rows(observations),
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
