from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from .models import ApeParityObservation, ApeParityReport, ApeParitySummaryRow


def build_summary_rows(
    observations: list[ApeParityObservation],
) -> list[ApeParitySummaryRow]:
    rows: list[ApeParitySummaryRow] = []
    for function_name in sorted({item.function_name for item in observations}):
        selected = [
            item for item in observations if item.function_name == function_name
        ]
        rows.append(
            ApeParitySummaryRow(
                function_name=function_name,
                case_count=len(selected),
                passed_case_count=sum(
                    1 for item in selected if item.status == "passed"
                ),
                failed_case_count=sum(
                    1 for item in selected if item.status == "failed"
                ),
                skipped_case_count=sum(
                    1 for item in selected if item.status == "skipped"
                ),
            )
        )
    return rows


def build_ape_parity_report(
    observations: list[ApeParityObservation],
) -> ApeParityReport:
    case_count = len(observations)
    passed_case_count = sum(1 for item in observations if item.status == "passed")
    failed_case_count = sum(1 for item in observations if item.status == "failed")
    skipped_case_count = sum(1 for item in observations if item.status == "skipped")
    return ApeParityReport(
        observations=observations,
        summary_rows=build_summary_rows(observations),
        case_count=case_count,
        passed_case_count=passed_case_count,
        failed_case_count=failed_case_count,
        skipped_case_count=skipped_case_count,
        all_passed=case_count > 0
        and passed_case_count == case_count
        and failed_case_count == 0
        and skipped_case_count == 0,
        limitations=[
            "The governed live `ape` parity registry is intentionally narrow until later rounds expand the shared fixture surface.",
            "This harness requires Rscript plus the `ape` and `jsonlite` R packages for live reference execution.",
        ],
    )


def write_ape_parity_summary_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` function summary."""
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


def write_ape_parity_observation_table(path: Path, report: ApeParityReport) -> Path:
    """Write one row per governed `ape` parity observation."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case_id",
                "fixture_kind",
                "fixture_id",
                "function_name",
                "python_function_name",
                "input_fixture",
                "tolerance",
                "r_version",
                "ape_version",
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
                    "fixture_kind": observation.fixture_kind,
                    "fixture_id": observation.fixture_id,
                    "function_name": observation.function_name,
                    "python_function_name": observation.python_function_name,
                    "input_fixture": str(observation.input_fixture),
                    "tolerance": format(observation.tolerance, ".12g"),
                    "r_version": observation.r_version or "",
                    "ape_version": observation.ape_version or "",
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
