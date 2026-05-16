from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import json
from pathlib import Path

from bijux_phylogenetics.errors import EngineWorkflowError

from .common import build_file_checksums, utc_now_text
from .workflows import EngineWorkflowReport


_CANONICAL_ENGINE_NAMES = {
    "mafft": "MAFFT",
    "trimal": "trimAl",
    "iqtree": "IQ-TREE",
    "fasttree": "FastTree",
    "mrbayes": "MrBayes",
    "beast": "BEAST",
}


def _parse_utc_timestamp(text: str) -> datetime:
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def _runtime_seconds(*, started_at_utc: str, ended_at_utc: str) -> float:
    return max(
        (_parse_utc_timestamp(ended_at_utc) - _parse_utc_timestamp(started_at_utc)).total_seconds(),
        0.0,
    )


def _canonical_engine_name(engine_name: str) -> str:
    return _CANONICAL_ENGINE_NAMES.get(engine_name.lower(), engine_name)


def _labeled_output_checksums(
    output_paths: dict[str, Path],
    *,
    stored_checksums: dict[str, str] | None = None,
) -> dict[str, str]:
    checksum_payload = (
        build_file_checksums(list(output_paths.values()))
        if stored_checksums is None or not stored_checksums
        else dict(stored_checksums)
    )
    checksums_by_path = {str(key): value for key, value in checksum_payload.items()}
    return {
        label: checksums_by_path[str(path)]
        for label, path in output_paths.items()
    }


@dataclass(slots=True)
class ExternalEngineValidationCase:
    engine_name: str
    validation_name: str
    validation_mode: str
    manifest_path: Path | None
    executable: str | None
    version_text: str | None
    command: list[str]
    exit_code: int | None
    runtime_seconds: float | None
    output_paths: dict[str, Path]
    output_checksums: dict[str, str]
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ExternalEngineValidationMatrixReport:
    generated_at_utc: str
    cases: list[ExternalEngineValidationCase]


def build_external_engine_validation_case(
    validation_name: str,
    report: EngineWorkflowReport,
    *,
    validation_mode: str = "workflow-run",
    notes: list[str] | None = None,
) -> ExternalEngineValidationCase:
    """Normalize one workflow report into a reviewer-facing validation row."""
    output_checksums = _labeled_output_checksums(
        dict(report.output_paths),
        stored_checksums=(
            None
            if not report.output_checksums
            else dict(report.output_checksums)
        ),
    )
    return ExternalEngineValidationCase(
        engine_name=_canonical_engine_name(report.engine_name),
        validation_name=validation_name,
        validation_mode=validation_mode,
        manifest_path=report.manifest_path,
        executable=report.run.executable,
        version_text=report.run.version.text,
        command=list(report.run.command),
        exit_code=report.run.exit_code,
        runtime_seconds=_runtime_seconds(
            started_at_utc=report.run.started_at_utc,
            ended_at_utc=report.run.ended_at_utc,
        ),
        output_paths=dict(report.output_paths),
        output_checksums=output_checksums,
        notes=[*report.notes, *(notes or [])],
    )


def build_beast_artifact_validation_case(
    validation_name: str,
    *,
    xml_path: Path,
    log_path: Path,
    tree_path: Path,
    burnin_fraction: float = 0.1,
) -> ExternalEngineValidationCase:
    """Build one validation row from governed real BEAST analysis artifacts."""
    from bijux_phylogenetics.bayesian.beast import (
        parse_beast_posterior_tree_samples,
        summarize_beast_analysis_xml,
        summarize_beast_log,
    )

    xml_report = summarize_beast_analysis_xml(xml_path)
    if not xml_report.valid:
        raise EngineWorkflowError(
            f"BEAST analysis XML failed validation for matrix case '{validation_name}': {xml_path}"
        )
    log_summary = summarize_beast_log(log_path, burnin_fraction=burnin_fraction)
    tree_report = parse_beast_posterior_tree_samples(
        tree_path,
        burnin_fraction=burnin_fraction,
    )
    output_paths = {
        "analysis_xml": xml_path,
        "posterior_log": log_path,
        "posterior_trees": tree_path,
    }
    return ExternalEngineValidationCase(
        engine_name="BEAST",
        validation_name=validation_name,
        validation_mode="fixture-parse",
        manifest_path=None,
        executable=None,
        version_text=xml_report.beast_version,
        command=[],
        exit_code=None,
        runtime_seconds=None,
        output_paths=output_paths,
        output_checksums=_labeled_output_checksums(output_paths),
        notes=[
            f"beast xml taxon count: {xml_report.taxon_count}",
            f"beast log kept rows after burn-in: {log_summary.kept_row_count}",
            f"beast posterior trees kept after burn-in: {tree_report.kept_tree_count}",
        ],
    )


def build_external_engine_validation_matrix(
    cases: list[ExternalEngineValidationCase],
) -> ExternalEngineValidationMatrixReport:
    """Assemble one ordered external-engine validation matrix."""
    return ExternalEngineValidationMatrixReport(
        generated_at_utc=utc_now_text(),
        cases=list(cases),
    )


def merge_external_engine_validation_matrices(
    reports: list[ExternalEngineValidationMatrixReport],
) -> ExternalEngineValidationMatrixReport:
    """Merge multiple ordered validation matrices into one reviewer-facing report."""
    merged_cases: list[ExternalEngineValidationCase] = []
    for report in reports:
        merged_cases.extend(report.cases)
    return ExternalEngineValidationMatrixReport(
        generated_at_utc=utc_now_text(),
        cases=merged_cases,
    )


def write_external_engine_validation_matrix(
    path: Path,
    report: ExternalEngineValidationMatrixReport,
) -> Path:
    """Persist one portable JSON validation matrix for reviewer inspection."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at_utc": report.generated_at_utc,
        "case_count": len(report.cases),
        "engine_names": sorted({case.engine_name for case in report.cases}),
        "cases": [
            {
                "engine_name": case.engine_name,
                "validation_name": case.validation_name,
                "validation_mode": case.validation_mode,
                "manifest_path": (
                    None if case.manifest_path is None else str(case.manifest_path)
                ),
                "executable": case.executable,
                "version_text": case.version_text,
                "command": case.command,
                "exit_code": case.exit_code,
                "runtime_seconds": case.runtime_seconds,
                "output_paths": {
                    key: str(value) for key, value in case.output_paths.items()
                },
                "output_checksums": dict(case.output_checksums),
                "notes": list(case.notes),
            }
            for case in report.cases
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path
