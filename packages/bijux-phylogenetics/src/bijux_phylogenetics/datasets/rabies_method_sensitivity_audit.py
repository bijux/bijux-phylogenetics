from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.engines.common import file_sha256
from bijux_phylogenetics.io.fasta import load_fasta_alignment

__all__ = [
    "RabiesMethodSensitivityReproducibilityAuditReport",
    "RabiesMethodSensitivityReproducibilityCheckRow",
    "RabiesMethodSensitivityVariantAuditRow",
    "audit_rabies_method_sensitivity_workflow_bundle",
    "write_rabies_method_sensitivity_reproducibility_audit_json",
    "write_rabies_method_sensitivity_reproducibility_checks_table",
    "write_rabies_method_sensitivity_variant_audit_table",
]

_WORKFLOW_MANIFEST_FILENAME = "rabies-method-sensitivity.manifest.json"
_REPORT_MANIFEST_FILENAME = (
    "report-artifacts/rabies-method-sensitivity-report.manifest.json"
)
_CONFIG_FILENAME = "workflow-config.resolved.json"
_PARALLEL_SUMMARY_FILENAME = "parallel-execution-summary.tsv"
_VARIANT_SUMMARY_FILENAME = "variant-summary.tsv"
_TASK_LOGS_DIRECTORY = "parallel-logs"
_VARIANTS_DIRECTORY = "variants"
_EXPECTED_VARIANT_FILENAMES = (
    "fasttree.nwk",
    "iqtree-support.nwk",
    "rooted-engine-comparison.tsv",
    "rooted-fasttree.nwk",
    "rooted-iqtree-support.nwk",
    "rooting-summary.tsv",
    "unrooted-comparison.tsv",
    "unrooted-conclusions.tsv",
    "unrooted-conflicting-clades.tsv",
    "unrooted-shared-clades.tsv",
    "unrooted-stability-summary.tsv",
    "unrooted-support-weighted-conflicts.tsv",
)


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityReproducibilityCheckRow:
    """One machine-readable pass/fail check within the bundle audit."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariantAuditRow:
    """One per-variant provenance and file-inventory summary."""

    variant_id: str
    status: str
    output_file_count: int
    output_byte_count: int
    output_digest: str
    missing_required_files: tuple[str, ...]
    unexpected_files: tuple[str, ...]
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityReproducibilityAuditReport:
    """One reviewer-facing reproducibility audit for the workflow bundle."""

    dataset_id: str
    bundle_root: Path
    workflow_manifest_path: Path
    report_manifest_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    all_passed: bool
    check_count: int
    failed_check_count: int
    variant_count: int
    failed_variant_count: int
    checks: tuple[RabiesMethodSensitivityReproducibilityCheckRow, ...]
    variants: tuple[RabiesMethodSensitivityVariantAuditRow, ...]


def audit_rabies_method_sensitivity_workflow_bundle(
    bundle_root: Path,
    *,
    sequences_path: Path,
    metadata_path: Path,
) -> RabiesMethodSensitivityReproducibilityAuditReport:
    """Audit one written bundle against current inputs, settings, and outputs."""
    bundle_root = bundle_root.resolve()
    workflow_manifest_path = bundle_root / _WORKFLOW_MANIFEST_FILENAME
    report_manifest_path = bundle_root / _REPORT_MANIFEST_FILENAME
    config_path = bundle_root / _CONFIG_FILENAME
    parallel_summary_path = bundle_root / _PARALLEL_SUMMARY_FILENAME
    variant_summary_path = bundle_root / _VARIANT_SUMMARY_FILENAME
    task_logs_root = bundle_root / _TASK_LOGS_DIRECTORY
    variants_root = bundle_root / _VARIANTS_DIRECTORY

    workflow_manifest = _load_json(workflow_manifest_path)
    report_manifest = _load_json(report_manifest_path)
    resolved_config = _load_json(config_path)
    parallel_rows = _read_tsv_rows(parallel_summary_path)
    variant_rows = _read_tsv_rows(variant_summary_path)

    checks: list[RabiesMethodSensitivityReproducibilityCheckRow] = []
    variant_audit_rows: list[RabiesMethodSensitivityVariantAuditRow] = []

    def add_check(
        check_id: str,
        *,
        surface: str,
        condition: bool,
        expected: object,
        observed: object,
        detail: str,
    ) -> None:
        checks.append(
            RabiesMethodSensitivityReproducibilityCheckRow(
                check_id=check_id,
                surface=surface,
                status="passed" if condition else "failed",
                expected="" if expected is None else str(expected),
                observed="" if observed is None else str(observed),
                detail=detail,
            )
        )

    input_paths = {
        "sequences.fasta": sequences_path,
        "metadata.csv": metadata_path,
    }
    recorded_input_checksums = {
        str(key): str(value)
        for key, value in dict(resolved_config.get("input_checksums", {})).items()
    }
    for filename, input_path in input_paths.items():
        recorded_checksum = recorded_input_checksums.get(filename)
        observed_checksum = file_sha256(input_path)
        add_check(
            f"input-checksum:{filename}",
            surface="input-checksum",
            condition=recorded_checksum == observed_checksum,
            expected=recorded_checksum,
            observed=observed_checksum,
            detail=f"{filename} matches the checksum recorded in workflow-config.resolved.json",
        )

    manifest_output_paths = {
        str(key): workflow_manifest_path.parent / Path(value)
        for key, value in dict(workflow_manifest.get("output_paths", {})).items()
    }
    manifest_output_checksums = {
        str(key): str(value)
        for key, value in dict(workflow_manifest.get("output_checksums", {})).items()
    }
    for key, output_path in sorted(manifest_output_paths.items()):
        if key in {"task_logs_root", "variants_root"}:
            add_check(
                f"workflow-manifest:{key}",
                surface="workflow-manifest",
                condition=output_path.is_dir(),
                expected="directory exists",
                observed="present" if output_path.is_dir() else "missing",
                detail=f"{key} resolves to an existing directory",
            )
            continue
        recorded_checksum = manifest_output_checksums.get(key)
        output_exists = output_path.is_file()
        observed_checksum = None if not output_exists else file_sha256(output_path)
        add_check(
            f"workflow-manifest:{key}",
            surface="workflow-manifest",
            condition=output_exists and recorded_checksum == observed_checksum,
            expected=recorded_checksum,
            observed=observed_checksum,
            detail=f"{key} matches the checksum recorded in {_WORKFLOW_MANIFEST_FILENAME}",
        )

    linked_artifacts = dict(report_manifest.get("linked_artifacts", {}))
    linked_artifact_count = int(report_manifest.get("linked_artifact_count", 0))
    add_check(
        "report-manifest:linked-artifact-count",
        surface="report-manifest",
        condition=linked_artifact_count == len(linked_artifacts),
        expected=linked_artifact_count,
        observed=len(linked_artifacts),
        detail="linked artifact count matches the linked_artifacts payload",
    )
    for key, payload in sorted(linked_artifacts.items()):
        artifact_path = (report_manifest_path.parent / str(payload["path"])).resolve()
        output_exists = artifact_path.is_file()
        observed_checksum = None if not output_exists else file_sha256(artifact_path)
        add_check(
            f"report-manifest:{key}",
            surface="report-manifest",
            condition=output_exists and str(payload["sha256"]) == observed_checksum,
            expected=payload["sha256"],
            observed=observed_checksum,
            detail=f"{key} matches the checksum recorded in the report manifest",
        )

    config_variants = {
        str(row["variant_id"]): row for row in list(resolved_config.get("variants", []))
    }
    manifest_task_records = {
        str(row["variant_id"]): row for row in list(workflow_manifest.get("task_records", []))
    }
    parallel_summary_rows = {
        str(row["variant_id"]): row for row in parallel_rows
    }
    variant_summary_rows = {
        str(row["variant_id"]): row for row in variant_rows
    }
    logged_variants = {
        path.stem: _parse_task_log(path) for path in sorted(task_logs_root.glob("*.log"))
    }

    config_variant_ids = sorted(config_variants)
    manifest_variant_ids = sorted(manifest_task_records)
    parallel_variant_ids = sorted(parallel_summary_rows)
    summary_variant_ids = sorted(variant_summary_rows)
    logged_variant_ids = sorted(logged_variants)
    written_variant_ids = sorted(
        path.name for path in variants_root.iterdir() if path.is_dir()
    )
    expected_variant_count = len(config_variant_ids)
    add_check(
        "variant-sets:manifest",
        surface="variant-sets",
        condition=config_variant_ids == manifest_variant_ids,
        expected=config_variant_ids,
        observed=manifest_variant_ids,
        detail="workflow manifest task records cover the configured variant ids",
    )
    add_check(
        "variant-sets:parallel-summary",
        surface="variant-sets",
        condition=config_variant_ids == parallel_variant_ids,
        expected=config_variant_ids,
        observed=parallel_variant_ids,
        detail="parallel summary rows cover the configured variant ids",
    )
    add_check(
        "variant-sets:variant-summary",
        surface="variant-sets",
        condition=config_variant_ids == summary_variant_ids,
        expected=config_variant_ids,
        observed=summary_variant_ids,
        detail="variant summary rows cover the configured variant ids",
    )
    add_check(
        "variant-sets:task-logs",
        surface="variant-sets",
        condition=config_variant_ids == logged_variant_ids,
        expected=config_variant_ids,
        observed=logged_variant_ids,
        detail="task log files cover the configured variant ids",
    )
    add_check(
        "variant-sets:variant-directories",
        surface="variant-sets",
        condition=config_variant_ids == written_variant_ids,
        expected=config_variant_ids,
        observed=written_variant_ids,
        detail="variant output directories cover the configured variant ids",
    )
    add_check(
        "variant-count:parallel-summary",
        surface="variant-count",
        condition=expected_variant_count == len(parallel_rows),
        expected=expected_variant_count,
        observed=len(parallel_rows),
        detail="parallel summary row count matches the configured variant count",
    )
    add_check(
        "variant-count:variant-summary",
        surface="variant-count",
        condition=expected_variant_count == len(variant_rows),
        expected=expected_variant_count,
        observed=len(variant_rows),
        detail="variant summary row count matches the configured variant count",
    )

    for variant_id in config_variant_ids:
        config_row = config_variants[variant_id]
        manifest_row = manifest_task_records.get(variant_id)
        parallel_row = parallel_summary_rows.get(variant_id)
        summary_row = variant_summary_rows.get(variant_id)
        log_row = logged_variants.get(variant_id)
        variant_root = variants_root / variant_id
        issues: list[str] = []

        if manifest_row is None:
            issues.append("workflow manifest task record is missing")
        if parallel_row is None:
            issues.append("parallel execution summary row is missing")
        if summary_row is None:
            issues.append("variant summary row is missing")
        if log_row is None:
            issues.append("task log is missing")

        expected_output_root = Path("variants", variant_id).as_posix()
        if manifest_row is not None and str(manifest_row.get("output_root")) != expected_output_root:
            issues.append("workflow manifest output_root differs from the expected variant directory")
        if parallel_row is not None and str(parallel_row.get("log_path")) != Path("parallel-logs", f"{variant_id}.log").as_posix():
            issues.append("parallel summary log_path differs from the expected task log path")
        if log_row is not None and str(log_row.get("output_root")) != expected_output_root:
            issues.append("task log output_root differs from the expected variant directory")
        if manifest_row is not None and str(manifest_row.get("status")) != "succeeded":
            issues.append("workflow manifest does not record the variant as succeeded")
        if parallel_row is not None and str(parallel_row.get("status")) != "succeeded":
            issues.append("parallel execution summary does not record the variant as succeeded")
        if log_row is not None and str(log_row.get("status")) != "succeeded":
            issues.append("task log does not record the variant as succeeded")

        for field_name in ("alignment_mode", "trimming_mode"):
            expected_value = str(config_row[field_name])
            if log_row is not None and str(log_row.get(field_name)) != expected_value:
                issues.append(f"task log {field_name} does not match the resolved config")
            if summary_row is not None and str(summary_row.get(field_name)) != expected_value:
                issues.append(f"variant summary {field_name} does not match the resolved config")

        expected_trim_gap_threshold = _format_float(float(config_row["trim_gap_threshold"]))
        if log_row is not None and str(log_row.get("trim_gap_threshold")) != expected_trim_gap_threshold:
            issues.append("task log trim_gap_threshold does not match the resolved config")
        if summary_row is not None and str(summary_row.get("trim_gap_threshold")) != expected_trim_gap_threshold:
            issues.append(
                "variant summary trim_gap_threshold does not match the resolved config"
            )

        variant_output_paths = sorted(
            path for path in variant_root.iterdir() if path.is_file()
        ) if variant_root.is_dir() else []
        variant_filenames = tuple(path.name for path in variant_output_paths)
        missing_required_files = tuple(
            name
            for name in (
                f"{variant_id}.aln",
                f"{variant_id}.trimmed.aln",
                *_EXPECTED_VARIANT_FILENAMES,
            )
            if not (variant_root / name).is_file()
        )
        unexpected_files = tuple(
            name
            for name in variant_filenames
            if name
            not in {
                f"{variant_id}.aln",
                f"{variant_id}.trimmed.aln",
                *_EXPECTED_VARIANT_FILENAMES,
            }
        )
        if not variant_root.is_dir():
            issues.append("variant output directory is missing")
        if missing_required_files:
            issues.append("one or more required variant output files are missing")
        if unexpected_files:
            issues.append("variant output directory contains unexpected files")

        alignment_path = variant_root / f"{variant_id}.aln"
        trimmed_alignment_path = variant_root / f"{variant_id}.trimmed.aln"
        if summary_row is not None and alignment_path.is_file():
            alignment_length = _alignment_length(alignment_path)
            if int(summary_row["alignment_length"]) != alignment_length:
                issues.append(
                    "variant summary alignment_length does not match the written alignment"
                )
        if summary_row is not None and trimmed_alignment_path.is_file():
            trimmed_alignment_length = _alignment_length(trimmed_alignment_path)
            if int(summary_row["trimmed_alignment_length"]) != trimmed_alignment_length:
                issues.append(
                    "variant summary trimmed_alignment_length does not match the written trimmed alignment"
                )

        output_byte_count = sum(path.stat().st_size for path in variant_output_paths)
        output_digest = _directory_digest(variant_root)
        variant_audit_rows.append(
            RabiesMethodSensitivityVariantAuditRow(
                variant_id=variant_id,
                status="passed" if not issues else "failed",
                output_file_count=len(variant_output_paths),
                output_byte_count=output_byte_count,
                output_digest=output_digest,
                missing_required_files=missing_required_files,
                unexpected_files=unexpected_files,
                issues=tuple(issues),
            )
        )

    failed_check_count = sum(1 for row in checks if row.status == "failed")
    failed_variant_count = sum(
        1 for row in variant_audit_rows if row.status == "failed"
    )
    return RabiesMethodSensitivityReproducibilityAuditReport(
        dataset_id=str(workflow_manifest["dataset_id"]),
        bundle_root=bundle_root,
        workflow_manifest_path=workflow_manifest_path,
        report_manifest_path=report_manifest_path,
        config_path=config_path,
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        all_passed=failed_check_count == 0 and failed_variant_count == 0,
        check_count=len(checks),
        failed_check_count=failed_check_count,
        variant_count=len(variant_audit_rows),
        failed_variant_count=failed_variant_count,
        checks=tuple(checks),
        variants=tuple(variant_audit_rows),
    )


def write_rabies_method_sensitivity_reproducibility_checks_table(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one tabular ledger of top-level audit checks."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_variant_audit_table(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one per-variant audit ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "status",
            "output_file_count",
            "output_byte_count",
            "output_digest",
            "missing_required_files",
            "unexpected_files",
            "issues",
        ),
        rows=[
            {
                **asdict(row),
                "missing_required_files": "; ".join(row.missing_required_files),
                "unexpected_files": "; ".join(row.unexpected_files),
                "issues": "; ".join(row.issues),
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_reproducibility_audit_json(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one machine-readable JSON summary for the bundle audit."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return path


def _alignment_length(path: Path) -> int:
    records = load_fasta_alignment(path)
    return len(records[0].sequence)


def _directory_digest(path: Path) -> str:
    if not path.is_dir():
        return ""
    lines = [
        f"{entry.relative_to(path).as_posix()}\t{file_sha256(entry)}"
        for entry in sorted(path.rglob("*"))
        if entry.is_file()
    ]
    payload = "\n".join(lines).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _format_float(value: float) -> str:
    return format(value, ".12g")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_task_log(path: Path) -> dict[str, str]:
    payload: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        payload[key.strip()] = value.strip()
    return payload


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ""
                    if value is None
                    else value
                    for key, value in row.items()
                }
            )
    return path
