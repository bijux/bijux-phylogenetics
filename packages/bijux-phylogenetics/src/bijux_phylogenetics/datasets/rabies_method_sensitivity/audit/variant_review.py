from __future__ import annotations

from pathlib import Path

from .contracts import RabiesMethodSensitivityVariantAuditRow
from .inventory import _EXPECTED_VARIANT_FILENAMES, RabiesMethodSensitivityAuditSnapshot
from .io import _alignment_length, _directory_digest, _format_float
from .review_context import RabiesMethodSensitivityAuditReviewContext


def build_rabies_method_sensitivity_variant_audit_rows(
    *,
    snapshot: RabiesMethodSensitivityAuditSnapshot,
    review_context: RabiesMethodSensitivityAuditReviewContext,
) -> list[RabiesMethodSensitivityVariantAuditRow]:
    variant_audit_rows: list[RabiesMethodSensitivityVariantAuditRow] = []
    for variant_id in review_context.config_variant_ids:
        config_row = review_context.config_variants[variant_id]
        manifest_row = review_context.manifest_task_records.get(variant_id)
        parallel_row = review_context.parallel_summary_rows.get(variant_id)
        summary_row = review_context.variant_summary_rows.get(variant_id)
        log_row = review_context.logged_variants.get(variant_id)
        variant_root = snapshot.variants_root / variant_id
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
        if (
            manifest_row is not None
            and str(manifest_row.get("output_root")) != expected_output_root
        ):
            issues.append(
                "workflow manifest output_root differs from the expected variant directory"
            )
        if (
            parallel_row is not None
            and str(parallel_row.get("log_path"))
            != Path("parallel-logs", f"{variant_id}.log").as_posix()
        ):
            issues.append(
                "parallel summary log_path differs from the expected task log path"
            )
        if (
            log_row is not None
            and str(log_row.get("output_root")) != expected_output_root
        ):
            issues.append(
                "task log output_root differs from the expected variant directory"
            )
        if manifest_row is not None and str(manifest_row.get("status")) != "succeeded":
            issues.append("workflow manifest does not record the variant as succeeded")
        if parallel_row is not None and str(parallel_row.get("status")) != "succeeded":
            issues.append(
                "parallel execution summary does not record the variant as succeeded"
            )
        if log_row is not None and str(log_row.get("status")) != "succeeded":
            issues.append("task log does not record the variant as succeeded")

        for field_name in ("alignment_mode", "trimming_mode"):
            expected_value = str(config_row[field_name])
            if log_row is not None and str(log_row.get(field_name)) != expected_value:
                issues.append(
                    f"task log {field_name} does not match the resolved config"
                )
            if (
                summary_row is not None
                and str(summary_row.get(field_name)) != expected_value
            ):
                issues.append(
                    f"variant summary {field_name} does not match the resolved config"
                )

        expected_trim_gap_threshold = _format_float(
            float(config_row["trim_gap_threshold"])
        )
        if (
            log_row is not None
            and str(log_row.get("trim_gap_threshold")) != expected_trim_gap_threshold
        ):
            issues.append(
                "task log trim_gap_threshold does not match the resolved config"
            )
        if (
            summary_row is not None
            and str(summary_row.get("trim_gap_threshold"))
            != expected_trim_gap_threshold
        ):
            issues.append(
                "variant summary trim_gap_threshold does not match the resolved config"
            )

        variant_output_paths = (
            sorted(path for path in variant_root.iterdir() if path.is_file())
            if variant_root.is_dir()
            else []
        )
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
    return variant_audit_rows
