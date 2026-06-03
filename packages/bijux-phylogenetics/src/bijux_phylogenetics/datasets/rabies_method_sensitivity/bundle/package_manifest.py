from __future__ import annotations

import hashlib
import json
from pathlib import Path

from ..models import RabiesMethodSensitivityPanelWorkflowReport


def _write_resolved_config(
    path: Path, report: RabiesMethodSensitivityPanelWorkflowReport
) -> Path:
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "label": report.dataset.label,
        "sequence_type": report.dataset.sequence_type,
        "workflow_prefix": report.dataset.workflow_prefix,
        "outgroup_taxa": list(report.dataset.outgroup_taxa),
        "iqtree_seed": report.iqtree_seed,
        "iqtree_threads": report.iqtree_threads,
        "bootstrap_replicates": report.bootstrap_replicates,
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "selected_variant_ids": [
            variant.variant_id for variant in report.dataset.variants
        ],
        "input_checksums": {
            "sequences.fasta": _sha256(report.dataset.sequences_path),
            "metadata.csv": _sha256(report.dataset.metadata_path),
        },
        "variants": [
            {
                "variant_id": variant.variant_id,
                "label": variant.label,
                "alignment_mode": variant.alignment_mode,
                "trimming_mode": variant.trimming_mode,
                "trim_gap_threshold": variant.trim_gap_threshold,
            }
            for variant in report.dataset.variants
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_manifest(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
) -> Path:
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "label": report.dataset.label,
        "report_kind": "rabies_method_sensitivity_workflow_bundle",
        "variant_count": len(report.variant_runs),
        "parallel_execution": {
            "execution_mode": report.execution_mode,
            "parallel_workers": report.parallel_workers,
            "requested_task_count": len(report.task_records),
            "completed_task_count": len(
                [task for task in report.task_records if task.status == "succeeded"]
            ),
            "failed_task_count": len(
                [task for task in report.task_records if task.status != "succeeded"]
            ),
        },
        "task_records": [
            {
                "variant_id": task.variant_id,
                "label": task.label,
                "status": task.status,
                "execution_mode": task.execution_mode,
                "log_path": Path("parallel-logs", task.log_path.name).as_posix(),
                "output_root": Path("variants", task.variant_id).as_posix(),
                "error_code": task.error_code,
                "error_message": task.error_message,
            }
            for task in report.task_records
        ],
        "output_paths": {
            key: value.name
            if value.parent == path.parent
            else value.relative_to(path.parent).as_posix()
            for key, value in bundle_paths.items()
        },
        "output_checksums": {
            key: _sha256(value)
            for key, value in bundle_paths.items()
            if value.is_file()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
