from __future__ import annotations

from collections.abc import Callable
import json
from pathlib import Path

from ..models import RabiesMethodSensitivityPanelWorkflowReport
from .shared import _relative_bundle_path


def _write_report_manifest(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
    sha256: Callable[[Path], str],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    linked_files = {
        key: value for key, value in bundle_paths.items() if value.is_file()
    }
    payload = {
        "dataset_id": report.dataset.dataset_id,
        "report_kind": "rabies_method_sensitivity_html_report",
        "variant_count": len(report.variant_runs),
        "parallel_workers": report.parallel_workers,
        "execution_mode": report.execution_mode,
        "stable_clade_count": len(report.stable_clade_rows),
        "changed_clade_count": len(report.changed_clade_rows),
        "linked_artifact_count": len(linked_files),
        "linked_artifacts": {
            key: {
                "path": _relative_bundle_path(path, value),
                "byte_count": value.stat().st_size,
                "sha256": sha256(value),
            }
            for key, value in linked_files.items()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path
