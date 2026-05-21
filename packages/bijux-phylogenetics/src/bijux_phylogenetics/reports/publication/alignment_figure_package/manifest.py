from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentForensicReport,
    AlignmentSummary,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist

from .contracts import AlignmentFigureAudit


def checksum(path: Path) -> str:
    """Return the stable checksum for one artifact path."""
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(payload: object) -> object:
    """Normalize dataclass-rich payloads into JSON-safe structures."""
    return json.loads(json.dumps(payload, default=str))


def build_pre_review_manifest(
    *,
    alignment_path: Path,
    artifact_paths: list[Path],
    reproducibility_manifest_path: Path,
    maximum_site_bins: int,
    window_size: int,
    step_size: int,
    panel_row_limit: int,
    summary: AlignmentSummary,
    forensic: AlignmentForensicReport,
    audit: AlignmentFigureAudit,
) -> dict[str, object]:
    """Build the manifest payload used before the reviewer checklist is written."""
    existing_artifact_paths = artifact_paths[:-1]
    return {
        "report_kind": "alignment_quality_figure_package",
        "input_path": str(alignment_path),
        "input_checksum": checksum(alignment_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {
            str(path): checksum(path) for path in existing_artifact_paths
        },
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "settings": {
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
        },
        "metrics": {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "quality_score": forensic.quality.quality_score,
            "publication_ready": audit.publication_ready,
            "heatmap_row_count": audit.heatmap_row_count,
            "heatmap_bin_count": audit.heatmap_bin_count,
            "plotted_window_count": audit.plotted_window_count,
            "plotted_sequence_count": audit.plotted_sequence_count,
        },
        "alignment_summary": json_ready(asdict(summary)),
        "alignment_quality": json_ready(asdict(forensic.quality)),
        "alignment_readiness": json_ready(asdict(forensic.readiness)),
        "alignment_low_information": json_ready(asdict(forensic.low_information)),
        "audit": json_ready(asdict(audit)),
    }


def build_machine_manifest(
    *,
    alignment_path: Path,
    artifact_paths: list[Path],
    reproducibility_manifest_path: Path,
    reproducibility_manifest: dict[str, object],
    maximum_site_bins: int,
    window_size: int,
    step_size: int,
    panel_row_limit: int,
    summary: AlignmentSummary,
    forensic: AlignmentForensicReport,
    audit: AlignmentFigureAudit,
) -> dict[str, object]:
    """Build the final machine manifest before the checklist is appended."""
    return {
        "report_kind": "alignment_quality_figure_package",
        "input_path": str(alignment_path),
        "input_checksum": checksum(alignment_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": checksum(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
        },
        "metrics": {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "quality_score": forensic.quality.quality_score,
            "publication_ready": audit.publication_ready,
            "heatmap_row_count": audit.heatmap_row_count,
            "heatmap_bin_count": audit.heatmap_bin_count,
            "plotted_window_count": audit.plotted_window_count,
            "plotted_sequence_count": audit.plotted_sequence_count,
        },
        "alignment_summary": json_ready(asdict(summary)),
        "alignment_quality": json_ready(asdict(forensic.quality)),
        "alignment_readiness": json_ready(asdict(forensic.readiness)),
        "alignment_low_information": json_ready(asdict(forensic.low_information)),
        "audit": json_ready(asdict(audit)),
    }


def attach_reviewer_audit_checklist(
    *,
    machine_manifest: dict[str, object],
    reviewer_audit_checklist_path: Path,
    reviewer_audit_checklist: ReviewerAuditChecklist,
) -> dict[str, object]:
    """Attach reviewer checklist outputs to the final machine manifest."""
    machine_manifest["output_paths"].append(str(reviewer_audit_checklist_path))
    machine_manifest["output_checksums"][str(reviewer_audit_checklist_path)] = checksum(
        reviewer_audit_checklist_path
    )
    machine_manifest["reviewer_audit_checklist_path"] = str(
        reviewer_audit_checklist_path
    )
    machine_manifest["reviewer_audit_checklist"] = json_ready(
        asdict(reviewer_audit_checklist)
    )
    return machine_manifest


def write_machine_manifest(path: Path, machine_manifest: dict[str, object]) -> Path:
    """Write the final machine manifest artifact."""
    path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
