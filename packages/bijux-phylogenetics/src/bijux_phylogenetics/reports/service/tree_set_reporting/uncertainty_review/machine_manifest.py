from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from ...ledger import sha256


def build_tree_uncertainty_manifest(
    *,
    title: str,
    tree_set_path: Path,
    out_path: Path,
    artifact_root: Path,
    summary,
    processing,
    budget_report,
    scaled_report_mode: bool,
    methods_summary_result,
    limitations: list[str],
    artifact_paths: dict[str, Path],
    core_sections: list[tuple[str, object]],
    supplemental_sections: list[tuple[str, object]],
) -> tuple[dict[str, object], list[tuple[str, str, str | None]], Path]:
    machine_manifest = {
        "report_kind": "tree-uncertainty",
        "title": title,
        "source_path": str(tree_set_path),
        "input_checksum": sha256(tree_set_path),
        "tree_count": summary.tree_count,
        "rooted_topology_count": summary.rooted_topology_count,
        "processing": asdict(processing),
        "budget": asdict(budget_report),
        "report_mode": "scaled-summary" if scaled_report_mode else "full-review",
        "artifact_root": str(artifact_root),
        "linked_artifact_count": len(artifact_paths) + 1,
        "methods_summary_path": artifact_paths["methods_summary"]
        .relative_to(out_path.parent)
        .as_posix(),
        "methods_summary_warning_count": methods_summary_result.warning_count,
        "limitations": limitations,
        "linked_artifacts": {
            name: {
                "path": path.relative_to(out_path.parent).as_posix(),
                "byte_count": path.stat().st_size,
            }
            for name, path in artifact_paths.items()
        },
        "sections": [name for name, _ in core_sections],
        "supplemental_sections": [name for name, _ in supplemental_sections],
    }
    artifact_links = [
        (
            name.replace("_", "-"),
            path.relative_to(out_path.parent).as_posix(),
            f"{path.stat().st_size} bytes",
        )
        for name, path in artifact_paths.items()
    ]
    artifact_manifest_path = artifact_root / "tree-uncertainty.manifest.json"
    machine_manifest["artifact_manifest_path"] = artifact_manifest_path.relative_to(
        out_path.parent
    ).as_posix()
    machine_manifest["linked_artifacts"]["tree_uncertainty_manifest"] = {
        "path": artifact_manifest_path.relative_to(out_path.parent).as_posix(),
        "byte_count": 0,
    }
    return machine_manifest, artifact_links, artifact_manifest_path
