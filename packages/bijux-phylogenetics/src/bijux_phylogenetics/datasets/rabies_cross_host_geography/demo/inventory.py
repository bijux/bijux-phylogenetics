from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    RabiesCrossHostGeographyPanelExportResult,
    RabiesCrossHostGeographyPanelWorkflowBundle,
)
from ..shared import _checksum


def _artifact_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if path.name.endswith(".manifest.json"):
        return "manifest"
    if suffix in {".html", ".htm"}:
        return "report"
    if suffix == ".md":
        return "markdown"
    if suffix == ".json":
        return "json"
    if suffix == ".tsv":
        return "table"
    if suffix == ".svg":
        return "figure"
    if suffix == ".log":
        return "log"
    if suffix in {".nwk", ".tree"}:
        return "tree"
    if suffix in {".aln", ".fasta"}:
        return "alignment"
    if suffix == ".csv":
        return "metadata"
    return "artifact"


def _workflow_artifact_section(relative_path: Path) -> str:
    parts = relative_path.parts
    if len(parts) >= 2 and parts[1] in {
        "bootstrap-review",
        "engine-artifacts",
        "biogeography",
        "comparative",
        "conclusion-stability",
    }:
        return parts[1]
    return "workflow"


def _relative_to_package_root(package_root: Path, path: Path) -> str:
    return path.relative_to(package_root).as_posix()


def _package_inventory_rows(
    *,
    output_root: Path,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    overview_path: Path,
    overview_html_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    dataset_paths = [
        dataset_export.readme_path,
        dataset_export.workflow_config_path,
        dataset_export.sequences_path,
        dataset_export.metadata_path,
        dataset_export.centroids_path,
        dataset_export.accession_table_path,
    ]
    for path in dataset_paths:
        rows.append(
            {
                "section": "dataset",
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    workflow_paths = sorted(
        path for path in workflow_bundle.output_root.rglob("*") if path.is_file()
    )
    for path in workflow_paths:
        rows.append(
            {
                "section": _workflow_artifact_section(
                    path.relative_to(workflow_bundle.output_root)
                ),
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    for path in (overview_path, overview_html_path):
        rows.append(
            {
                "section": "package",
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    return rows


def _write_package_artifact_inventory(
    path: Path,
    *,
    output_root: Path,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    overview_path: Path,
    overview_html_path: Path,
) -> tuple[Path, list[dict[str, str]]]:
    rows = _package_inventory_rows(
        output_root=output_root,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
    )
    return (
        write_taxon_rows(
            path,
            columns=["section", "kind", "relative_path", "sha256", "size_bytes"],
            rows=rows,
        ),
        rows,
    )
