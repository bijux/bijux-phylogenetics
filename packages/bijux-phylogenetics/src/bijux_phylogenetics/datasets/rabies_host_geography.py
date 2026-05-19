# ruff: noqa: F401
from __future__ import annotations

import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from .rabies_cross_host_geography.bundle import (
    write_rabies_cross_host_geography_panel_workflow_bundle as _write_workflow_bundle,
)
from .rabies_cross_host_geography.config import (
    export_rabies_cross_host_geography_panel_dataset as _export_dataset,
    load_rabies_cross_host_geography_panel_dataset as _load_dataset,
)
from .rabies_cross_host_geography.models import (
    RabiesComparativeBranchRepair,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelDemoResult,
    RabiesCrossHostGeographyPanelExportResult,
    RabiesCrossHostGeographyPanelWorkflowBundle,
    RabiesCrossHostGeographyPanelWorkflowConfig,
    RabiesCrossHostGeographyPanelWorkflowReport,
    RabiesScientificFindingRow,
    RabiesWorkflowConfigAuditRow,
)
from .rabies_cross_host_geography.package import (
    _build_flagship_answer_summary,
    _write_demo_overview_html,
    _write_demo_package_manifest,
    _write_overview,
    _write_package_artifact_inventory,
    _write_package_reproducibility_checklist,
)
from .rabies_cross_host_geography.workflow import (
    run_rabies_cross_host_geography_panel_workflow as _run_workflow,
)

load_rabies_cross_host_geography_panel_dataset = _load_dataset
export_rabies_cross_host_geography_panel_dataset = _export_dataset
run_rabies_cross_host_geography_panel_workflow = _run_workflow
write_rabies_cross_host_geography_panel_workflow_bundle = _write_workflow_bundle


def run_rabies_cross_host_geography_panel_demo(
    output_root: Path,
    *,
    config_path: Path | None = None,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
) -> RabiesCrossHostGeographyPanelDemoResult:
    """Materialize the packaged integrated rabies dataset and rerun the full workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_rabies_cross_host_geography_panel_dataset(config_path)
    dataset_export = export_rabies_cross_host_geography_panel_dataset(
        output_root / "dataset",
        config_path=config_path,
    )
    with TemporaryDirectory(prefix="rabies-cross-host-geography-") as temporary_root:
        workflow_report = run_rabies_cross_host_geography_panel_workflow(
            Path(temporary_root),
            config_path=config_path,
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    short_answer = _build_flagship_answer_summary(workflow_bundle)
    artifact_inventory_path = output_root / "rabies-cross-host-geography-artifacts.tsv"
    reproducibility_checklist_path = (
        output_root / "rabies-cross-host-geography-reproducibility-checklist.tsv"
    )
    overview_path = _write_overview(
        output_root / "overview.md",
        dataset=dataset,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
    )
    overview_html_path = _write_demo_overview_html(
        output_root / "rabies-cross-host-geography-overview.html",
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
    )
    artifact_inventory_path, artifact_inventory_rows = _write_package_artifact_inventory(
        artifact_inventory_path,
        output_root=output_root,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
    )
    reproducibility_checklist_path, checklist_rows = (
        _write_package_reproducibility_checklist(
            reproducibility_checklist_path,
            workflow_bundle=workflow_bundle,
            inventory_rows=artifact_inventory_rows,
            artifact_inventory_path=artifact_inventory_path,
        )
    )
    package_manifest_path = _write_demo_package_manifest(
        output_root / "rabies-cross-host-geography-package.manifest.json",
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        artifact_inventory_rows=artifact_inventory_rows,
        reproducibility_checklist_path=reproducibility_checklist_path,
        checklist_rows=checklist_rows,
    )
    return RabiesCrossHostGeographyPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
        package_manifest_path=package_manifest_path,
    )
