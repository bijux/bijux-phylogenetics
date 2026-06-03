# ruff: noqa: F401
from __future__ import annotations

from pathlib import Path

from .rabies_cross_host_geography import (
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
from .rabies_cross_host_geography import (
    export_rabies_cross_host_geography_panel_dataset as _export_dataset,
)
from .rabies_cross_host_geography import (
    load_rabies_cross_host_geography_panel_dataset as _load_dataset,
)
from .rabies_cross_host_geography import (
    run_rabies_cross_host_geography_panel_workflow as _run_workflow,
)
from .rabies_cross_host_geography import (
    write_rabies_cross_host_geography_panel_workflow_bundle as _write_workflow_bundle,
)
from .rabies_cross_host_geography.demo.builder import (
    _materialize_rabies_cross_host_geography_panel_demo,
)
from .rabies_cross_host_geography.demo.inventory import (
    _write_package_artifact_inventory,
)
from .rabies_cross_host_geography.demo.manifest import (
    _write_demo_package_manifest,
)
from .rabies_cross_host_geography.demo.overview import (
    _build_flagship_answer_summary,
    _write_overview,
)
from .rabies_cross_host_geography.demo.presentation import (
    _write_demo_overview_html,
)
from .rabies_cross_host_geography.demo.reproducibility import (
    _write_package_reproducibility_checklist,
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
    return _materialize_rabies_cross_host_geography_panel_demo(
        output_root,
        config_path=config_path,
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=iqtree_seed,
        iqtree_threads=iqtree_threads,
        bootstrap_replicates=bootstrap_replicates,
        load_dataset=load_rabies_cross_host_geography_panel_dataset,
        export_dataset=export_rabies_cross_host_geography_panel_dataset,
        run_workflow=run_rabies_cross_host_geography_panel_workflow,
        write_workflow_bundle=write_rabies_cross_host_geography_panel_workflow_bundle,
    )
