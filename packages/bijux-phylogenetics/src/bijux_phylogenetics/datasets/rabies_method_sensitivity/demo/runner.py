from __future__ import annotations

from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from ..bundle import write_rabies_method_sensitivity_panel_workflow_bundle
from ..config import (
    export_rabies_method_sensitivity_panel_dataset,
    load_rabies_method_sensitivity_panel_dataset,
)
from ..models import RabiesMethodSensitivityPanelDemoResult
from ..workflow import run_rabies_method_sensitivity_panel_workflow
from .overview import _write_overview


def run_rabies_method_sensitivity_panel_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
    parallel_workers: int | None = None,
    variant_ids: tuple[str, ...] | None = None,
) -> RabiesMethodSensitivityPanelDemoResult:
    """Materialize the packaged dataset and rerun the governed sensitivity workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_rabies_method_sensitivity_panel_dataset()
    dataset_export = export_rabies_method_sensitivity_panel_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="rabies-method-sensitivity-") as temporary_root:
        workflow_report = run_rabies_method_sensitivity_panel_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
            parallel_workers=parallel_workers,
            variant_ids=variant_ids,
        )
        workflow_bundle = write_rabies_method_sensitivity_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return RabiesMethodSensitivityPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )
