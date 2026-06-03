from __future__ import annotations

from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from .bundle import write_influenza_a_ha_reference_workflow_bundle
from .export import export_influenza_a_ha_reference_dataset
from .models import (
    InfluenzaAHAReferenceDataset,
    InfluenzaAHAReferenceDemoResult,
    InfluenzaAHAReferenceWorkflowBundle,
)
from .panel import (
    BOOTSTRAP_REPLICATES,
    IQTREE_SEED,
    IQTREE_THREADS,
    load_influenza_a_ha_reference_dataset,
)
from .workflow import run_influenza_a_ha_reference_workflow


def run_influenza_a_ha_reference_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = IQTREE_SEED,
    iqtree_threads: int = IQTREE_THREADS,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> InfluenzaAHAReferenceDemoResult:
    """Materialize the packaged viral dataset and rerun the governed workflow bundle."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_influenza_a_ha_reference_dataset()
    dataset_export = export_influenza_a_ha_reference_dataset(output_root / "dataset")
    with TemporaryDirectory(prefix="influenza-a-ha-reference-") as temporary_root:
        workflow_report = run_influenza_a_ha_reference_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_influenza_a_ha_reference_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md",
        dataset,
        workflow_bundle,
    )
    return InfluenzaAHAReferenceDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_overview(
    path: Path,
    dataset: InfluenzaAHAReferenceDataset,
    workflow_bundle: InfluenzaAHAReferenceWorkflowBundle,
) -> Path:
    lines = [
        "# Influenza A HA Reference Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- sequence type: `{dataset.sequence_type}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.summary_path.name}`",
        f"- final supported tree: `{workflow_bundle.tree_path.name}`",
        f"- selected model table: `{workflow_bundle.model_table_path.name}`",
        f"- branch support table: `{workflow_bundle.support_table_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
