from __future__ import annotations

from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from .bundle import write_pleistocene_bear_cytb_fragment_workflow_bundle
from .export import export_pleistocene_bear_cytb_fragment_dataset
from .models import (
    PleistoceneBearCytbFragmentDataset,
    PleistoceneBearCytbFragmentDemoResult,
    PleistoceneBearCytbFragmentWorkflowBundle,
)
from .panel import BOOTSTRAP_REPLICATES, IQTREE_SEED, IQTREE_THREADS
from .workflow import run_pleistocene_bear_cytb_fragment_workflow


def run_pleistocene_bear_cytb_fragment_demo(
    output_root: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    iqtree_seed: int = IQTREE_SEED,
    iqtree_threads: int = IQTREE_THREADS,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
) -> PleistoceneBearCytbFragmentDemoResult:
    """Materialize the packaged degraded bear dataset and rerun the governed workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = export_pleistocene_bear_cytb_fragment_dataset(output_root / "dataset")
    with TemporaryDirectory(prefix="pleistocene-bear-workflow-") as temporary_root:
        workflow_report = run_pleistocene_bear_cytb_fragment_workflow(
            Path(temporary_root),
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_pleistocene_bear_cytb_fragment_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md",
        workflow_report.dataset,
        workflow_bundle,
    )
    return PleistoceneBearCytbFragmentDemoResult(
        output_root=output_root,
        dataset=workflow_report.dataset,
        dataset_export=dataset,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_overview(
    path: Path,
    dataset: PleistoceneBearCytbFragmentDataset,
    workflow_bundle: PleistoceneBearCytbFragmentWorkflowBundle,
) -> Path:
    lines = [
        "# Pleistocene Bear CYTB Fragment Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- degraded sequence count: `{len(dataset.degraded_sequence_ids)}`",
        f"- selected model: `{workflow_bundle.selected_model}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.summary_path.name}`",
        f"- missingness effects: `{workflow_bundle.missingness_effects_path.name}`",
        f"- final supported tree: `{workflow_bundle.tree_path.name}`",
        f"- cleaned alignment: `{workflow_bundle.cleaned_alignment_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
