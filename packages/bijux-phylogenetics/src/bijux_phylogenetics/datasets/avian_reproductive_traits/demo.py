from __future__ import annotations

from pathlib import Path
import shutil

from .bundle import write_avian_reproductive_trait_workflow_bundle
from .export import export_avian_reproductive_trait_dataset
from .models import (
    AvianReproductiveTraitDemoResult,
    AvianReproductiveTraitWorkflowBundle,
)
from .workflow import run_avian_reproductive_trait_workflow


def run_avian_reproductive_trait_demo(
    output_root: Path,
) -> AvianReproductiveTraitDemoResult:
    """Materialize the packaged avian dataset and its governed workflow outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_avian_reproductive_trait_workflow()
    dataset_export = export_avian_reproductive_trait_dataset(output_root / "dataset")
    workflow_bundle = write_avian_reproductive_trait_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_demo_overview(output_root / "overview.md", workflow_bundle)
    return AvianReproductiveTraitDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_demo_overview(
    path: Path,
    bundle: AvianReproductiveTraitWorkflowBundle,
) -> Path:
    lines = [
        "# Avian Reproductive Trait Demo",
        "",
        "This demo materializes the packaged bird dataset and regenerates the",
        "governed comparative workflow outputs that ship with the package.",
        "",
        "## Workflow Outputs",
        "",
        f"- workflow summary: `{bundle.summary_path}`",
        f"- PGLS lambda profile: `{bundle.pgls_lambda_profile_path}`",
        f"- Brownian summary: `{bundle.brownian_summary_path}`",
        f"- OU summary: `{bundle.ou_summary_path}`",
        f"- signal summary: `{bundle.signal_summary_path}`",
        f"- continuous ancestral summary: `{bundle.continuous_ancestral_summary_path}`",
        f"- discrete ancestral summary: `{bundle.discrete_ancestral_summary_path}`",
        f"- clade trait summary: `{bundle.clade_summary_path}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
