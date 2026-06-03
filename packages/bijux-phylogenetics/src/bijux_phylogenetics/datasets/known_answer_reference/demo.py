from __future__ import annotations

from pathlib import Path
import shutil

from .bundle import write_known_answer_reference_workflow_bundle
from .export import export_known_answer_reference_dataset
from .models import (
    KnownAnswerReferenceDemoResult,
    KnownAnswerReferenceWorkflowBundle,
    KnownAnswerReferenceWorkflowReport,
)
from .policy import format_number
from .workflow import run_known_answer_reference_workflow


def run_known_answer_reference_demo(
    output_root: Path,
) -> KnownAnswerReferenceDemoResult:
    """Materialize the packaged simulation dataset and rerun the recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_known_answer_reference_workflow()
    dataset_export = export_known_answer_reference_dataset(output_root / "dataset")
    workflow_bundle = write_known_answer_reference_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_overview(
        output_root / "overview.md", report, workflow_bundle
    )
    return KnownAnswerReferenceDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_overview(
    path: Path,
    report: KnownAnswerReferenceWorkflowReport,
    bundle: KnownAnswerReferenceWorkflowBundle,
) -> Path:
    lines = [
        "# Known-Answer Simulation Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- taxon count: `{report.dataset.taxon_count}`",
        f"- alignment length: `{report.dataset.sequence_length}`",
        f"- distance recovery preserves rooted topology: `{str(bundle.rooted_topology_equal).lower()}`",
        f"- distance recovery preserves unrooted topology: `{str(bundle.same_unrooted_topology).lower()}`",
        f"- continuous internal-node mean absolute error: `{format_number(bundle.continuous_internal_node_mean_absolute_error)}`",
        f"- discrete internal-node accuracy: `{format_number(bundle.discrete_internal_node_accuracy)}`",
        f"- host internal-node accuracy: `{format_number(bundle.host_internal_node_accuracy)}`",
        f"- host event accuracy: `{format_number(bundle.host_event_accuracy)}`",
        f"- geographic internal-node accuracy: `{format_number(bundle.geographic_internal_node_accuracy)}`",
        f"- geographic event accuracy: `{format_number(bundle.geographic_event_accuracy)}`",
        f"- threshold passes: `{bundle.threshold_pass_count}/{bundle.threshold_row_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- recovered distance tree: `{bundle.distance_tree_path.name}`",
        f"- tree recovery ledger: `{bundle.tree_recovery_path.name}`",
        f"- parameter recovery ledger: `{bundle.parameter_recovery_path.name}`",
        f"- OU fit summary: `{bundle.ou_fit_summary_path.name}`",
        f"- continuous node recovery ledger: `{bundle.continuous_node_recovery_path.name}`",
        f"- discrete node recovery ledger: `{bundle.discrete_node_recovery_path.name}`",
        f"- host event recovery ledger: `{bundle.host_event_recovery_path.name}`",
        f"- geographic event recovery ledger: `{bundle.geographic_event_recovery_path.name}`",
        f"- threshold evaluation ledger: `{bundle.threshold_evaluation_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
