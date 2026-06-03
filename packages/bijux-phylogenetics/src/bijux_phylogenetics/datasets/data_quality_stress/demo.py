from __future__ import annotations

from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from .bundle import write_catarrhine_data_quality_stress_panel_workflow_bundle
from .export import export_catarrhine_data_quality_stress_panel_dataset
from .models import (
    CatarrhineDataQualityStressPanelDataset,
    CatarrhineDataQualityStressPanelDemoResult,
    CatarrhineDataQualityStressPanelWorkflowBundle,
)
from .panel import load_catarrhine_data_quality_stress_panel_dataset
from .workflow import run_catarrhine_data_quality_stress_panel_workflow


def run_catarrhine_data_quality_stress_panel_demo(
    output_root: Path,
) -> CatarrhineDataQualityStressPanelDemoResult:
    """Materialize the packaged stress dataset and rerun the governed cleanup workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    dataset_export = export_catarrhine_data_quality_stress_panel_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="catarrhine-data-quality-stress-") as temporary_root:
        workflow_report = run_catarrhine_data_quality_stress_panel_workflow(
            Path(temporary_root)
        )
        workflow_bundle = write_catarrhine_data_quality_stress_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return CatarrhineDataQualityStressPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_overview(
    path: Path,
    dataset: CatarrhineDataQualityStressPanelDataset,
    workflow_bundle: CatarrhineDataQualityStressPanelWorkflowBundle,
) -> Path:
    lines = [
        "# Catarrhine Data Quality Stress Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- raw taxon count: `{workflow_bundle.raw_taxon_count}`",
        f"- cleaned taxon count: `{workflow_bundle.cleaned_taxon_count}`",
        f"- duplicate sequence identifiers: `{workflow_bundle.duplicate_sequence_identifier_count}`",
        f"- illegal FASTA characters: `{workflow_bundle.illegal_character_count}`",
        f"- empty FASTA records: `{workflow_bundle.empty_sequence_count}`",
        f"- raw-sequence length outliers: `{workflow_bundle.raw_sequence_length_outlier_count}`",
        f"- coding frame errors: `{workflow_bundle.coding_frame_error_count}`",
        f"- coding internal stop codons: `{workflow_bundle.coding_internal_stop_count}`",
        f"- duplicate trait taxa: `{workflow_bundle.duplicate_trait_taxon_count}`",
        f"- sequence outliers: `{workflow_bundle.sequence_outlier_count}`",
        f"- repaired branch count: `{workflow_bundle.repaired_branch_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- raw sequence findings: `{workflow_bundle.raw_sequence_findings_path.name}`",
        f"- raw sequence repair ledger: `{workflow_bundle.raw_sequence_repair_path.name}`",
        f"- repaired sequence input: `{workflow_bundle.repaired_sequence_input_path.name}`",
        f"- repaired sequence validation: `{workflow_bundle.repaired_sequence_validation_path.name}`",
        f"- coding sequence exclusions: `{workflow_bundle.coding_sequence_exclusions_path.name}`",
        f"- prepared coding sequences: `{workflow_bundle.prepared_coding_sequences_path.name}`",
        f"- raw trait linkage mismatch: `{workflow_bundle.raw_trait_linkage_path.name}`",
        f"- trait duplicates: `{workflow_bundle.trait_duplicates_path.name}`",
        f"- trait missing values: `{workflow_bundle.trait_missing_values_path.name}`",
        f"- sequence outliers: `{workflow_bundle.sequence_outliers_path.name}`",
        f"- tree issues: `{workflow_bundle.tree_issues_path.name}`",
        f"- repair actions: `{workflow_bundle.repair_actions_path.name}`",
        f"- cleaned traits: `{workflow_bundle.cleaned_traits_path.name}`",
        f"- cleaned alignment: `{workflow_bundle.cleaned_alignment_path.name}`",
        f"- cleaned tree: `{workflow_bundle.cleaned_tree_path.name}`",
        f"- cleaned validation: `{workflow_bundle.cleaned_validation_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
