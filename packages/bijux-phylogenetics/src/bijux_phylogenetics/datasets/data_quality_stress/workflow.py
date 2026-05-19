from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from .bundle import write_catarrhine_data_quality_stress_panel_workflow_bundle
from .cleanup import build_catarrhine_data_quality_stress_panel_workflow_report
from .export import export_catarrhine_data_quality_stress_panel_dataset
from bijux_phylogenetics.core.alignment import (
    AlignmentRecord,
    CodingSequencePreparationReport,
    FastaInputValidationReport,
    FastaRepairReport,
    SequenceCompositionOutlier,
    SequenceLengthOutlier,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.pruning import drop_tree_taxa
from bijux_phylogenetics.core.traits import (
    TraitLinkageReport,
    TraitValidationReport,
    link_tree_to_traits,
    validate_traits_table,
)
from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.diagnostics.validation import (
    TreeInspectionReport,
    TreeValidationReport,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError
from bijux_phylogenetics.io.fasta import (
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.cleaning import detect_composition_outlier_sequences
from bijux_phylogenetics.io.fasta.coding import prepare_coding_sequences_for_alignment
from bijux_phylogenetics.io.fasta.records import (
    detect_sequence_length_outliers,
    repair_fasta_input,
    validate_fasta_input,
)
from bijux_phylogenetics.io.newick import write_newick
from .panel import (
    load_catarrhine_data_quality_stress_panel_dataset as load_packaged_catarrhine_data_quality_stress_panel_dataset,
)
from .traits import (
    detect_missing_traits,
    load_permissive_trait_rows,
    resolve_duplicate_traits,
    selected_trait_rows,
)

_DATASET_ID = "catarrhine_data_quality_stress_panel"
_DATASET_LABEL = "Catarrhine data quality stress panel"
_RAW_ALIGNMENT_NAME = "alignment.fasta"
_RAW_SEQUENCE_INPUT_NAME = "sequence-input.fasta"
_RAW_CODING_SEQUENCE_NAME = "coding-sequences.fasta"
_RAW_TREE_NAME = "tree.nwk"
_RAW_TRAITS_NAME = "traits.csv"
_RAW_TRAIT_MISMATCH_NAME = "traits-mismatch.csv"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_REQUIRED_TRAITS = ("body_mass_g", "gestation_days")
_TREE_BRANCH_FLOOR = 1e-6


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelDataset:
    """Packaged stress panel with explicit sequence, tree, and trait defects."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    raw_alignment_path: Path
    raw_tree_path: Path
    raw_traits_path: Path
    raw_sequence_input_path: Path
    raw_coding_sequences_path: Path
    raw_trait_mismatch_path: Path
    reference_output_root: Path
    taxon_count: int
    raw_trait_row_count: int
    required_traits: tuple[str, ...]
    sequence_type: str
    source_summary: str


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelExportResult:
    """Materialized copy of the packaged stress panel."""

    output_root: Path
    readme_path: Path
    raw_alignment_path: Path
    raw_tree_path: Path
    raw_traits_path: Path
    raw_sequence_input_path: Path
    raw_coding_sequences_path: Path
    raw_trait_mismatch_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class TraitDuplicateResolution:
    """Deterministic duplicate-taxon choice for one raw trait surface."""

    taxon: str
    occurrence_count: int
    selected_row_number: int
    selected_non_missing_field_count: int
    discarded_row_numbers: list[int]
    selected_reason: str


@dataclass(slots=True)
class TraitMissingObservation:
    """One missing trait observation in the raw traits table."""

    taxon: str
    row_number: int
    trait: str
    required_for_analysis: bool
    action: str


@dataclass(slots=True)
class DataQualityRepairAction:
    """One explicit action taken to move from raw inputs to the cleaned subset."""

    action_kind: str
    affected_taxa: list[str]
    affected_nodes: list[str]
    reason: str
    result: str


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit and cleanup report for the packaged stress panel."""

    dataset: CatarrhineDataQualityStressPanelDataset
    raw_sequence_input_validation: FastaInputValidationReport
    raw_sequence_input_repair: FastaRepairReport
    raw_sequence_length_outliers: list[SequenceLengthOutlier]
    repaired_sequence_input_validation: FastaInputValidationReport
    coding_sequence_preparation: CodingSequencePreparationReport
    raw_trait_mismatch_linkage: TraitLinkageReport
    raw_trait_mismatch_error: str | None
    raw_alignment_validation: FastaInputValidationReport
    sequence_outliers: list[SequenceCompositionOutlier]
    raw_tree_inspection: TreeInspectionReport
    raw_tree_validation: TreeValidationReport
    trait_duplicates: list[TraitDuplicateResolution]
    missing_traits: list[TraitMissingObservation]
    cleaned_trait_validation: TraitValidationReport
    cleaned_tree_validation: TreeValidationReport
    cleaned_linkage: TraitLinkageReport
    cleaned_alignment_validation: FastaInputValidationReport
    cleaned_alignment_records: list[AlignmentRecord]
    repaired_sequence_input_path: Path
    prepared_coding_sequences_path: Path
    cleaned_tree_path: Path
    cleaned_traits_path: Path
    cleaned_alignment_path: Path
    cleaned_taxa: list[str]
    dropped_taxa: list[str]
    repair_actions: list[DataQualityRepairAction]
    repaired_branch_nodes: list[str]


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelWorkflowBundle:
    """Written workflow outputs for the packaged stress panel."""

    output_root: Path
    raw_taxon_count: int
    cleaned_taxon_count: int
    duplicate_sequence_identifier_count: int
    illegal_character_count: int
    empty_sequence_count: int
    raw_sequence_length_outlier_count: int
    duplicate_trait_taxon_count: int
    missing_trait_value_count: int
    sequence_outlier_count: int
    tree_zero_length_branch_count: int
    tree_negative_branch_count: int
    tree_long_branch_outlier_count: int
    coding_frame_error_count: int
    coding_internal_stop_count: int
    raw_trait_missing_from_traits_count: int
    raw_trait_extra_taxon_count: int
    dropped_taxon_count: int
    repaired_branch_count: int
    workflow_summary_path: Path
    raw_sequence_findings_path: Path
    raw_sequence_repair_path: Path
    repaired_sequence_input_path: Path
    repaired_sequence_validation_path: Path
    coding_sequence_exclusions_path: Path
    prepared_coding_sequences_path: Path
    raw_trait_linkage_path: Path
    trait_duplicates_path: Path
    trait_missing_values_path: Path
    sequence_outliers_path: Path
    tree_issues_path: Path
    repair_actions_path: Path
    cleaned_traits_path: Path
    cleaned_alignment_path: Path
    cleaned_tree_path: Path
    cleaned_linkage_path: Path
    cleaned_validation_path: Path


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelDemoResult:
    """Dataset export plus written workflow bundle for the public demo."""

    output_root: Path
    dataset: CatarrhineDataQualityStressPanelDataset
    dataset_export: CatarrhineDataQualityStressPanelExportResult
    workflow_bundle: CatarrhineDataQualityStressPanelWorkflowBundle
    overview_path: Path


def load_catarrhine_data_quality_stress_panel_dataset() -> (
    CatarrhineDataQualityStressPanelDataset
):
    """Expose the packaged catarrhine stress panel as a first-class dataset surface."""
    return load_packaged_catarrhine_data_quality_stress_panel_dataset()


def run_catarrhine_data_quality_stress_panel_workflow(
    out_dir: Path,
) -> CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit the raw stress panel and write one cleaned comparative subset."""
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    return build_catarrhine_data_quality_stress_panel_workflow_report(dataset, out_dir)


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
def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "stress"
        / _DATASET_ID
    )
