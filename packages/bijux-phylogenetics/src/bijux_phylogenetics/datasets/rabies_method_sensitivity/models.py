from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology import TreeComparisonReport
from bijux_phylogenetics.phylo.topology import TreeRootingReport
from bijux_phylogenetics.engines.inference import (
    InferenceComparisonWorkflowReport,
)
from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport

__all__ = [
    "RabiesMethodSensitivityCladeRow",
    "RabiesMethodSensitivityConclusionRow",
    "RabiesMethodSensitivityPanelDataset",
    "RabiesMethodSensitivityPanelDemoResult",
    "RabiesMethodSensitivityPanelExportResult",
    "RabiesMethodSensitivityPanelWorkflowBundle",
    "RabiesMethodSensitivityPanelWorkflowReport",
    "RabiesMethodSensitivityPreprocessingComparisonRow",
    "RabiesMethodSensitivityTaskRecord",
    "RabiesMethodSensitivityVariant",
    "RabiesMethodSensitivityVariantRun",
]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariant:
    """One declared method combination in the packaged sensitivity matrix."""

    variant_id: str
    label: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDataset:
    """Packaged rabies dataset for method-sensitivity workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_type: str
    workflow_prefix: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    parallel_workers: int
    source_accessions: tuple[str, ...]
    variants: tuple[RabiesMethodSensitivityVariant, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelExportResult:
    """Materialized copy of the packaged rabies method-sensitivity dataset."""

    output_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class RabiesMethodSensitivityVariantRun:
    """One executed variant in the governed method-sensitivity matrix."""

    config: RabiesMethodSensitivityVariant
    alignment_workflow: EngineWorkflowReport
    trimming_workflow: EngineWorkflowReport
    inference_comparison: InferenceComparisonWorkflowReport
    rooted_fasttree_path: Path
    rooted_iqtree_path: Path
    fasttree_rooting: TreeRootingReport
    iqtree_rooting: TreeRootingReport
    rooted_engine_comparison: TreeComparisonReport
    rooted_engine_comparison_table_path: Path
    alignment_length: int
    trimmed_alignment_length: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityTaskRecord:
    """One isolated variant execution record within the workflow batch."""

    variant_id: str
    label: str
    execution_mode: str
    status: str
    log_path: Path
    output_root: Path
    error_code: str | None
    error_message: str | None


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityPreprocessingComparisonRow:
    """One rooted IQ-TREE comparison across two preprocessing variants."""

    left_variant_id: str
    right_variant_id: str
    comparison_axis: str
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    same_taxa_different_rooting: bool


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityCladeRow:
    """One aggregated stable or changed clade-level conclusion across variants."""

    split_id: str
    conclusion_class: str
    evidence_class: str
    occurrence_count: int
    variant_count: int
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityConclusionRow:
    """One high-level biological or analytical conclusion from the workflow."""

    conclusion_id: str
    method_axis: str
    stability_status: str
    claim: str
    evidence: str
    caution: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelWorkflowReport:
    """One governed method-sensitivity workflow run over the packaged rabies panel."""

    dataset: RabiesMethodSensitivityPanelDataset
    execution_record_path: Path
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    parallel_workers: int
    execution_mode: str
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...]
    variant_runs: tuple[RabiesMethodSensitivityVariantRun, ...]
    preprocessing_comparison_rows: tuple[
        RabiesMethodSensitivityPreprocessingComparisonRow, ...
    ]
    stable_clade_rows: tuple[RabiesMethodSensitivityCladeRow, ...]
    changed_clade_rows: tuple[RabiesMethodSensitivityCladeRow, ...]
    conclusion_rows: tuple[RabiesMethodSensitivityConclusionRow, ...]


@dataclass(slots=True)
class RabiesMethodSensitivityPanelWorkflowBundle:
    """Written reviewer-facing outputs for the rabies method-sensitivity workflow."""

    output_root: Path
    variant_count: int
    stable_clade_count: int
    changed_clade_count: int
    preprocessing_change_pair_count: int
    rooted_engine_change_variant_count: int
    serious_conflict_variant_count: int
    execution_record_path: Path
    parallel_workers: int
    execution_mode: str
    workflow_summary_path: Path
    variant_summary_path: Path
    parallel_summary_path: Path
    preprocessing_comparison_path: Path
    stable_clades_path: Path
    changed_clades_path: Path
    conclusion_summary_path: Path
    config_path: Path
    manifest_path: Path
    report_manifest_path: Path
    slurm_job_plan_path: Path
    slurm_assumptions_path: Path
    slurm_summary_path: Path
    slurm_array_partitions_path: Path
    slurm_array_members_path: Path
    slurm_array_strategy_path: Path
    slurm_array_scripts_root: Path
    slurm_job_evidence_root: Path
    slurm_job_evidence_index_path: Path
    slurm_job_evidence_summary_path: Path
    slurm_storage_categories_path: Path
    slurm_storage_variants_path: Path
    slurm_storage_summary_path: Path
    slurm_storage_report_path: Path
    slurm_output_explosion_checks_path: Path
    slurm_output_explosion_variants_path: Path
    slurm_output_explosion_summary_path: Path
    slurm_output_explosion_report_path: Path
    slurm_tree_retention_checks_path: Path
    slurm_tree_retention_files_path: Path
    slurm_tree_retention_summary_path: Path
    slurm_tree_retention_report_path: Path
    slurm_merge_checks_path: Path
    slurm_merge_variants_path: Path
    slurm_merge_summary_path: Path
    slurm_merge_report_path: Path
    slurm_job_count: int
    slurm_total_estimated_core_hours: float
    slurm_maximum_estimated_memory_mib: int
    slurm_maximum_estimated_wallclock_minutes: int
    slurm_total_estimated_scratch_mib: int
    slurm_total_estimated_output_mib: int
    slurm_array_partition_count: int
    slurm_array_script_count: int
    slurm_array_largest_partition_size: int
    slurm_job_evidence_file_count: int
    slurm_job_evidence_total_runtime_seconds: float
    slurm_job_evidence_total_output_byte_count: int
    slurm_storage_total_estimated_mib: int
    slurm_storage_output_byte_count: int
    slurm_storage_log_byte_count: int
    slurm_storage_tree_byte_count: int
    slurm_storage_posterior_sample_byte_count: int
    slurm_storage_report_byte_count: int
    slurm_storage_largest_variant_id: str
    slurm_output_explosion_status: str
    slurm_output_explosion_global_issue_count: int
    slurm_output_explosion_warning_variant_count: int
    slurm_output_explosion_high_risk_variant_count: int
    slurm_tree_retention_status: str
    slurm_tree_set_file_count: int
    slurm_tree_posterior_sample_file_count: int
    slurm_tree_thinning_recommended_file_count: int
    slurm_tree_thinning_required_file_count: int
    slurm_tree_compression_recommended_file_count: int
    slurm_tree_compression_required_file_count: int
    slurm_merge_status: str
    slurm_merge_ready: bool
    slurm_mergeable_variant_count: int
    slurm_merge_failed_check_count: int
    slurm_output_freshness_path: Path
    slurm_output_freshness_checks_path: Path
    slurm_output_freshness_summary_path: Path
    slurm_job_status_path: Path
    slurm_partition_status_path: Path
    slurm_workflow_status_path: Path
    slurm_failure_recovery_jobs_path: Path
    slurm_failure_recovery_partitions_path: Path
    slurm_failure_recovery_summary_path: Path
    slurm_failure_recovery_report_path: Path
    slurm_output_freshness_check_count: int
    slurm_output_freshness_failed_check_count: int
    slurm_fresh_output_job_count: int
    slurm_stale_output_job_count: int
    slurm_completed_job_count: int
    slurm_failed_job_count: int
    slurm_pending_job_count: int
    slurm_stale_job_count: int
    slurm_failure_recovery_status: str
    slurm_failure_recovery_rerunnable_job_count: int
    slurm_failure_recovery_blocked_job_count: int
    slurm_failure_recovery_partition_count: int
    reproducibility_checks_path: Path
    reproducibility_variant_audit_path: Path
    reproducibility_audit_path: Path
    reproducibility_passed: bool
    reproducibility_check_count: int
    reproducibility_failed_check_count: int
    reproducibility_failed_variant_count: int
    report_path: Path
    report_linked_artifact_count: int
    report_html_size_bytes: int
    report_linked_artifact_bytes: int
    report_total_output_bytes: int
    task_logs_root: Path
    variants_root: Path


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDemoResult:
    """Dataset export plus workflow outputs for the public method-sensitivity demo."""

    output_root: Path
    dataset: RabiesMethodSensitivityPanelDataset
    dataset_export: RabiesMethodSensitivityPanelExportResult
    workflow_bundle: RabiesMethodSensitivityPanelWorkflowBundle
    overview_path: Path
