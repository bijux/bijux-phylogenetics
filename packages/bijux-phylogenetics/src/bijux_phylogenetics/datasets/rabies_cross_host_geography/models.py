from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
    from bijux_phylogenetics.ancestral.tree_set import (
        DiscreteAncestralTreeSetReport,
    )
    from bijux_phylogenetics.biogeography import (
        BiogeographyReportPackageResult,
    )
    from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
        PGLSCategoricalContrastReport,
    )
    from bijux_phylogenetics.comparative.pgls.posterior_tree import (
        PosteriorTreePGLSReport,
    )
    from bijux_phylogenetics.comparative.reporting import ComparativeMethodReport
    from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
        RabiesMethodSensitivityPanelWorkflowReport,
    )
    from bijux_phylogenetics.diagnostics.conclusion_stability import (
        ConclusionStabilityReport,
    )
    from bijux_phylogenetics.ecology import HostSwitchingReport
    from bijux_phylogenetics.engines.inference import FastaToTreeWorkflowReport
    from bijux_phylogenetics.phylo.alignment import (
        AlignmentQualityReport,
        SequenceQualityRankingReport,
    )
    from bijux_phylogenetics.phylo.topology import TreeRootingReport

_DATASET_ID = "rabies_cross_host_geography_panel"
_DATASET_LABEL = "Rabies cross-host geography panel"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_PREFIX = "rabies-cross-host-geography-panel"
_HOST_TRAIT = "host_group"
_GEOGRAPHY_TRAIT = "region_group"
_HOST_MODEL = "ard"
_GEOGRAPHY_MODEL = "ard"
_IQTREE_SEED = 1
_IQTREE_THREADS = 1
_BOOTSTRAP_REPLICATES = 1000
_WORKFLOW_TIMEOUT_SECONDS = 300.0
_MAX_BOOTSTRAP_TREE_COUNT = 1500
_MAX_REPORT_TABLE_ROWS = 25
_MEMORY_WARNING_THRESHOLD_BYTES = 67108864
_OUTGROUP_TAXA = ("bat_chile_rv108",)
_SOURCE_ACCESSIONS = (
    "MG458305",
    "MG458304",
    "PV641713",
    "PX845689",
    "OQ693985",
    "PX845683",
    "PX845681",
    "PX845678",
    "PX845676",
)
_WORKFLOW_CONFIG_NAME = "workflow-config.json"
_CLADE_METADATA_COLUMNS = ("host_species", "host_group", "country", "region_group")
_COMPARATIVE_FORMULA = "region_longitude ~ host_group"
_COMPARATIVE_RESPONSE = "region_longitude"
_COMPARATIVE_BRANCH_LENGTH_FLOOR = 1e-6
_BOOTSTRAP_CONSENSUS_THRESHOLD = 0.5
_BOOTSTRAP_ROBUST_SUPPORT_THRESHOLD = 0.9
_ALIGNMENT_MODE = "auto"
_TRIMMING_MODE = "gap-threshold"
_TRIM_GAP_THRESHOLD = 0.1
_FLAGSHIP_QUESTION = (
    "Do the host-associated rabies lineages in this compact panel occupy one "
    "distinct geographic regime while retaining one coherent phylogenetic signal?"
)


@dataclass(frozen=True, slots=True)
class RabiesCrossHostGeographyPanelWorkflowConfig:
    """One config-driven scientific workflow definition for the packaged rabies panel."""

    config_path: Path
    dataset_id: str
    label: str
    sequences_path: Path
    metadata_path: Path
    centroids_path: Path
    sequence_type: str
    workflow_prefix: str
    host_trait: str
    geography_trait: str
    host_model: str
    geography_model: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    timeout_seconds: float | None
    max_bootstrap_tree_count: int | None
    max_report_table_rows: int | None
    memory_warning_threshold_bytes: int | None
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float
    bootstrap_consensus_threshold: float
    bootstrap_robust_support_threshold: float
    clade_metadata_columns: tuple[str, ...]
    comparative_formula: str
    comparative_response: str
    comparative_branch_length_floor: float


@dataclass(frozen=True, slots=True)
class RabiesComparativeBranchRepair:
    """One explicit branch-length adjustment needed before comparative fitting."""

    node_label: str
    original_branch_length: float
    repaired_branch_length: float
    reason: str


@dataclass(slots=True)
class RabiesWorkflowConfigAuditRow:
    """One governed validation check for the packaged rabies workflow config."""

    check_id: str
    status: str
    observed_value: str
    detail: str


@dataclass(slots=True)
class RabiesScientificFindingRow:
    """One reviewer-facing scientific finding from the rabies workflow bundle."""

    finding_id: str
    question: str
    claim: str
    evidence: str
    caution: str
    source_artifact: str


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelDataset:
    """Packaged rabies panel for one complete host and geography workflow."""

    dataset_id: str
    label: str
    dataset_root: Path
    workflow_config_path: Path
    sequences_path: Path
    metadata_path: Path
    centroids_path: Path
    reference_output_root: Path
    sequence_count: int
    sequence_type: str
    workflow_prefix: str
    host_trait: str
    geography_trait: str
    host_model: str
    geography_model: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    timeout_seconds: float | None
    max_bootstrap_tree_count: int | None
    max_report_table_rows: int | None
    memory_warning_threshold_bytes: int | None
    outgroup_taxa: tuple[str, ...]
    observed_host_group_count: int
    observed_region_group_count: int
    clade_metadata_columns: tuple[str, ...]
    comparative_formula: str
    comparative_response: str
    comparative_branch_length_floor: float
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelExportResult:
    """Materialized copy of the packaged rabies integrated dataset."""

    output_root: Path
    readme_path: Path
    workflow_config_path: Path
    sequences_path: Path
    metadata_path: Path
    centroids_path: Path
    accession_table_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelWorkflowReport:
    """One full raw-sequence-to-result workflow run over the packaged rabies panel."""

    dataset: RabiesCrossHostGeographyPanelDataset
    config: RabiesCrossHostGeographyPanelWorkflowConfig
    config_audit_rows: list[RabiesWorkflowConfigAuditRow]
    fasta_to_tree: FastaToTreeWorkflowReport
    rooted_tree_path: Path
    rooting_report: TreeRootingReport
    host_switching: HostSwitchingReport
    biogeography_report: BiogeographyReportPackageResult
    aligned_quality: AlignmentQualityReport
    trimmed_quality: AlignmentQualityReport
    trimmed_sequence_ranking: SequenceQualityRankingReport
    comparative_traits_rows: list[dict[str, str]]
    comparative_tree_path: Path
    comparative_branch_repairs: list[RabiesComparativeBranchRepair]
    comparative_report: ComparativeMethodReport
    comparative_categorical_contrasts: PGLSCategoricalContrastReport
    method_sensitivity_report: RabiesMethodSensitivityPanelWorkflowReport
    rooted_bootstrap_tree_set_path: Path
    comparative_bootstrap_tree_set_path: Path
    host_ancestral_report: DiscreteAncestralReport
    host_ancestral_tree_set_report: DiscreteAncestralTreeSetReport
    geography_ancestral_report: DiscreteAncestralReport
    geography_ancestral_tree_set_report: DiscreteAncestralTreeSetReport
    comparative_posterior_tree_report: PosteriorTreePGLSReport
    conclusion_stability_report: ConclusionStabilityReport


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelWorkflowBundle:
    """Written integrated workflow outputs for the packaged rabies panel."""

    output_root: Path
    selected_model: str
    sequence_type: str
    inferred_sequence_type: str
    input_repair_applied: bool
    aligned_quality_score: float
    trimmed_quality_score: float
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    clade_row_count: int
    bootstrap_tree_count: int
    bootstrap_topology_count: int
    bootstrap_unstable_branch_count: int
    bootstrap_consensus_rooted_rf_distance: int
    bootstrap_consensus_same_unrooted_topology: bool
    bootstrap_consensus_high_support_conflict_count: int
    bootstrap_consensus_branch_score_distance: float | None
    rooted_outgroup_taxa: tuple[str, ...]
    root_host: str
    root_host_confidence: float
    host_switch_count: int
    certain_host_switch_count: int
    uncertain_host_switch_count: int
    root_region: str
    root_region_probability: float
    changed_region_branch_count: int
    migration_event_count: int
    strongly_supported_migration_event_count: int
    comparative_selected_model: str
    comparative_response: str
    comparative_formula: str
    comparative_pgls_lambda: float
    comparative_pgls_r_squared: float
    comparative_branch_repair_count: int
    conclusion_stable_count: int
    conclusion_weak_count: int
    conclusion_unstable_count: int
    timeout_seconds: float | None
    max_bootstrap_tree_count: int | None
    max_report_table_rows: int | None
    memory_warning_threshold_bytes: int | None
    workflow_runtime_seconds: float
    bootstrap_review_runtime_seconds: float
    bootstrap_review_peak_memory_bytes: int
    budget_warning_count: int
    config_check_count: int
    scientific_finding_count: int
    workflow_summary_path: Path
    resource_observations_path: Path
    config_audit_path: Path
    resolved_config_path: Path
    input_validation_path: Path
    alignment_quality_path: Path
    alignment_sequence_ranking_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    rooting_report_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path
    clade_table_path: Path
    bootstrap_output_root: Path
    bootstrap_summary_path: Path
    bootstrap_consensus_tree_path: Path
    bootstrap_clade_frequencies_path: Path
    bootstrap_unstable_branches_path: Path
    bootstrap_unstable_clades_path: Path
    bootstrap_distance_matrix_path: Path
    bootstrap_topology_clusters_path: Path
    bootstrap_tree_comparison_summary_path: Path
    bootstrap_tree_comparison_table_path: Path
    bootstrap_tree_comparison_report_path: Path
    host_switch_summary_path: Path
    host_state_nodes_path: Path
    host_switch_branches_path: Path
    host_switch_counts_path: Path
    host_switch_fits_path: Path
    host_switch_unsupported_path: Path
    host_switch_exclusions_path: Path
    biogeography_output_root: Path
    biogeography_report_path: Path
    biogeography_tree_figure_path: Path
    biogeography_map_path: Path
    comparative_traits_path: Path
    comparative_tree_path: Path
    comparative_repairs_path: Path
    comparative_output_root: Path
    comparative_report_path: Path
    comparative_summary_path: Path
    comparative_coefficients_path: Path
    comparative_residuals_path: Path
    comparative_signal_path: Path
    comparative_model_comparison_path: Path
    comparative_interpretation_path: Path
    comparative_audit_path: Path
    comparative_contrasts_path: Path
    comparative_model_matrix_path: Path
    comparative_categorical_contrasts_path: Path
    comparative_lambda_profile_path: Path
    comparative_manifest_path: Path
    conclusion_stability_output_root: Path
    conclusion_stability_summary_path: Path
    key_clade_stability_path: Path
    support_value_stability_path: Path
    ancestral_state_stability_path: Path
    comparative_coefficient_stability_path: Path
    conclusion_stability_report_path: Path
    scientific_findings_path: Path
    final_report_path: Path
    final_manifest_path: Path


@dataclass(slots=True)
class RabiesCrossHostGeographyPanelDemoResult:
    """Dataset export plus integrated workflow outputs for the public rabies demo."""

    output_root: Path
    dataset: RabiesCrossHostGeographyPanelDataset
    dataset_export: RabiesCrossHostGeographyPanelExportResult
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle
    overview_path: Path
    overview_html_path: Path
    artifact_inventory_path: Path
    reproducibility_checklist_path: Path
    package_manifest_path: Path
