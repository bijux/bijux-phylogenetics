from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from Bio import Phylo

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.tree_set import (
    DiscreteAncestralTreeSetReport,
    summarize_discrete_ancestral_tree_set,
)
from bijux_phylogenetics.biogeography import (
    BiogeographyReportPackageResult,
    build_biogeography_report_package,
)
from bijux_phylogenetics.trees import (
    CladeTableReport,
    CladeTableRow,
    extract_tree_clades,
    write_clade_table,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSResult,
    run_pgls,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls_categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls_lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.posterior_tree_pgls import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.report_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeCoefficientTableRow,
    ComparativeInterpretationRow,
    ComparativeResidualTableRow,
    ComparativeSignalTableRow,
    summarize_comparative_analysis,
    summarize_comparative_audit,
    summarize_comparative_coefficients,
    summarize_comparative_interpretation,
    summarize_comparative_residuals,
    summarize_comparative_signal,
    write_comparative_audit_table,
    write_comparative_coefficient_table,
    write_comparative_contrast_table,
    write_comparative_interpretation_table,
    write_comparative_model_comparison_table,
    write_comparative_residual_table,
    write_comparative_signal_table,
    write_comparative_summary_table,
)
from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    build_comparative_method_report,
)
from bijux_phylogenetics.compare.reports import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.core.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    RabiesMethodSensitivityPanelWorkflowReport,
    run_rabies_method_sensitivity_panel_workflow,
)
from bijux_phylogenetics.diagnostics.conclusion_stability import (
    ConclusionStabilityReport,
    build_ancestral_state_stability_rows,
    build_comparative_coefficient_stability_rows,
    build_conclusion_stability_report,
    build_key_clade_stability_rows,
    build_support_value_stability_rows,
    write_ancestral_state_stability_table,
    write_comparative_coefficient_stability_table,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
    write_key_clade_stability_table,
    write_support_value_stability_table,
)
from bijux_phylogenetics.engines.fasta_to_tree import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.host_association import (
    HostSwitchingReport,
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
)
from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.fasta._shared import load_permissive_fasta_records
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import validate_fasta_input
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.trees import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    compute_clade_frequency_table,
    write_bootstrap_tree_set_artifacts,
)

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


def load_rabies_cross_host_geography_panel_dataset(
    config_path: Path | None = None,
) -> RabiesCrossHostGeographyPanelDataset:
    """Expose the packaged rabies host-and-geography panel as one owned surface."""
    resolved_config = _load_workflow_config(config_path)
    _raise_for_failed_config_audit(_build_workflow_config_audit_rows(resolved_config))
    dataset_root = resolved_config.config_path.parent
    validation = validate_fasta_input(
        resolved_config.sequences_path,
        sequence_type=resolved_config.sequence_type,
    )
    observed_host_groups, observed_region_groups = _read_observed_groups(
        resolved_config.metadata_path,
        host_trait=resolved_config.host_trait,
        geography_trait=resolved_config.geography_trait,
    )
    return RabiesCrossHostGeographyPanelDataset(
        dataset_id=resolved_config.dataset_id,
        label=resolved_config.label,
        dataset_root=dataset_root,
        workflow_config_path=resolved_config.config_path,
        sequences_path=resolved_config.sequences_path,
        metadata_path=resolved_config.metadata_path,
        centroids_path=resolved_config.centroids_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=resolved_config.sequence_type,
        workflow_prefix=resolved_config.workflow_prefix,
        host_trait=resolved_config.host_trait,
        geography_trait=resolved_config.geography_trait,
        host_model=resolved_config.host_model,
        geography_model=resolved_config.geography_model,
        iqtree_seed=resolved_config.iqtree_seed,
        iqtree_threads=resolved_config.iqtree_threads,
        bootstrap_replicates=resolved_config.bootstrap_replicates,
        timeout_seconds=resolved_config.timeout_seconds,
        max_bootstrap_tree_count=resolved_config.max_bootstrap_tree_count,
        max_report_table_rows=resolved_config.max_report_table_rows,
        memory_warning_threshold_bytes=resolved_config.memory_warning_threshold_bytes,
        outgroup_taxa=resolved_config.outgroup_taxa,
        observed_host_group_count=len(observed_host_groups),
        observed_region_group_count=len(observed_region_groups),
        clade_metadata_columns=resolved_config.clade_metadata_columns,
        comparative_formula=resolved_config.comparative_formula,
        comparative_response=resolved_config.comparative_response,
        comparative_branch_length_floor=resolved_config.comparative_branch_length_floor,
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Real rabies virus nucleoprotein sequences paired with grouped host "
            "and macroregion metadata so one governed workflow can rerun tree "
            "inference, host switching, geography review, bootstrap topology "
            "summary, clade extraction, and one comparative model from raw "
            "sequence inputs."
        ),
    )


def export_rabies_cross_host_geography_panel_dataset(
    destination: Path,
    *,
    config_path: Path | None = None,
) -> RabiesCrossHostGeographyPanelExportResult:
    """Copy the packaged integrated rabies dataset and stable expected outputs."""
    dataset = load_rabies_cross_host_geography_panel_dataset(config_path)
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md", destination / "README.md"
    )
    workflow_config_path = shutil.copy2(
        dataset.workflow_config_path, destination / _WORKFLOW_CONFIG_NAME
    )
    sequences_path = shutil.copy2(
        dataset.sequences_path, destination / "sequences.fasta"
    )
    metadata_path = shutil.copy2(dataset.metadata_path, destination / "metadata.csv")
    centroids_path = shutil.copy2(
        dataset.centroids_path, destination / "region-centroids.csv"
    )
    accession_table_path = _write_source_accession_table(
        destination / "source-accessions.tsv",
        dataset=dataset,
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return RabiesCrossHostGeographyPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        workflow_config_path=Path(workflow_config_path),
        sequences_path=Path(sequences_path),
        metadata_path=Path(metadata_path),
        centroids_path=Path(centroids_path),
        accession_table_path=accession_table_path,
        expected_output_root=expected_output_root,
    )


def _write_source_accession_table(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
) -> Path:
    metadata_rows = _read_metadata_rows(dataset.metadata_path)
    accession_index = {
        str(row["accession"]): row for row in metadata_rows if row.get("accession")
    }
    ordered_rows = []
    for accession in dataset.source_accessions:
        row = accession_index[accession]
        ordered_rows.append(
            {
                "accession": accession,
                "accession_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "taxon": str(row["taxon"]),
                "isolate": str(row["isolate"]),
                "host_species": str(row["host_species"]),
                "host_group": str(row["host_group"]),
                "country": str(row["country"]),
                "region_group": str(row["region_group"]),
                "collection_date": str(row.get("collection_date", "")),
            }
        )
    return write_taxon_rows(
        path,
        columns=[
            "accession",
            "accession_url",
            "taxon",
            "isolate",
            "host_species",
            "host_group",
            "country",
            "region_group",
            "collection_date",
        ],
        rows=ordered_rows,
    )


def _read_metadata_rows(metadata_path: Path) -> list[dict[str, str]]:
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def run_rabies_cross_host_geography_panel_workflow(
    out_dir: Path,
    *,
    config_path: Path | None = None,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
) -> RabiesCrossHostGeographyPanelWorkflowReport:
    """Run the full integrated rabies workflow from sequences and metadata."""
    dataset = load_rabies_cross_host_geography_panel_dataset(config_path)
    config = _load_workflow_config(config_path)
    config_audit_rows = _build_workflow_config_audit_rows(config)
    _raise_for_failed_config_audit(config_audit_rows)
    workflow = run_fasta_to_tree_workflow(
        dataset.sequences_path,
        out_dir=out_dir,
        prefix=dataset.workflow_prefix,
        sequence_type=dataset.sequence_type,
        mafft_executable=mafft_executable,
        alignment_mode=config.alignment_mode,
        trimal_executable=trimal_executable,
        trimming_mode=config.trimming_mode,
        iqtree_executable=iqtree_executable,
        iqtree_seed=config.iqtree_seed if iqtree_seed is None else iqtree_seed,
        iqtree_threads=(
            config.iqtree_threads if iqtree_threads is None else iqtree_threads
        ),
        trim_gap_threshold=config.trim_gap_threshold,
        bootstrap_replicates=(
            config.bootstrap_replicates
            if bootstrap_replicates is None
            else bootstrap_replicates
        ),
        timeout_seconds=config.timeout_seconds,
    )
    rooted_tree, rooting_report = root_tree_on_outgroup(
        workflow.output_paths["tree"],
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    rooted_tree_path = out_dir / f"{dataset.workflow_prefix}.rooted.tree"
    write_newick(rooted_tree_path, rooted_tree)
    host_switching = summarize_host_switching(
        rooted_tree_path,
        dataset.metadata_path,
        trait=dataset.host_trait,
        taxon_column="taxon",
        model=dataset.host_model,
    )
    biogeography_report = build_biogeography_report_package(
        tree_path=rooted_tree_path,
        traits_path=dataset.metadata_path,
        centroids_path=dataset.centroids_path,
        trait=dataset.geography_trait,
        out_dir=out_dir / "biogeography-report",
        taxon_column="taxon",
        model=dataset.geography_model,
        region_column="region",
        latitude_column="latitude",
        longitude_column="longitude",
    )
    aligned_quality = build_alignment_quality_report(workflow.output_paths["alignment"])
    trimmed_quality = build_alignment_quality_report(
        workflow.output_paths["trimmed_alignment"]
    )
    trimmed_sequence_ranking = build_sequence_quality_ranking(
        workflow.output_paths["trimmed_alignment"]
    )
    comparative_traits_rows = _build_comparative_trait_rows(
        metadata_path=dataset.metadata_path,
        centroids_path=dataset.centroids_path,
        host_trait=dataset.host_trait,
        geography_trait=dataset.geography_trait,
    )
    comparative_traits_path = (
        out_dir / f"{dataset.workflow_prefix}.comparative-traits.tsv"
    )
    write_taxon_rows(
        comparative_traits_path,
        columns=list(comparative_traits_rows[0].keys()),
        rows=comparative_traits_rows,
    )
    comparative_tree_path, comparative_branch_repairs = _write_comparative_tree(
        rooted_tree_path,
        out_path=out_dir / f"{dataset.workflow_prefix}.comparative.tree",
        branch_length_floor=config.comparative_branch_length_floor,
    )
    comparative_report = build_comparative_method_report(
        comparative_tree_path,
        comparative_traits_path,
        formula=config.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )
    comparative_categorical_contrasts = summarize_pgls_categorical_contrasts(
        comparative_tree_path,
        comparative_traits_path,
        formula=config.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )
    method_sensitivity_report = run_rabies_method_sensitivity_panel_workflow(
        out_dir / "method-sensitivity-review",
        mafft_executable=mafft_executable,
        trimal_executable=trimal_executable,
        iqtree_executable=iqtree_executable,
        fasttree_executable=fasttree_executable,
        iqtree_seed=config.iqtree_seed if iqtree_seed is None else iqtree_seed,
        iqtree_threads=(
            config.iqtree_threads if iqtree_threads is None else iqtree_threads
        ),
        bootstrap_replicates=(
            config.bootstrap_replicates
            if bootstrap_replicates is None
            else bootstrap_replicates
        ),
    )
    rooted_bootstrap_tree_set_path = _write_rooted_tree_set_on_outgroup(
        workflow.bootstrap_workflow.output_paths["bootstrap_trees"],
        out_path=out_dir / f"{dataset.workflow_prefix}.rooted-bootstrap.trees",
        outgroup_taxa=list(dataset.outgroup_taxa),
    )
    comparative_bootstrap_tree_set_path = _write_comparative_tree_set(
        rooted_bootstrap_tree_set_path,
        out_path=out_dir / f"{dataset.workflow_prefix}.comparative-bootstrap.trees",
        reference_tree_path=comparative_tree_path,
        branch_length_floor=config.comparative_branch_length_floor,
    )
    host_ancestral_report = reconstruct_discrete_ancestral_states(
        rooted_tree_path,
        dataset.metadata_path,
        trait=dataset.host_trait,
        taxon_column="taxon",
        model=dataset.host_model,
    )
    host_ancestral_tree_set_report = summarize_discrete_ancestral_tree_set(
        rooted_bootstrap_tree_set_path,
        dataset.metadata_path,
        trait=dataset.host_trait,
        taxon_column="taxon",
        model=_canonicalize_discrete_tree_set_model(dataset.host_model),
    )
    geography_ancestral_report = reconstruct_discrete_ancestral_states(
        rooted_tree_path,
        dataset.metadata_path,
        trait=dataset.geography_trait,
        taxon_column="taxon",
        model=dataset.geography_model,
    )
    geography_ancestral_tree_set_report = summarize_discrete_ancestral_tree_set(
        rooted_bootstrap_tree_set_path,
        dataset.metadata_path,
        trait=dataset.geography_trait,
        taxon_column="taxon",
        model=_canonicalize_discrete_tree_set_model(dataset.geography_model),
    )
    comparative_posterior_tree_report = run_posterior_tree_pgls(
        comparative_bootstrap_tree_set_path,
        comparative_traits_path,
        formula=config.comparative_formula,
        taxon_column="taxon",
        lambda_value="estimate",
    )
    conclusion_stability_report = _build_conclusion_stability_report(
        rooted_tree_path=rooted_tree_path,
        rooted_bootstrap_tree_set_path=rooted_bootstrap_tree_set_path,
        comparative_traits_path=comparative_traits_path,
        comparative_tree_path=comparative_tree_path,
        comparative_report=comparative_report,
        comparative_posterior_tree_report=comparative_posterior_tree_report,
        method_sensitivity_report=method_sensitivity_report,
        metadata_path=dataset.metadata_path,
        host_trait=dataset.host_trait,
        host_model=dataset.host_model,
        geography_trait=dataset.geography_trait,
        geography_model=dataset.geography_model,
        comparative_formula=config.comparative_formula,
        comparative_branch_length_floor=config.comparative_branch_length_floor,
        host_ancestral_report=host_ancestral_report,
        host_ancestral_tree_set_report=host_ancestral_tree_set_report,
        geography_ancestral_report=geography_ancestral_report,
        geography_ancestral_tree_set_report=geography_ancestral_tree_set_report,
    )
    return RabiesCrossHostGeographyPanelWorkflowReport(
        dataset=dataset,
        config=config,
        config_audit_rows=config_audit_rows,
        fasta_to_tree=workflow,
        rooted_tree_path=rooted_tree_path,
        rooting_report=rooting_report,
        host_switching=host_switching,
        biogeography_report=biogeography_report,
        aligned_quality=aligned_quality,
        trimmed_quality=trimmed_quality,
        trimmed_sequence_ranking=trimmed_sequence_ranking,
        comparative_traits_rows=comparative_traits_rows,
        comparative_tree_path=comparative_tree_path,
        comparative_branch_repairs=comparative_branch_repairs,
        comparative_report=comparative_report,
        comparative_categorical_contrasts=comparative_categorical_contrasts,
        method_sensitivity_report=method_sensitivity_report,
        rooted_bootstrap_tree_set_path=rooted_bootstrap_tree_set_path,
        comparative_bootstrap_tree_set_path=comparative_bootstrap_tree_set_path,
        host_ancestral_report=host_ancestral_report,
        host_ancestral_tree_set_report=host_ancestral_tree_set_report,
        geography_ancestral_report=geography_ancestral_report,
        geography_ancestral_tree_set_report=geography_ancestral_tree_set_report,
        comparative_posterior_tree_report=comparative_posterior_tree_report,
        conclusion_stability_report=conclusion_stability_report,
    )


def write_rabies_cross_host_geography_panel_workflow_bundle(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
) -> RabiesCrossHostGeographyPanelWorkflowBundle:
    """Write the complete integrated workflow bundle for the packaged rabies panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow = report.fasta_to_tree
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary

    config_audit_path = _write_workflow_config_audit_table(
        output_root / "workflow-config-audit.tsv",
        report.config_audit_rows,
    )
    resolved_config_path = _write_resolved_workflow_config(
        output_root / "workflow-config.resolved.json",
        report.config,
    )
    input_validation_path = _write_input_validation_table(
        output_root / "input-validation.tsv",
        workflow=workflow,
    )
    alignment_quality_path = _write_alignment_quality_table(
        output_root / "alignment-quality.tsv",
        aligned=report.aligned_quality,
        trimmed=report.trimmed_quality,
    )
    alignment_sequence_ranking_path = _write_sequence_ranking_table(
        output_root / "alignment-sequence-ranking.tsv",
        report.trimmed_sequence_ranking,
    )

    alignment_path = _copy_output(
        workflow.output_paths["alignment"],
        output_root / workflow.output_paths["alignment"].name,
    )
    trimmed_alignment_path = _copy_output(
        workflow.output_paths["trimmed_alignment"],
        output_root / workflow.output_paths["trimmed_alignment"].name,
    )
    tree_path = _copy_output(
        report.rooted_tree_path,
        output_root / report.rooted_tree_path.name,
    )
    stable_rooting_report = replace(
        report.rooting_report, tree_path=Path(tree_path.name)
    )
    rooting_report_path = write_tree_rooting_report(
        output_root / f"{report.dataset.workflow_prefix}.rooting.tsv",
        stable_rooting_report,
    )
    model_table_path = _copy_output(
        workflow.output_paths["model_table"],
        output_root / workflow.output_paths["model_table"].name,
    )
    support_table_path = _copy_output(
        workflow.output_paths["support_table"],
        output_root / workflow.output_paths["support_table"].name,
    )
    log_path = _copy_output(
        workflow.output_paths["log"],
        output_root / workflow.output_paths["log"].name,
    )
    manifest_path = _copy_output(
        workflow.manifest_path,
        output_root / workflow.manifest_path.name,
    )
    engine_artifact_root = (
        output_root / "engine-artifacts" / report.dataset.workflow_prefix
    )
    shutil.copytree(workflow.engine_artifact_dir, engine_artifact_root)

    clade_report = extract_tree_clades(
        report.rooted_tree_path,
        metadata_path=report.dataset.metadata_path,
        taxon_column="taxon",
        metadata_columns=list(report.dataset.clade_metadata_columns),
    )
    stable_clade_report = _stabilize_clade_report(
        clade_report,
        stable_source_path=Path(tree_path.name),
    )
    clade_table_path = write_clade_table(
        output_root / "clade-table.tsv",
        stable_clade_report,
    )

    bootstrap_output_root = output_root / "bootstrap-review"
    bootstrap_artifacts = write_bootstrap_tree_set_artifacts(
        report.fasta_to_tree.bootstrap_workflow.output_paths["bootstrap_trees"],
        out_dir=bootstrap_output_root,
        prefix="bootstrap-review",
        consensus_threshold=report.config.bootstrap_consensus_threshold,
        robust_support_threshold=report.config.bootstrap_robust_support_threshold,
        max_tree_count=report.config.max_bootstrap_tree_count,
        memory_warning_threshold_bytes=report.config.memory_warning_threshold_bytes,
    )
    bootstrap_summary_path = _write_stable_bootstrap_summary_table(
        bootstrap_output_root / "bootstrap-review.summary.tsv",
        bootstrap_artifacts.summary_report,
    )
    bootstrap_tree_comparison_report = build_tree_comparison_report(
        tree_path,
        bootstrap_artifacts.output_paths["consensus_tree"],
        out_path=bootstrap_output_root
        / "rooted-tree-vs-bootstrap-consensus.report.html",
    )
    _stabilize_bundle_report_paths(
        bootstrap_tree_comparison_report.output_path,
        output_root=output_root,
    )
    bootstrap_tree_comparison_table_path = write_tree_comparison_table(
        bootstrap_output_root / "rooted-tree-vs-bootstrap-consensus.comparison.tsv",
        tree_path,
        bootstrap_artifacts.output_paths["consensus_tree"],
    )
    bootstrap_tree_comparison_summary_path = _write_bootstrap_tree_comparison_summary(
        bootstrap_output_root / "rooted-tree-vs-bootstrap-consensus.summary.tsv",
        bootstrap_tree_comparison_report,
    )

    host_switch_summary_path = write_host_switch_summary_table(
        output_root / "host-switch-summary.tsv",
        report.host_switching,
    )
    host_state_nodes_path = write_host_state_node_table(
        output_root / "host-state-nodes.tsv",
        report.host_switching,
    )
    host_switch_branches_path = write_host_switch_branch_table(
        output_root / "host-switch-branches.tsv",
        report.host_switching,
    )
    host_switch_counts_path = write_host_switch_count_table(
        output_root / "host-switch-counts.tsv",
        report.host_switching,
    )
    host_switch_fits_path = write_host_switch_fit_table(
        output_root / "host-switch-fits.tsv",
        report.host_switching,
    )
    host_switch_unsupported_path = write_unsupported_host_switch_claim_table(
        output_root / "host-switch-unsupported.tsv",
        report.host_switching,
    )
    host_switch_exclusions_path = write_host_switch_exclusion_table(
        output_root / "host-switch-exclusions.tsv",
        report.host_switching,
    )

    biogeography_output_root = output_root / "biogeography"
    shutil.copytree(report.biogeography_report.output_dir, biogeography_output_root)
    biogeography_report_path = biogeography_output_root / "biogeography-report.html"
    biogeography_tree_figure_path = (
        biogeography_output_root / "ancestral-region-tree.svg"
    )
    biogeography_map_path = biogeography_output_root / "geographic-region-map.html"

    comparative_traits_path = write_taxon_rows(
        output_root / "comparative-traits.tsv",
        columns=list(report.comparative_traits_rows[0].keys()),
        rows=report.comparative_traits_rows,
    )
    comparative_tree_path = _copy_output(
        report.comparative_tree_path,
        output_root / "comparative-tree.nwk",
    )
    comparative_repairs_path = _write_comparative_branch_repairs_table(
        output_root / "comparative-tree-adjustments.tsv",
        report.comparative_branch_repairs,
    )
    comparative_output_root = output_root / "comparative"
    comparative_output_root.mkdir(parents=True, exist_ok=True)
    comparative_report = report.comparative_report
    comparative_summary_row = summarize_comparative_analysis(comparative_report)
    comparative_coefficient_rows = summarize_comparative_coefficients(
        comparative_report
    )
    comparative_residual_rows = summarize_comparative_residuals(comparative_report)
    comparative_signal_row = summarize_comparative_signal(comparative_report)
    comparative_interpretation_rows = summarize_comparative_interpretation(
        comparative_report
    )
    comparative_audit_rows = summarize_comparative_audit(comparative_report)
    comparative_report_path = _write_comparative_report(
        comparative_output_root / "comparative-report.html",
        summary_row=comparative_summary_row,
        coefficient_rows=comparative_coefficient_rows,
        residual_rows=comparative_residual_rows,
        signal_row=comparative_signal_row,
        interpretation_rows=comparative_interpretation_rows,
        branch_repairs=report.comparative_branch_repairs,
    )
    comparative_summary_path = write_comparative_summary_table(
        comparative_output_root / "comparative-summary.tsv",
        comparative_summary_row,
    )
    comparative_coefficients_path = write_comparative_coefficient_table(
        comparative_output_root / "coefficient-table.tsv",
        comparative_coefficient_rows,
    )
    comparative_residuals_path = write_comparative_residual_table(
        comparative_output_root / "residual-summary.tsv",
        comparative_residual_rows,
    )
    comparative_signal_path = write_comparative_signal_table(
        comparative_output_root / "signal-summary.tsv",
        comparative_signal_row,
    )
    comparative_model_comparison_path = write_comparative_model_comparison_table(
        comparative_output_root / "model-comparison.tsv",
        comparative_report,
    )
    comparative_interpretation_path = write_comparative_interpretation_table(
        comparative_output_root / "interpretation-table.tsv",
        comparative_interpretation_rows,
    )
    comparative_audit_path = write_comparative_audit_table(
        comparative_output_root / "audit-table.tsv",
        comparative_audit_rows,
    )
    comparative_contrasts_path = write_comparative_contrast_table(
        comparative_output_root / "contrast-table.tsv",
        comparative_report,
    )
    comparative_model_matrix_path = comparative_output_root / "model-matrix.tsv"
    write_pgls_model_matrix_table(
        comparative_model_matrix_path,
        comparative_report.snapshot.pgls_inputs.model_matrix,
    )
    comparative_categorical_contrasts_path = write_pgls_categorical_contrast_table(
        comparative_output_root / "categorical-contrasts.tsv",
        report.comparative_categorical_contrasts,
    )
    comparative_lambda_profile_path = write_pgls_lambda_profile_table(
        comparative_output_root / "lambda-profile.tsv",
        comparative_report.snapshot.pgls_model.lambda_fit,
    )
    comparative_manifest_path = _write_comparative_manifest(
        comparative_output_root / "comparative.manifest.json",
        comparative_summary_row=comparative_summary_row,
        branch_repairs=report.comparative_branch_repairs,
        output_paths={
            "comparative_report": comparative_report_path,
            "summary_table": comparative_summary_path,
            "coefficient_table": comparative_coefficients_path,
            "residual_table": comparative_residuals_path,
            "signal_table": comparative_signal_path,
            "model_comparison_table": comparative_model_comparison_path,
            "interpretation_table": comparative_interpretation_path,
            "audit_table": comparative_audit_path,
            "contrast_table": comparative_contrasts_path,
            "model_matrix_table": comparative_model_matrix_path,
            "categorical_contrast_table": comparative_categorical_contrasts_path,
            "lambda_profile_table": comparative_lambda_profile_path,
        },
    )
    conclusion_stability_output_root = output_root / "conclusion-stability"
    conclusion_stability_output_root.mkdir(parents=True, exist_ok=True)
    conclusion_stability_summary_path = write_conclusion_stability_summary_table(
        conclusion_stability_output_root / "conclusion-stability-summary.tsv",
        report.conclusion_stability_report,
    )
    key_clade_stability_path = write_key_clade_stability_table(
        conclusion_stability_output_root / "key-clade-stability.tsv",
        report.conclusion_stability_report.key_clade_rows,
    )
    support_value_stability_path = write_support_value_stability_table(
        conclusion_stability_output_root / "support-value-stability.tsv",
        report.conclusion_stability_report.support_value_rows,
    )
    ancestral_state_stability_path = write_ancestral_state_stability_table(
        conclusion_stability_output_root / "ancestral-state-stability.tsv",
        report.conclusion_stability_report.ancestral_state_rows,
    )
    comparative_coefficient_stability_path = (
        write_comparative_coefficient_stability_table(
            conclusion_stability_output_root / "comparative-coefficient-stability.tsv",
            report.conclusion_stability_report.comparative_coefficient_rows,
        )
    )
    conclusion_stability_report_path = write_conclusion_stability_report_html(
        conclusion_stability_output_root / "conclusion-stability-report.html",
        report.conclusion_stability_report,
    )
    scientific_finding_rows = _build_scientific_finding_rows(
        report=report,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        comparative_summary_row=comparative_summary_row,
        comparative_interpretation_rows=comparative_interpretation_rows,
    )
    scientific_findings_path = _write_scientific_findings_table(
        output_root / "scientific-findings.tsv",
        scientific_finding_rows,
    )

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report=report,
        clade_row_count=len(stable_clade_report.rows),
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        comparative_summary_row=comparative_summary_row,
        scientific_finding_count=len(scientific_finding_rows),
    )
    resource_observations_path = _write_resource_observation_table(
        output_root / "resource-observations.tsv",
        report=report,
        bootstrap_artifacts=bootstrap_artifacts,
    )
    final_report_path = _write_integrated_report(
        output_root / "rabies-cross-host-geography-report.html",
        report=report,
        workflow_summary_path=workflow_summary_path,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        clade_row_count=len(stable_clade_report.rows),
        comparative_summary_row=comparative_summary_row,
        comparative_interpretation_rows=comparative_interpretation_rows,
        comparative_branch_repair_count=len(report.comparative_branch_repairs),
        scientific_finding_rows=scientific_finding_rows,
        max_report_table_rows=report.config.max_report_table_rows,
    )
    final_manifest_path = _write_manifest(
        output_root / "rabies-cross-host-geography.manifest.json",
        report=report,
        comparative_summary_row=comparative_summary_row,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        clade_row_count=len(stable_clade_report.rows),
        scientific_finding_count=len(scientific_finding_rows),
        bundle_paths={
            "workflow_summary": workflow_summary_path,
            "config_audit": config_audit_path,
            "resolved_config": resolved_config_path,
            "input_validation": input_validation_path,
            "alignment_quality": alignment_quality_path,
            "alignment_sequence_ranking": alignment_sequence_ranking_path,
            "alignment": alignment_path,
            "trimmed_alignment": trimmed_alignment_path,
            "rooted_tree": tree_path,
            "rooting_report": rooting_report_path,
            "model_table": model_table_path,
            "support_table": support_table_path,
            "clade_table": clade_table_path,
            "bootstrap_summary": bootstrap_summary_path,
            "bootstrap_consensus_tree": bootstrap_artifacts.output_paths[
                "consensus_tree"
            ],
            "bootstrap_clade_frequencies": bootstrap_artifacts.output_paths[
                "clade_frequencies"
            ],
            "bootstrap_unstable_branches": bootstrap_artifacts.output_paths[
                "unstable_branches"
            ],
            "bootstrap_unstable_clades": bootstrap_artifacts.output_paths[
                "unstable_clades"
            ],
            "bootstrap_distance_matrix": bootstrap_artifacts.output_paths[
                "distance_matrix"
            ],
            "bootstrap_topology_clusters": bootstrap_artifacts.output_paths[
                "topology_clusters"
            ],
            "bootstrap_tree_comparison_summary": bootstrap_tree_comparison_summary_path,
            "bootstrap_tree_comparison_table": bootstrap_tree_comparison_table_path,
            "bootstrap_tree_comparison_report": (
                bootstrap_tree_comparison_report.output_path
            ),
            "host_switch_summary": host_switch_summary_path,
            "host_state_nodes": host_state_nodes_path,
            "host_switch_branches": host_switch_branches_path,
            "host_switch_counts": host_switch_counts_path,
            "host_switch_fits": host_switch_fits_path,
            "host_switch_unsupported": host_switch_unsupported_path,
            "host_switch_exclusions": host_switch_exclusions_path,
            "biogeography_report": biogeography_report_path,
            "biogeography_tree_figure": biogeography_tree_figure_path,
            "biogeography_map": biogeography_map_path,
            "comparative_traits": comparative_traits_path,
            "comparative_tree": comparative_tree_path,
            "comparative_repairs": comparative_repairs_path,
            "comparative_report": comparative_report_path,
            "comparative_summary": comparative_summary_path,
            "comparative_coefficients": comparative_coefficients_path,
            "comparative_residuals": comparative_residuals_path,
            "comparative_signal": comparative_signal_path,
            "comparative_model_comparison": comparative_model_comparison_path,
            "comparative_interpretation": comparative_interpretation_path,
            "comparative_audit": comparative_audit_path,
            "comparative_contrasts": comparative_contrasts_path,
            "comparative_model_matrix": comparative_model_matrix_path,
            "comparative_categorical_contrasts": comparative_categorical_contrasts_path,
            "comparative_lambda_profile": comparative_lambda_profile_path,
            "comparative_manifest": comparative_manifest_path,
            "conclusion_stability_summary": conclusion_stability_summary_path,
            "key_clade_stability": key_clade_stability_path,
            "support_value_stability": support_value_stability_path,
            "ancestral_state_stability": ancestral_state_stability_path,
            "comparative_coefficient_stability": (
                comparative_coefficient_stability_path
            ),
            "conclusion_stability_report": conclusion_stability_report_path,
            "scientific_findings": scientific_findings_path,
            "final_report": final_report_path,
        },
    )

    return RabiesCrossHostGeographyPanelWorkflowBundle(
        output_root=output_root,
        selected_model=workflow.selected_model,
        sequence_type=report.dataset.sequence_type,
        inferred_sequence_type=workflow.sequence_type,
        input_repair_applied=workflow.input_repair is not None,
        aligned_quality_score=report.aligned_quality.quality_score,
        trimmed_quality_score=report.trimmed_quality.quality_score,
        minimum_support=workflow.support_summary.minimum_support,
        maximum_support=workflow.support_summary.maximum_support,
        median_support=workflow.support_summary.median_support,
        weakly_supported_clade_count=workflow.support_summary.weakly_supported_clade_count,
        clade_row_count=len(stable_clade_report.rows),
        bootstrap_tree_count=bootstrap_artifacts.summary_report.tree_count,
        bootstrap_topology_count=(
            bootstrap_artifacts.summary_report.diversity.rooted_topology_count
        ),
        bootstrap_unstable_branch_count=(
            bootstrap_artifacts.summary_report.unstable_branch_count
        ),
        bootstrap_consensus_rooted_rf_distance=(
            bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
        ),
        bootstrap_consensus_same_unrooted_topology=(
            bootstrap_tree_comparison_report.topology.same_unrooted_topology
        ),
        bootstrap_consensus_high_support_conflict_count=len(
            [
                row
                for row in bootstrap_tree_comparison_report.support.conflicting_clades
                if row.conflict_classification == "high_support_conflict"
            ]
        ),
        bootstrap_consensus_branch_score_distance=(
            bootstrap_tree_comparison_report.branch_lengths.branch_score.branch_score_distance
        ),
        rooted_outgroup_taxa=tuple(report.rooting_report.rooted_outgroup_taxa),
        root_host=host_summary.root_host,
        root_host_confidence=host_summary.root_confidence,
        host_switch_count=host_summary.host_switch_count,
        certain_host_switch_count=host_summary.certain_host_switch_count,
        uncertain_host_switch_count=host_summary.uncertain_host_switch_count,
        root_region=geography_summary.root_region,
        root_region_probability=geography_summary.root_region_probability,
        changed_region_branch_count=geography_summary.changed_branch_count,
        migration_event_count=migration_summary.event_count,
        strongly_supported_migration_event_count=(
            migration_summary.strongly_supported_event_count
        ),
        comparative_selected_model=comparative_summary_row.selected_model,
        comparative_response=comparative_summary_row.response,
        comparative_formula=comparative_summary_row.formula,
        comparative_pgls_lambda=comparative_summary_row.pgls_lambda,
        comparative_pgls_r_squared=comparative_summary_row.pgls_r_squared,
        comparative_branch_repair_count=len(report.comparative_branch_repairs),
        conclusion_stable_count=report.conclusion_stability_report.summary.stable_count,
        conclusion_weak_count=report.conclusion_stability_report.summary.weak_count,
        conclusion_unstable_count=(
            report.conclusion_stability_report.summary.unstable_count
        ),
        timeout_seconds=report.config.timeout_seconds,
        max_bootstrap_tree_count=report.config.max_bootstrap_tree_count,
        max_report_table_rows=report.config.max_report_table_rows,
        memory_warning_threshold_bytes=report.config.memory_warning_threshold_bytes,
        workflow_runtime_seconds=report.fasta_to_tree.runtime_seconds,
        bootstrap_review_runtime_seconds=(
            bootstrap_artifacts.summary_report.processing.runtime_seconds
        ),
        bootstrap_review_peak_memory_bytes=(
            bootstrap_artifacts.summary_report.processing.peak_memory_bytes
        ),
        budget_warning_count=len(bootstrap_artifacts.budget_report.warning_messages),
        config_check_count=len(report.config_audit_rows),
        scientific_finding_count=len(scientific_finding_rows),
        workflow_summary_path=workflow_summary_path,
        resource_observations_path=resource_observations_path,
        config_audit_path=config_audit_path,
        resolved_config_path=resolved_config_path,
        input_validation_path=input_validation_path,
        alignment_quality_path=alignment_quality_path,
        alignment_sequence_ranking_path=alignment_sequence_ranking_path,
        alignment_path=alignment_path,
        trimmed_alignment_path=trimmed_alignment_path,
        tree_path=tree_path,
        rooting_report_path=rooting_report_path,
        model_table_path=model_table_path,
        support_table_path=support_table_path,
        log_path=log_path,
        manifest_path=manifest_path,
        engine_artifact_root=engine_artifact_root,
        clade_table_path=clade_table_path,
        bootstrap_output_root=bootstrap_output_root,
        bootstrap_summary_path=bootstrap_summary_path,
        bootstrap_consensus_tree_path=bootstrap_artifacts.output_paths[
            "consensus_tree"
        ],
        bootstrap_clade_frequencies_path=bootstrap_artifacts.output_paths[
            "clade_frequencies"
        ],
        bootstrap_unstable_branches_path=bootstrap_artifacts.output_paths[
            "unstable_branches"
        ],
        bootstrap_unstable_clades_path=bootstrap_artifacts.output_paths[
            "unstable_clades"
        ],
        bootstrap_distance_matrix_path=bootstrap_artifacts.output_paths[
            "distance_matrix"
        ],
        bootstrap_topology_clusters_path=bootstrap_artifacts.output_paths[
            "topology_clusters"
        ],
        bootstrap_tree_comparison_summary_path=(bootstrap_tree_comparison_summary_path),
        bootstrap_tree_comparison_table_path=bootstrap_tree_comparison_table_path,
        bootstrap_tree_comparison_report_path=(
            bootstrap_tree_comparison_report.output_path
        ),
        host_switch_summary_path=host_switch_summary_path,
        host_state_nodes_path=host_state_nodes_path,
        host_switch_branches_path=host_switch_branches_path,
        host_switch_counts_path=host_switch_counts_path,
        host_switch_fits_path=host_switch_fits_path,
        host_switch_unsupported_path=host_switch_unsupported_path,
        host_switch_exclusions_path=host_switch_exclusions_path,
        biogeography_output_root=biogeography_output_root,
        biogeography_report_path=biogeography_report_path,
        biogeography_tree_figure_path=biogeography_tree_figure_path,
        biogeography_map_path=biogeography_map_path,
        comparative_traits_path=comparative_traits_path,
        comparative_tree_path=comparative_tree_path,
        comparative_repairs_path=comparative_repairs_path,
        comparative_output_root=comparative_output_root,
        comparative_report_path=comparative_report_path,
        comparative_summary_path=comparative_summary_path,
        comparative_coefficients_path=comparative_coefficients_path,
        comparative_residuals_path=comparative_residuals_path,
        comparative_signal_path=comparative_signal_path,
        comparative_model_comparison_path=comparative_model_comparison_path,
        comparative_interpretation_path=comparative_interpretation_path,
        comparative_audit_path=comparative_audit_path,
        comparative_contrasts_path=comparative_contrasts_path,
        comparative_model_matrix_path=comparative_model_matrix_path,
        comparative_categorical_contrasts_path=comparative_categorical_contrasts_path,
        comparative_lambda_profile_path=comparative_lambda_profile_path,
        comparative_manifest_path=comparative_manifest_path,
        conclusion_stability_output_root=conclusion_stability_output_root,
        conclusion_stability_summary_path=conclusion_stability_summary_path,
        key_clade_stability_path=key_clade_stability_path,
        support_value_stability_path=support_value_stability_path,
        ancestral_state_stability_path=ancestral_state_stability_path,
        comparative_coefficient_stability_path=(comparative_coefficient_stability_path),
        conclusion_stability_report_path=conclusion_stability_report_path,
        scientific_findings_path=scientific_findings_path,
        final_report_path=final_report_path,
        final_manifest_path=final_manifest_path,
    )


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
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_rabies_cross_host_geography_panel_dataset(config_path)
    dataset_export = export_rabies_cross_host_geography_panel_dataset(
        output_root / "dataset",
        config_path=config_path,
    )
    with TemporaryDirectory(prefix="rabies-cross-host-geography-") as temporary_root:
        workflow_report = run_rabies_cross_host_geography_panel_workflow(
            Path(temporary_root),
            config_path=config_path,
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
        workflow_bundle = write_rabies_cross_host_geography_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    short_answer = _build_flagship_answer_summary(workflow_bundle)
    artifact_inventory_path = output_root / "rabies-cross-host-geography-artifacts.tsv"
    reproducibility_checklist_path = (
        output_root / "rabies-cross-host-geography-reproducibility-checklist.tsv"
    )
    overview_path = _write_overview(
        output_root / "overview.md",
        dataset=dataset,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
    )
    overview_html_path = _write_demo_overview_html(
        output_root / "rabies-cross-host-geography-overview.html",
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
    )
    artifact_inventory_path, artifact_inventory_rows = _write_package_artifact_inventory(
        artifact_inventory_path,
        output_root=output_root,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
    )
    reproducibility_checklist_path, checklist_rows = (
        _write_package_reproducibility_checklist(
            reproducibility_checklist_path,
            workflow_bundle=workflow_bundle,
            inventory_rows=artifact_inventory_rows,
            artifact_inventory_path=artifact_inventory_path,
        )
    )
    package_manifest_path = _write_demo_package_manifest(
        output_root / "rabies-cross-host-geography-package.manifest.json",
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        config=workflow_report.config,
        short_answer=short_answer,
        artifact_inventory_path=artifact_inventory_path,
        artifact_inventory_rows=artifact_inventory_rows,
        reproducibility_checklist_path=reproducibility_checklist_path,
        checklist_rows=checklist_rows,
    )
    return RabiesCrossHostGeographyPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
        artifact_inventory_path=artifact_inventory_path,
        reproducibility_checklist_path=reproducibility_checklist_path,
        package_manifest_path=package_manifest_path,
    )


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "pathogens"
        / _DATASET_ID
    )


def _default_workflow_config_path() -> Path:
    return _resource_root() / _WORKFLOW_CONFIG_NAME


def _load_workflow_config(
    config_path: Path | None,
) -> RabiesCrossHostGeographyPanelWorkflowConfig:
    resolved_path = (
        _default_workflow_config_path()
        if config_path is None
        else Path(config_path).expanduser().resolve()
    )
    payload = json.loads(resolved_path.read_text(encoding="utf-8"))
    dataset_root = resolved_path.parent
    dataset_id = payload.get("dataset_id", _DATASET_ID)
    if dataset_id != _DATASET_ID:
        raise ValueError(
            f"workflow config dataset_id must be '{_DATASET_ID}', got '{dataset_id}'"
        )
    return RabiesCrossHostGeographyPanelWorkflowConfig(
        config_path=resolved_path,
        dataset_id=dataset_id,
        label=payload.get("label", _DATASET_LABEL),
        sequences_path=dataset_root / payload.get("sequences_path", "sequences.fasta"),
        metadata_path=dataset_root / payload.get("metadata_path", "metadata.csv"),
        centroids_path=dataset_root
        / payload.get("centroids_path", "region-centroids.csv"),
        sequence_type=payload.get("sequence_type", _SEQUENCE_TYPE),
        workflow_prefix=payload.get("workflow_prefix", _WORKFLOW_PREFIX),
        host_trait=payload.get("host_trait", _HOST_TRAIT),
        geography_trait=payload.get("geography_trait", _GEOGRAPHY_TRAIT),
        host_model=payload.get("host_model", _HOST_MODEL),
        geography_model=payload.get("geography_model", _GEOGRAPHY_MODEL),
        outgroup_taxa=tuple(payload.get("outgroup_taxa", list(_OUTGROUP_TAXA))),
        iqtree_seed=int(payload.get("iqtree_seed", _IQTREE_SEED)),
        iqtree_threads=int(payload.get("iqtree_threads", _IQTREE_THREADS)),
        bootstrap_replicates=int(
            payload.get("bootstrap_replicates", _BOOTSTRAP_REPLICATES)
        ),
        timeout_seconds=(
            None
            if payload.get("timeout_seconds", _WORKFLOW_TIMEOUT_SECONDS) is None
            else float(payload.get("timeout_seconds", _WORKFLOW_TIMEOUT_SECONDS))
        ),
        max_bootstrap_tree_count=(
            None
            if payload.get("max_bootstrap_tree_count", _MAX_BOOTSTRAP_TREE_COUNT)
            is None
            else int(payload.get("max_bootstrap_tree_count", _MAX_BOOTSTRAP_TREE_COUNT))
        ),
        max_report_table_rows=(
            None
            if payload.get("max_report_table_rows", _MAX_REPORT_TABLE_ROWS) is None
            else int(payload.get("max_report_table_rows", _MAX_REPORT_TABLE_ROWS))
        ),
        memory_warning_threshold_bytes=(
            None
            if payload.get(
                "memory_warning_threshold_bytes",
                _MEMORY_WARNING_THRESHOLD_BYTES,
            )
            is None
            else int(
                payload.get(
                    "memory_warning_threshold_bytes",
                    _MEMORY_WARNING_THRESHOLD_BYTES,
                )
            )
        ),
        alignment_mode=payload.get("alignment_mode", _ALIGNMENT_MODE),
        trimming_mode=payload.get("trimming_mode", _TRIMMING_MODE),
        trim_gap_threshold=float(
            payload.get("trim_gap_threshold", _TRIM_GAP_THRESHOLD)
        ),
        bootstrap_consensus_threshold=float(
            payload.get(
                "bootstrap_consensus_threshold",
                _BOOTSTRAP_CONSENSUS_THRESHOLD,
            )
        ),
        bootstrap_robust_support_threshold=float(
            payload.get(
                "bootstrap_robust_support_threshold",
                _BOOTSTRAP_ROBUST_SUPPORT_THRESHOLD,
            )
        ),
        clade_metadata_columns=tuple(
            payload.get("clade_metadata_columns", list(_CLADE_METADATA_COLUMNS))
        ),
        comparative_formula=payload.get("comparative_formula", _COMPARATIVE_FORMULA),
        comparative_response=payload.get("comparative_response", _COMPARATIVE_RESPONSE),
        comparative_branch_length_floor=float(
            payload.get(
                "comparative_branch_length_floor",
                _COMPARATIVE_BRANCH_LENGTH_FLOOR,
            )
        ),
    )


def _build_workflow_config_audit_rows(
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
) -> list[RabiesWorkflowConfigAuditRow]:
    rows: list[RabiesWorkflowConfigAuditRow] = []
    input_files = (
        ("sequences_path", config.sequences_path),
        ("metadata_path", config.metadata_path),
        ("centroids_path", config.centroids_path),
    )
    missing_input_paths: list[Path] = []
    for check_id, path in input_files:
        exists = path.is_file()
        rows.append(
            RabiesWorkflowConfigAuditRow(
                check_id=check_id,
                status="pass" if exists else "fail",
                observed_value=path.name,
                detail="input file is present"
                if exists
                else "configured input file is missing",
            )
        )
        if not exists:
            missing_input_paths.append(path)
    if missing_input_paths:
        return rows

    records = load_permissive_fasta_records(config.sequences_path)
    sequence_ids = sorted(
        {record.identifier.strip() for record in records if record.identifier.strip()}
    )
    sequence_id_set = set(sequence_ids)
    metadata_rows: list[dict[str, str]] = []
    metadata_columns: list[str] = []
    with config.metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        metadata_columns = [] if reader.fieldnames is None else list(reader.fieldnames)
        metadata_rows = list(reader)
    required_metadata_columns = [
        "taxon",
        config.host_trait,
        config.geography_trait,
        *config.clade_metadata_columns,
    ]
    missing_metadata_columns = sorted(
        {
            column
            for column in required_metadata_columns
            if column not in set(metadata_columns)
        }
    )
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="metadata_required_columns",
            status="pass" if not missing_metadata_columns else "fail",
            observed_value=str(
                len(required_metadata_columns) - len(missing_metadata_columns)
            ),
            detail=(
                "metadata exposes the required workflow columns"
                if not missing_metadata_columns
                else "missing metadata columns: " + ", ".join(missing_metadata_columns)
            ),
        )
    )
    if missing_metadata_columns:
        return rows

    metadata_taxa = sorted(
        {row["taxon"].strip() for row in metadata_rows if row["taxon"].strip()}
    )
    metadata_taxon_set = set(metadata_taxa)
    missing_metadata_taxa = sorted(sequence_id_set - metadata_taxon_set)
    missing_sequence_taxa = sorted(metadata_taxon_set - sequence_id_set)
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="taxon_crosswalk",
            status=(
                "pass"
                if not missing_metadata_taxa and not missing_sequence_taxa
                else "fail"
            ),
            observed_value=str(len(metadata_taxa)),
            detail=(
                "metadata taxa match the FASTA identifiers"
                if not missing_metadata_taxa and not missing_sequence_taxa
                else (
                    "sequence-only taxa: "
                    + (", ".join(missing_metadata_taxa) or "none")
                    + "; metadata-only taxa: "
                    + (", ".join(missing_sequence_taxa) or "none")
                )
            ),
        )
    )

    outgroup_taxa = sorted(config.outgroup_taxa)
    missing_outgroup_taxa = sorted(set(outgroup_taxa) - sequence_id_set)
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="outgroup_taxa",
            status="pass" if not missing_outgroup_taxa else "fail",
            observed_value="|".join(outgroup_taxa),
            detail=(
                "all outgroup taxa are present in the FASTA panel"
                if not missing_outgroup_taxa
                else "missing outgroup taxa: " + ", ".join(missing_outgroup_taxa)
            ),
        )
    )

    centroid_rows: list[dict[str, str]] = []
    with config.centroids_path.open("r", encoding="utf-8", newline="") as handle:
        centroid_rows = list(csv.DictReader(handle))
    centroid_region_set = {
        row["region"].strip() for row in centroid_rows if row["region"].strip()
    }
    metadata_region_set = {
        row[config.geography_trait].strip()
        for row in metadata_rows
        if row[config.geography_trait].strip()
    }
    missing_centroid_regions = sorted(metadata_region_set - centroid_region_set)
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="centroid_region_coverage",
            status="pass" if not missing_centroid_regions else "fail",
            observed_value=str(len(metadata_region_set)),
            detail=(
                "each grouped geography state has one centroid row"
                if not missing_centroid_regions
                else "missing centroid rows for: " + ", ".join(missing_centroid_regions)
            ),
        )
    )

    comparative_columns = {
        "taxon",
        "host_group",
        "region_group",
        "region_latitude",
        "region_longitude",
    }
    response_supported = config.comparative_response in comparative_columns
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="comparative_response_column",
            status="pass" if response_supported else "fail",
            observed_value=config.comparative_response,
            detail=(
                "comparative response is present in the derived trait table"
                if response_supported
                else "expected one of: " + ", ".join(sorted(comparative_columns))
            ),
        )
    )
    timeout_valid = config.timeout_seconds is None or config.timeout_seconds > 0.0
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="timeout_seconds",
            status="pass" if timeout_valid else "fail",
            observed_value=(
                ""
                if config.timeout_seconds is None
                else _format_number(config.timeout_seconds)
            ),
            detail=(
                "workflow timeout budget is positive"
                if timeout_valid
                else "timeout_seconds must be greater than zero when configured"
            ),
        )
    )
    max_tree_count_valid = (
        config.max_bootstrap_tree_count is None or config.max_bootstrap_tree_count >= 1
    )
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="max_bootstrap_tree_count",
            status="pass" if max_tree_count_valid else "fail",
            observed_value=(
                ""
                if config.max_bootstrap_tree_count is None
                else str(config.max_bootstrap_tree_count)
            ),
            detail=(
                "bootstrap summary tree budget is positive"
                if max_tree_count_valid
                else "max_bootstrap_tree_count must be at least 1 when configured"
            ),
        )
    )
    max_report_rows_valid = (
        config.max_report_table_rows is None or config.max_report_table_rows >= 1
    )
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="max_report_table_rows",
            status="pass" if max_report_rows_valid else "fail",
            observed_value=(
                ""
                if config.max_report_table_rows is None
                else str(config.max_report_table_rows)
            ),
            detail=(
                "review table row budget is positive"
                if max_report_rows_valid
                else "max_report_table_rows must be at least 1 when configured"
            ),
        )
    )
    memory_threshold_valid = (
        config.memory_warning_threshold_bytes is None
        or config.memory_warning_threshold_bytes >= 1
    )
    rows.append(
        RabiesWorkflowConfigAuditRow(
            check_id="memory_warning_threshold_bytes",
            status="pass" if memory_threshold_valid else "fail",
            observed_value=(
                ""
                if config.memory_warning_threshold_bytes is None
                else str(config.memory_warning_threshold_bytes)
            ),
            detail=(
                "memory warning threshold is positive"
                if memory_threshold_valid
                else "memory_warning_threshold_bytes must be at least 1 when configured"
            ),
        )
    )
    return rows


def _raise_for_failed_config_audit(rows: list[RabiesWorkflowConfigAuditRow]) -> None:
    failures = [row for row in rows if row.status == "fail"]
    if not failures:
        return
    details = "; ".join(f"{row.check_id}: {row.detail}" for row in failures)
    raise ValueError(f"rabies workflow config failed validation: {details}")


def _read_observed_groups(
    metadata_path: Path,
    *,
    host_trait: str,
    geography_trait: str,
) -> tuple[set[str], set[str]]:
    host_groups: set[str] = set()
    region_groups: set[str] = set()
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            host_group = row.get(host_trait, "").strip()
            region_group = row.get(geography_trait, "").strip()
            if host_group:
                host_groups.add(host_group)
            if region_group:
                region_groups.add(region_group)
    return host_groups, region_groups


def _build_comparative_trait_rows(
    *,
    metadata_path: Path,
    centroids_path: Path,
    host_trait: str,
    geography_trait: str,
) -> list[dict[str, str]]:
    centroids_by_region: dict[str, dict[str, str]] = {}
    with centroids_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            centroids_by_region[row["region"].strip()] = row
    rows: list[dict[str, str]] = []
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            region = row[geography_trait].strip()
            centroid = centroids_by_region[region]
            rows.append(
                {
                    "taxon": row["taxon"].strip(),
                    "host_group": row[host_trait].strip(),
                    "region_group": region,
                    "region_latitude": centroid["latitude"].strip(),
                    "region_longitude": centroid["longitude"].strip(),
                }
            )
    return rows


def _write_comparative_tree(
    rooted_tree_path: Path,
    *,
    out_path: Path,
    branch_length_floor: float,
) -> tuple[Path, list[RabiesComparativeBranchRepair]]:
    tree = load_tree(rooted_tree_path)
    repairs = _apply_branch_length_floor(tree.root, floor=branch_length_floor)
    write_newick(out_path, tree)
    return out_path, repairs


def _write_comparative_tree_set(
    rooted_tree_set_path: Path,
    *,
    out_path: Path,
    reference_tree_path: Path,
    branch_length_floor: float,
) -> Path:
    reference_length_lookup = _build_branch_length_lookup(
        load_tree(reference_tree_path)
    )
    rooted_trees = _load_tree_set_trees(rooted_tree_set_path)
    adjusted_trees = []
    for tree in rooted_trees:
        _overlay_branch_lengths_from_reference(
            tree.root,
            reference_length_lookup=reference_length_lookup,
            floor=branch_length_floor,
        )
        _apply_branch_length_floor(tree.root, floor=branch_length_floor)
        adjusted_trees.append(tree)
    return write_tree_set(out_path, adjusted_trees)


def _write_rooted_tree_set_on_outgroup(
    tree_set_path: Path,
    *,
    out_path: Path,
    outgroup_taxa: list[str],
) -> Path:
    rooted_trees = [
        _root_tree_on_outgroup_from_tree(tree, outgroup_taxa=outgroup_taxa)
        for tree in _load_tree_set_trees(tree_set_path)
    ]
    return write_tree_set(out_path, rooted_trees)


def _load_tree_set_trees(path: Path) -> list[PhyloTree]:
    source_format = "newick"
    return [
        tree_from_biophylo(tree, source_format=source_format)
        for tree in Phylo.parse(path, source_format)
    ]


def _root_tree_on_outgroup_from_tree(
    tree: PhyloTree, *, outgroup_taxa: list[str]
) -> PhyloTree:
    biophylo_tree = tree_to_biophylo(tree)
    matched = [
        next(biophylo_tree.find_clades(name=taxon), None) for taxon in outgroup_taxa
    ]
    if not any(clade is not None for clade in matched):
        raise ValueError(
            "none of the requested outgroup taxa were found while rooting a tree set"
        )
    biophylo_tree.root_with_outgroup(*[clade for clade in matched if clade is not None])
    return tree_from_biophylo(biophylo_tree, source_format="newick")


def _build_branch_length_lookup(tree: PhyloTree) -> dict[str, float]:
    lookup: dict[str, float] = {}
    for node in tree.iter_nodes():
        if node is tree.root or node.branch_length is None:
            continue
        lookup[_node_signature(node)] = node.branch_length
    return lookup


def _overlay_branch_lengths_from_reference(
    node: TreeNode,
    *,
    reference_length_lookup: dict[str, float],
    floor: float,
) -> None:
    for child in node.children:
        signature = _node_signature(child)
        child.branch_length = reference_length_lookup.get(signature, floor)
        _overlay_branch_lengths_from_reference(
            child,
            reference_length_lookup=reference_length_lookup,
            floor=floor,
        )


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    return "|".join(sorted(_descendant_taxa(node)))


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return taxa


def _build_conclusion_stability_report(
    *,
    rooted_tree_path: Path,
    rooted_bootstrap_tree_set_path: Path,
    comparative_traits_path: Path,
    comparative_tree_path: Path,
    comparative_report: ComparativeMethodReport,
    comparative_posterior_tree_report: PosteriorTreePGLSReport,
    method_sensitivity_report: RabiesMethodSensitivityPanelWorkflowReport,
    metadata_path: Path,
    host_trait: str,
    host_model: str,
    geography_trait: str,
    geography_model: str,
    comparative_formula: str,
    comparative_branch_length_floor: float,
    host_ancestral_report: DiscreteAncestralReport,
    host_ancestral_tree_set_report: DiscreteAncestralTreeSetReport,
    geography_ancestral_report: DiscreteAncestralReport,
    geography_ancestral_tree_set_report: DiscreteAncestralTreeSetReport,
) -> ConclusionStabilityReport:
    baseline_clades = extract_tree_clades(rooted_tree_path)
    bootstrap_frequencies = compute_clade_frequency_table(
        rooted_bootstrap_tree_set_path
    )
    method_tree_paths = [
        variant.rooted_iqtree_path for variant in method_sensitivity_report.variant_runs
    ] + [
        variant.rooted_fasttree_path
        for variant in method_sensitivity_report.variant_runs
    ]
    method_clade_reports = [extract_tree_clades(path) for path in method_tree_paths]
    key_clade_rows = build_key_clade_stability_rows(
        baseline_clades=baseline_clades,
        bootstrap_frequencies=bootstrap_frequencies,
        method_clade_reports=method_clade_reports,
    )
    support_value_rows = build_support_value_stability_rows(
        baseline_clades=baseline_clades,
        bootstrap_frequencies=bootstrap_frequencies,
        method_clade_reports=method_clade_reports,
    )
    host_method_reports = [
        reconstruct_discrete_ancestral_states(
            path,
            metadata_path,
            trait=host_trait,
            taxon_column="taxon",
            model=host_model,
        )
        for path in method_tree_paths
    ]
    geography_method_reports = [
        reconstruct_discrete_ancestral_states(
            path,
            metadata_path,
            trait=geography_trait,
            taxon_column="taxon",
            model=geography_model,
        )
        for path in method_tree_paths
    ]
    ancestral_state_rows = build_ancestral_state_stability_rows(
        baseline_report=host_ancestral_report,
        bootstrap_report=host_ancestral_tree_set_report,
        method_reports=host_method_reports,
    ) + build_ancestral_state_stability_rows(
        baseline_report=geography_ancestral_report,
        bootstrap_report=geography_ancestral_tree_set_report,
        method_reports=geography_method_reports,
    )
    comparative_method_results = [
        _run_comparative_pgls_on_tree(
            tree_path=path,
            comparative_traits_path=comparative_traits_path,
            formula=comparative_formula,
            branch_length_floor=comparative_branch_length_floor,
        )
        for path in method_tree_paths
    ]
    comparative_coefficient_rows = build_comparative_coefficient_stability_rows(
        baseline_result=comparative_report.snapshot.pgls_model,
        bootstrap_report=comparative_posterior_tree_report,
        method_results=comparative_method_results,
    )
    return build_conclusion_stability_report(
        key_clade_rows=key_clade_rows,
        support_value_rows=support_value_rows,
        ancestral_state_rows=ancestral_state_rows,
        comparative_coefficient_rows=comparative_coefficient_rows,
    )


def _canonicalize_discrete_tree_set_model(model: str) -> str:
    normalized = model.strip().lower()
    alias_map = {
        "er": "equal-rates",
        "equal-rates": "equal-rates",
        "sym": "symmetric",
        "symmetric": "symmetric",
        "ard": "all-rates-different",
        "all-rates-different": "all-rates-different",
        "fitch": "fitch",
    }
    return alias_map.get(normalized, model)


def _run_comparative_pgls_on_tree(
    *,
    tree_path: Path,
    comparative_traits_path: Path,
    formula: str,
    branch_length_floor: float,
) -> PGLSResult:
    with TemporaryDirectory(prefix="bijux-rabies-comparative-tree-") as temporary_root:
        adjusted_tree_path, _repairs = _write_comparative_tree(
            tree_path,
            out_path=Path(temporary_root) / "comparative-tree.nwk",
            branch_length_floor=branch_length_floor,
        )
        return run_pgls(
            adjusted_tree_path,
            comparative_traits_path,
            formula=formula,
            taxon_column="taxon",
            lambda_value="estimate",
        )


def _apply_branch_length_floor(
    root: TreeNode,
    *,
    floor: float,
) -> list[RabiesComparativeBranchRepair]:
    repairs: list[RabiesComparativeBranchRepair] = []

    def visit(node: TreeNode, *, is_root: bool) -> None:
        if not is_root and node.branch_length is not None and node.branch_length <= 0.0:
            repairs.append(
                RabiesComparativeBranchRepair(
                    node_label=node.name or "<internal>",
                    original_branch_length=float(node.branch_length),
                    repaired_branch_length=floor,
                    reason=(
                        "comparative branch-length methods require strictly positive "
                        "nonroot branch lengths"
                    ),
                )
            )
            node.branch_length = floor
        for child in node.children:
            visit(child, is_root=False)

    visit(root, is_root=True)
    return repairs


def _stabilize_clade_report(
    report: CladeTableReport,
    *,
    stable_source_path: Path,
) -> CladeTableReport:
    return CladeTableReport(
        path=stable_source_path,
        source_format=report.source_format,
        tree_count=report.tree_count,
        metadata_path=report.metadata_path,
        taxon_column=report.taxon_column,
        metadata_columns=list(report.metadata_columns),
        rows=[
            CladeTableRow(
                source_path=stable_source_path,
                tree_index=row.tree_index,
                node_kind=row.node_kind,
                clade_id=row.clade_id,
                node_label=row.node_label,
                taxon_count=row.taxon_count,
                taxa=list(row.taxa),
                support=row.support,
                support_fraction=row.support_fraction,
                branch_length=row.branch_length,
                root_depth=row.root_depth,
                descendant_tip_depth_min=row.descendant_tip_depth_min,
                descendant_tip_depth_max=row.descendant_tip_depth_max,
                node_age=row.node_age,
                metadata=list(row.metadata),
            )
            for row in report.rows
        ],
    )


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _write_workflow_config_audit_table(
    path: Path,
    rows: list[RabiesWorkflowConfigAuditRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["check_id", "status", "observed_value", "detail"],
        rows=[
            {
                "check_id": row.check_id,
                "status": row.status,
                "observed_value": row.observed_value,
                "detail": row.detail,
            }
            for row in rows
        ],
    )


def _write_resolved_workflow_config(
    path: Path,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
) -> Path:
    payload = {
        "report_kind": "rabies_cross_host_geography_workflow_config",
        "dataset_id": config.dataset_id,
        "label": config.label,
        "source_config": config.config_path.name,
        "input_files": {
            "sequences_path": {
                "path": config.sequences_path.name,
                "sha256": _checksum(config.sequences_path),
            },
            "metadata_path": {
                "path": config.metadata_path.name,
                "sha256": _checksum(config.metadata_path),
            },
            "centroids_path": {
                "path": config.centroids_path.name,
                "sha256": _checksum(config.centroids_path),
            },
        },
        "workflow": {
            "sequence_type": config.sequence_type,
            "workflow_prefix": config.workflow_prefix,
            "host_trait": config.host_trait,
            "geography_trait": config.geography_trait,
            "host_model": config.host_model,
            "geography_model": config.geography_model,
            "outgroup_taxa": list(config.outgroup_taxa),
            "iqtree_seed": config.iqtree_seed,
            "iqtree_threads": config.iqtree_threads,
            "bootstrap_replicates": config.bootstrap_replicates,
            "timeout_seconds": config.timeout_seconds,
            "max_bootstrap_tree_count": config.max_bootstrap_tree_count,
            "max_report_table_rows": config.max_report_table_rows,
            "memory_warning_threshold_bytes": (config.memory_warning_threshold_bytes),
            "alignment_mode": config.alignment_mode,
            "trimming_mode": config.trimming_mode,
            "trim_gap_threshold": config.trim_gap_threshold,
            "bootstrap_consensus_threshold": config.bootstrap_consensus_threshold,
            "bootstrap_robust_support_threshold": (
                config.bootstrap_robust_support_threshold
            ),
            "clade_metadata_columns": list(config.clade_metadata_columns),
            "comparative_formula": config.comparative_formula,
            "comparative_response": config.comparative_response,
            "comparative_branch_length_floor": (config.comparative_branch_length_floor),
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_bootstrap_tree_comparison_summary(
    path: Path,
    comparison_report: ComparisonReportBuildResult,
) -> Path:
    high_support_conflict_count = len(
        [
            row
            for row in comparison_report.support.conflicting_clades
            if row.conflict_classification == "high_support_conflict"
        ]
    )
    row = {
        "left_tree": comparison_report.topology.left_path.name,
        "right_tree": comparison_report.topology.right_path.name,
        "shared_taxon_count": str(len(comparison_report.topology.shared_taxa)),
        "rooted_rf_distance": str(
            comparison_report.topology.rooted_robinson_foulds_distance
        ),
        "rooted_normalized_rf": _format_number(
            comparison_report.topology.rooted_normalized_robinson_foulds
        ),
        "topology_equal": (
            "true" if comparison_report.topology.topology_equal else "false"
        ),
        "same_unrooted_topology": (
            "true" if comparison_report.topology.same_unrooted_topology else "false"
        ),
        "same_taxa_different_rooting": (
            "true"
            if comparison_report.topology.same_taxa_different_rooting
            else "false"
        ),
        "same_topology_different_branch_lengths": (
            "true"
            if comparison_report.topology.same_topology_different_branch_lengths
            else "false"
        ),
        "support_conflict_count": str(
            len(comparison_report.support.conflicting_clades)
        ),
        "high_support_conflict_count": str(high_support_conflict_count),
        "branch_score_distance": _format_number(
            comparison_report.branch_lengths.branch_score.branch_score_distance
        ),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_stable_bootstrap_summary_table(
    path: Path,
    report: BootstrapTreeSetSummaryReport,
) -> Path:
    row = {
        "tree_count": str(report.tree_count),
        "shared_taxon_count": str(len(report.shared_taxa)),
        "rooted_topology_count": str(report.diversity.rooted_topology_count),
        "dominant_topology_frequency": _format_number(
            report.diversity.dominant_topology_frequency
        ),
        "effective_topology_count": _format_number(
            report.diversity.effective_topology_count
        ),
        "mean_robinson_foulds_distance": _format_number(
            report.diversity.mean_robinson_foulds_distance
        ),
        "mean_normalized_robinson_foulds_distance": _format_number(
            report.diversity.mean_normalized_robinson_foulds_distance
        ),
        "consensus_threshold": _format_number(report.consensus_threshold),
        "robust_support_threshold": _format_number(report.robust_support_threshold),
        "unstable_branch_count": str(report.unstable_branch_count),
        "warning_count": str(len(report.warnings)),
        "consensus_newick": report.consensus.consensus_newick,
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _stabilize_bundle_report_paths(path: Path, *, output_root: Path) -> Path:
    text = path.read_text(encoding="utf-8")
    normalized_root = output_root.as_posix().rstrip("/") + "/"
    path.write_text(text.replace(normalized_root, ""), encoding="utf-8")
    return path


def _build_scientific_finding_rows(
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
) -> list[RabiesScientificFindingRow]:
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_question = (
        "Does the bootstrap consensus preserve the rooted ML conclusion?"
    )
    if bootstrap_tree_comparison_report.topology.topology_equal:
        bootstrap_claim = "The bootstrap consensus preserves the rooted ML topology on the shared taxon set."
    else:
        bootstrap_claim = "The bootstrap consensus differs from the rooted ML topology after support-driven summarization."
    comparative_claim = next(
        (
            row.claim
            for row in comparative_interpretation_rows
            if row.topic == "coefficient" and row.claim
        ),
        "The comparative layer did not expose one stable host-associated longitude shift.",
    )
    return [
        RabiesScientificFindingRow(
            finding_id="root_host_state",
            question="What host state anchors the rooted rabies panel?",
            claim=f"The rooted tree places the ancestral host state in {host_summary.root_host}.",
            evidence=(
                f"root host confidence {_format_number(host_summary.root_confidence)} "
                f"with outgroup {','.join(report.dataset.outgroup_taxa)}"
            ),
            caution=(
                "The panel is compact and grouped by broad host classes rather than species-level host states."
            ),
            source_artifact=report.dataset.workflow_prefix + ".rooting.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="root_region_state",
            question="What geographic regime anchors the rooted rabies panel?",
            claim=(
                f"The rooted tree places the ancestral region in {geography_summary.root_region}."
            ),
            evidence=(
                f"root region probability {_format_number(geography_summary.root_region_probability)} "
                f"across {geography_summary.changed_branch_count} changed branches"
            ),
            caution=(
                "Grouped macroregions simplify the raw locality labels so the result should be treated as regional rather than site-level history."
            ),
            source_artifact="biogeography/summary.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="host_switching",
            question="How much host-switching signal appears in the rooted tree?",
            claim=(
                f"The host reconstruction inferred {host_summary.host_switch_count} host-switch branch changes."
            ),
            evidence=(
                f"certain changes {host_summary.certain_host_switch_count}; "
                f"uncertain changes {host_summary.uncertain_host_switch_count}"
            ),
            caution=(
                "Branch-wise host changes depend on the grouped host coding and should not be over-read as one exhaustive host-jump catalogue."
            ),
            source_artifact="host-switch-summary.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="bootstrap_consensus",
            question=bootstrap_question,
            claim=bootstrap_claim,
            evidence=(
                f"rooted RF distance {bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance}; "
                f"high-support conflicts "
                f"{len([row for row in bootstrap_tree_comparison_report.support.conflicting_clades if row.conflict_classification == 'high_support_conflict'])}"
            ),
            caution=(
                "Consensus trees can collapse low-support branches, so exact rooted agreement is stricter than shared major clades."
            ),
            source_artifact=(
                "bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv"
            ),
        ),
        RabiesScientificFindingRow(
            finding_id="comparative_longitude",
            question=(
                "Do host-associated rabies lineages occupy one distinct longitudinal regime in this panel?"
            ),
            claim=comparative_claim,
            evidence=(
                f"selected model {comparative_summary_row.selected_model}; "
                f"PGLS lambda {_format_number(comparative_summary_row.pgls_lambda)}; "
                f"r-squared {_format_number(comparative_summary_row.pgls_r_squared)}"
            ),
            caution=(
                "The comparative claim is associational, uses only nine taxa, and retains residual-diagnostic cautions."
            ),
            source_artifact="comparative/interpretation-table.tsv",
        ),
        RabiesScientificFindingRow(
            finding_id="migration_events",
            question="How much regional movement is implied by the geographic reconstruction?",
            claim=(
                f"The biogeography layer inferred {migration_summary.event_count} migration events across the rooted tree."
            ),
            evidence=(
                f"strongly supported migration events {migration_summary.strongly_supported_event_count}"
            ),
            caution=(
                "Event counts summarize transitions over grouped regions and do not replace one dated dispersal analysis."
            ),
            source_artifact="biogeography/event-table.tsv",
        ),
    ]


def _write_scientific_findings_table(
    path: Path,
    rows: list[RabiesScientificFindingRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "finding_id",
            "question",
            "claim",
            "evidence",
            "caution",
            "source_artifact",
        ],
        rows=[
            {
                "finding_id": row.finding_id,
                "question": row.question,
                "claim": row.claim,
                "evidence": row.evidence,
                "caution": row.caution,
                "source_artifact": row.source_artifact,
            }
            for row in rows
        ],
    )


def _write_input_validation_table(
    path: Path,
    *,
    workflow: FastaToTreeWorkflowReport,
) -> Path:
    validation = (
        workflow.input_validation
        if workflow.repaired_input_validation is None
        else workflow.repaired_input_validation
    )
    sequence_type_report = validation.sequence_type_report
    row = {
        "sequence_count": str(validation.summary.sequence_count),
        "detected_type": sequence_type_report.detected_type or "",
        "selected_type": sequence_type_report.selected_type or "",
        "confidence": sequence_type_report.confidence or "",
        "repair_required": "true"
        if (
            workflow.input_validation.duplicate_identifiers
            or workflow.input_validation.illegal_characters
            or workflow.input_validation.empty_sequences
        )
        else "false",
        "repair_applied": "true" if workflow.input_repair is not None else "false",
        "duplicate_identifier_count": str(
            len(workflow.input_validation.duplicate_identifiers)
        ),
        "illegal_character_count": str(
            len(workflow.input_validation.illegal_characters)
        ),
        "empty_sequence_count": str(len(workflow.input_validation.empty_sequences)),
        "warning_count": str(len(validation.warnings)),
        "warnings": " | ".join(validation.warnings),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_alignment_quality_table(
    path: Path,
    *,
    aligned: AlignmentQualityReport,
    trimmed: AlignmentQualityReport,
) -> Path:
    rows = []
    for stage, report in (("aligned", aligned), ("trimmed", trimmed)):
        rows.append(
            {
                "stage": stage,
                "sequence_count": str(report.sequence_count),
                "alignment_length": str(report.alignment_length),
                "missing_data_fraction": _format_number(report.missing_data_fraction),
                "gap_fraction": _format_number(report.gap_fraction),
                "ambiguity_fraction": _format_number(report.ambiguity_fraction),
                "variable_site_count": str(report.variable_site_count),
                "parsimony_informative_site_count": str(
                    report.parsimony_informative_site_count
                ),
                "quality_score": _format_number(report.quality_score),
                "suspicious_alignment": (
                    "true" if report.suspicious_alignment else "false"
                ),
                "suspicious_reasons": " | ".join(report.suspicious_reasons),
            }
        )
    return write_taxon_rows(path, columns=list(rows[0].keys()), rows=rows)


def _write_sequence_ranking_table(
    path: Path,
    report: SequenceQualityRankingReport,
) -> Path:
    rows = [
        {
            "identifier": row.identifier,
            "rank": str(row.rank),
            "score": _format_number(row.score),
            "missing_fraction": _format_number(row.missing_fraction),
            "gap_fraction": _format_number(row.gap_fraction),
            "ambiguity_fraction": _format_number(row.ambiguity_fraction),
            "composition_outlier": "true" if row.composition_outlier else "false",
            "duplicate_status": row.duplicate_status,
            "note": row.note,
        }
        for row in report.rows
    ]
    return write_taxon_rows(path, columns=list(rows[0].keys()), rows=rows)


def _write_comparative_branch_repairs_table(
    path: Path,
    rows: list[RabiesComparativeBranchRepair],
) -> Path:
    if not rows:
        return write_taxon_rows(
            path,
            columns=[
                "node_label",
                "original_branch_length",
                "repaired_branch_length",
                "reason",
            ],
            rows=[],
        )
    return write_taxon_rows(
        path,
        columns=[
            "node_label",
            "original_branch_length",
            "repaired_branch_length",
            "reason",
        ],
        rows=[
            {
                "node_label": row.node_label,
                "original_branch_length": _format_number(row.original_branch_length),
                "repaired_branch_length": _format_number(row.repaired_branch_length),
                "reason": row.reason,
            }
            for row in rows
        ],
    )


def _write_workflow_summary_table(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    clade_row_count: int,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    scientific_finding_count: int,
) -> Path:
    support = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_summary = bootstrap_artifacts.summary_report
    row = {
        "dataset_id": report.dataset.dataset_id,
        "sequence_count": str(report.dataset.sequence_count),
        "sequence_type": report.dataset.sequence_type,
        "inferred_sequence_type": report.fasta_to_tree.sequence_type,
        "selected_model": report.fasta_to_tree.selected_model,
        "input_repair_applied": (
            "true" if report.fasta_to_tree.input_repair is not None else "false"
        ),
        "aligned_quality_score": _format_number(report.aligned_quality.quality_score),
        "trimmed_quality_score": _format_number(report.trimmed_quality.quality_score),
        "minimum_support": _format_number(support.minimum_support),
        "maximum_support": _format_number(support.maximum_support),
        "median_support": _format_number(support.median_support),
        "weakly_supported_clade_count": str(support.weakly_supported_clade_count),
        "clade_row_count": str(clade_row_count),
        "bootstrap_tree_count": str(bootstrap_summary.tree_count),
        "bootstrap_topology_count": str(
            bootstrap_summary.diversity.rooted_topology_count
        ),
        "bootstrap_unstable_branch_count": str(bootstrap_summary.unstable_branch_count),
        "bootstrap_consensus_rooted_rf_distance": str(
            bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
        ),
        "bootstrap_consensus_same_unrooted_topology": (
            "true"
            if bootstrap_tree_comparison_report.topology.same_unrooted_topology
            else "false"
        ),
        "bootstrap_consensus_high_support_conflict_count": str(
            len(
                [
                    row
                    for row in bootstrap_tree_comparison_report.support.conflicting_clades
                    if row.conflict_classification == "high_support_conflict"
                ]
            )
        ),
        "outgroup_taxa": ",".join(report.dataset.outgroup_taxa),
        "root_host": host_summary.root_host,
        "root_host_confidence": _format_number(host_summary.root_confidence),
        "host_switch_count": str(host_summary.host_switch_count),
        "certain_host_switch_count": str(host_summary.certain_host_switch_count),
        "uncertain_host_switch_count": str(host_summary.uncertain_host_switch_count),
        "root_region": geography_summary.root_region,
        "root_region_probability": _format_number(
            geography_summary.root_region_probability
        ),
        "changed_region_branch_count": str(geography_summary.changed_branch_count),
        "migration_event_count": str(migration_summary.event_count),
        "strongly_supported_migration_event_count": str(
            migration_summary.strongly_supported_event_count
        ),
        "comparative_response": comparative_summary_row.response,
        "comparative_formula": comparative_summary_row.formula,
        "comparative_selected_model": comparative_summary_row.selected_model,
        "comparative_pgls_lambda": _format_number(comparative_summary_row.pgls_lambda),
        "comparative_pgls_r_squared": _format_number(
            comparative_summary_row.pgls_r_squared
        ),
        "comparative_branch_repair_count": str(len(report.comparative_branch_repairs)),
        "conclusion_stable_count": str(
            report.conclusion_stability_report.summary.stable_count
        ),
        "conclusion_weak_count": str(
            report.conclusion_stability_report.summary.weak_count
        ),
        "conclusion_unstable_count": str(
            report.conclusion_stability_report.summary.unstable_count
        ),
        "timeout_seconds": _format_number(report.config.timeout_seconds),
        "max_bootstrap_tree_count": (
            ""
            if report.config.max_bootstrap_tree_count is None
            else str(report.config.max_bootstrap_tree_count)
        ),
        "max_report_table_rows": (
            ""
            if report.config.max_report_table_rows is None
            else str(report.config.max_report_table_rows)
        ),
        "memory_warning_threshold_bytes": (
            ""
            if report.config.memory_warning_threshold_bytes is None
            else str(report.config.memory_warning_threshold_bytes)
        ),
        "budget_warning_count": str(
            len(bootstrap_artifacts.budget_report.warning_messages)
        ),
        "config_check_count": str(len(report.config_audit_rows)),
        "scientific_finding_count": str(scientific_finding_count),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_resource_observation_table(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
) -> Path:
    row = {
        "dataset_id": report.dataset.dataset_id,
        "timeout_seconds": _format_number(report.config.timeout_seconds),
        "workflow_runtime_seconds": _format_number(
            report.fasta_to_tree.runtime_seconds
        ),
        "bootstrap_review_runtime_seconds": _format_number(
            bootstrap_artifacts.summary_report.processing.runtime_seconds
        ),
        "bootstrap_review_peak_memory_bytes": str(
            bootstrap_artifacts.summary_report.processing.peak_memory_bytes
        ),
        "budget_warning_count": str(
            len(bootstrap_artifacts.budget_report.warning_messages)
        ),
        "budget_warnings": " | ".join(
            bootstrap_artifacts.budget_report.warning_messages
        ),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _artifact_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if path.name.endswith(".manifest.json"):
        return "manifest"
    if suffix in {".html", ".htm"}:
        return "report"
    if suffix == ".md":
        return "markdown"
    if suffix == ".json":
        return "json"
    if suffix == ".tsv":
        return "table"
    if suffix == ".svg":
        return "figure"
    if suffix == ".log":
        return "log"
    if suffix in {".nwk", ".tree"}:
        return "tree"
    if suffix in {".aln", ".fasta"}:
        return "alignment"
    if suffix == ".csv":
        return "metadata"
    return "artifact"


def _workflow_artifact_section(relative_path: Path) -> str:
    parts = relative_path.parts
    if len(parts) >= 2 and parts[1] in {
        "bootstrap-review",
        "engine-artifacts",
        "biogeography",
        "comparative",
        "conclusion-stability",
    }:
        return parts[1]
    return "workflow"


def _relative_to_package_root(package_root: Path, path: Path) -> str:
    return path.relative_to(package_root).as_posix()


def _package_inventory_rows(
    *,
    output_root: Path,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    overview_path: Path,
    overview_html_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    dataset_paths = [
        dataset_export.readme_path,
        dataset_export.workflow_config_path,
        dataset_export.sequences_path,
        dataset_export.metadata_path,
        dataset_export.centroids_path,
        dataset_export.accession_table_path,
    ]
    for path in dataset_paths:
        rows.append(
            {
                "section": "dataset",
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    workflow_paths = sorted(
        path
        for path in workflow_bundle.output_root.rglob("*")
        if path.is_file()
    )
    for path in workflow_paths:
        rows.append(
            {
                "section": _workflow_artifact_section(
                    path.relative_to(workflow_bundle.output_root)
                ),
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    for path in (overview_path, overview_html_path):
        rows.append(
            {
                "section": "package",
                "kind": _artifact_kind(path),
                "relative_path": _relative_to_package_root(output_root, path),
                "sha256": _checksum(path),
                "size_bytes": str(path.stat().st_size),
            }
        )
    return rows


def _write_package_artifact_inventory(
    path: Path,
    *,
    output_root: Path,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    overview_path: Path,
    overview_html_path: Path,
) -> tuple[Path, list[dict[str, str]]]:
    rows = _package_inventory_rows(
        output_root=output_root,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
        overview_html_path=overview_html_path,
    )
    return (
        write_taxon_rows(
            path,
            columns=["section", "kind", "relative_path", "sha256", "size_bytes"],
            rows=rows,
        ),
        rows,
    )


def _has_package_artifact(
    inventory_rows: list[dict[str, str]],
    relative_path: str,
) -> bool:
    return any(row["relative_path"] == relative_path for row in inventory_rows)


def _package_inventory_counts(
    inventory_rows: list[dict[str, str]],
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in inventory_rows:
        section = row["section"]
        counts[section] = counts.get(section, 0) + 1
    return counts


def _write_package_reproducibility_checklist(
    path: Path,
    *,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    inventory_rows: list[dict[str, str]],
    artifact_inventory_path: Path,
) -> tuple[Path, list[dict[str, str]]]:
    rows = [
        {
            "section": "inputs",
            "check_id": "dataset-inputs-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        "dataset/sequences.fasta",
                        "dataset/metadata.csv",
                        "dataset/region-centroids.csv",
                        "dataset/source-accessions.tsv",
                        "dataset/workflow-config.json",
                    )
                )
                else "blocked"
            ),
            "summary": "raw sequences, metadata, centroids, accession ledger, and workflow config are exported together",
            "evidence": "dataset export includes all package-level study inputs required to rerun the workflow",
            "artifact_path": "dataset/workflow-config.json",
        },
        {
            "section": "tree-inference",
            "check_id": "tree-inference-evidence-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        f"workflow/{workflow_bundle.tree_path.name}",
                        f"workflow/{workflow_bundle.model_table_path.name}",
                        f"workflow/{workflow_bundle.support_table_path.name}",
                        f"workflow/{workflow_bundle.manifest_path.name}",
                        f"workflow/{workflow_bundle.log_path.name}",
                    )
                )
                else "blocked"
            ),
            "summary": "tree inference outputs, manifest, and engine log are preserved",
            "evidence": "the package retains the rooted tree plus model, support, manifest, and log artifacts for the sequence-to-tree workflow",
            "artifact_path": f"workflow/{workflow_bundle.manifest_path.name}",
        },
        {
            "section": "uncertainty",
            "check_id": "bootstrap-uncertainty-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        "workflow/bootstrap-review/bootstrap-review.summary.tsv",
                        "workflow/bootstrap-review/bootstrap-consensus.tree",
                        "workflow/bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv",
                        "workflow/bootstrap-review/rooted-tree-vs-bootstrap-consensus.report.html",
                    )
                )
                else "blocked"
            ),
            "summary": "bootstrap consensus and topology-conflict review artifacts are retained",
            "evidence": (
                "bootstrap review exported "
                f"{workflow_bundle.bootstrap_tree_count} trees over "
                f"{workflow_bundle.bootstrap_topology_count} rooted topologies"
            ),
            "artifact_path": "workflow/bootstrap-review/bootstrap-review.summary.tsv",
        },
        {
            "section": "analysis",
            "check_id": "downstream-analysis-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        f"workflow/{workflow_bundle.host_switch_summary_path.name}",
                        f"workflow/biogeography/{workflow_bundle.biogeography_report_path.name}",
                        f"workflow/comparative/{workflow_bundle.comparative_report_path.name}",
                        (
                            "workflow/conclusion-stability/"
                            f"{workflow_bundle.conclusion_stability_report_path.name}"
                        ),
                        f"workflow/{workflow_bundle.scientific_findings_path.name}",
                    )
                )
                else "blocked"
            ),
            "summary": "host-switching, biogeography, comparative, stability, and findings surfaces are preserved together",
            "evidence": (
                "the package retains one integrated downstream evidence chain from "
                "rooted tree to host/geography interpretation and stability review"
            ),
            "artifact_path": f"workflow/{workflow_bundle.scientific_findings_path.name}",
        },
        {
            "section": "package",
            "check_id": "package-navigation-exported",
            "status": (
                "pass"
                if all(
                    _has_package_artifact(inventory_rows, relative_path)
                    for relative_path in (
                        "overview.md",
                        "rabies-cross-host-geography-overview.html",
                    )
                )
                else "blocked"
            ),
            "summary": "reviewer overview surfaces are included with the package",
            "evidence": (
                "the package includes one markdown overview, one reviewer HTML overview, "
                f"and one artifact inventory at {artifact_inventory_path.name}"
            ),
            "artifact_path": "rabies-cross-host-geography-overview.html",
        },
        {
            "section": "limitations",
            "check_id": "interpretation-risks-surfaced",
            "status": (
                "risk"
                if (
                    workflow_bundle.budget_warning_count > 0
                    or workflow_bundle.conclusion_weak_count > 0
                    or workflow_bundle.conclusion_unstable_count > 0
                )
                else "pass"
            ),
            "summary": "the package records interpretation limits and stability caveats",
            "evidence": (
                "budget warnings="
                f"{workflow_bundle.budget_warning_count}; weak conclusions="
                f"{workflow_bundle.conclusion_weak_count}; unstable conclusions="
                f"{workflow_bundle.conclusion_unstable_count}"
            ),
            "artifact_path": (
                "workflow/conclusion-stability/"
                f"{workflow_bundle.conclusion_stability_report_path.name}"
            ),
        },
    ]
    return (
        write_taxon_rows(
            path,
            columns=[
                "section",
                "check_id",
                "status",
                "summary",
                "evidence",
                "artifact_path",
            ],
            rows=rows,
        ),
        rows,
    )


def _write_overview(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    reproducibility_checklist_path: Path,
) -> Path:
    lines = [
        "# Rabies Cross-Host Geography Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- sequence count: `{dataset.sequence_count}`",
        f"- workflow config: `{config.config_path.name}`",
        f"- biological question: {_FLAGSHIP_QUESTION}",
        f"- short answer: {short_answer}",
        f"- host workflow trait: `{dataset.host_trait}`",
        f"- geography workflow trait: `{dataset.geography_trait}`",
        f"- comparative formula: `{workflow_bundle.comparative_formula}`",
        "",
        "Generated outputs:",
        "",
        "- source accession ledger: `dataset/source-accessions.tsv`",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- resource observations: `{workflow_bundle.resource_observations_path.name}`",
        f"- clade table: `{workflow_bundle.clade_table_path.name}`",
        f"- bootstrap review: `bootstrap-review/{workflow_bundle.bootstrap_summary_path.name}`",
        (
            "- rooted-versus-consensus comparison: "
            f"`bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}`"
        ),
        f"- comparative report: `comparative/{workflow_bundle.comparative_report_path.name}`",
        (
            "- conclusion stability report: "
            f"`conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}`"
        ),
        f"- final report: `{workflow_bundle.final_report_path.name}`",
        f"- package artifact inventory: `{artifact_inventory_path.name}`",
        f"- package reproducibility checklist: `{reproducibility_checklist_path.name}`",
        "- package overview html: `rabies-cross-host-geography-overview.html`",
        "- package manifest: `rabies-cross-host-geography-package.manifest.json`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _build_flagship_answer_summary(
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
) -> str:
    return (
        "The rooted panel remains anchored in "
        f"{workflow_bundle.root_host} and {workflow_bundle.root_region}, and "
        "`host_group[canid]` shows a nominally supported positive longitude "
        "association under the selected comparative model, but the inference "
        "remains cautionary because the panel is intentionally compact."
    )


def _write_demo_overview_html(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    reproducibility_checklist_path: Path,
) -> Path:
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Rabies Cross-Host Geography Package Overview</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #eef5ef 100%); color: #173024; }",
            "    main { max-width: 1040px; margin: 0 auto; padding: 28px; }",
            "    h1, h2 { margin: 0 0 10px; }",
            "    p { line-height: 1.6; }",
            "    .panel { background: rgba(255,255,255,0.9); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 18px; margin-top: 18px; box-shadow: 0 14px 36px rgba(23,48,36,0.08); }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin-top: 18px; }",
            "    .card { background: rgba(255,255,255,0.9); border: 1px solid rgba(23,48,36,0.12); border-radius: 16px; padding: 16px; }",
            "    .label { color: #5f7469; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .value { display: block; font-size: 21px; margin-top: 6px; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    code { background: rgba(23,48,36,0.06); padding: 0 4px; border-radius: 4px; }",
            "    a { color: #23523b; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Rabies Cross-Host Geography Package</h1>",
            "  <p>This public flagship workflow starts from one real rabies nucleoprotein FASTA panel plus grouped host and geographic metadata, reruns the full owned sequence-to-tree workflow, and then connects the rooted tree to host-switching, geographic transition, migration, and comparative review surfaces.</p>",
            '  <section class="panel">',
            "    <h2>Biological Question</h2>",
            f"    <p>{escape(_FLAGSHIP_QUESTION)}</p>",
            '    <h2 style="margin-top: 16px;">Short Answer</h2>',
            f"    <p>{escape(short_answer)}</p>",
            "  </section>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">dataset id</span><span class="value">{escape(dataset.dataset_id)}</span></div>',
            f'    <div class="card"><span class="label">sequence count</span><span class="value">{dataset.sequence_count}</span></div>',
            f'    <div class="card"><span class="label">selected model</span><span class="value">{escape(workflow_bundle.selected_model)}</span></div>',
            f'    <div class="card"><span class="label">comparative model</span><span class="value">{escape(workflow_bundle.comparative_selected_model)}</span></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Reproducibility Surface</h2>",
            _html_list(
                [
                    f'workflow config: <a href="dataset/{dataset_export.workflow_config_path.name}">dataset/{dataset_export.workflow_config_path.name}</a>',
                    f'source accession ledger: <a href="dataset/{dataset_export.accession_table_path.name}">dataset/{dataset_export.accession_table_path.name}</a>',
                    f'final workflow manifest: <a href="workflow/{workflow_bundle.final_manifest_path.name}">workflow/{workflow_bundle.final_manifest_path.name}</a>',
                    f'package artifact inventory: <a href="{artifact_inventory_path.name}">{artifact_inventory_path.name}</a>',
                    f'package reproducibility checklist: <a href="{reproducibility_checklist_path.name}">{reproducibility_checklist_path.name}</a>',
                    f'package manifest: <a href="{path.name.replace("-overview.html", "-package.manifest.json")}">{path.name.replace("-overview.html", "-package.manifest.json")}</a>',
                ]
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Key Outputs</h2>",
            _html_list(
                [
                    f'<a href="workflow/{workflow_bundle.final_report_path.name}">workflow/{workflow_bundle.final_report_path.name}</a>',
                    f'<a href="workflow/{workflow_bundle.workflow_summary_path.name}">workflow/{workflow_bundle.workflow_summary_path.name}</a>',
                    f'<a href="workflow/bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}">workflow/bootstrap-review/{workflow_bundle.bootstrap_tree_comparison_summary_path.name}</a>',
                    f'<a href="workflow/comparative/{workflow_bundle.comparative_report_path.name}">workflow/comparative/{workflow_bundle.comparative_report_path.name}</a>',
                    f'<a href="workflow/conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}">workflow/conclusion-stability/{workflow_bundle.conclusion_stability_report_path.name}</a>',
                    f'<a href="workflow/{workflow_bundle.scientific_findings_path.name}">workflow/{workflow_bundle.scientific_findings_path.name}</a>',
                    f'<a href="dataset/{dataset_export.workflow_config_path.name}">dataset/{dataset_export.workflow_config_path.name}</a>',
                ]
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _write_demo_package_manifest(
    path: Path,
    *,
    dataset: RabiesCrossHostGeographyPanelDataset,
    dataset_export: RabiesCrossHostGeographyPanelExportResult,
    workflow_bundle: RabiesCrossHostGeographyPanelWorkflowBundle,
    config: RabiesCrossHostGeographyPanelWorkflowConfig,
    short_answer: str,
    artifact_inventory_path: Path,
    artifact_inventory_rows: list[dict[str, str]],
    reproducibility_checklist_path: Path,
    checklist_rows: list[dict[str, str]],
) -> Path:
    inventory_counts = _package_inventory_counts(artifact_inventory_rows)
    blocked_check_count = len(
        [row for row in checklist_rows if row["status"] == "blocked"]
    )
    risk_check_count = len([row for row in checklist_rows if row["status"] == "risk"])
    payload = {
        "report_kind": "rabies_cross_host_geography_package",
        "dataset_id": dataset.dataset_id,
        "label": dataset.label,
        "biological_question": _FLAGSHIP_QUESTION,
        "short_answer": short_answer,
        "package_files": {
            "overview_markdown": {
                "path": "overview.md",
                "checksum": _checksum(path.parent / "overview.md"),
            },
            "overview_html": {
                "path": "rabies-cross-host-geography-overview.html",
                "checksum": _checksum(
                    path.parent / "rabies-cross-host-geography-overview.html"
                ),
            },
            "artifact_inventory": {
                "path": artifact_inventory_path.name,
                "checksum": _checksum(artifact_inventory_path),
                "artifact_count": len(artifact_inventory_rows),
                "section_counts": inventory_counts,
            },
            "reproducibility_checklist": {
                "path": reproducibility_checklist_path.name,
                "checksum": _checksum(reproducibility_checklist_path),
                "item_count": len(checklist_rows),
                "blocked_count": blocked_check_count,
                "risk_count": risk_check_count,
            },
        },
        "config": {
            "path": f"dataset/{dataset_export.workflow_config_path.name}",
            "checksum": _checksum(dataset_export.workflow_config_path),
            "workflow_prefix": config.workflow_prefix,
            "alignment_mode": config.alignment_mode,
            "trimming_mode": config.trimming_mode,
            "trim_gap_threshold": config.trim_gap_threshold,
            "bootstrap_consensus_threshold": config.bootstrap_consensus_threshold,
            "bootstrap_robust_support_threshold": (
                config.bootstrap_robust_support_threshold
            ),
            "comparative_formula": config.comparative_formula,
            "comparative_response": config.comparative_response,
            "comparative_branch_length_floor": (
                config.comparative_branch_length_floor
            ),
            "timeout_seconds": config.timeout_seconds,
            "max_bootstrap_tree_count": config.max_bootstrap_tree_count,
            "max_report_table_rows": config.max_report_table_rows,
            "memory_warning_threshold_bytes": config.memory_warning_threshold_bytes,
        },
        "dataset_files": {
            "readme": {
                "path": f"dataset/{dataset_export.readme_path.name}",
                "checksum": _checksum(dataset_export.readme_path),
            },
            "sequences": {
                "path": f"dataset/{dataset_export.sequences_path.name}",
                "checksum": _checksum(dataset_export.sequences_path),
            },
            "metadata": {
                "path": f"dataset/{dataset_export.metadata_path.name}",
                "checksum": _checksum(dataset_export.metadata_path),
            },
            "centroids": {
                "path": f"dataset/{dataset_export.centroids_path.name}",
                "checksum": _checksum(dataset_export.centroids_path),
            },
            "source_accessions": {
                "path": f"dataset/{dataset_export.accession_table_path.name}",
                "checksum": _checksum(dataset_export.accession_table_path),
            },
        },
        "workflow_files": {
            "final_report": {
                "path": f"workflow/{workflow_bundle.final_report_path.name}",
                "checksum": _checksum(workflow_bundle.final_report_path),
            },
            "workflow_log": {
                "path": f"workflow/{workflow_bundle.log_path.name}",
                "checksum": _checksum(workflow_bundle.log_path),
            },
            "workflow_summary": {
                "path": f"workflow/{workflow_bundle.workflow_summary_path.name}",
                "checksum": _checksum(workflow_bundle.workflow_summary_path),
            },
            "resource_observations": {
                "path": f"workflow/{workflow_bundle.resource_observations_path.name}",
                "checksum": _checksum(workflow_bundle.resource_observations_path),
            },
            "final_manifest": {
                "path": f"workflow/{workflow_bundle.final_manifest_path.name}",
                "checksum": _checksum(workflow_bundle.final_manifest_path),
            },
            "rooted_tree": {
                "path": f"workflow/{workflow_bundle.tree_path.name}",
                "checksum": _checksum(workflow_bundle.tree_path),
            },
            "rooting_report": {
                "path": f"workflow/{workflow_bundle.rooting_report_path.name}",
                "checksum": _checksum(workflow_bundle.rooting_report_path),
            },
            "model_table": {
                "path": f"workflow/{workflow_bundle.model_table_path.name}",
                "checksum": _checksum(workflow_bundle.model_table_path),
            },
            "support_table": {
                "path": f"workflow/{workflow_bundle.support_table_path.name}",
                "checksum": _checksum(workflow_bundle.support_table_path),
            },
            "bootstrap_summary": {
                "path": (
                    "workflow/bootstrap-review/"
                    f"{workflow_bundle.bootstrap_summary_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.bootstrap_summary_path),
            },
            "bootstrap_tree_comparison_summary": {
                "path": (
                    "workflow/bootstrap-review/"
                    f"{workflow_bundle.bootstrap_tree_comparison_summary_path.name}"
                ),
                "checksum": _checksum(
                    workflow_bundle.bootstrap_tree_comparison_summary_path
                ),
            },
            "host_switch_summary": {
                "path": f"workflow/{workflow_bundle.host_switch_summary_path.name}",
                "checksum": _checksum(workflow_bundle.host_switch_summary_path),
            },
            "biogeography_report": {
                "path": (
                    "workflow/biogeography/"
                    f"{workflow_bundle.biogeography_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.biogeography_report_path),
            },
            "comparative_report": {
                "path": (
                    "workflow/comparative/"
                    f"{workflow_bundle.comparative_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.comparative_report_path),
            },
            "conclusion_stability_report": {
                "path": (
                    "workflow/conclusion-stability/"
                    f"{workflow_bundle.conclusion_stability_report_path.name}"
                ),
                "checksum": _checksum(workflow_bundle.conclusion_stability_report_path),
            },
            "scientific_findings": {
                "path": f"workflow/{workflow_bundle.scientific_findings_path.name}",
                "checksum": _checksum(workflow_bundle.scientific_findings_path),
            },
        },
        "metrics": {
            "sequence_count": dataset.sequence_count,
            "selected_model": workflow_bundle.selected_model,
            "root_host": workflow_bundle.root_host,
            "root_region": workflow_bundle.root_region,
            "bootstrap_tree_count": workflow_bundle.bootstrap_tree_count,
            "workflow_runtime_seconds": workflow_bundle.workflow_runtime_seconds,
            "bootstrap_review_runtime_seconds": (
                workflow_bundle.bootstrap_review_runtime_seconds
            ),
            "bootstrap_review_peak_memory_bytes": (
                workflow_bundle.bootstrap_review_peak_memory_bytes
            ),
            "budget_warning_count": workflow_bundle.budget_warning_count,
            "host_switch_count": workflow_bundle.host_switch_count,
            "migration_event_count": workflow_bundle.migration_event_count,
            "comparative_selected_model": workflow_bundle.comparative_selected_model,
            "comparative_pgls_lambda": workflow_bundle.comparative_pgls_lambda,
            "comparative_pgls_r_squared": workflow_bundle.comparative_pgls_r_squared,
            "comparative_branch_repair_count": (
                workflow_bundle.comparative_branch_repair_count
            ),
            "conclusion_stable_count": workflow_bundle.conclusion_stable_count,
            "conclusion_weak_count": workflow_bundle.conclusion_weak_count,
            "conclusion_unstable_count": workflow_bundle.conclusion_unstable_count,
            "scientific_finding_count": workflow_bundle.scientific_finding_count,
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _write_comparative_report(
    path: Path,
    *,
    summary_row: ComparativeAnalysisSummaryRow,
    coefficient_rows: list[ComparativeCoefficientTableRow],
    residual_rows: list[ComparativeResidualTableRow],
    signal_row: ComparativeSignalTableRow,
    interpretation_rows: list[ComparativeInterpretationRow],
    branch_repairs: list[RabiesComparativeBranchRepair],
) -> Path:
    key_claim = next(
        (
            row.claim
            for row in interpretation_rows
            if row.topic == "coefficient" and "nominally supported" in row.claim
        ),
        "no coefficient reached nominal support",
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Rabies Comparative Report</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #f3f7ef 100%); color: #173024; }",
            "    main { max-width: 1040px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { margin: 0 0 10px; }",
            "    p { line-height: 1.55; }",
            "    .panel { background: rgba(255,255,255,0.88); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 18px; margin-top: 18px; box-shadow: 0 14px 36px rgba(23,48,36,0.08); }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin-top: 18px; }",
            "    .card { background: rgba(255,255,255,0.88); border: 1px solid rgba(23,48,36,0.12); border-radius: 18px; padding: 16px; }",
            "    .label { color: #5f7469; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .value { display: block; font-size: 22px; margin-top: 6px; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(23,48,36,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Rabies Comparative Report</h1>",
            "  <p>This comparative section asks whether the host-associated lineages in the rabies demonstration tree are associated with a consistent eastward geographic placement when geography is summarized as regional longitude. The result is interpretive evidence, not causal proof, and it inherits the small-panel limits of this dataset.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">formula</span><span class="value">{escape(summary_row.formula)}</span></div>',
            f'    <div class="card"><span class="label">analysis taxa</span><span class="value">{summary_row.analysis_taxa}</span></div>',
            f'    <div class="card"><span class="label">selected trait model</span><span class="value">{escape(summary_row.selected_model)}</span></div>',
            f'    <div class="card"><span class="label">pgls r-squared</span><span class="value">{_format_number(summary_row.pgls_r_squared)}</span></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Question and Answer</h2>",
            _html_list(
                [
                    "Question: does host association coincide with a consistent longitudinal shift in this rabies panel?",
                    f"Answer: {key_claim}.",
                    (
                        f"Phylogenetic signal remains strong for the response trait "
                        f"(Blomberg's K {_format_number(signal_row.blombergs_k)}, "
                        f"Pagel's lambda {_format_number(signal_row.pagels_lambda)})."
                    ),
                    (
                        "Interpret the coefficient evidence cautiously because the "
                        "residual diagnostics retain review warnings and the sample "
                        "is intentionally compact."
                    ),
                ]
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Coefficient Summary</h2>",
            _table(
                headers=[
                    "term",
                    "estimate",
                    "standard_error",
                    "p_value",
                    "significant",
                ],
                rows=[
                    [
                        row.term,
                        _format_number(row.estimate),
                        _format_number(row.standard_error),
                        _format_number(row.p_value),
                        "true" if row.significant else "false",
                    ]
                    for row in coefficient_rows
                ],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Residual Diagnostics</h2>",
            _table(
                headers=[
                    "analysis",
                    "residual_variance",
                    "max_abs_standardized_residual",
                    "phylogenetic_residual_lambda",
                    "warnings",
                ],
                rows=[
                    [
                        row.analysis,
                        _format_number(row.residual_variance),
                        _format_number(row.max_abs_standardized_residual),
                        _format_number(row.phylogenetic_residual_lambda),
                        "; ".join(row.warnings),
                    ]
                    for row in residual_rows
                ],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Comparative Tree Adjustments</h2>",
            _html_list(
                [
                    "The comparative fit uses the rooted demonstration tree after flooring any nonpositive nonroot branch lengths to a tiny positive value.",
                    f"Adjusted branch count: {len(branch_repairs)}",
                ]
            ),
            _table(
                headers=[
                    "node_label",
                    "original_branch_length",
                    "repaired_branch_length",
                    "reason",
                ],
                rows=[
                    [
                        row.node_label,
                        _format_number(row.original_branch_length),
                        _format_number(row.repaired_branch_length),
                        row.reason,
                    ]
                    for row in branch_repairs
                ]
                or [["", "", "", "no branch-length repair was required"]],
            ),
            "  </section>",
            '  <section class="panel">',
            "    <h2>Interpretation Ledger</h2>",
            _table(
                headers=["topic", "claim", "evidence", "caution"],
                rows=[
                    [row.topic, row.claim, row.evidence, row.caution]
                    for row in interpretation_rows
                ],
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _write_comparative_manifest(
    path: Path,
    *,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    branch_repairs: list[RabiesComparativeBranchRepair],
    output_paths: dict[str, Path],
) -> Path:
    payload = {
        "report_kind": "rabies_cross_host_geography_comparative_bundle",
        "metrics": {
            "response": comparative_summary_row.response,
            "formula": comparative_summary_row.formula,
            "analysis_taxa": comparative_summary_row.analysis_taxa,
            "selected_model": comparative_summary_row.selected_model,
            "pgls_lambda": comparative_summary_row.pgls_lambda,
            "pgls_r_squared": comparative_summary_row.pgls_r_squared,
            "branch_repair_count": len(branch_repairs),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in output_paths.items()
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_manifest(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    clade_row_count: int,
    scientific_finding_count: int,
    bundle_paths: dict[str, Path],
) -> Path:
    manifest = {
        "report_kind": "rabies_cross_host_geography_workflow_bundle",
        "dataset_id": report.dataset.dataset_id,
        "input_checksums": {
            "workflow-config.json": _checksum(report.dataset.workflow_config_path),
            "sequences.fasta": _checksum(report.dataset.sequences_path),
            "metadata.csv": _checksum(report.dataset.metadata_path),
            "region-centroids.csv": _checksum(report.dataset.centroids_path),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in bundle_paths.items()
        },
        "metrics": {
            "sequence_count": report.dataset.sequence_count,
            "selected_model": report.fasta_to_tree.selected_model,
            "minimum_support": report.fasta_to_tree.support_summary.minimum_support,
            "maximum_support": report.fasta_to_tree.support_summary.maximum_support,
            "host_switch_count": report.host_switching.summary.host_switch_count,
            "migration_event_count": report.biogeography_report.event_report.summary.event_count,
            "root_host": report.host_switching.summary.root_host,
            "root_region": report.biogeography_report.state_report.summary.root_region,
            "clade_row_count": clade_row_count,
            "bootstrap_tree_count": bootstrap_artifacts.summary_report.tree_count,
            "bootstrap_unstable_branch_count": (
                bootstrap_artifacts.summary_report.unstable_branch_count
            ),
            "bootstrap_consensus_rooted_rf_distance": (
                bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
            ),
            "bootstrap_consensus_same_unrooted_topology": (
                bootstrap_tree_comparison_report.topology.same_unrooted_topology
            ),
            "bootstrap_consensus_high_support_conflict_count": len(
                [
                    row
                    for row in bootstrap_tree_comparison_report.support.conflicting_clades
                    if row.conflict_classification == "high_support_conflict"
                ]
            ),
            "comparative_selected_model": comparative_summary_row.selected_model,
            "comparative_pgls_lambda": comparative_summary_row.pgls_lambda,
            "comparative_pgls_r_squared": comparative_summary_row.pgls_r_squared,
            "conclusion_stable_count": (
                report.conclusion_stability_report.summary.stable_count
            ),
            "conclusion_weak_count": (
                report.conclusion_stability_report.summary.weak_count
            ),
            "conclusion_unstable_count": (
                report.conclusion_stability_report.summary.unstable_count
            ),
            "config_check_count": len(report.config_audit_rows),
            "scientific_finding_count": scientific_finding_count,
        },
    }
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _write_integrated_report(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    workflow_summary_path: Path,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    clade_row_count: int,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
    comparative_branch_repair_count: int,
    scientific_finding_rows: list[RabiesScientificFindingRow],
    max_report_table_rows: int | None,
) -> Path:
    support_summary = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_summary = bootstrap_artifacts.summary_report
    core_question = (
        "Do the host-associated rabies lineages in this compact panel occupy one "
        "distinct geographic regime while retaining one coherent phylogenetic signal?"
    )
    core_answer = next(
        (
            row.claim
            for row in comparative_interpretation_rows
            if row.topic == "coefficient" and "nominally supported" in row.claim
        ),
        "the comparative layer did not recover a nominally supported host effect",
    )
    html = "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Rabies Host and Geography Workflow</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f1ea 0%, #e7efe7 100%); color: #163222; }",
            "    main { max-width: 1360px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 34px; }",
            "    h2 { margin: 0 0 10px; font-size: 22px; }",
            "    p { line-height: 1.55; }",
            "    .cards { display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 14px; margin: 18px 0 24px; }",
            "    .card, .panel { background: rgba(255,255,255,0.86); border: 1px solid rgba(22,50,34,0.12); border-radius: 18px; padding: 18px; box-shadow: 0 16px 42px rgba(22,50,34,0.08); }",
            "    .label { color: #5b7466; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .card strong { display: block; font-size: 21px; margin-top: 6px; }",
            "    .grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }",
            "    .full { grid-column: 1 / -1; }",
            "    .figure-shell { overflow: auto; }",
            "    .figure-shell img { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(22,50,34,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #16543a; }",
            "    iframe { width: 100%; min-height: 760px; border: 1px solid rgba(22,50,34,0.12); border-radius: 14px; background: white; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Rabies Host and Geography Workflow</h1>",
            "  <p>Complete end-to-end review for one real rabies nucleoprotein panel. The workflow starts from raw sequences plus combined host and geography metadata, validates the FASTA surface, aligns and trims the panel, infers a bootstrap-supported maximum-likelihood tree, roots that tree on one explicit outgroup, summarizes bootstrap topology uncertainty, extracts clades, reconstructs host and geographic histories, and fits one comparative model over a derived geographic trait.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">sequences</span><strong>{report.dataset.sequence_count}</strong></div>',
            f'    <div class="card"><span class="label">selected model</span><strong>{escape(report.fasta_to_tree.selected_model)}</strong></div>',
            f'    <div class="card"><span class="label">aligned quality</span><strong>{_format_number(report.aligned_quality.quality_score)}</strong></div>',
            f'    <div class="card"><span class="label">trimmed quality</span><strong>{_format_number(report.trimmed_quality.quality_score)}</strong></div>',
            f'    <div class="card"><span class="label">root host</span><strong>{escape(host_summary.root_host)}</strong></div>',
            f'    <div class="card"><span class="label">root region</span><strong>{escape(geography_summary.root_region)}</strong></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Scientific Question</h2>",
            f"    <p>{escape(core_question)}</p>",
            f"    <p><strong>Working answer:</strong> {escape(core_answer)}. The comparative layer selects {escape(comparative_summary_row.selected_model)} as the better continuous-trait surface, but the residual diagnostics remain cautionary and the panel is intentionally small.</p>",
            _html_list(
                [
                    f"FASTA validation resolved the raw sequence type as {report.fasta_to_tree.sequence_type}.",
                    f"Bootstrap support spans {_support_range_text(support_summary.minimum_support, support_summary.maximum_support)} across the final rooted tree.",
                    f"Host reconstruction inferred {host_summary.host_switch_count} host-switch branches, with {host_summary.certain_host_switch_count} certain and {host_summary.uncertain_host_switch_count} uncertain changes.",
                    f"Geographic reconstruction inferred {migration_summary.event_count} migration events across {geography_summary.changed_branch_count} changed branches.",
                    f"Bootstrap replicate review retained {bootstrap_summary.tree_count} trees across {bootstrap_summary.diversity.rooted_topology_count} rooted topologies.",
                    (
                        "The rooted ML tree versus bootstrap consensus comparison "
                        f"returned rooted RF distance {bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance}."
                    ),
                    f"The clade table contains {clade_row_count} node rows and the comparative tree required {comparative_branch_repair_count} explicit branch-length repair(s).",
                    (
                        "Bootstrap review emitted budget warnings: "
                        + "; ".join(bootstrap_artifacts.budget_report.warning_messages)
                    )
                    if bootstrap_artifacts.budget_report.warning_messages
                    else (
                        "Configured workflow budgets covered the bootstrap review "
                        "without tree-count failure or peak-memory warning."
                    ),
                ]
            ),
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Sequence-to-Tree Outputs</h2>",
            _html_list(
                [
                    "input validation: input-validation.tsv",
                    "alignment quality: alignment-quality.tsv",
                    "alignment sequence ranking: alignment-sequence-ranking.tsv",
                    "alignment: rabies-cross-host-geography-panel.aln",
                    "trimmed alignment: rabies-cross-host-geography-panel.trimmed.aln",
                    "rooted tree: rabies-cross-host-geography-panel.rooted.tree",
                    "support table: rabies-cross-host-geography-panel.support.tsv",
                    "workflow summary: workflow-summary.tsv",
                    "resource observations: resource-observations.tsv",
                ]
            ),
            _support_table(
                report.fasta_to_tree,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Bootstrap and Clade Review</h2>",
            _html_list(
                [
                    f"bootstrap tree count: {bootstrap_summary.tree_count}",
                    f"rooted topology count: {bootstrap_summary.diversity.rooted_topology_count}",
                    f"unstable branch count: {bootstrap_summary.unstable_branch_count}",
                    f"clade row count: {clade_row_count}",
                    (
                        "see bootstrap-review/ for consensus, clade frequencies, "
                        "instability, distances, topology clusters, and rooted-tree comparison"
                    ),
                ]
            ),
            _table(
                headers=[
                    "tree_count",
                    "rooted_topology_count",
                    "dominant_topology_frequency",
                    "effective_topology_count",
                    "unstable_branch_count",
                ],
                rows=[
                    [
                        str(bootstrap_summary.tree_count),
                        str(bootstrap_summary.diversity.rooted_topology_count),
                        _format_number(
                            bootstrap_summary.diversity.dominant_topology_frequency
                        ),
                        _format_number(
                            bootstrap_summary.diversity.effective_topology_count
                        ),
                        str(bootstrap_summary.unstable_branch_count),
                    ]
                ],
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Host Switching</h2>",
            _html_list(
                [
                    f"workflow trait: {report.dataset.host_trait}",
                    f"root host confidence: {_format_number(host_summary.root_confidence)}",
                    f"host-switch rows: {len(report.host_switching.count_rows)}",
                    "see host-switch-summary.tsv, host-state-nodes.tsv, host-switch-branches.tsv, and host-switch-counts.tsv",
                ]
            ),
            _host_count_table(
                report.host_switching,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Comparative Layer</h2>",
            _html_list(
                [
                    f"formula: {comparative_summary_row.formula}",
                    f"selected trait model: {comparative_summary_row.selected_model}",
                    f"PGLS lambda: {_format_number(comparative_summary_row.pgls_lambda)}",
                    f"PGLS r-squared: {_format_number(comparative_summary_row.pgls_r_squared)}",
                    "see comparative/ for coefficients, model comparison, diagnostics, signal summary, and interpretation tables",
                ]
            ),
            _table(
                headers=["topic", "claim", "evidence"],
                rows=[
                    [row.topic, row.claim, row.evidence]
                    for row in comparative_interpretation_rows[:5]
                ],
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Biogeography</h2>",
            '      <p>The bundle includes the detailed biogeography package at <a href="biogeography/biogeography-report.html">biogeography/biogeography-report.html</a> together with the ancestral-region tree SVG and the self-contained geographic map.</p>',
            '      <div class="grid">',
            '        <div class="panel">',
            "          <h2>Ancestral-Region Tree</h2>",
            '          <div class="figure-shell">',
            '            <img src="biogeography/ancestral-region-tree.svg" alt="Ancestral region tree">',
            "          </div>",
            "        </div>",
            '        <div class="panel">',
            "          <h2>Geographic Map</h2>",
            '          <iframe src="biogeography/geographic-region-map.html" title="Geographic region map"></iframe>',
            "        </div>",
            "      </div>",
            _migration_event_table(
                report.biogeography_report,
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            '    <section class="panel full">',
            "      <h2>Scientific Findings Ledger</h2>",
            _table(
                headers=[
                    "finding_id",
                    "question",
                    "claim",
                    "evidence",
                    "caution",
                    "source_artifact",
                ],
                rows=[
                    [
                        row.finding_id,
                        row.question,
                        row.claim,
                        row.evidence,
                        row.caution,
                        row.source_artifact,
                    ]
                    for row in scientific_finding_rows
                ],
                max_rows=max_report_table_rows,
            ),
            "    </section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Key Files</h2>",
            _html_list(
                [
                    f'<a href="{workflow_summary_path.name}">{workflow_summary_path.name}</a>',
                    '<a href="resource-observations.tsv">resource-observations.tsv</a>',
                    '<a href="workflow-config-audit.tsv">workflow-config-audit.tsv</a>',
                    '<a href="clade-table.tsv">clade-table.tsv</a>',
                    '<a href="bootstrap-review/bootstrap-review.summary.tsv">bootstrap-review/bootstrap-review.summary.tsv</a>',
                    '<a href="bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv">bootstrap-review/rooted-tree-vs-bootstrap-consensus.summary.tsv</a>',
                    '<a href="comparative/comparative-report.html">comparative/comparative-report.html</a>',
                    '<a href="comparative/interpretation-table.tsv">comparative/interpretation-table.tsv</a>',
                    '<a href="scientific-findings.tsv">scientific-findings.tsv</a>',
                    '<a href="host-switch-summary.tsv">host-switch-summary.tsv</a>',
                    '<a href="biogeography/event-table.tsv">biogeography/event-table.tsv</a>',
                    '<a href="rabies-cross-host-geography.manifest.json">rabies-cross-host-geography.manifest.json</a>',
                ]
            ),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    path.write_text(html + "\n", encoding="utf-8")
    return path


def _support_table(
    report: FastaToTreeWorkflowReport,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=["node", "descendant_taxa", "support", "support_fraction"],
        rows=[
            [
                row.node,
                ", ".join(row.descendant_taxa),
                _format_number(row.support),
                _format_number(row.support_fraction),
            ]
            for row in report.support_rows
        ],
        max_rows=max_rows,
    )


def _host_count_table(
    report: HostSwitchingReport,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=[
            "transition",
            "certain_switch_count",
            "uncertain_switch_count",
            "total_switch_count",
        ],
        rows=[
            [
                row.transition,
                str(row.certain_switch_count),
                str(row.uncertain_switch_count),
                str(row.total_switch_count),
            ]
            for row in report.count_rows
        ],
        max_rows=max_rows,
    )


def _migration_event_table(
    report: BiogeographyReportPackageResult,
    *,
    max_rows: int | None = None,
) -> str:
    return _table(
        headers=[
            "branch_id",
            "source_region",
            "target_region",
            "support",
            "midpoint_depth",
        ],
        rows=[
            [
                row.branch_id,
                row.source_region,
                row.target_region,
                _format_number(row.support),
                _format_number(row.midpoint_depth),
            ]
            for row in report.event_report.event_rows
        ],
        max_rows=max_rows,
    )


def _html_list(items: list[str]) -> str:
    return "<ul>" + "".join(f"<li>{item}</li>" for item in items) + "</ul>"


def _table(
    headers: list[str],
    rows: list[list[str]],
    *,
    max_rows: int | None = None,
) -> str:
    rendered_rows = rows if max_rows is None else rows[:max_rows]
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rendered_rows
    )
    truncation_note = ""
    if max_rows is not None and len(rows) > max_rows:
        truncation_note = (
            f"<p><em>Showing the first {max_rows} of {len(rows)} rows. "
            "Use the linked TSV artifacts for the full table.</em></p>"
        )
    return (
        truncation_note
        + f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    )


def _support_range_text(
    minimum_support: float | None,
    maximum_support: float | None,
) -> str:
    if minimum_support is None or maximum_support is None:
        return "not available"
    return f"{_format_number(minimum_support)}-{_format_number(maximum_support)}"


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
