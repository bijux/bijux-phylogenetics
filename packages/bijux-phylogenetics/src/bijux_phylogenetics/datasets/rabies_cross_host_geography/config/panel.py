# ruff: noqa: F401, F403, F405
from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from hashlib import sha256
from html import escape
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
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    PGLSCategoricalContrastReport,
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    PosteriorTreePGLSReport,
    run_posterior_tree_pgls,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
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
from bijux_phylogenetics.compare.presentation import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.phylo.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.topology import (
    TreeRootingReport,
    root_tree_on_outgroup,
    write_tree_rooting_report,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
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
from bijux_phylogenetics.engines.inference import (
    FastaToTreeWorkflowReport,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.ecology import (
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

from ..models import (
    _SOURCE_ACCESSIONS,
    _WORKFLOW_CONFIG_NAME,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelExportResult,
)
from .audit import (
    _build_workflow_config_audit_rows,
    _raise_for_failed_config_audit,
    _read_observed_groups,
)
from .workflow_config import _load_workflow_config


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


