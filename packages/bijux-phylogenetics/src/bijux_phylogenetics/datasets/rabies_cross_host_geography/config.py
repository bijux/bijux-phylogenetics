# ruff: noqa: F401, F403, F405
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
from bijux_phylogenetics.engines.fasta_to_tree import (
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
from bijux_phylogenetics.io.fasta import load_permissive_fasta_records
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

from .models import (
    _ALIGNMENT_MODE,
    _BOOTSTRAP_CONSENSUS_THRESHOLD,
    _BOOTSTRAP_REPLICATES,
    _BOOTSTRAP_ROBUST_SUPPORT_THRESHOLD,
    _CLADE_METADATA_COLUMNS,
    _COMPARATIVE_BRANCH_LENGTH_FLOOR,
    _COMPARATIVE_FORMULA,
    _COMPARATIVE_RESPONSE,
    _DATASET_ID,
    _DATASET_LABEL,
    _GEOGRAPHY_MODEL,
    _GEOGRAPHY_TRAIT,
    _HOST_MODEL,
    _HOST_TRAIT,
    _IQTREE_SEED,
    _IQTREE_THREADS,
    _MAX_BOOTSTRAP_TREE_COUNT,
    _MAX_REPORT_TABLE_ROWS,
    _MEMORY_WARNING_THRESHOLD_BYTES,
    _OUTGROUP_TAXA,
    _SEQUENCE_TYPE,
    _SOURCE_ACCESSIONS,
    _TRIMMING_MODE,
    _TRIM_GAP_THRESHOLD,
    _WORKFLOW_CONFIG_NAME,
    _WORKFLOW_PREFIX,
    _WORKFLOW_TIMEOUT_SECONDS,
    RabiesCrossHostGeographyPanelDataset,
    RabiesCrossHostGeographyPanelExportResult,
    RabiesCrossHostGeographyPanelWorkflowConfig,
    RabiesWorkflowConfigAuditRow,
)
from .shared import _format_number

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



def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent.parent
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
