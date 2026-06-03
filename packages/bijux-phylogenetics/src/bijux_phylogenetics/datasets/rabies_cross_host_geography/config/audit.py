from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.fasta import load_permissive_fasta_records

from ..models import (
    RabiesCrossHostGeographyPanelWorkflowConfig,
    RabiesWorkflowConfigAuditRow,
)
from ..shared import _format_number


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
