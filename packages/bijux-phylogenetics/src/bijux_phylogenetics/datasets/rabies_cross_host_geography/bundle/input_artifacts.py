from __future__ import annotations

import json
from pathlib import Path
import shutil

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.engines.inference import FastaToTreeWorkflowReport
from bijux_phylogenetics.phylo.alignment import (
    AlignmentQualityReport,
    SequenceQualityRankingReport,
)

from ..models import (
    RabiesComparativeBranchRepair,
    RabiesCrossHostGeographyPanelWorkflowConfig,
    RabiesWorkflowConfigAuditRow,
)
from ..shared import _checksum, _format_number


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
            "memory_warning_threshold_bytes": config.memory_warning_threshold_bytes,
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
            "comparative_branch_length_floor": config.comparative_branch_length_floor,
        },
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


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
