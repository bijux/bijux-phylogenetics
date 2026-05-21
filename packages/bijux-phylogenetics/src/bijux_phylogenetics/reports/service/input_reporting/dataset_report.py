from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.core.dataset import (
    audit_dataset_inputs,
    summarize_dataset_readiness,
)
from bijux_phylogenetics.datasets.study_inputs import detect_missing_trait_values
from bijux_phylogenetics.diagnostics.validation import (
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import build_input_ledger, serialize_input_ledger
from ..linkage import annotate_tree_against_table
from ..models import ReportBuildResult
from ..summary import build_machine_manifest, report_summary_and_limitations


def render_dataset_report(
    *,
    tree_path: Path,
    metadata_path: Path,
    out_path: Path,
    traits_path: Path | None = None,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> ReportBuildResult:
    """Build the explicit tree plus table dataset report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    metadata_linkage = annotate_tree_against_table(tree_path, metadata_path)
    traits_linkage = (
        annotate_tree_against_table(tree_path, traits_path)
        if traits_path is not None
        else None
    )
    trait_missing_values = (
        detect_missing_trait_values(traits_path) if traits_path is not None else None
    )
    dataset_readiness = (
        summarize_dataset_readiness(tree_path, metadata_path, traits_path)
        if traits_path is not None
        else None
    )
    dataset_audit = (
        audit_dataset_inputs(
            tree_path,
            metadata_path,
            traits_path,
            alignment_path=alignment_path,
            tip_dates_path=tip_dates_path,
            calibration_path=calibration_path,
        )
        if traits_path is not None
        else None
    )
    title = "Bijux Dataset Report"
    reviewer_summary, limitations = report_summary_and_limitations(
        report_kind="dataset",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        dataset_audit=dataset_audit,
    )
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("tree-validation", asdict(validation)),
        section("tree-inspection", asdict(inspection)),
        section("tree-forensic", asdict(forensic)),
        section("metadata-linkage", asdict(metadata_linkage)),
    ]
    if traits_linkage is not None:
        sections.append(section("traits-linkage", asdict(traits_linkage)))
    if trait_missing_values is not None:
        sections.append(section("trait-missing-values", asdict(trait_missing_values)))
    if dataset_readiness is not None:
        sections.append(section("dataset-readiness", asdict(dataset_readiness)))
    if dataset_audit is not None:
        sections.append(section("dataset-audit", asdict(dataset_audit)))
        sections.append(
            section("dataset-findings", [asdict(row) for row in dataset_audit.findings])
        )
        sections.append(
            section(
                "dataset-analysis-decisions",
                [asdict(row) for row in dataset_audit.analysis_decisions],
            )
        )
        sections.append(
            section(
                "dataset-readiness-levels",
                [asdict(row) for row in dataset_audit.readiness_levels],
            )
        )
        sections.append(section("dataset-crosswalk", asdict(dataset_audit.crosswalk)))
        sections.append(
            section("dataset-completeness", asdict(dataset_audit.completeness_matrix))
        )
        sections.append(
            section("dataset-exclusions", asdict(dataset_audit.exclusion_table))
        )
        sections.append(
            section("dataset-mismatches", asdict(dataset_audit.mismatch_report))
        )
        sections.append(section("dataset-risk-score", asdict(dataset_audit.risk_score)))
        sections.append(
            section("dataset-minimal-fix-plan", asdict(dataset_audit.minimal_fix_plan))
        )
        sections.append(
            section(
                "dataset-reviewer-checklist", asdict(dataset_audit.reviewer_checklist)
            )
        )
        sections.append(
            section("dataset-ordering", asdict(dataset_audit.ordering_audit))
        )
        sections.append(
            section(
                "dataset-pruning", [asdict(row) for row in dataset_audit.pruning_steps]
            )
        )
        sections.append(
            section(
                "dataset-group-imbalance",
                [asdict(row) for row in dataset_audit.group_imbalance_warnings],
            )
        )
    input_ledger_entries: list[tuple[Path, str, list[str]]] = [
        (
            tree_path,
            "tree",
            ["tree_validation", "tree_inspection", "tree_forensic", "dataset_audit"],
        ),
        (metadata_path, "metadata", ["metadata_linkage", "dataset_audit"]),
    ]
    if traits_path is not None:
        input_ledger_entries.append(
            (
                traits_path,
                "traits",
                ["traits_linkage", "trait_missing_values", "dataset_audit"],
            )
        )
    if alignment_path is not None:
        input_ledger_entries.append(
            (alignment_path, "alignment", ["alignment_forensic", "dataset_audit"])
        )
    if tip_dates_path is not None:
        input_ledger_entries.append((tip_dates_path, "tip_dates", ["dataset_audit"]))
    if calibration_path is not None:
        input_ledger_entries.append(
            (calibration_path, "calibrations", ["dataset_audit"])
        )
    input_ledger = build_input_ledger(input_ledger_entries)
    sections.append(
        section("dataset-input-ledger", serialize_input_ledger(input_ledger))
    )
    sections.append(section("limitations", limitations))
    input_paths = [tree_path, metadata_path]
    if traits_path is not None:
        input_paths.append(traits_path)
    if alignment_path is not None:
        input_paths.append(alignment_path)
    if tip_dates_path is not None:
        input_paths.append(tip_dates_path)
    if calibration_path is not None:
        input_paths.append(calibration_path)
    machine_manifest = build_machine_manifest(
        report_kind="dataset",
        title=title,
        input_paths=input_paths,
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    machine_manifest["input_ledger"] = serialize_input_ledger(input_ledger)
    machine_manifest_path = write_machine_manifest(
        report_sidecar_path(out_path), machine_manifest
    )
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return ReportBuildResult(
        output_path=out_path,
        machine_manifest_path=machine_manifest_path,
        report_kind="dataset",
        title=title,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        metadata_linkage=metadata_linkage,
        traits_linkage=traits_linkage,
        trait_missing_values=trait_missing_values,
        alignment=None,
        alignment_quality=None,
        alignment_forensic=None,
        alignment_low_information=None,
        alignment_duplicate_policy=None,
        alignment_ambiguous_columns=None,
        alignment_sequence_ranking=None,
        alignment_coding=None,
        alignment_identity_matrix=None,
        alignment_linkage=None,
        dataset_readiness=dataset_readiness,
        dataset_audit=dataset_audit,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )
