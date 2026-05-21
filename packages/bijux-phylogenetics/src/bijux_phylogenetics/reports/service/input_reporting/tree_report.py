from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    forensic_tree_path,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.render.html import write_html_report

from ..artifacts import report_sidecar_path, section, write_machine_manifest
from ..ledger import build_input_ledger, serialize_input_ledger
from ..models import ReportBuildResult
from ..summary import build_machine_manifest, report_summary_and_limitations


def render_tree_report(*, tree_path: Path, out_path: Path) -> ReportBuildResult:
    """Build the explicit single-tree report contract."""
    validation = validate_tree_path(tree_path)
    inspection = inspect_tree_path(tree_path)
    forensic = forensic_tree_path(tree_path)
    title = "Bijux Tree Report"
    reviewer_summary, limitations = report_summary_and_limitations(
        report_kind="tree",
        validation=validation,
        inspection=inspection,
        forensic=forensic,
    )
    sections = [
        section("reviewer-summary", reviewer_summary),
        section("tree-validation", asdict(validation)),
        section("tree-inspection", asdict(inspection)),
        section("tree-forensic", asdict(forensic)),
        section("limitations", limitations),
    ]
    machine_manifest = build_machine_manifest(
        report_kind="tree",
        title=title,
        input_paths=[tree_path],
        sections=sections,
        inspection=inspection,
    )
    machine_manifest["reviewer_summary"] = reviewer_summary
    machine_manifest["limitations"] = limitations
    input_ledger = build_input_ledger(
        [(tree_path, "tree", ["tree_validation", "tree_inspection", "tree_forensic"])]
    )
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
        report_kind="tree",
        title=title,
        validation=validation,
        inspection=inspection,
        forensic=forensic,
        metadata_linkage=None,
        traits_linkage=None,
        trait_missing_values=None,
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
        dataset_readiness=None,
        dataset_audit=None,
        input_ledger=input_ledger,
        machine_manifest=machine_manifest,
    )
