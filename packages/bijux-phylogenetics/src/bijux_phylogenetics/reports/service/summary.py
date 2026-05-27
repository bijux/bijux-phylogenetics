from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.dataset import DatasetAuditReport
from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
)
from bijux_phylogenetics.phylo.alignment import AlignmentForensicReport

from .ledger import sha256


def build_machine_manifest(
    *,
    report_kind: str,
    title: str,
    input_paths: list[Path],
    sections: list[tuple[str, str]],
    inspection: TreeInspectionReport,
) -> dict[str, object]:
    return {
        "report_kind": report_kind,
        "title": title,
        "input_paths": [str(path) for path in input_paths],
        "input_checksums": {str(path): sha256(path) for path in input_paths},
        "sections": [name for name, _ in sections],
        "metrics": {
            "tip_count": inspection.tip_count,
            "node_count": inspection.node_count,
            "clade_count": inspection.clade_count,
        },
    }


def report_summary_and_limitations(
    *,
    report_kind: str,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
    dataset_audit: DatasetAuditReport | None = None,
    alignment_forensic: AlignmentForensicReport | None = None,
) -> tuple[list[str], list[str]]:
    summary = [
        f"tree validity decision: {validation.validity_decision}",
        f"tree quality score: {inspection.tree_quality_score}",
        (
            "tree is currently safe for publication-oriented use"
            if forensic.safe_for_publication
            else "tree still carries publication-facing risks that should be reviewed"
        ),
    ]
    limitations = list(forensic.warnings)
    if report_kind == "dataset" and dataset_audit is not None:
        summary.append(
            f"dataset readiness decision: {dataset_audit.readiness_decision}"
        )
        summary.append(
            f"blocked analyses: {len(dataset_audit.blocked_analyses)}, risky analyses: "
            f"{sum(1 for row in dataset_audit.analysis_decisions if row.decision == 'risky')}"
        )
        limitations.extend(
            finding.message
            for finding in dataset_audit.findings
            if finding.severity in {"warning", "blocker"}
        )
    if report_kind == "phylo-inputs" and alignment_forensic is not None:
        safe_methods = sum(
            1
            for flag in (
                alignment_forensic.safe_for_distance_analysis,
                alignment_forensic.safe_for_maximum_likelihood,
                alignment_forensic.safe_for_bayesian_inference,
                alignment_forensic.safe_for_coding_analysis,
                alignment_forensic.safe_for_publication,
            )
            if flag
        )
        summary.append(f"alignment safe-for flags passed: {safe_methods}/5")
        summary.append(
            "alignment suspicious diagnostics: flagged"
            if alignment_forensic.quality.suspicious_alignment
            else "alignment suspicious diagnostics: clear"
        )
        limitations.extend(alignment_forensic.limitations)
    if inspection.mixed_support_scales:
        limitations.append(
            "support labels originate from mixed scales and should be interpreted only after normalization audit"
        )
    return summary, sorted(set(limitations))


def distance_method_limitations() -> list[str]:
    return [
        "distance methods collapse site-by-site sequence evidence into pairwise summaries before tree building",
        "different evolutionary histories can yield similar pairwise distances, so topology is not uniquely identified by the matrix alone",
        "UPGMA additionally assumes an ultrametric clock-like process and can misplace taxa when rates vary across lineages",
        "Neighbor-Joining is often useful for quick structure, but it is still a summary approximation rather than a full likelihood inference",
        "BIONJ remains a distance-summary method, but it still approximates full-sequence evidence and can stabilize noisy reductions relative to classic Neighbor-Joining by using variance-aware joining",
    ]
