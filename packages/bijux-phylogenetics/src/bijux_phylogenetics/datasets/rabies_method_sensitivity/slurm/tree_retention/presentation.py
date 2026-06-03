from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import RabiesMethodSensitivitySlurmTreeRetentionReport


def write_rabies_method_sensitivity_slurm_tree_retention_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the reviewer-facing tree-retention policy report."""
    return write_html_report(
        title="Rabies Slurm Tree Retention Report",
        sections=[
            (
                "policy-summary",
                "\n".join(
                    [
                        f"overall_policy_status: {report.overall_policy_status}",
                        f"variant_count: {report.variant_count}",
                        f"file_count: {report.file_count}",
                        f"tree_set_file_count: {report.tree_set_file_count}",
                        (
                            "thinning_required_file_count: "
                            f"{report.thinning_required_file_count}"
                        ),
                        (
                            "compression_required_file_count: "
                            f"{report.compression_required_file_count}"
                        ),
                    ]
                ),
            ),
            (
                "global-issues",
                "none"
                if report.global_issue_count == 0
                else "\n".join(report.global_issues),
            ),
            (
                "required-actions",
                "none"
                if report.overall_policy_status != "required"
                else "\n".join(
                    f"{row.relative_path}: {'; '.join(row.issues)}"
                    for row in report.files
                    if row.thinning_policy == "thin_required"
                    or row.compression_policy == "compress_required"
                ),
            ),
            (
                "recommended-actions",
                "none"
                if report.overall_policy_status == "no_action"
                else "\n".join(
                    f"{row.relative_path}: {'; '.join(row.issues)}"
                    for row in report.files
                    if row.issue_count > 0
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "overall_policy_status": report.overall_policy_status,
            "variant_count": report.variant_count,
            "file_count": report.file_count,
            "tree_set_file_count": report.tree_set_file_count,
            "posterior_sample_file_count": report.posterior_sample_file_count,
            "thinning_recommended_file_count": report.thinning_recommended_file_count,
            "thinning_required_file_count": report.thinning_required_file_count,
            "compression_recommended_file_count": report.compression_recommended_file_count,
            "compression_required_file_count": report.compression_required_file_count,
            "total_tree_count": report.total_tree_count,
            "total_tree_byte_count": report.total_tree_byte_count,
            "largest_tree_set_path": report.largest_tree_set_path,
            "largest_tree_set_tree_count": report.largest_tree_set_tree_count,
        },
        summary_metrics=[
            ("overall status", report.overall_policy_status),
            ("tree-set files", report.tree_set_file_count),
            ("posterior files", report.posterior_sample_file_count),
            ("thinning required", report.thinning_required_file_count),
            ("compression required", report.compression_required_file_count),
            ("total tree count", report.total_tree_count),
            ("total tree bytes", report.total_tree_byte_count),
        ],
        artifact_links=[
            ("tree retention checks", "slurm-tree-retention-checks.tsv", None),
            ("tree retention files", "slurm-tree-retention-files.tsv", None),
            ("tree retention summary", "slurm-tree-retention-policy.json", None),
        ],
    )
