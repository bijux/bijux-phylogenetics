from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    RabiesMethodSensitivitySlurmMergeCheckRow,
    RabiesMethodSensitivitySlurmMergeReport,
    RabiesMethodSensitivitySlurmMergeVariantRow,
)
from .inputs import load_slurm_merge_inputs
from .shared import (
    _write_tsv,
)

__all__ = [
    "RabiesMethodSensitivitySlurmMergeCheckRow",
    "RabiesMethodSensitivitySlurmMergeReport",
    "RabiesMethodSensitivitySlurmMergeVariantRow",
    "build_rabies_method_sensitivity_slurm_merge_report",
    "write_rabies_method_sensitivity_slurm_merge_checks_table",
    "write_rabies_method_sensitivity_slurm_merge_html_report",
    "write_rabies_method_sensitivity_slurm_merge_summary_json",
    "write_rabies_method_sensitivity_slurm_merge_variants_table",
]

def build_rabies_method_sensitivity_slurm_merge_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmMergeReport:
    """Assess whether distributed batch outputs merge into one coherent result."""
    loaded_inputs = load_slurm_merge_inputs(bundle_root)
    bundle_root = loaded_inputs.bundle_root
    config = loaded_inputs.config
    checks = list(loaded_inputs.checks)
    workflow_summary = loaded_inputs.workflow_summary
    configured_variant_ids = loaded_inputs.configured_variant_ids
    variant_summary_rows = loaded_inputs.variant_summary_rows
    preprocessing_rows = loaded_inputs.preprocessing_rows
    stable_clade_rows = loaded_inputs.stable_clade_rows
    changed_clade_rows = loaded_inputs.changed_clade_rows
    conclusion_rows = loaded_inputs.conclusion_rows
    variant_summary_by_variant = loaded_inputs.variant_summary_by_variant
    job_status_by_variant = loaded_inputs.job_status_by_variant
    freshness_by_variant = loaded_inputs.freshness_by_variant
    job_evidence_by_variant = loaded_inputs.job_evidence_by_variant

    def add_check(
        check_id: str,
        *,
        surface: str,
        condition: bool,
        expected: object,
        observed: object,
        detail: str,
    ) -> None:
        checks.append(
            RabiesMethodSensitivitySlurmMergeCheckRow(
                check_id=check_id,
                surface=surface,
                status="passed" if condition else "failed",
                expected="" if expected is None else str(expected),
                observed="" if observed is None else str(observed),
                detail=detail,
            )
        )

    stable_clade_count = len(stable_clade_rows)
    changed_clade_count = len(changed_clade_rows)
    preprocessing_change_pair_count = sum(
        1
        for row in preprocessing_rows
        if int(row["robinson_foulds_distance"]) > 0
        or str(row["same_taxa_different_rooting"]) == "true"
    )
    serious_conflict_variant_count = sum(
        1
        for row in variant_summary_rows
        if int(row["serious_conflict_count"]) > 0
    )
    rooted_engine_change_variant_count = sum(
        1
        for row in variant_summary_rows
        if int(row["rooted_engine_rf_distance"]) > 0
        or str(row["rooted_engine_same_taxa_different_rooting"]) == "true"
    )
    maximum_serious_conflict_count = max(
        int(row["serious_conflict_count"]) for row in variant_summary_rows
    )
    selected_models = tuple(
        sorted({str(row["selected_model"]) for row in variant_summary_rows})
    )

    add_check(
        "workflow-summary:variant-count",
        surface="workflow-summary",
        condition=int(workflow_summary["variant_count"]) == len(variant_summary_rows),
        expected=workflow_summary["variant_count"],
        observed=len(variant_summary_rows),
        detail="workflow summary variant_count matches the merged variant summary rows",
    )
    add_check(
        "workflow-summary:stable-clades",
        surface="workflow-summary",
        condition=int(workflow_summary["stable_clade_count"]) == stable_clade_count,
        expected=workflow_summary["stable_clade_count"],
        observed=stable_clade_count,
        detail="workflow summary stable_clade_count matches the merged stable-clade ledger",
    )
    add_check(
        "workflow-summary:changed-clades",
        surface="workflow-summary",
        condition=int(workflow_summary["changed_clade_count"]) == changed_clade_count,
        expected=workflow_summary["changed_clade_count"],
        observed=changed_clade_count,
        detail="workflow summary changed_clade_count matches the merged changed-clade ledger",
    )
    add_check(
        "workflow-summary:preprocessing-change-pairs",
        surface="workflow-summary",
        condition=int(workflow_summary["preprocessing_change_pair_count"])
        == preprocessing_change_pair_count,
        expected=workflow_summary["preprocessing_change_pair_count"],
        observed=preprocessing_change_pair_count,
        detail="workflow summary preprocessing_change_pair_count matches the merged preprocessing comparisons",
    )
    add_check(
        "workflow-summary:rooted-engine-change-variants",
        surface="workflow-summary",
        condition=int(workflow_summary["rooted_engine_change_variant_count"])
        == rooted_engine_change_variant_count,
        expected=workflow_summary["rooted_engine_change_variant_count"],
        observed=rooted_engine_change_variant_count,
        detail="workflow summary rooted_engine_change_variant_count matches merged variant summaries",
    )
    add_check(
        "workflow-summary:serious-conflict-variants",
        surface="workflow-summary",
        condition=int(workflow_summary["serious_conflict_variant_count"])
        == serious_conflict_variant_count,
        expected=workflow_summary["serious_conflict_variant_count"],
        observed=serious_conflict_variant_count,
        detail="workflow summary serious_conflict_variant_count matches merged variant summaries",
    )
    add_check(
        "workflow-summary:maximum-serious-conflicts",
        surface="workflow-summary",
        condition=int(workflow_summary["maximum_serious_conflict_count"])
        == maximum_serious_conflict_count,
        expected=workflow_summary["maximum_serious_conflict_count"],
        observed=maximum_serious_conflict_count,
        detail="workflow summary maximum_serious_conflict_count matches merged variant summaries",
    )

    variant_rows: list[RabiesMethodSensitivitySlurmMergeVariantRow] = []
    for variant_id in configured_variant_ids:
        variant_summary_row = variant_summary_by_variant.get(variant_id)
        job_status_row = job_status_by_variant.get(variant_id)
        freshness_row = freshness_by_variant.get(variant_id)
        job_evidence_row = job_evidence_by_variant.get(variant_id)
        issues: list[str] = []
        evidence_status = "present"
        if job_status_row is None:
            issues.append("job-status row is missing")
        elif str(job_status_row["status"]) != "completed":
            issues.append("job-status does not mark the variant as completed")
        if freshness_row is None:
            issues.append("output-freshness row is missing")
        elif str(freshness_row["freshness_status"]) != "fresh":
            issues.append("output-freshness does not mark the variant as fresh")
        if variant_summary_row is None:
            issues.append("variant summary row is missing")
        if job_evidence_row is None:
            evidence_status = "missing"
            issues.append("job-evidence row is missing")
        else:
            evidence_json_path = bundle_root / str(job_evidence_row["evidence_json_path"])
            evidence_html_path = bundle_root / str(job_evidence_row["evidence_html_path"])
            if not evidence_json_path.is_file():
                evidence_status = "missing"
                issues.append("job evidence JSON is missing")
            if not evidence_html_path.is_file():
                evidence_status = "missing"
                issues.append("job evidence HTML is missing")
        merge_status = "merged" if not issues else "blocked"
        variant_rows.append(
            RabiesMethodSensitivitySlurmMergeVariantRow(
                variant_id=variant_id,
                merge_status=merge_status,
                job_status=(
                    "missing"
                    if job_status_row is None
                    else str(job_status_row["status"])
                ),
                output_freshness_status=(
                    "missing"
                    if freshness_row is None
                    else str(freshness_row["freshness_status"])
                ),
                evidence_status=evidence_status,
                included_in_merge=not issues,
                selected_model=(
                    ""
                    if variant_summary_row is None
                    else str(variant_summary_row["selected_model"])
                ),
                serious_conflict_count=(
                    0
                    if variant_summary_row is None
                    else int(variant_summary_row["serious_conflict_count"])
                ),
                rooted_engine_rf_distance=(
                    0
                    if variant_summary_row is None
                    else int(variant_summary_row["rooted_engine_rf_distance"])
                ),
                rooted_engine_same_taxa_different_rooting=(
                    False
                    if variant_summary_row is None
                    else str(
                        variant_summary_row["rooted_engine_same_taxa_different_rooting"]
                    )
                    == "true"
                ),
                issue_count=len(issues),
                issues=tuple(issues),
                evidence_json_path=(
                    ""
                    if job_evidence_row is None
                    else str(job_evidence_row["evidence_json_path"])
                ),
                evidence_html_path=(
                    ""
                    if job_evidence_row is None
                    else str(job_evidence_row["evidence_html_path"])
                ),
            )
        )

    merged_variant_count = sum(1 for row in variant_rows if row.included_in_merge)
    add_check(
        "merge-readiness:variant-count",
        surface="merge-readiness",
        condition=merged_variant_count == len(configured_variant_ids),
        expected=len(configured_variant_ids),
        observed=merged_variant_count,
        detail="every configured variant must be completed, fresh, and evidenced before the global merge is ready",
    )
    failed_check_count = sum(1 for row in checks if row.status == "failed")
    failed_variant_count = sum(1 for row in variant_rows if not row.included_in_merge)
    merge_ready = failed_check_count == 0 and failed_variant_count == 0
    merge_status = "merge-ready" if merge_ready else "merge-blocked"
    return RabiesMethodSensitivitySlurmMergeReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        merge_status=merge_status,
        merge_ready=merge_ready,
        expected_variant_count=len(configured_variant_ids),
        merged_variant_count=merged_variant_count,
        mergeable_variant_count=merged_variant_count,
        failed_variant_count=failed_variant_count,
        failed_check_count=failed_check_count,
        check_count=len(checks),
        stable_clade_count=stable_clade_count,
        changed_clade_count=changed_clade_count,
        preprocessing_comparison_count=len(preprocessing_rows),
        conclusion_count=len(conclusion_rows),
        serious_conflict_variant_count=serious_conflict_variant_count,
        rooted_engine_change_variant_count=rooted_engine_change_variant_count,
        maximum_serious_conflict_count=maximum_serious_conflict_count,
        selected_models=selected_models,
        checks=tuple(checks),
        variants=tuple(variant_rows),
    )


def write_rabies_method_sensitivity_slurm_merge_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write one check-level merge ledger."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_slurm_merge_variants_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write one per-variant merge-decision ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "merge_status",
            "job_status",
            "output_freshness_status",
            "evidence_status",
            "included_in_merge",
            "selected_model",
            "serious_conflict_count",
            "rooted_engine_rf_distance",
            "rooted_engine_same_taxa_different_rooting",
            "issue_count",
            "issues",
            "evidence_json_path",
            "evidence_html_path",
        ),
        rows=[
            {
                "variant_id": row.variant_id,
                "merge_status": row.merge_status,
                "job_status": row.job_status,
                "output_freshness_status": row.output_freshness_status,
                "evidence_status": row.evidence_status,
                "included_in_merge": str(row.included_in_merge).lower(),
                "selected_model": row.selected_model,
                "serious_conflict_count": row.serious_conflict_count,
                "rooted_engine_rf_distance": row.rooted_engine_rf_distance,
                "rooted_engine_same_taxa_different_rooting": str(
                    row.rooted_engine_same_taxa_different_rooting
                ).lower(),
                "issue_count": row.issue_count,
                "issues": " | ".join(row.issues),
                "evidence_json_path": row.evidence_json_path,
                "evidence_html_path": row.evidence_html_path,
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_slurm_merge_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write the structured merge summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_merge_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write the reviewer-facing global merge report."""
    sections = [
        (
            "merge-status",
            "\n".join(
                [
                    f"merge_status: {report.merge_status}",
                    f"merge_ready: {str(report.merge_ready).lower()}",
                    f"expected_variant_count: {report.expected_variant_count}",
                    f"merged_variant_count: {report.merged_variant_count}",
                    f"failed_variant_count: {report.failed_variant_count}",
                ]
            ),
        ),
        (
            "merged-science-summary",
            "\n".join(
                [
                    f"stable_clade_count: {report.stable_clade_count}",
                    f"changed_clade_count: {report.changed_clade_count}",
                    f"preprocessing_comparison_count: {report.preprocessing_comparison_count}",
                    f"conclusion_count: {report.conclusion_count}",
                    f"serious_conflict_variant_count: {report.serious_conflict_variant_count}",
                    f"rooted_engine_change_variant_count: {report.rooted_engine_change_variant_count}",
                    f"maximum_serious_conflict_count: {report.maximum_serious_conflict_count}",
                    f"selected_models: {', '.join(report.selected_models)}",
                ]
            ),
        ),
        (
            "failed-checks",
            "none"
            if report.failed_check_count == 0
            else "\n".join(
                f"{row.check_id}: {row.detail}"
                for row in report.checks
                if row.status == "failed"
            ),
        ),
        (
            "blocked-variants",
            "none"
            if report.failed_variant_count == 0
            else "\n".join(
                f"{row.variant_id}: {'; '.join(row.issues)}"
                for row in report.variants
                if not row.included_in_merge
            ),
        ),
    ]
    return write_html_report(
        title="Rabies Slurm Merge Report",
        sections=sections,
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "merge_status": report.merge_status,
            "merge_ready": report.merge_ready,
            "expected_variant_count": report.expected_variant_count,
            "merged_variant_count": report.merged_variant_count,
            "failed_variant_count": report.failed_variant_count,
            "failed_check_count": report.failed_check_count,
            "selected_models": list(report.selected_models),
        },
        summary_metrics=[
            ("merge status", report.merge_status),
            ("merge ready", str(report.merge_ready).lower()),
            ("merged variants", report.merged_variant_count),
            ("failed variants", report.failed_variant_count),
            ("failed checks", report.failed_check_count),
            ("stable clades", report.stable_clade_count),
            ("changed clades", report.changed_clade_count),
        ],
        artifact_links=[
            ("merge checks", "slurm-merge-checks.tsv", None),
            ("merge variants", "slurm-merge-variants.tsv", None),
            ("merge summary", "slurm-merge-report.json", None),
        ],
    )
