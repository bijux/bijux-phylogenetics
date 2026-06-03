from __future__ import annotations

from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmMergeVariantRow


def count_preprocessing_change_pairs(
    preprocessing_rows: list[dict[str, str]],
) -> int:
    return sum(
        1
        for row in preprocessing_rows
        if int(row["robinson_foulds_distance"]) > 0
        or str(row["same_taxa_different_rooting"]) == "true"
    )


def count_serious_conflict_variants(
    variant_summary_rows: list[dict[str, str]],
) -> int:
    return sum(
        1 for row in variant_summary_rows if int(row["serious_conflict_count"]) > 0
    )


def count_rooted_engine_change_variants(
    variant_summary_rows: list[dict[str, str]],
) -> int:
    return sum(
        1
        for row in variant_summary_rows
        if int(row["rooted_engine_rf_distance"]) > 0
        or str(row["rooted_engine_same_taxa_different_rooting"]) == "true"
    )


def maximum_serious_conflicts(
    variant_summary_rows: list[dict[str, str]],
) -> int:
    return max(int(row["serious_conflict_count"]) for row in variant_summary_rows)


def collect_selected_models(
    variant_summary_rows: list[dict[str, str]],
) -> tuple[str, ...]:
    return tuple(sorted({str(row["selected_model"]) for row in variant_summary_rows}))


def build_merge_variant_rows(
    *,
    bundle_root: Path,
    configured_variant_ids: list[str],
    variant_summary_by_variant: dict[str, dict[str, str]],
    job_status_by_variant: dict[str, dict[str, str]],
    freshness_by_variant: dict[str, dict[str, str]],
    job_evidence_by_variant: dict[str, dict[str, str]],
) -> list[RabiesMethodSensitivitySlurmMergeVariantRow]:
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
            evidence_json_path = bundle_root / str(
                job_evidence_row["evidence_json_path"]
            )
            evidence_html_path = bundle_root / str(
                job_evidence_row["evidence_html_path"]
            )
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
    return variant_rows
