from __future__ import annotations

from pathlib import Path

from .contracts import (
    RabiesMethodSensitivitySlurmMergeCheckRow,
    RabiesMethodSensitivitySlurmMergeReport,
)
from .inputs import load_slurm_merge_inputs
from .merge_policy import (
    build_merge_variant_rows,
    collect_selected_models,
    count_preprocessing_change_pairs,
    count_rooted_engine_change_variants,
    count_serious_conflict_variants,
    maximum_serious_conflicts,
)


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
    preprocessing_change_pair_count = count_preprocessing_change_pairs(
        preprocessing_rows
    )
    serious_conflict_variant_count = count_serious_conflict_variants(
        variant_summary_rows
    )
    rooted_engine_change_variant_count = count_rooted_engine_change_variants(
        variant_summary_rows
    )
    maximum_serious_conflict_count = maximum_serious_conflicts(variant_summary_rows)
    selected_models = collect_selected_models(variant_summary_rows)

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

    variant_rows = build_merge_variant_rows(
        bundle_root=bundle_root,
        configured_variant_ids=configured_variant_ids,
        variant_summary_by_variant=variant_summary_by_variant,
        job_status_by_variant=job_status_by_variant,
        freshness_by_variant=freshness_by_variant,
        job_evidence_by_variant=job_evidence_by_variant,
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
