from __future__ import annotations

from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmTreeRetentionReport
from .file_policy import build_tree_retention_file_row, iter_tree_relative_paths
from .inputs import _add_check, load_tree_retention_inputs


def build_rabies_method_sensitivity_slurm_tree_retention_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmTreeRetentionReport:
    """Derive safe thinning and compression policies for retained tree-bearing files."""
    loaded_inputs = load_tree_retention_inputs(bundle_root)
    bundle_root = loaded_inputs.bundle_root
    config = loaded_inputs.config
    configured_variant_ids = loaded_inputs.configured_variant_ids
    output_explosion_summary = loaded_inputs.output_explosion_summary
    storage_category_by_id = loaded_inputs.storage_category_by_id
    checks = list(loaded_inputs.checks)

    file_rows = tuple(
        build_tree_retention_file_row(
            bundle_root=bundle_root, relative_path=relative_path
        )
        for relative_path in iter_tree_relative_paths(bundle_root)
    )
    observed_variant_ids = sorted({row.variant_id for row in file_rows})
    _add_check(
        checks,
        "tree-files:variant-coverage",
        surface="tree-files",
        condition=observed_variant_ids == configured_variant_ids,
        expected=configured_variant_ids,
        observed=observed_variant_ids,
        detail="tree-bearing file rows cover every configured rabies workflow variant",
    )

    tree_artifact_file_count = sum(
        1 for row in file_rows if row.artifact_scope == "tree_artifact"
    )
    posterior_sample_file_count = sum(
        1 for row in file_rows if row.artifact_scope == "posterior_sample"
    )
    total_tree_byte_count = sum(row.byte_count for row in file_rows)
    total_tree_count = sum(row.tree_count for row in file_rows)
    tree_set_file_count = sum(1 for row in file_rows if row.tree_count > 1)
    thinning_recommended_file_count = sum(
        1 for row in file_rows if row.thinning_policy == "thin_recommended"
    )
    thinning_required_file_count = sum(
        1 for row in file_rows if row.thinning_policy == "thin_required"
    )
    compression_recommended_file_count = sum(
        1 for row in file_rows if row.compression_policy == "compress_recommended"
    )
    compression_required_file_count = sum(
        1 for row in file_rows if row.compression_policy == "compress_required"
    )
    _add_check(
        checks,
        "storage:tree-file-count",
        surface="storage",
        condition=int(storage_category_by_id["trees"]["total_file_count"])
        == tree_artifact_file_count,
        expected=storage_category_by_id["trees"]["total_file_count"],
        observed=tree_artifact_file_count,
        detail="storage tree file count matches the inspected tree-artifact files",
    )
    _add_check(
        checks,
        "storage:tree-byte-count",
        surface="storage",
        condition=int(storage_category_by_id["trees"]["total_byte_count"])
        == sum(
            row.byte_count for row in file_rows if row.artifact_scope == "tree_artifact"
        ),
        expected=storage_category_by_id["trees"]["total_byte_count"],
        observed=sum(
            row.byte_count for row in file_rows if row.artifact_scope == "tree_artifact"
        ),
        detail="storage tree byte count matches the inspected tree-artifact files",
    )
    _add_check(
        checks,
        "storage:posterior-file-count",
        surface="storage",
        condition=int(storage_category_by_id["posterior_samples"]["total_file_count"])
        == posterior_sample_file_count,
        expected=storage_category_by_id["posterior_samples"]["total_file_count"],
        observed=posterior_sample_file_count,
        detail="posterior-sample file count matches the inspected tree-bearing sample files",
    )
    _add_check(
        checks,
        "storage:posterior-byte-count",
        surface="storage",
        condition=int(storage_category_by_id["posterior_samples"]["total_byte_count"])
        == sum(
            row.byte_count
            for row in file_rows
            if row.artifact_scope == "posterior_sample"
        ),
        expected=storage_category_by_id["posterior_samples"]["total_byte_count"],
        observed=sum(
            row.byte_count
            for row in file_rows
            if row.artifact_scope == "posterior_sample"
        ),
        detail="posterior-sample byte count matches the inspected tree-bearing sample files",
    )
    _add_check(
        checks,
        "output-explosion:tree-bytes",
        surface="output-explosion",
        condition=int(output_explosion_summary["total_tree_byte_count"])
        == sum(
            row.byte_count for row in file_rows if row.artifact_scope == "tree_artifact"
        ),
        expected=output_explosion_summary["total_tree_byte_count"],
        observed=sum(
            row.byte_count for row in file_rows if row.artifact_scope == "tree_artifact"
        ),
        detail="output-explosion tree bytes match the tree-retention inspection surface",
    )
    _add_check(
        checks,
        "output-explosion:posterior-bytes",
        surface="output-explosion",
        condition=int(output_explosion_summary["total_posterior_sample_byte_count"])
        == sum(
            row.byte_count
            for row in file_rows
            if row.artifact_scope == "posterior_sample"
        ),
        expected=output_explosion_summary["total_posterior_sample_byte_count"],
        observed=sum(
            row.byte_count
            for row in file_rows
            if row.artifact_scope == "posterior_sample"
        ),
        detail="output-explosion posterior bytes match the tree-retention inspection surface",
    )

    global_issues: list[str] = []
    overall_policy_status = "no_action"
    if thinning_required_file_count > 0 or compression_required_file_count > 0:
        overall_policy_status = "required"
    elif thinning_recommended_file_count > 0 or compression_recommended_file_count > 0:
        overall_policy_status = "recommended"
    if tree_set_file_count == 0:
        global_issues.append(
            "no retained multi-tree artifacts were found, so thinning is not currently applicable"
        )
    if posterior_sample_file_count == 0:
        global_issues.append(
            "no retained posterior tree samples were found, so compression remains a forward-looking policy"
        )
    _add_check(
        checks,
        "policy-summary:overall-status",
        surface="policy-summary",
        condition=overall_policy_status
        == _derive_overall_policy_status(
            failed_check_count=sum(1 for row in checks if row.status == "failed"),
            thinning_recommended_file_count=thinning_recommended_file_count,
            thinning_required_file_count=thinning_required_file_count,
            compression_recommended_file_count=compression_recommended_file_count,
            compression_required_file_count=compression_required_file_count,
        ),
        expected=_derive_overall_policy_status(
            failed_check_count=sum(1 for row in checks if row.status == "failed"),
            thinning_recommended_file_count=thinning_recommended_file_count,
            thinning_required_file_count=thinning_required_file_count,
            compression_recommended_file_count=compression_recommended_file_count,
            compression_required_file_count=compression_required_file_count,
        ),
        observed=overall_policy_status,
        detail="overall tree-retention status matches the file-level policy counts",
    )
    failed_check_count = sum(1 for row in checks if row.status == "failed")
    if failed_check_count > 0:
        overall_policy_status = "required"

    largest_tree_set_row = max(
        file_rows,
        key=lambda row: (row.tree_count, row.byte_count, row.relative_path),
    )
    return RabiesMethodSensitivitySlurmTreeRetentionReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        overall_policy_status=overall_policy_status,
        variant_count=len(configured_variant_ids),
        file_count=len(file_rows),
        check_count=len(checks),
        failed_check_count=failed_check_count,
        tree_artifact_file_count=tree_artifact_file_count,
        tree_set_file_count=tree_set_file_count,
        posterior_sample_file_count=posterior_sample_file_count,
        thinning_recommended_file_count=thinning_recommended_file_count,
        thinning_required_file_count=thinning_required_file_count,
        compression_recommended_file_count=compression_recommended_file_count,
        compression_required_file_count=compression_required_file_count,
        total_tree_count=total_tree_count,
        total_tree_byte_count=total_tree_byte_count,
        largest_tree_set_path=largest_tree_set_row.relative_path,
        largest_tree_set_tree_count=largest_tree_set_row.tree_count,
        global_issue_count=len(global_issues),
        global_issues=tuple(global_issues),
        checks=tuple(checks),
        files=file_rows,
    )


def _derive_overall_policy_status(
    *,
    failed_check_count: int,
    thinning_recommended_file_count: int,
    thinning_required_file_count: int,
    compression_recommended_file_count: int,
    compression_required_file_count: int,
) -> str:
    if (
        failed_check_count > 0
        or thinning_required_file_count > 0
        or compression_required_file_count > 0
    ):
        return "required"
    if thinning_recommended_file_count > 0 or compression_recommended_file_count > 0:
        return "recommended"
    return "no_action"
