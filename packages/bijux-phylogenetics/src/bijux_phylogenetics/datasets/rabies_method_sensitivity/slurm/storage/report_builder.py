from __future__ import annotations

from pathlib import Path

from .contracts import (
    RabiesMethodSensitivitySlurmStorageAssumptionRow,
    RabiesMethodSensitivitySlurmStorageCategoryRow,
    RabiesMethodSensitivitySlurmStorageReport,
    RabiesMethodSensitivitySlurmStorageVariantRow,
)
from .inventory import scan_storage_inventory
from .shared import (
    _CATEGORY_DETAILS,
    _CATEGORY_LABELS,
    _CONFIG_FILENAME,
    _STORAGE_CATEGORIES,
    _load_json,
    _to_mib,
)


def build_rabies_method_sensitivity_slurm_storage_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmStorageReport:
    """Estimate retained storage by category before scaling this workflow further."""
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    configured_variant_ids = [
        str(row["variant_id"]) for row in list(config.get("variants", []))
    ]
    variant_totals, shared_totals = scan_storage_inventory(
        bundle_root=bundle_root,
        configured_variant_ids=configured_variant_ids,
    )

    variant_rows = tuple(
        _build_variant_row(variant_id=variant_id, totals=variant_totals[variant_id])
        for variant_id in configured_variant_ids
    )
    category_rows = tuple(
        _build_category_row(
            category_id=category_id,
            variant_totals=variant_totals,
            shared_totals=shared_totals,
        )
        for category_id in _STORAGE_CATEGORIES
    )
    largest_variant = max(
        variant_rows,
        key=lambda row: (row.total_byte_count, row.variant_id),
    )
    assumptions = (
        RabiesMethodSensitivitySlurmStorageAssumptionRow(
            assumption_id="observed-retained-bytes",
            parameter="total_byte_count",
            value="sum of bytes currently written into the governed workflow bundle",
            rationale=(
                "The storage estimate is anchored to real retained outputs from the "
                "current governed workflow run instead of an invented storage budget."
            ),
        ),
        RabiesMethodSensitivitySlurmStorageAssumptionRow(
            assumption_id="canonical-log-boundary",
            parameter="log_byte_count",
            value="only parallel-logs/*.log count as canonical logs",
            rationale=(
                "Copied task logs inside per-job evidence packages duplicate the same "
                "execution story, so they remain reviewer artifacts instead of being "
                "counted twice as canonical logging storage."
            ),
        ),
        RabiesMethodSensitivitySlurmStorageAssumptionRow(
            assumption_id="tree-artifact-separation",
            parameter="tree_byte_count",
            value="Newick and other retained tree files are separated from general outputs",
            rationale=(
                "Large workflow storage pressure often comes from retained tree "
                "artifacts, so they need their own reviewer-visible category."
            ),
        ),
        RabiesMethodSensitivitySlurmStorageAssumptionRow(
            assumption_id="explicit-zero-posterior-category",
            parameter="posterior_sample_byte_count",
            value="0 when the governed workflow emits no posterior chains or tree sets",
            rationale=(
                "A zero posterior-sample footprint should be explicit and reviewable "
                "rather than silently folded into a generic no-data category."
            ),
        ),
        RabiesMethodSensitivitySlurmStorageAssumptionRow(
            assumption_id="shared-report-overhead",
            parameter="workflow_shared_byte_count",
            value="workflow-wide reports, manifests, and ledgers are counted once outside per-variant rows",
            rationale=(
                "The retained storage plan needs to distinguish variant-scaled bytes "
                "from one-time workflow-wide reporting overhead."
            ),
        ),
    )
    output_byte_count = sum(
        row.total_byte_count for row in category_rows if row.category_id == "outputs"
    )
    log_byte_count = sum(
        row.total_byte_count for row in category_rows if row.category_id == "logs"
    )
    tree_byte_count = sum(
        row.total_byte_count for row in category_rows if row.category_id == "trees"
    )
    posterior_sample_byte_count = sum(
        row.total_byte_count
        for row in category_rows
        if row.category_id == "posterior_samples"
    )
    report_byte_count = sum(
        row.total_byte_count for row in category_rows if row.category_id == "reports"
    )
    total_file_count = sum(row.total_file_count for row in category_rows)
    total_byte_count = sum(row.total_byte_count for row in category_rows)
    variant_scoped_file_count = sum(row.total_file_count for row in variant_rows)
    variant_scoped_byte_count = sum(row.total_byte_count for row in variant_rows)
    workflow_shared_file_count = sum(row.workflow_file_count for row in category_rows)
    workflow_shared_byte_count = sum(row.workflow_byte_count for row in category_rows)
    return RabiesMethodSensitivitySlurmStorageReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        variant_count=len(configured_variant_ids),
        total_file_count=total_file_count,
        variant_scoped_file_count=variant_scoped_file_count,
        workflow_shared_file_count=workflow_shared_file_count,
        total_byte_count=total_byte_count,
        total_estimated_storage_mib=_to_mib(total_byte_count),
        variant_scoped_byte_count=variant_scoped_byte_count,
        workflow_shared_byte_count=workflow_shared_byte_count,
        output_byte_count=output_byte_count,
        log_byte_count=log_byte_count,
        tree_byte_count=tree_byte_count,
        posterior_sample_byte_count=posterior_sample_byte_count,
        report_byte_count=report_byte_count,
        largest_variant_id=largest_variant.variant_id,
        largest_variant_total_byte_count=largest_variant.total_byte_count,
        assumptions=assumptions,
        categories=category_rows,
        variants=variant_rows,
    )


def _build_category_row(
    *,
    category_id: str,
    variant_totals: dict[str, dict[str, int]],
    shared_totals: dict[str, int],
) -> RabiesMethodSensitivitySlurmStorageCategoryRow:
    variant_file_count = sum(
        totals[f"{category_id}_file_count"] for totals in variant_totals.values()
    )
    variant_byte_count = sum(
        totals[f"{category_id}_byte_count"] for totals in variant_totals.values()
    )
    workflow_file_count = shared_totals[f"{category_id}_file_count"]
    workflow_byte_count = shared_totals[f"{category_id}_byte_count"]
    total_byte_count = variant_byte_count + workflow_byte_count
    return RabiesMethodSensitivitySlurmStorageCategoryRow(
        category_id=category_id,
        category_label=_CATEGORY_LABELS[category_id],
        variant_file_count=variant_file_count,
        workflow_file_count=workflow_file_count,
        total_file_count=variant_file_count + workflow_file_count,
        variant_byte_count=variant_byte_count,
        workflow_byte_count=workflow_byte_count,
        total_byte_count=total_byte_count,
        estimated_storage_mib=_to_mib(total_byte_count),
        detail=_CATEGORY_DETAILS[category_id],
    )


def _build_variant_row(
    *,
    variant_id: str,
    totals: dict[str, int],
) -> RabiesMethodSensitivitySlurmStorageVariantRow:
    return RabiesMethodSensitivitySlurmStorageVariantRow(
        variant_id=variant_id,
        output_file_count=totals["outputs_file_count"],
        log_file_count=totals["logs_file_count"],
        tree_file_count=totals["trees_file_count"],
        posterior_sample_file_count=totals["posterior_samples_file_count"],
        report_file_count=totals["reports_file_count"],
        total_file_count=totals["total_file_count"],
        output_byte_count=totals["outputs_byte_count"],
        log_byte_count=totals["logs_byte_count"],
        tree_byte_count=totals["trees_byte_count"],
        posterior_sample_byte_count=totals["posterior_samples_byte_count"],
        report_byte_count=totals["reports_byte_count"],
        total_byte_count=totals["total_byte_count"],
        estimated_storage_mib=_to_mib(totals["total_byte_count"]),
    )
