from __future__ import annotations

from dataclasses import asdict
import csv
import json
from math import ceil
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    RabiesMethodSensitivitySlurmStorageAssumptionRow,
    RabiesMethodSensitivitySlurmStorageCategoryRow,
    RabiesMethodSensitivitySlurmStorageReport,
    RabiesMethodSensitivitySlurmStorageVariantRow,
)

__all__ = [
    "RabiesMethodSensitivitySlurmStorageAssumptionRow",
    "RabiesMethodSensitivitySlurmStorageCategoryRow",
    "RabiesMethodSensitivitySlurmStorageReport",
    "RabiesMethodSensitivitySlurmStorageVariantRow",
    "build_rabies_method_sensitivity_slurm_storage_report",
    "write_rabies_method_sensitivity_slurm_storage_categories_table",
    "write_rabies_method_sensitivity_slurm_storage_html_report",
    "write_rabies_method_sensitivity_slurm_storage_summary_json",
    "write_rabies_method_sensitivity_slurm_storage_variants_table",
]

_CONFIG_FILENAME = "workflow-config.resolved.json"
_MEBIBYTE = 1024 * 1024
_STORAGE_CATEGORIES = (
    "outputs",
    "logs",
    "trees",
    "posterior_samples",
    "reports",
)
_CATEGORY_LABELS = {
    "outputs": "workflow outputs",
    "logs": "canonical logs",
    "trees": "tree artifacts",
    "posterior_samples": "posterior samples",
    "reports": "review artifacts",
}
_CATEGORY_DETAILS = {
    "outputs": (
        "Retained scientific outputs such as alignments, comparison ledgers, and "
        "variant-scoped tables that are not tree sets or reviewer reports."
    ),
    "logs": (
        "Canonical orchestration logs under parallel-logs. Copied task logs inside "
        "per-job evidence packages are counted under review artifacts to avoid "
        "double-counting the same execution record as canonical storage."
    ),
    "trees": (
        "Tree products such as rooted or unrooted Newick outputs and other retained "
        "tree-topology artifacts."
    ),
    "posterior_samples": (
        "Bayesian posterior chains or sampled tree sets. This governed workflow does "
        "not currently emit them, but the category is kept explicit so zero remains "
        "an audited, reviewer-visible statement."
    ),
    "reports": (
        "Workflow-wide summaries, manifests, HTML reports, Slurm ledgers, evidence "
        "packages, and other reviewer-facing accountability artifacts."
    ),
}
def build_rabies_method_sensitivity_slurm_storage_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmStorageReport:
    """Estimate retained storage by category before scaling this workflow further."""
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    configured_variant_ids = [
        str(row["variant_id"]) for row in list(config.get("variants", []))
    ]
    variant_totals = {
        variant_id: _empty_category_totals() for variant_id in configured_variant_ids
    }
    shared_totals = _empty_category_totals()

    for path in sorted(bundle_root.rglob("*")):
        if not path.is_file():
            continue
        relative_path = path.relative_to(bundle_root)
        category_id, variant_id = _classify_storage_path(relative_path)
        byte_count = path.stat().st_size
        if variant_id is None:
            totals = shared_totals
        else:
            totals = variant_totals.setdefault(variant_id, _empty_category_totals())
        file_key = f"{category_id}_file_count"
        byte_key = f"{category_id}_byte_count"
        totals[file_key] += 1
        totals[byte_key] += byte_count
        totals["total_file_count"] += 1
        totals["total_byte_count"] += byte_count

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
    output_byte_count = sum(row.total_byte_count for row in category_rows if row.category_id == "outputs")
    log_byte_count = sum(row.total_byte_count for row in category_rows if row.category_id == "logs")
    tree_byte_count = sum(row.total_byte_count for row in category_rows if row.category_id == "trees")
    posterior_sample_byte_count = sum(
        row.total_byte_count for row in category_rows if row.category_id == "posterior_samples"
    )
    report_byte_count = sum(row.total_byte_count for row in category_rows if row.category_id == "reports")
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


def write_rabies_method_sensitivity_slurm_storage_categories_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write the category-level retained-storage estimate."""
    return _write_tsv(
        path,
        fieldnames=(
            "category_id",
            "category_label",
            "variant_file_count",
            "workflow_file_count",
            "total_file_count",
            "variant_byte_count",
            "workflow_byte_count",
            "total_byte_count",
            "estimated_storage_mib",
            "detail",
        ),
        rows=[asdict(row) for row in report.categories],
    )


def write_rabies_method_sensitivity_slurm_storage_variants_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write one per-variant retained-storage ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "output_file_count",
            "log_file_count",
            "tree_file_count",
            "posterior_sample_file_count",
            "report_file_count",
            "total_file_count",
            "output_byte_count",
            "log_byte_count",
            "tree_byte_count",
            "posterior_sample_byte_count",
            "report_byte_count",
            "total_byte_count",
            "estimated_storage_mib",
        ),
        rows=[asdict(row) for row in report.variants],
    )


def write_rabies_method_sensitivity_slurm_storage_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write the machine-readable retained-storage summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_storage_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write the reviewer-facing retained-storage estimate."""
    category_lines = [
        (
            f"{row.category_id}: {row.total_file_count} files, "
            f"{row.total_byte_count} bytes, {row.estimated_storage_mib} MiB"
        )
        for row in report.categories
    ]
    largest_variant = next(
        row for row in report.variants if row.variant_id == report.largest_variant_id
    )
    return write_html_report(
        title="Rabies Slurm Storage Report",
        sections=[
            (
                "storage-summary",
                "\n".join(
                    [
                        f"variant_count: {report.variant_count}",
                        f"total_file_count: {report.total_file_count}",
                        f"total_byte_count: {report.total_byte_count}",
                        f"total_estimated_storage_mib: {report.total_estimated_storage_mib}",
                        (
                            "variant_scoped_byte_count: "
                            f"{report.variant_scoped_byte_count}"
                        ),
                        (
                            "workflow_shared_byte_count: "
                            f"{report.workflow_shared_byte_count}"
                        ),
                    ]
                ),
            ),
            ("storage-categories", "\n".join(category_lines)),
            (
                "largest-variant",
                "\n".join(
                    [
                        f"variant_id: {largest_variant.variant_id}",
                        f"total_file_count: {largest_variant.total_file_count}",
                        f"total_byte_count: {largest_variant.total_byte_count}",
                        (
                            "estimated_storage_mib: "
                            f"{largest_variant.estimated_storage_mib}"
                        ),
                    ]
                ),
            ),
            (
                "posterior-samples",
                (
                    "this governed workflow writes no posterior chains or posterior "
                    "tree sets, so posterior_samples remains an explicit zero-valued "
                    "storage category"
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "variant_count": report.variant_count,
            "total_byte_count": report.total_byte_count,
            "total_estimated_storage_mib": report.total_estimated_storage_mib,
            "largest_variant_id": report.largest_variant_id,
            "largest_variant_total_byte_count": report.largest_variant_total_byte_count,
            "output_byte_count": report.output_byte_count,
            "log_byte_count": report.log_byte_count,
            "tree_byte_count": report.tree_byte_count,
            "posterior_sample_byte_count": report.posterior_sample_byte_count,
            "report_byte_count": report.report_byte_count,
        },
        summary_metrics=[
            ("estimated storage MiB", report.total_estimated_storage_mib),
            ("total files", report.total_file_count),
            ("workflow outputs bytes", report.output_byte_count),
            ("log bytes", report.log_byte_count),
            ("tree bytes", report.tree_byte_count),
            ("posterior bytes", report.posterior_sample_byte_count),
            ("report bytes", report.report_byte_count),
        ],
        artifact_links=[
            ("storage categories", "slurm-storage-categories.tsv", None),
            ("storage variants", "slurm-storage-variants.tsv", None),
            ("storage summary", "slurm-storage-report.json", None),
        ],
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


def _classify_storage_path(relative_path: Path) -> tuple[str, str | None]:
    parts = relative_path.parts
    if not parts:
        raise ValueError("relative_path must not be empty")
    top_level = parts[0]
    if top_level == "parallel-logs":
        return "logs", relative_path.stem
    if top_level == "variants" and len(parts) >= 3:
        return _classify_variant_file(relative_path.name), parts[1]
    if top_level == "slurm-job-evidence" and len(parts) >= 3:
        return "reports", parts[1]
    if top_level in {"report-artifacts", "slurm-arrays"}:
        return "reports", None
    return "reports", None


def _classify_variant_file(filename: str) -> str:
    normalized = filename.lower()
    if normalized.endswith((".nwk", ".tree", ".treefile", ".contree", ".ufboot")):
        return "trees"
    if _looks_like_posterior_sample(normalized):
        return "posterior_samples"
    return "outputs"


def _looks_like_posterior_sample(filename: str) -> bool:
    if filename.endswith((".trees", ".state", ".trace", ".p", ".t")):
        return True
    return any(
        token in filename for token in ("posterior", "mcmc", "beast", "mrbayes")
    )


def _empty_category_totals() -> dict[str, int]:
    return {
        "outputs_file_count": 0,
        "logs_file_count": 0,
        "trees_file_count": 0,
        "posterior_samples_file_count": 0,
        "reports_file_count": 0,
        "outputs_byte_count": 0,
        "logs_byte_count": 0,
        "trees_byte_count": 0,
        "posterior_samples_byte_count": 0,
        "reports_byte_count": 0,
        "total_file_count": 0,
        "total_byte_count": 0,
    }


def _to_mib(byte_count: int) -> int:
    return max(0, ceil(byte_count / _MEBIBYTE))


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
