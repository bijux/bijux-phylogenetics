from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

__all__ = [
    "RabiesMethodSensitivitySlurmOutputExplosionCheckRow",
    "RabiesMethodSensitivitySlurmOutputExplosionReport",
    "RabiesMethodSensitivitySlurmOutputExplosionVariantRow",
    "build_rabies_method_sensitivity_slurm_output_explosion_report",
    "write_rabies_method_sensitivity_slurm_output_explosion_checks_table",
    "write_rabies_method_sensitivity_slurm_output_explosion_html_report",
    "write_rabies_method_sensitivity_slurm_output_explosion_summary_json",
    "write_rabies_method_sensitivity_slurm_output_explosion_variants_table",
]

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_JOB_PLAN_FILENAME = "slurm-job-plan.tsv"
_SLURM_STORAGE_CATEGORIES_FILENAME = "slurm-storage-categories.tsv"
_SLURM_STORAGE_VARIANTS_FILENAME = "slurm-storage-variants.tsv"
_SLURM_STORAGE_SUMMARY_FILENAME = "slurm-storage-report.json"
_CATEGORY_IDS = ["logs", "outputs", "posterior_samples", "reports", "trees"]
_MEBIBYTE = 1024 * 1024
_OUTPUT_WARNING_MIB = 128
_OUTPUT_HIGH_MIB = 512
_STORAGE_WARNING_MIB = 256
_STORAGE_HIGH_MIB = 1024
_TREE_WARNING_BYTES = 64 * _MEBIBYTE
_TREE_HIGH_BYTES = 256 * _MEBIBYTE
_TREE_WARNING_FILES = 64
_TREE_HIGH_FILES = 512
_POSTERIOR_WARNING_BYTES = 64 * _MEBIBYTE
_POSTERIOR_HIGH_BYTES = 512 * _MEBIBYTE
_POSTERIOR_WARNING_FILES = 64
_POSTERIOR_HIGH_FILES = 512
_REPORT_WARNING_BYTES = 128 * _MEBIBYTE
_REPORT_HIGH_BYTES = 512 * _MEBIBYTE
_TOTAL_OUTPUT_WARNING_MIB = 512
_TOTAL_OUTPUT_HIGH_MIB = 2048
_TOTAL_STORAGE_WARNING_MIB = 1024
_TOTAL_STORAGE_HIGH_MIB = 4096
_DOMINANT_VARIANT_OUTPUT_SHARE = 0.65
_DOMINANT_VARIANT_OUTPUT_MIB = 128
_BOOTSTRAP_WARNING_REPLICATES = 1000
_BOOTSTRAP_HIGH_REPLICATES = 5000


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputExplosionCheckRow:
    """One machine-readable consistency check behind the explosion-risk report."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputExplosionVariantRow:
    """One per-variant retained-output explosion assessment."""

    variant_id: str
    risk_status: str
    estimated_output_mib: int
    estimated_storage_mib: int
    tree_file_count: int
    tree_byte_count: int
    posterior_sample_file_count: int
    posterior_sample_byte_count: int
    report_byte_count: int
    output_share: float
    issue_count: int
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputExplosionReport:
    """One retained-output explosion-risk summary over the governed rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    bootstrap_replicates: int
    overall_risk_status: str
    variant_count: int
    check_count: int
    failed_check_count: int
    low_risk_variant_count: int
    warning_variant_count: int
    high_risk_variant_count: int
    global_issue_count: int
    total_estimated_output_mib: int
    total_estimated_storage_mib: int
    total_tree_byte_count: int
    total_tree_file_count: int
    total_posterior_sample_byte_count: int
    total_posterior_sample_file_count: int
    total_report_byte_count: int
    largest_variant_id: str
    largest_variant_output_share: float
    global_issues: tuple[str, ...]
    checks: tuple[RabiesMethodSensitivitySlurmOutputExplosionCheckRow, ...]
    variants: tuple[RabiesMethodSensitivitySlurmOutputExplosionVariantRow, ...]


def build_rabies_method_sensitivity_slurm_output_explosion_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmOutputExplosionReport:
    """Assess whether retained workflow outputs are trending toward an explosion risk."""
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    job_plan_rows = _read_tsv_rows(bundle_root / _SLURM_JOB_PLAN_FILENAME)
    storage_category_rows = _read_tsv_rows(bundle_root / _SLURM_STORAGE_CATEGORIES_FILENAME)
    storage_variant_rows = _read_tsv_rows(bundle_root / _SLURM_STORAGE_VARIANTS_FILENAME)
    storage_summary = _load_json(bundle_root / _SLURM_STORAGE_SUMMARY_FILENAME)

    checks: list[RabiesMethodSensitivitySlurmOutputExplosionCheckRow] = []

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
            RabiesMethodSensitivitySlurmOutputExplosionCheckRow(
                check_id=check_id,
                surface=surface,
                status="passed" if condition else "failed",
                expected="" if expected is None else str(expected),
                observed="" if observed is None else str(observed),
                detail=detail,
            )
        )

    configured_variant_ids = sorted(
        str(row["variant_id"]) for row in list(config.get("variants", []))
    )
    job_plan_by_variant = {
        str(row["variant_id"]): row for row in job_plan_rows
    }
    storage_variant_by_variant = {
        str(row["variant_id"]): row for row in storage_variant_rows
    }
    storage_category_by_id = {
        str(row["category_id"]): row for row in storage_category_rows
    }
    add_check(
        "job-plan:variant-coverage",
        surface="job-plan",
        condition=sorted(job_plan_by_variant) == configured_variant_ids,
        expected=configured_variant_ids,
        observed=sorted(job_plan_by_variant),
        detail="job-plan rows cover the configured variant ids",
    )
    add_check(
        "storage:variant-coverage",
        surface="storage",
        condition=sorted(storage_variant_by_variant) == configured_variant_ids,
        expected=configured_variant_ids,
        observed=sorted(storage_variant_by_variant),
        detail="storage variant rows cover the configured variant ids",
    )
    add_check(
        "storage:category-coverage",
        surface="storage",
        condition=sorted(storage_category_by_id) == _CATEGORY_IDS,
        expected=_CATEGORY_IDS,
        observed=sorted(storage_category_by_id),
        detail="storage category rows cover the explicit retained-output categories",
    )
    total_estimated_output_mib = sum(
        int(row["estimated_output_mib"]) for row in job_plan_rows
    )
    add_check(
        "job-plan:total-output-mib",
        surface="job-plan",
        condition=int(storage_summary["variant_count"]) == len(configured_variant_ids),
        expected=len(configured_variant_ids),
        observed=storage_summary["variant_count"],
        detail="storage summary variant_count matches the configured variant count",
    )
    add_check(
        "storage:total-storage-mib",
        surface="storage",
        condition=int(storage_summary["total_estimated_storage_mib"])
        >= int(storage_summary["total_byte_count"]) // _MEBIBYTE,
        expected="storage summary rounds total retained bytes up into MiB",
        observed=storage_summary["total_estimated_storage_mib"],
        detail="storage summary total_estimated_storage_mib is a rounded-up retained-byte estimate",
    )
    add_check(
        "storage:largest-variant",
        surface="storage",
        condition=str(storage_summary["largest_variant_id"])
        == max(
            storage_variant_rows,
            key=lambda row: (int(row["total_byte_count"]), str(row["variant_id"])),
        )["variant_id"],
        expected=storage_summary["largest_variant_id"],
        observed=max(
            storage_variant_rows,
            key=lambda row: (int(row["total_byte_count"]), str(row["variant_id"])),
        )["variant_id"],
        detail="storage summary largest_variant_id matches the written storage-variant rows",
    )

    variant_rows: list[RabiesMethodSensitivitySlurmOutputExplosionVariantRow] = []
    total_output_denom = max(1, total_estimated_output_mib)
    for variant_id in configured_variant_ids:
        plan_row = job_plan_by_variant.get(variant_id)
        storage_row = storage_variant_by_variant.get(variant_id)
        if plan_row is None or storage_row is None:
            issues = ["planning or storage row is missing"]
            variant_rows.append(
                RabiesMethodSensitivitySlurmOutputExplosionVariantRow(
                    variant_id=variant_id,
                    risk_status="high",
                    estimated_output_mib=0 if plan_row is None else int(plan_row["estimated_output_mib"]),
                    estimated_storage_mib=0 if storage_row is None else int(storage_row["estimated_storage_mib"]),
                    tree_file_count=0 if storage_row is None else int(storage_row["tree_file_count"]),
                    tree_byte_count=0 if storage_row is None else int(storage_row["tree_byte_count"]),
                    posterior_sample_file_count=0 if storage_row is None else int(storage_row["posterior_sample_file_count"]),
                    posterior_sample_byte_count=0 if storage_row is None else int(storage_row["posterior_sample_byte_count"]),
                    report_byte_count=0 if storage_row is None else int(storage_row["report_byte_count"]),
                    output_share=0.0,
                    issue_count=1,
                    issues=tuple(issues),
                )
            )
            continue

        estimated_output_mib = int(plan_row["estimated_output_mib"])
        estimated_storage_mib = int(storage_row["estimated_storage_mib"])
        tree_file_count = int(storage_row["tree_file_count"])
        tree_byte_count = int(storage_row["tree_byte_count"])
        posterior_sample_file_count = int(storage_row["posterior_sample_file_count"])
        posterior_sample_byte_count = int(storage_row["posterior_sample_byte_count"])
        report_byte_count = int(storage_row["report_byte_count"])
        output_share = round(estimated_output_mib / total_output_denom, 4)
        issues: list[str] = []
        severity = "low"

        severity = _raise_variant_risk(
            severity,
            issues,
            value=estimated_output_mib,
            warning_threshold=_OUTPUT_WARNING_MIB,
            high_threshold=_OUTPUT_HIGH_MIB,
            warning_detail="estimated retained output MiB is large for one variant",
            high_detail="estimated retained output MiB is very large for one variant",
        )
        severity = _raise_variant_risk(
            severity,
            issues,
            value=estimated_storage_mib,
            warning_threshold=_STORAGE_WARNING_MIB,
            high_threshold=_STORAGE_HIGH_MIB,
            warning_detail="estimated retained storage MiB is large for one variant",
            high_detail="estimated retained storage MiB is very large for one variant",
        )
        severity = _raise_variant_risk(
            severity,
            issues,
            value=tree_byte_count,
            warning_threshold=_TREE_WARNING_BYTES,
            high_threshold=_TREE_HIGH_BYTES,
            warning_detail="tree artifact bytes are large enough to pressure retained storage",
            high_detail="tree artifact bytes are large enough to dominate retained storage",
        )
        severity = _raise_variant_risk(
            severity,
            issues,
            value=tree_file_count,
            warning_threshold=_TREE_WARNING_FILES,
            high_threshold=_TREE_HIGH_FILES,
            warning_detail="tree artifact file counts suggest tree-set growth pressure",
            high_detail="tree artifact file counts suggest a large retained tree-set burden",
        )
        severity = _raise_variant_risk(
            severity,
            issues,
            value=posterior_sample_byte_count,
            warning_threshold=_POSTERIOR_WARNING_BYTES,
            high_threshold=_POSTERIOR_HIGH_BYTES,
            warning_detail="posterior sample bytes are large enough to pressure retained storage",
            high_detail="posterior sample bytes are large enough to dominate retained storage",
        )
        severity = _raise_variant_risk(
            severity,
            issues,
            value=posterior_sample_file_count,
            warning_threshold=_POSTERIOR_WARNING_FILES,
            high_threshold=_POSTERIOR_HIGH_FILES,
            warning_detail="posterior sample file counts suggest chain or tree-set growth pressure",
            high_detail="posterior sample file counts suggest a very large retained chain or tree-set burden",
        )
        severity = _raise_variant_risk(
            severity,
            issues,
            value=report_byte_count,
            warning_threshold=_REPORT_WARNING_BYTES,
            high_threshold=_REPORT_HIGH_BYTES,
            warning_detail="review artifact bytes are becoming large for one variant",
            high_detail="review artifact bytes are dominating retained storage for one variant",
        )
        if (
            output_share >= _DOMINANT_VARIANT_OUTPUT_SHARE
            and estimated_output_mib >= _DOMINANT_VARIANT_OUTPUT_MIB
        ):
            severity = _escalate_risk(
                severity,
                "warning",
            )
            issues.append("one variant dominates most retained output MiB")
        variant_rows.append(
            RabiesMethodSensitivitySlurmOutputExplosionVariantRow(
                variant_id=variant_id,
                risk_status=severity,
                estimated_output_mib=estimated_output_mib,
                estimated_storage_mib=estimated_storage_mib,
                tree_file_count=tree_file_count,
                tree_byte_count=tree_byte_count,
                posterior_sample_file_count=posterior_sample_file_count,
                posterior_sample_byte_count=posterior_sample_byte_count,
                report_byte_count=report_byte_count,
                output_share=output_share,
                issue_count=len(issues),
                issues=tuple(issues),
            )
        )

    storage_trees_row = storage_category_by_id["trees"]
    storage_posterior_row = storage_category_by_id["posterior_samples"]
    storage_reports_row = storage_category_by_id["reports"]
    global_issues: list[str] = []
    global_severity = "low"
    global_severity = _raise_variant_risk(
        global_severity,
        global_issues,
        value=total_estimated_output_mib,
        warning_threshold=_TOTAL_OUTPUT_WARNING_MIB,
        high_threshold=_TOTAL_OUTPUT_HIGH_MIB,
        warning_detail="total retained output MiB suggests scaling pressure",
        high_detail="total retained output MiB suggests a severe retained-output burden",
    )
    global_severity = _raise_variant_risk(
        global_severity,
        global_issues,
        value=int(storage_summary["total_estimated_storage_mib"]),
        warning_threshold=_TOTAL_STORAGE_WARNING_MIB,
        high_threshold=_TOTAL_STORAGE_HIGH_MIB,
        warning_detail="total retained storage MiB suggests scaling pressure",
        high_detail="total retained storage MiB suggests severe retained-storage pressure",
    )
    global_severity = _raise_variant_risk(
        global_severity,
        global_issues,
        value=int(storage_trees_row["total_byte_count"]),
        warning_threshold=_TREE_WARNING_BYTES,
        high_threshold=_TREE_HIGH_BYTES,
        warning_detail="total tree artifact bytes are large enough to pressure retained storage",
        high_detail="total tree artifact bytes are dominating retained storage",
    )
    global_severity = _raise_variant_risk(
        global_severity,
        global_issues,
        value=int(storage_posterior_row["total_byte_count"]),
        warning_threshold=_POSTERIOR_WARNING_BYTES,
        high_threshold=_POSTERIOR_HIGH_BYTES,
        warning_detail="posterior sample bytes are becoming large enough to pressure retained storage",
        high_detail="posterior sample bytes are dominating retained storage",
    )
    global_severity = _raise_variant_risk(
        global_severity,
        global_issues,
        value=int(storage_reports_row["total_byte_count"]),
        warning_threshold=_REPORT_WARNING_BYTES,
        high_threshold=_REPORT_HIGH_BYTES,
        warning_detail="review artifact bytes are becoming large enough to pressure retained storage",
        high_detail="review artifact bytes are dominating retained storage",
    )
    bootstrap_replicates = int(config["bootstrap_replicates"])
    total_tree_file_count = int(storage_trees_row["total_file_count"])
    if (
        bootstrap_replicates >= _BOOTSTRAP_HIGH_REPLICATES
        and total_tree_file_count >= _TREE_WARNING_FILES
    ):
        global_severity = _escalate_risk(global_severity, "high")
        global_issues.append(
            "bootstrap replicates and retained tree-file counts together suggest a severe tree-output explosion risk"
        )
    elif (
        bootstrap_replicates >= _BOOTSTRAP_WARNING_REPLICATES
        and total_tree_file_count >= _TREE_WARNING_FILES
    ):
        global_severity = _escalate_risk(global_severity, "warning")
        global_issues.append(
            "bootstrap replicates and retained tree-file counts together suggest tree-output growth pressure"
        )

    low_risk_variant_count = sum(1 for row in variant_rows if row.risk_status == "low")
    warning_variant_count = sum(
        1 for row in variant_rows if row.risk_status == "warning"
    )
    high_risk_variant_count = sum(
        1 for row in variant_rows if row.risk_status == "high"
    )
    failed_check_count = sum(1 for row in checks if row.status == "failed")
    overall_risk_status = global_severity
    if high_risk_variant_count > 0:
        overall_risk_status = "high"
    elif overall_risk_status != "high" and warning_variant_count > 0:
        overall_risk_status = "warning"
    if failed_check_count > 0:
        overall_risk_status = "high"

    largest_variant = max(
        variant_rows,
        key=lambda row: (row.output_share, row.variant_id),
    )
    add_check(
        "risk-summary:variant-counts",
        surface="risk-summary",
        condition=low_risk_variant_count + warning_variant_count + high_risk_variant_count
        == len(variant_rows),
        expected=len(variant_rows),
        observed=low_risk_variant_count + warning_variant_count + high_risk_variant_count,
        detail="risk-status counts cover every configured variant exactly once",
    )
    add_check(
        "risk-summary:overall-status",
        surface="risk-summary",
        condition=overall_risk_status
        == _derive_overall_risk_status(
            failed_check_count=failed_check_count,
            warning_variant_count=warning_variant_count,
            high_risk_variant_count=high_risk_variant_count,
            global_severity=global_severity,
        ),
        expected=_derive_overall_risk_status(
            failed_check_count=failed_check_count,
            warning_variant_count=warning_variant_count,
            high_risk_variant_count=high_risk_variant_count,
            global_severity=global_severity,
        ),
        observed=overall_risk_status,
        detail="overall explosion-risk status matches the global and per-variant risk counts",
    )

    return RabiesMethodSensitivitySlurmOutputExplosionReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        bootstrap_replicates=bootstrap_replicates,
        overall_risk_status=overall_risk_status,
        variant_count=len(configured_variant_ids),
        check_count=len(checks),
        failed_check_count=failed_check_count,
        low_risk_variant_count=low_risk_variant_count,
        warning_variant_count=warning_variant_count,
        high_risk_variant_count=high_risk_variant_count,
        global_issue_count=len(global_issues),
        total_estimated_output_mib=total_estimated_output_mib,
        total_estimated_storage_mib=int(storage_summary["total_estimated_storage_mib"]),
        total_tree_byte_count=int(storage_trees_row["total_byte_count"]),
        total_tree_file_count=int(storage_trees_row["total_file_count"]),
        total_posterior_sample_byte_count=int(storage_posterior_row["total_byte_count"]),
        total_posterior_sample_file_count=int(storage_posterior_row["total_file_count"]),
        total_report_byte_count=int(storage_reports_row["total_byte_count"]),
        largest_variant_id=largest_variant.variant_id,
        largest_variant_output_share=largest_variant.output_share,
        global_issues=tuple(global_issues),
        checks=tuple(checks),
        variants=tuple(variant_rows),
    )


def write_rabies_method_sensitivity_slurm_output_explosion_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the check-level output explosion ledger."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_slurm_output_explosion_variants_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the per-variant output explosion assessment."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "risk_status",
            "estimated_output_mib",
            "estimated_storage_mib",
            "tree_file_count",
            "tree_byte_count",
            "posterior_sample_file_count",
            "posterior_sample_byte_count",
            "report_byte_count",
            "output_share",
            "issue_count",
            "issues",
        ),
        rows=[
            {
                "variant_id": row.variant_id,
                "risk_status": row.risk_status,
                "estimated_output_mib": row.estimated_output_mib,
                "estimated_storage_mib": row.estimated_storage_mib,
                "tree_file_count": row.tree_file_count,
                "tree_byte_count": row.tree_byte_count,
                "posterior_sample_file_count": row.posterior_sample_file_count,
                "posterior_sample_byte_count": row.posterior_sample_byte_count,
                "report_byte_count": row.report_byte_count,
                "output_share": _format_float(row.output_share),
                "issue_count": row.issue_count,
                "issues": " | ".join(row.issues),
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_slurm_output_explosion_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the machine-readable output explosion report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_output_explosion_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the reviewer-facing output explosion report."""
    return write_html_report(
        title="Rabies Slurm Output Explosion Report",
        sections=[
            (
                "risk-summary",
                "\n".join(
                    [
                        f"overall_risk_status: {report.overall_risk_status}",
                        f"bootstrap_replicates: {report.bootstrap_replicates}",
                        f"variant_count: {report.variant_count}",
                        (
                            "total_estimated_output_mib: "
                            f"{report.total_estimated_output_mib}"
                        ),
                        (
                            "total_estimated_storage_mib: "
                            f"{report.total_estimated_storage_mib}"
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
                "high-risk-variants",
                "none"
                if report.high_risk_variant_count == 0
                else "\n".join(
                    f"{row.variant_id}: {'; '.join(row.issues)}"
                    for row in report.variants
                    if row.risk_status == "high"
                ),
            ),
            (
                "warning-variants",
                "none"
                if report.warning_variant_count == 0
                else "\n".join(
                    f"{row.variant_id}: {'; '.join(row.issues)}"
                    for row in report.variants
                    if row.risk_status == "warning"
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "overall_risk_status": report.overall_risk_status,
            "bootstrap_replicates": report.bootstrap_replicates,
            "variant_count": report.variant_count,
            "warning_variant_count": report.warning_variant_count,
            "high_risk_variant_count": report.high_risk_variant_count,
            "total_estimated_output_mib": report.total_estimated_output_mib,
            "total_estimated_storage_mib": report.total_estimated_storage_mib,
            "total_tree_byte_count": report.total_tree_byte_count,
            "total_posterior_sample_byte_count": report.total_posterior_sample_byte_count,
            "total_report_byte_count": report.total_report_byte_count,
            "largest_variant_id": report.largest_variant_id,
            "largest_variant_output_share": report.largest_variant_output_share,
        },
        summary_metrics=[
            ("overall risk", report.overall_risk_status),
            ("warning variants", report.warning_variant_count),
            ("high-risk variants", report.high_risk_variant_count),
            ("global issues", report.global_issue_count),
            ("estimated output MiB", report.total_estimated_output_mib),
            ("estimated storage MiB", report.total_estimated_storage_mib),
            ("posterior bytes", report.total_posterior_sample_byte_count),
        ],
        artifact_links=[
            ("output explosion checks", "slurm-output-explosion-checks.tsv", None),
            ("output explosion variants", "slurm-output-explosion-variants.tsv", None),
            ("output explosion summary", "slurm-output-explosion-report.json", None),
        ],
    )


def _derive_overall_risk_status(
    *,
    failed_check_count: int,
    warning_variant_count: int,
    high_risk_variant_count: int,
    global_severity: str,
) -> str:
    if failed_check_count > 0 or high_risk_variant_count > 0 or global_severity == "high":
        return "high"
    if warning_variant_count > 0 or global_severity == "warning":
        return "warning"
    return "low"


def _raise_variant_risk(
    current_severity: str,
    issues: list[str],
    *,
    value: int,
    warning_threshold: int,
    high_threshold: int,
    warning_detail: str,
    high_detail: str,
) -> str:
    if value >= high_threshold:
        issues.append(high_detail)
        return _escalate_risk(current_severity, "high")
    if value >= warning_threshold:
        issues.append(warning_detail)
        return _escalate_risk(current_severity, "warning")
    return current_severity


def _escalate_risk(current_severity: str, next_severity: str) -> str:
    order = {"low": 0, "warning": 1, "high": 2}
    return next_severity if order[next_severity] > order[current_severity] else current_severity


def _format_float(value: float) -> str:
    return f"{value:.4f}".rstrip("0").rstrip(".")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


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
