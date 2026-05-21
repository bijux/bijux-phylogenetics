from __future__ import annotations

from dataclasses import asdict
import csv
import gzip
import json
from math import ceil
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from .contracts import (
    RabiesMethodSensitivitySlurmTreeRetentionCheckRow,
    RabiesMethodSensitivitySlurmTreeRetentionFileRow,
    RabiesMethodSensitivitySlurmTreeRetentionReport,
)

__all__ = [
    "RabiesMethodSensitivitySlurmTreeRetentionCheckRow",
    "RabiesMethodSensitivitySlurmTreeRetentionFileRow",
    "RabiesMethodSensitivitySlurmTreeRetentionReport",
    "build_rabies_method_sensitivity_slurm_tree_retention_report",
    "write_rabies_method_sensitivity_slurm_tree_retention_checks_table",
    "write_rabies_method_sensitivity_slurm_tree_retention_files_table",
    "write_rabies_method_sensitivity_slurm_tree_retention_html_report",
    "write_rabies_method_sensitivity_slurm_tree_retention_summary_json",
]

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_STORAGE_CATEGORIES_FILENAME = "slurm-storage-categories.tsv"
_SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME = "slurm-output-explosion-report.json"
_MEBIBYTE = 1024 * 1024
_TREE_FILE_SUFFIXES = {
    ".boottrees",
    ".contree",
    ".newick",
    ".nwk",
    ".phy",
    ".tree",
    ".treefile",
    ".trees",
    ".tre",
}
_COMPRESSED_SUFFIXES = {".bz2", ".gz", ".xz"}
_THINNING_RECOMMENDED_TREE_COUNT = 2_000
_THINNING_REQUIRED_TREE_COUNT = 10_000
_THINNING_TARGET_TREE_COUNT = 1_000
_COMPRESSION_RECOMMENDED_TREE_COUNT = 2_000
_COMPRESSION_REQUIRED_TREE_COUNT = 10_000
_COMPRESSION_RECOMMENDED_BYTES = 4 * _MEBIBYTE
_COMPRESSION_REQUIRED_BYTES = 64 * _MEBIBYTE

def build_rabies_method_sensitivity_slurm_tree_retention_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmTreeRetentionReport:
    """Derive safe thinning and compression policies for retained tree-bearing files."""
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    storage_category_rows = _read_tsv_rows(bundle_root / _SLURM_STORAGE_CATEGORIES_FILENAME)
    output_explosion_summary = _load_json(
        bundle_root / _SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME
    )
    checks: list[RabiesMethodSensitivitySlurmTreeRetentionCheckRow] = []

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
            RabiesMethodSensitivitySlurmTreeRetentionCheckRow(
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
    storage_category_by_id = {
        str(row["category_id"]): row for row in storage_category_rows
    }
    add_check(
        "storage:category-coverage",
        surface="storage",
        condition={"logs", "outputs", "posterior_samples", "reports", "trees"}
        == set(storage_category_by_id),
        expected=sorted(["logs", "outputs", "posterior_samples", "reports", "trees"]),
        observed=sorted(storage_category_by_id),
        detail="tree-retention policy reads the explicit storage-category surface",
    )

    file_rows = tuple(
        _build_file_row(bundle_root=bundle_root, relative_path=relative_path)
        for relative_path in _iter_tree_relative_paths(bundle_root)
    )
    observed_variant_ids = sorted({row.variant_id for row in file_rows})
    add_check(
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
    add_check(
        "storage:tree-file-count",
        surface="storage",
        condition=int(storage_category_by_id["trees"]["total_file_count"])
        == tree_artifact_file_count,
        expected=storage_category_by_id["trees"]["total_file_count"],
        observed=tree_artifact_file_count,
        detail="storage tree file count matches the inspected tree-artifact files",
    )
    add_check(
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
    add_check(
        "storage:posterior-file-count",
        surface="storage",
        condition=int(storage_category_by_id["posterior_samples"]["total_file_count"])
        == posterior_sample_file_count,
        expected=storage_category_by_id["posterior_samples"]["total_file_count"],
        observed=posterior_sample_file_count,
        detail="posterior-sample file count matches the inspected tree-bearing sample files",
    )
    add_check(
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
    add_check(
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
    add_check(
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
    elif (
        thinning_recommended_file_count > 0
        or compression_recommended_file_count > 0
    ):
        overall_policy_status = "recommended"
    if tree_set_file_count == 0:
        global_issues.append(
            "no retained multi-tree artifacts were found, so thinning is not currently applicable"
        )
    if posterior_sample_file_count == 0:
        global_issues.append(
            "no retained posterior tree samples were found, so compression remains a forward-looking policy"
        )
    add_check(
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


def write_rabies_method_sensitivity_slurm_tree_retention_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the check-level tree-retention ledger."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_slurm_tree_retention_files_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the per-file tree-retention recommendations."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "relative_path",
            "artifact_scope",
            "tree_count",
            "byte_count",
            "thinning_policy",
            "thinning_interval",
            "retained_tree_count",
            "compression_policy",
            "recommended_suffix",
            "issue_count",
            "issues",
        ),
        rows=[
            {
                "variant_id": row.variant_id,
                "relative_path": row.relative_path,
                "artifact_scope": row.artifact_scope,
                "tree_count": row.tree_count,
                "byte_count": row.byte_count,
                "thinning_policy": row.thinning_policy,
                "thinning_interval": row.thinning_interval,
                "retained_tree_count": row.retained_tree_count,
                "compression_policy": row.compression_policy,
                "recommended_suffix": row.recommended_suffix,
                "issue_count": row.issue_count,
                "issues": " | ".join(row.issues),
            }
            for row in report.files
        ],
    )


def write_rabies_method_sensitivity_slurm_tree_retention_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the machine-readable tree-retention policy report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


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


def _iter_tree_relative_paths(bundle_root: Path) -> tuple[Path, ...]:
    return tuple(
        relative_path
        for relative_path in sorted(
            path.relative_to(bundle_root)
            for path in bundle_root.rglob("*")
            if path.is_file() and _is_tree_bearing_file(path.relative_to(bundle_root))
        )
    )


def _build_file_row(
    *,
    bundle_root: Path,
    relative_path: Path,
) -> RabiesMethodSensitivitySlurmTreeRetentionFileRow:
    path = bundle_root / relative_path
    byte_count = path.stat().st_size
    tree_count = _count_trees(path)
    artifact_scope = _classify_artifact_scope(relative_path)
    thinning_policy, thinning_interval, retained_tree_count = _derive_thinning_policy(
        tree_count=tree_count
    )
    compression_policy, recommended_suffix = _derive_compression_policy(
        relative_path=relative_path,
        tree_count=tree_count,
        byte_count=byte_count,
    )
    issues: list[str] = []
    if thinning_policy == "thin_recommended":
        issues.append(
            "tree-set size is large enough that interval thinning is recommended"
        )
    elif thinning_policy == "thin_required":
        issues.append(
            "tree-set size is large enough that interval thinning is required"
        )
    if compression_policy == "compress_recommended":
        issues.append("tree-set size is large enough that gzip compression is recommended")
    elif compression_policy == "compress_required":
        issues.append("tree-set size is large enough that gzip compression is required")
    if tree_count <= 1 and artifact_scope == "tree_artifact":
        issues.append("single-tree artifact should be kept in full without thinning")
    return RabiesMethodSensitivitySlurmTreeRetentionFileRow(
        variant_id=_variant_id_from_relative_path(relative_path),
        relative_path=relative_path.as_posix(),
        artifact_scope=artifact_scope,
        tree_count=tree_count,
        byte_count=byte_count,
        thinning_policy=thinning_policy,
        thinning_interval=thinning_interval,
        retained_tree_count=retained_tree_count,
        compression_policy=compression_policy,
        recommended_suffix=recommended_suffix,
        issue_count=len(issues),
        issues=tuple(issues),
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


def _derive_thinning_policy(*, tree_count: int) -> tuple[str, int, int]:
    if tree_count <= 1:
        return ("not_applicable", 1, tree_count)
    if tree_count >= _THINNING_REQUIRED_TREE_COUNT:
        interval = max(2, ceil(tree_count / _THINNING_TARGET_TREE_COUNT))
        return ("thin_required", interval, ceil(tree_count / interval))
    if tree_count >= _THINNING_RECOMMENDED_TREE_COUNT:
        interval = max(2, ceil(tree_count / _THINNING_TARGET_TREE_COUNT))
        return ("thin_recommended", interval, ceil(tree_count / interval))
    return ("keep_full", 1, tree_count)


def _derive_compression_policy(
    *,
    relative_path: Path,
    tree_count: int,
    byte_count: int,
) -> tuple[str, str]:
    if any(suffix in _COMPRESSED_SUFFIXES for suffix in relative_path.suffixes):
        return ("already_compressed", "")
    if (
        tree_count >= _COMPRESSION_REQUIRED_TREE_COUNT
        or byte_count >= _COMPRESSION_REQUIRED_BYTES
    ):
        return ("compress_required", ".gz")
    if (
        tree_count >= _COMPRESSION_RECOMMENDED_TREE_COUNT
        or byte_count >= _COMPRESSION_RECOMMENDED_BYTES
    ):
        return ("compress_recommended", ".gz")
    return ("keep_plain", "")


def _classify_artifact_scope(relative_path: Path) -> str:
    name = relative_path.name.lower()
    path_text = relative_path.as_posix().lower()
    if any(
        token in path_text
        for token in ("posterior", "boottrees", "bootstrap", "treeset", "samples")
    ) or name.endswith(".trees"):
        return "posterior_sample"
    return "tree_artifact"


def _variant_id_from_relative_path(relative_path: Path) -> str:
    parts = relative_path.parts
    if len(parts) >= 3 and parts[0] == "variants":
        return parts[1]
    if len(parts) >= 3 and parts[0] == "slurm-job-evidence":
        return parts[1]
    return "workflow_shared"


def _count_trees(path: Path) -> int:
    if ".gz" in path.suffixes:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return max(1, handle.read().count(";"))
    return max(1, path.read_text(encoding="utf-8").count(";"))


def _is_tree_bearing_file(relative_path: Path) -> bool:
    if relative_path.parts and relative_path.parts[0] != "variants":
        return False
    if len(relative_path.parts) < 3:
        return False
    suffixes = relative_path.suffixes
    if not suffixes:
        return False
    terminal_suffix = suffixes[-1].lower()
    if terminal_suffix in _COMPRESSED_SUFFIXES and len(suffixes) >= 2:
        terminal_suffix = suffixes[-2].lower()
    return terminal_suffix in _TREE_FILE_SUFFIXES


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
