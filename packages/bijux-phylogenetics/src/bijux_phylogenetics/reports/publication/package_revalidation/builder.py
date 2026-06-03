from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..support import (
    SUPPORTED_PUBLICATION_PACKAGE_KIND,
    checksum,
    entry_checksum,
    entry_path,
    entry_size,
    ignored_package_prefixes,
    manifest_file_entries,
    mapping,
    read_manifest,
    read_tsv_rows,
    section_counts,
    text,
)
from .contracts import (
    PublicationPackageRevalidationArtifactRow,
    PublicationPackageRevalidationCheckRow,
    PublicationPackageRevalidationResult,
)
from .inventory import artifact_row, unexpected_files
from .presentation import write_html_report
from .revalidation_policy import (
    check_row,
    status,
)
from .revalidation_policy import (
    overall_revalidation_status as derive_overall_revalidation_status,
)

_ARTIFACT_COLUMNS = [
    "artifact_scope",
    "section",
    "kind",
    "relative_path",
    "status",
    "expected_sha256",
    "observed_sha256",
    "expected_size_bytes",
    "observed_size_bytes",
    "detail",
]
_CHECK_COLUMNS = [
    "section",
    "check_id",
    "status",
    "summary",
    "evidence",
    "artifact_path",
]


def write_publication_package_revalidation_report(
    output_root: Path,
    manifest_path: Path,
) -> PublicationPackageRevalidationResult:
    """Revalidate one stored publication package against its recorded manifest."""

    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_path.resolve()
    package_root = manifest_path.parent.resolve()
    manifest = read_manifest(manifest_path)
    report_kind = text(manifest.get("report_kind"))
    artifact_rows: list[PublicationPackageRevalidationArtifactRow] = []
    check_rows: list[PublicationPackageRevalidationCheckRow] = []

    supported = report_kind == SUPPORTED_PUBLICATION_PACKAGE_KIND
    check_rows.append(
        check_row(
            section="manifest",
            check_id="supported-package-kind",
            status=status(blocked=not supported),
            summary=(
                "package manifest uses the governed rabies study package contract"
                if supported
                else "package manifest is not on the governed rabies study package contract"
            ),
            evidence=f"report_kind={report_kind or 'missing'}",
            artifact_path=manifest_path.name,
        )
    )

    package_files = mapping(manifest, "package_files")
    inventory_entry = mapping(package_files, "artifact_inventory")
    inventory_relative_path = entry_path(inventory_entry)
    inventory_path = (
        package_root / inventory_relative_path
        if inventory_relative_path
        else package_root
    )
    inventory_row = artifact_row(
        artifact_scope="package_control",
        section="package",
        relative_path=inventory_relative_path or manifest_path.name,
        expected_sha256=entry_checksum(inventory_entry),
        expected_size_bytes=None,
        observed_path=inventory_path,
    )
    artifact_rows.append(inventory_row)

    inventory_rows: list[dict[str, str]] = []
    inventory_parse_error: str | None = None
    if inventory_path.exists():
        try:
            inventory_rows = read_tsv_rows(inventory_path)
        except (csv.Error, UnicodeDecodeError) as error:
            inventory_parse_error = str(error)
    else:
        inventory_parse_error = "artifact inventory is missing"

    inventory_count_matches = (
        inventory_entry.get("artifact_count") == len(inventory_rows)
        if inventory_rows
        else False
    )
    manifest_section_counts = inventory_entry.get("section_counts")
    inventory_sections_match = (
        isinstance(manifest_section_counts, dict)
        and {str(key): int(value) for key, value in manifest_section_counts.items()}
        == section_counts(inventory_rows)
        if inventory_rows
        else False
    )
    inventory_blocked = (
        inventory_row.status == "blocked"
        or inventory_parse_error is not None
        or not inventory_count_matches
        or not inventory_sections_match
    )
    inventory_evidence_parts = [
        f"stored checksum match={str(inventory_row.status == 'pass').lower()}",
    ]
    if inventory_parse_error is not None:
        inventory_evidence_parts.append(inventory_parse_error)
    else:
        inventory_evidence_parts.append(
            f"stored artifact_count={inventory_entry.get('artifact_count', '')}; observed rows={len(inventory_rows)}"
        )
        inventory_evidence_parts.append(
            f"stored section counts match={str(inventory_sections_match).lower()}"
        )
    check_rows.append(
        check_row(
            section="package",
            check_id="artifact-inventory-preserved",
            status=status(blocked=inventory_blocked),
            summary="artifact inventory still matches the stored package-control record",
            evidence="; ".join(inventory_evidence_parts),
            artifact_path=inventory_relative_path,
        )
    )

    checklist_entry = mapping(package_files, "reproducibility_checklist")
    checklist_relative_path = entry_path(checklist_entry)
    checklist_path = (
        package_root / checklist_relative_path
        if checklist_relative_path
        else package_root
    )
    checklist_row = artifact_row(
        artifact_scope="package_control",
        section="package",
        relative_path=checklist_relative_path or manifest_path.name,
        expected_sha256=entry_checksum(checklist_entry),
        expected_size_bytes=None,
        observed_path=checklist_path,
    )
    artifact_rows.append(checklist_row)

    checklist_rows_payload: list[dict[str, str]] = []
    checklist_parse_error: str | None = None
    if checklist_path.exists():
        try:
            checklist_rows_payload = read_tsv_rows(checklist_path)
        except (csv.Error, UnicodeDecodeError) as error:
            checklist_parse_error = str(error)
    else:
        checklist_parse_error = "reproducibility checklist is missing"
    checklist_item_count_matches = (
        checklist_entry.get("item_count") == len(checklist_rows_payload)
        if checklist_rows_payload
        else False
    )
    checklist_blocked = (
        checklist_row.status == "blocked"
        or checklist_parse_error is not None
        or not checklist_item_count_matches
    )
    checklist_evidence_parts = [
        f"stored checksum match={str(checklist_row.status == 'pass').lower()}",
    ]
    if checklist_parse_error is not None:
        checklist_evidence_parts.append(checklist_parse_error)
    else:
        checklist_evidence_parts.append(
            f"stored item_count={checklist_entry.get('item_count', '')}; observed rows={len(checklist_rows_payload)}"
        )
        checklist_evidence_parts.append(
            f"stored blocked_count={checklist_entry.get('blocked_count', '')}; stored risk_count={checklist_entry.get('risk_count', '')}"
        )
    check_rows.append(
        check_row(
            section="package",
            check_id="reproducibility-checklist-preserved",
            status=status(blocked=checklist_blocked),
            summary="reproducibility checklist still matches the stored package-control record",
            evidence="; ".join(checklist_evidence_parts),
            artifact_path=checklist_relative_path,
        )
    )

    inventory_artifact_rows: list[PublicationPackageRevalidationArtifactRow] = []
    for row in inventory_rows:
        relative_path = row.get("relative_path", "").strip()
        if not relative_path:
            continue
        inventory_artifact_rows.append(
            artifact_row(
                artifact_scope="inventory",
                section=row.get("section", "").strip() or "artifact",
                relative_path=relative_path,
                expected_sha256=row.get("sha256", "").strip() or None,
                expected_size_bytes=entry_size(
                    {"size_bytes": row.get("size_bytes", "")}
                ),
                observed_path=package_root / relative_path,
            )
        )
    artifact_rows.extend(inventory_artifact_rows)

    inventory_missing_count = sum(
        1 for row in inventory_artifact_rows if row.observed_sha256 is None
    )
    inventory_checksum_mismatch_count = sum(
        1
        for row in inventory_artifact_rows
        if row.expected_sha256 is not None
        and row.observed_sha256 is not None
        and row.expected_sha256 != row.observed_sha256
    )
    inventory_size_mismatch_count = sum(
        1
        for row in inventory_artifact_rows
        if row.expected_size_bytes is not None
        and row.observed_size_bytes is not None
        and row.expected_size_bytes != row.observed_size_bytes
    )
    inventory_artifact_blocked = (
        inventory_parse_error is not None
        or inventory_missing_count > 0
        or inventory_checksum_mismatch_count > 0
        or inventory_size_mismatch_count > 0
    )
    check_rows.append(
        check_row(
            section="artifacts",
            check_id="inventory-listed-artifacts-match",
            status=status(blocked=inventory_artifact_blocked),
            summary="all inventory-listed package inputs and outputs still match the stored package inventory",
            evidence=(
                f"missing={inventory_missing_count}; checksum_mismatches={inventory_checksum_mismatch_count}; "
                f"size_mismatches={inventory_size_mismatch_count}"
            ),
            artifact_path=inventory_relative_path,
        )
    )

    manifest_file_rows: list[PublicationPackageRevalidationArtifactRow] = []
    manifest_file_mismatch_count = 0
    for (
        block_name,
        _entry_name,
        relative_path,
        expected_sha256,
    ) in manifest_file_entries(manifest):
        row = artifact_row(
            artifact_scope="manifest_registry",
            section=block_name,
            relative_path=relative_path,
            expected_sha256=expected_sha256,
            expected_size_bytes=None,
            observed_path=package_root / relative_path,
        )
        if relative_path in {
            inventory_relative_path,
            checklist_relative_path,
        }:
            continue
        manifest_file_rows.append(row)
        if row.status == "blocked":
            manifest_file_mismatch_count += 1
    manifest_row = PublicationPackageRevalidationArtifactRow(
        artifact_scope="manifest",
        section="package",
        kind="manifest",
        relative_path=manifest_path.relative_to(package_root).as_posix(),
        status="pass",
        expected_sha256=None,
        observed_sha256=checksum(manifest_path),
        expected_size_bytes=None,
        observed_size_bytes=manifest_path.stat().st_size,
        detail="revalidation used this manifest as the stored package trust root",
    )
    artifact_rows.extend(manifest_file_rows)
    artifact_rows.append(manifest_row)
    check_rows.append(
        check_row(
            section="manifest",
            check_id="manifest-declared-files-match",
            status=status(blocked=manifest_file_mismatch_count > 0),
            summary="manifest-declared package files still match their stored checksums",
            evidence=f"blocked manifest-declared files={manifest_file_mismatch_count}",
            artifact_path=manifest_path.name,
        )
    )

    expected_relative_paths = {
        row.relative_path for row in artifact_rows if row.relative_path
    }
    unexpected_relative_paths = unexpected_files(
        package_root=package_root,
        expected_relative_paths=expected_relative_paths,
        output_root=output_root,
        ignored_prefixes=ignored_package_prefixes(report_kind),
    )
    check_rows.append(
        check_row(
            section="package",
            check_id="unexpected-package-files",
            status=status(risk=bool(unexpected_relative_paths)),
            summary="package root does not contain undeclared extra files",
            evidence=(
                "no unexpected files detected"
                if not unexpected_relative_paths
                else f"unexpected files: {' | '.join(unexpected_relative_paths[:10])}"
            ),
            artifact_path=manifest_path.name,
        )
    )

    matched_artifact_count = sum(1 for row in artifact_rows if row.status == "pass")
    missing_artifact_count = sum(
        1 for row in artifact_rows if row.observed_sha256 is None
    )
    checksum_mismatch_count = sum(
        1
        for row in artifact_rows
        if row.expected_sha256 is not None
        and row.observed_sha256 is not None
        and row.expected_sha256 != row.observed_sha256
    )
    size_mismatch_count = sum(
        1
        for row in artifact_rows
        if row.expected_size_bytes is not None
        and row.observed_size_bytes is not None
        and row.expected_size_bytes != row.observed_size_bytes
    )
    blocked_check_count = sum(1 for row in check_rows if row.status == "blocked")
    risk_check_count = sum(1 for row in check_rows if row.status == "risk")
    all_original_artifacts_match = blocked_check_count == 0
    overall_revalidation_status = derive_overall_revalidation_status(
        blocked_check_count=blocked_check_count,
        risk_check_count=risk_check_count,
    )

    artifact_table_path = write_taxon_rows(
        output_root / "publication-package-revalidation-artifacts.tsv",
        columns=_ARTIFACT_COLUMNS,
        rows=[asdict(row) for row in artifact_rows],
    )
    check_table_path = write_taxon_rows(
        output_root / "publication-package-revalidation-checks.tsv",
        columns=_CHECK_COLUMNS,
        rows=[asdict(row) for row in check_rows],
    )
    summary_path = output_root / "publication-package-revalidation-summary.json"
    report_path = output_root / "publication-package-revalidation-report.html"

    result = PublicationPackageRevalidationResult(
        output_root=output_root,
        manifest_path=manifest_path,
        package_root=package_root,
        report_kind=report_kind,
        artifact_table_path=artifact_table_path,
        check_table_path=check_table_path,
        summary_path=summary_path,
        report_path=report_path,
        artifact_rows=artifact_rows,
        check_rows=check_rows,
        matched_artifact_count=matched_artifact_count,
        missing_artifact_count=missing_artifact_count,
        checksum_mismatch_count=checksum_mismatch_count,
        size_mismatch_count=size_mismatch_count,
        unexpected_file_count=len(unexpected_relative_paths),
        blocked_check_count=blocked_check_count,
        risk_check_count=risk_check_count,
        all_original_artifacts_match=all_original_artifacts_match,
        overall_revalidation_status=overall_revalidation_status,
    )
    summary_payload = {
        "report_kind": report_kind,
        "manifest_path": str(manifest_path),
        "package_root": str(package_root),
        "matched_artifact_count": matched_artifact_count,
        "missing_artifact_count": missing_artifact_count,
        "checksum_mismatch_count": checksum_mismatch_count,
        "size_mismatch_count": size_mismatch_count,
        "unexpected_file_count": len(unexpected_relative_paths),
        "blocked_check_count": blocked_check_count,
        "risk_check_count": risk_check_count,
        "all_original_artifacts_match": all_original_artifacts_match,
        "overall_revalidation_status": overall_revalidation_status,
    }
    summary_path.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_html_report(report_path, result=result, manifest_path=manifest_path)
    return result
