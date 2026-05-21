from __future__ import annotations

from .contracts import PublicationPackageComparisonCheckRow


def status(*, blocked: bool = False, risk: bool = False) -> str:
    """Return the governed comparison status for one reviewer-facing check."""
    if blocked:
        return "blocked"
    if risk:
        return "risk"
    return "pass"


def check_row(
    *,
    section: str,
    check_id: str,
    status: str,
    summary: str,
    evidence: str,
    left_artifact_path: str,
    right_artifact_path: str,
) -> PublicationPackageComparisonCheckRow:
    """Build one reviewer-facing comparison check row."""
    return PublicationPackageComparisonCheckRow(
        section=section,
        check_id=check_id,
        status=status,
        summary=summary,
        evidence=evidence,
        left_artifact_path=left_artifact_path,
        right_artifact_path=right_artifact_path,
    )


def config_differences(
    left_manifest: dict[str, object],
    right_manifest: dict[str, object],
    *,
    mapping,
) -> dict[str, tuple[object, object]]:
    """Build the governed config-difference index across both manifests."""
    left_config = mapping(left_manifest, "config")
    right_config = mapping(right_manifest, "config")
    differences: dict[str, tuple[object, object]] = {}
    for key in sorted(set(left_config) | set(right_config)):
        if key in {"path", "checksum"}:
            continue
        if left_config.get(key) != right_config.get(key):
            differences[key] = (left_config.get(key), right_config.get(key))
    return differences


def finding_difference_count(
    left_rows: dict[str, dict[str, str]],
    right_rows: dict[str, dict[str, str]],
) -> int:
    """Count biological finding differences across both package versions."""
    count = 0
    for finding_id in sorted(set(left_rows) | set(right_rows)):
        left_row = left_rows.get(finding_id)
        right_row = right_rows.get(finding_id)
        if left_row is None or right_row is None:
            count += 1
            continue
        if {key: value for key, value in left_row.items() if key != "finding_id"} != {
            key: value for key, value in right_row.items() if key != "finding_id"
        }:
            count += 1
    return count
