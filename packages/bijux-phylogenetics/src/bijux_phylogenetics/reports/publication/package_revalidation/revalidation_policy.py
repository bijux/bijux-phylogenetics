from __future__ import annotations

from .contracts import PublicationPackageRevalidationCheckRow


def status(*, blocked: bool = False, risk: bool = False) -> str:
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
    artifact_path: str,
) -> PublicationPackageRevalidationCheckRow:
    return PublicationPackageRevalidationCheckRow(
        section=section,
        check_id=check_id,
        status=status,
        summary=summary,
        evidence=evidence,
        artifact_path=artifact_path,
    )


def overall_revalidation_status(
    *,
    blocked_check_count: int,
    risk_check_count: int,
) -> str:
    if blocked_check_count > 0:
        return "blocked"
    if risk_check_count > 0:
        return "risk"
    return "pass"
