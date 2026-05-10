from __future__ import annotations

import json
from pathlib import Path


COVERAGE_GAPS_JSON = "coverage-gaps.json"
COVERAGE_GAPS_MARKDOWN = "coverage-gaps.md"


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def build_evidence_coverage_gap_report(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    book_root = repo_root / "evidence-book"
    debts_path = book_root / "index" / "scientific-debt-register.json"
    debt_payload = _load_json(debts_path)
    gap_entries = [
        entry
        for entry in debt_payload.get("debts", [])
        if isinstance(entry, dict) and entry.get("debt_kind") == "coverage_gap"
    ]
    family_gaps: list[dict[str, object]] = []
    return {
        "schema_version": 1,
        "coverage_gap_count": len(gap_entries),
        "family_gap_count": len(family_gaps),
        "coverage_gaps": gap_entries,
        "family_gaps": family_gaps,
    }


def render_evidence_coverage_gap_report(payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Coverage Gaps",
        "",
        f"- coverage gap debts: `{payload['coverage_gap_count']}`",
        f"- family gaps: `{payload['family_gap_count']}`",
        "",
        "## Family Gaps",
        "",
    ]
    for gap in payload["family_gaps"]:
        lines.append(
            f"- `{gap['study_id']}/{gap['family_id']}` — `{gap['coverage_status']}`"
        )
        for detail in gap["known_gaps"]:
            lines.append(f"  Detail: {detail}")
    lines.extend(
        [
            "",
            "## Debt Register Gaps",
            "",
        ]
    )
    for gap in payload["coverage_gaps"]:
        lines.append(f"- `{gap['debt_id']}`")
        if gap.get("detail"):
            lines.append(f"  Detail: {gap['detail']}")
    lines.append("")
    return "\n".join(lines)


def encode_evidence_coverage_gap_report(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"
