from __future__ import annotations

import json


REVIEWER_SUMMARY_JSON = "results/reviewer-summary.json"
REVIEWER_SUMMARY_MARKDOWN = "results/reviewer-summary.md"


def build_bundle_reviewer_summary(
    *,
    study_manifest: dict[str, object],
    bundle_manifest: dict[str, object],
    claims_payload: dict[str, object] | list[object] | None,
) -> dict[str, object]:
    claim_rows = []
    if isinstance(claims_payload, dict):
        rows = claims_payload.get("claims", [])
        if isinstance(rows, list):
            claim_rows = [row for row in rows if isinstance(row, dict)]
    elif isinstance(claims_payload, list):
        claim_rows = [row for row in claims_payload if isinstance(row, dict)]
    claim_titles = [
        str(row.get("claim_title", row.get("claim_id", "claim"))) for row in claim_rows
    ]
    fragment_ids = sorted(
        {
            str(fragment_id)
            for row in claim_rows
            for fragment_id in row.get("source_fragments", [])
            if isinstance(fragment_id, str)
        }
    )
    limitations = [
        str(item)
        for item in bundle_manifest.get("limitations", [])
        if isinstance(item, str)
    ]
    source_locators = [
        str(entry["locator"])
        for entry in bundle_manifest.get("source_basis", [])
        if isinstance(entry, dict) and isinstance(entry.get("locator"), str)
    ]
    payload = {
        "schema_version": 1,
        "study_id": study_manifest["study_id"],
        "study_title": study_manifest["study_title"],
        "evidence_id": bundle_manifest["evidence_id"],
        "evidence_title": bundle_manifest["evidence_title"],
        "owner_package": bundle_manifest["owner_package"],
        "comparison_mode": bundle_manifest["comparison_mode"],
        "bundle_verdict_status": bundle_manifest["verdict"]["status"],
        "bundle_verdict_summary": bundle_manifest["verdict"]["summary"],
        "claim_count": len(claim_rows),
        "claim_titles": claim_titles,
        "claim_tags": bundle_manifest["claim_tags"],
        "source_locators": source_locators,
        "source_fragment_ids": fragment_ids,
        "limitations": limitations,
        "reviewer_summary_lines": [
            str(bundle_manifest["summary"]),
            (
                f"Comparison mode: {bundle_manifest['comparison_mode']}."
                f" Bundle verdict: {bundle_manifest['verdict']['status']}."
            ),
            (
                "Claims stay explicit and reviewer-readable through governed claim rows, "
                "portable source locators, and tracked limitations."
            ),
        ],
    }
    return payload


def render_bundle_reviewer_summary(payload: dict[str, object]) -> str:
    lines = [
        f"# {payload['evidence_title']} Reviewer Summary",
        "",
        f"- study: `{payload['study_id']}`",
        f"- evidence id: `{payload['evidence_id']}`",
        f"- comparison mode: `{payload['comparison_mode']}`",
        f"- verdict: `{payload['bundle_verdict_status']}`",
        f"- claims: `{payload['claim_count']}`",
        "",
        "## Summary",
        "",
    ]
    lines.extend(f"- {line}" for line in payload["reviewer_summary_lines"])
    lines.extend(
        [
            "",
            "## Claims",
            "",
        ]
    )
    if payload["claim_titles"]:
        lines.extend(f"- {title}" for title in payload["claim_titles"])
    else:
        lines.append("- No governed claim rows were recorded for this bundle.")
    lines.extend(
        [
            "",
            "## Source Fragments",
            "",
        ]
    )
    if payload["source_fragment_ids"]:
        lines.extend(
            f"- `{fragment_id}`" for fragment_id in payload["source_fragment_ids"]
        )
    else:
        lines.append("- No source fragments were recorded for this bundle.")
    lines.extend(
        [
            "",
            "## Limitations",
            "",
        ]
    )
    if payload["limitations"]:
        lines.extend(f"- {item}" for item in payload["limitations"])
    else:
        lines.append("- No additional limitations were recorded for this bundle.")
    lines.append("")
    return "\n".join(lines)


def encode_bundle_reviewer_summary(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"
