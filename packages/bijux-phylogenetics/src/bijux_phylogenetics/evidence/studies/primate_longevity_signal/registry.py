from __future__ import annotations

import json
from pathlib import Path

from .definitions import (
    EVIDENCE_ID,
    FAMILY_DEFINITIONS,
    FRAGMENT_CLASSIFICATIONS,
    SOURCE_LOCATOR,
    STUDY_ID,
)


def build_primate_source_fragment_map(repo_root: Path) -> dict[str, object]:
    comparison = _comparison_payload(repo_root)
    supporting_bundles = _component_evidence_by_fragment(repo_root)
    fragments: list[dict[str, object]] = []
    for row in comparison["comparisons"]:
        block_id = row["block_id"]
        classification = FRAGMENT_CLASSIFICATIONS[block_id]
        fragments.append(
            {
                "fragment_id": block_id,
                "fragment_title": row["title"],
                "script_line_spec": row["script_lines"],
                "script_line_spans": _line_spans(row["script_lines"]),
                "script_locators": _line_locators(row["script_lines"]),
                "evidence_id": EVIDENCE_ID,
                "concept_family": classification["concept_family"],
                "claim_ids": classification["claim_ids"],
                "parity_expectation": classification["parity_expectation"],
                "scope": classification["scope"],
                "block_status": row["status"],
                "supporting_evidence_ids": supporting_bundles.get(block_id, []),
                "review_note": row["note"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "source_locator": SOURCE_LOCATOR,
        "fragment_count": len(fragments),
        "fragments": fragments,
    }


def build_primate_claim_registry(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    fragment_map["fragments"]
    from .bundles import build_primate_summary_bundle_claims

    claims = build_primate_summary_bundle_claims(repo_root)["claims"]
    claims.extend(_component_claim_rows(repo_root))
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "claim_count": len(claims),
        "claims": claims,
    }


def build_primate_family_index(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    claim_registry = build_primate_claim_registry(repo_root)
    fragments = fragment_map["fragments"]
    claims = claim_registry["claims"]
    families: list[dict[str, object]] = []
    for family_id, definition in FAMILY_DEFINITIONS.items():
        family_fragments = [
            fragment
            for fragment in fragments
            if fragment["concept_family"] == family_id
        ]
        fragment_ids = {fragment["fragment_id"] for fragment in family_fragments}
        family_claim_ids = sorted(
            {
                claim_id
                for fragment in family_fragments
                for claim_id in fragment["claim_ids"]
            }
            | {
                claim["claim_id"]
                for claim in claims
                if fragment_ids & set(claim.get("source_fragments", []))
            }
        )
        supporting_evidence_ids = sorted(
            {
                evidence_id
                for fragment in family_fragments
                for evidence_id in fragment.get("supporting_evidence_ids", [])
            }
        )
        claim_titles = {
            claim["claim_id"]: claim["claim_title"]
            for claim in claims
            if claim["claim_id"] in family_claim_ids
        }
        families.append(
            {
                "family_id": family_id,
                "family_title": definition["title"],
                "summary": definition["summary"],
                "evidence_ids": [EVIDENCE_ID, *supporting_evidence_ids],
                "fragment_ids": sorted(fragment_ids),
                "claim_ids": family_claim_ids,
                "claim_titles": claim_titles,
                "fragment_count": len(family_fragments),
                "family_verdict": _family_verdict(
                    [fragment["block_status"] for fragment in family_fragments]
                ),
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "source_family_id": "pcm1-plots-signal",
        "source_family_title": "PCM1 plots signal evidence family",
        "evidence_id": EVIDENCE_ID,
        "family_count": len(families),
        "families": families,
    }


def _study_root(repo_root: Path) -> Path:
    return Path(repo_root) / "evidence-book" / "studies" / STUDY_ID


def _bundle_root(repo_root: Path) -> Path:
    return _study_root(repo_root) / EVIDENCE_ID


def _comparison_payload(repo_root: Path) -> dict[str, object]:
    bundle_root = _bundle_root(repo_root)
    path = bundle_root / "comparison.json"
    if not path.exists():
        path = bundle_root / "results" / "comparison.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _component_evidence_by_fragment(repo_root: Path) -> dict[str, list[str]]:
    from ..primate_pcm1_component_bundles import build_primate_pcm1_component_bundles

    supporting_bundles: dict[str, set[str]] = {}
    for evidence_id, payload in build_primate_pcm1_component_bundles(repo_root).items():
        for fragment_id in payload["manifest"]["source_fragments"]:
            supporting_bundles.setdefault(fragment_id, set()).add(evidence_id)
    return {
        fragment_id: sorted(evidence_ids)
        for fragment_id, evidence_ids in supporting_bundles.items()
    }


def _component_claim_rows(repo_root: Path) -> list[dict[str, object]]:
    from ..primate_pcm1_component_bundles import build_primate_pcm1_component_bundles

    rows: list[dict[str, object]] = []
    for evidence_id, payload in sorted(
        build_primate_pcm1_component_bundles(repo_root).items()
    ):
        claim = payload["claims"]["claims"][0]
        rows.append(
            {
                "claim_id": claim["claim_id"],
                "claim_title": claim["claim_title"],
                "summary": claim["summary"],
                "verdict": claim["verdict"],
                "evidence_ids": [evidence_id],
                "source_fragments": claim["source_fragments"],
            }
        )
    return rows


def _line_spans(spec: str) -> list[dict[str, int]]:
    spans: list[dict[str, int]] = []
    for raw_part in spec.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            start_text, end_text = part.split("-", maxsplit=1)
            start = int(start_text)
            end = int(end_text)
        else:
            start = int(part)
            end = int(part)
        spans.append({"start_line": start, "end_line": end})
    return spans


def _line_locators(spec: str) -> list[str]:
    locators: list[str] = []
    for span in _line_spans(spec):
        if span["start_line"] == span["end_line"]:
            locators.append(f"{SOURCE_LOCATOR}#L{span['start_line']}")
        else:
            locators.append(
                f"{SOURCE_LOCATOR}#L{span['start_line']}-L{span['end_line']}"
            )
    return locators


def _family_verdict(fragment_statuses: list[str]) -> str:
    if all(
        status in {"workflow_only", "artifact_only", "plot_only"}
        for status in fragment_statuses
    ):
        return "not_comparable"
    if "verified_with_tolerance" in fragment_statuses:
        return "matched_with_tolerance"
    return "matched"
