from __future__ import annotations

from datetime import date
import json
from pathlib import Path

from .definitions import CLAIM_DEFINITIONS, EVIDENCE_ID, STUDY_ID
from .registry import _bundle_root, build_primate_source_fragment_map


def build_primate_summary_bundle_claims(repo_root: Path) -> dict[str, object]:
    fragment_map = build_primate_source_fragment_map(repo_root)
    fragments = fragment_map["fragments"]
    claims: list[dict[str, object]] = []
    for claim_id, definition in CLAIM_DEFINITIONS.items():
        claim_fragments = [
            fragment["fragment_id"]
            for fragment in fragments
            if claim_id in fragment["claim_ids"]
        ]
        claims.append(
            {
                "claim_id": claim_id,
                "claim_title": definition["claim_title"],
                "summary": definition["summary"],
                "verdict": definition["verdict"],
                "evidence_ids": [EVIDENCE_ID],
                "source_fragments": claim_fragments,
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "evidence_id": EVIDENCE_ID,
        "claim_count": len(claims),
        "claims": claims,
    }


def refresh_primate_summary_bundle(repo_root: Path) -> dict[str, object]:
    bundle_root = _bundle_root(Path(repo_root))
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("expected summary bundle manifest to be a JSON object")
    claims_payload = build_primate_summary_bundle_claims(Path(repo_root))
    manifest["claim_ids"] = [
        claim["claim_id"]
        for claim in claims_payload["claims"]
        if isinstance(claim, dict) and isinstance(claim.get("claim_id"), str)
    ]
    freshness = manifest.get("freshness")
    if isinstance(freshness, dict):
        freshness["last_generated_on"] = date.today().isoformat()
    return {"manifest": manifest, "claims": claims_payload}
