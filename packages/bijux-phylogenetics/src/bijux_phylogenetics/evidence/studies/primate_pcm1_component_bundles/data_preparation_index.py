from __future__ import annotations

from pathlib import Path

from .bundle_builder import build_primate_pcm1_component_bundles
from .definitions import DATA_PREPARATION_BUNDLE_IDS, STUDY_ID


def build_primate_data_preparation_bundle_index(repo_root: Path) -> dict[str, object]:
    bundles = build_primate_pcm1_component_bundles(repo_root)
    entries: list[dict[str, object]] = []
    for evidence_id in DATA_PREPARATION_BUNDLE_IDS:
        manifest = bundles[evidence_id]["manifest"]
        entries.append(
            {
                "evidence_id": evidence_id,
                "title": manifest["evidence_title"],
                "claim_id": manifest["claim_ids"][0],
                "relative_path": (
                    Path("studies") / "primate-longevity-signal" / evidence_id
                ).as_posix(),
                "claim_tags": manifest["claim_tags"],
                "source_fragments": manifest["source_fragments"],
            }
        )
    return {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "bundle_count": len(entries),
        "bundle_ids": DATA_PREPARATION_BUNDLE_IDS,
        "bundles": entries,
        "review_rule": (
            "Downstream lambda and ancestral-state evidence must rely on these preprocessing bundles instead of hiding data preparation inside later model comparisons."
        ),
    }


def render_primate_data_preparation_bundle_index_markdown(
    payload: dict[str, object],
) -> str:
    lines = [
        "# Primate Data-Preparation Parity Bundles",
        "",
        payload["review_rule"],
        "",
        f"Bundles: `{payload['bundle_count']}`",
        "",
    ]
    for entry in payload["bundles"]:
        lines.append(
            f"- `{entry['evidence_id']}` — {entry['title']} "
            f"(`{', '.join(entry['source_fragments'])}`)"
        )
    lines.append("")
    return "\n".join(lines)
