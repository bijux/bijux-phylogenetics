from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.closure import (
    build_closure_criteria,
    build_evidence_maturity_scorecard,
)
from bijux_phylogenetics.evidence.coverage import build_evidence_coverage_gap_report
from bijux_phylogenetics.evidence.freshness import build_evidence_freshness_report
from bijux_phylogenetics.evidence.integrity import build_evidence_integrity_report
from bijux_phylogenetics.evidence.reviewer import (
    build_bundle_reviewer_summary,
    render_bundle_reviewer_summary,
)
from bijux_phylogenetics.evidence.workbench import (
    DOCS_EVIDENCE_OVERVIEW,
    render_docs_evidence_overview,
)


REPO_ROOT = Path(__file__).resolve().parents[3]


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise AssertionError(f"expected JSON object at {path}")
    return payload


def test_repository_reviewer_summaries_match_generated_payloads() -> None:
    studies_root = REPO_ROOT / "evidence-book" / "studies"

    for study_root in sorted(path for path in studies_root.iterdir() if path.is_dir()):
        study_manifest = _load_json(study_root / "study.json")
        for bundle_root in sorted(study_root.glob("evidence-*")):
            bundle_manifest = _load_json(bundle_root / "manifest.json")
            claims_path = bundle_root / "claims.json"
            claim_verdicts_path = bundle_root / "claim_verdicts.json"
            if claims_path.exists():
                claims_payload: dict[str, object] | list[object] | None = _load_json(
                    claims_path
                )
            elif claim_verdicts_path.exists():
                claims_payload = json.loads(
                    claim_verdicts_path.read_text(encoding="utf-8")
                )
            else:
                claims_payload = None
            reviewer_payload = build_bundle_reviewer_summary(
                study_manifest=study_manifest,
                bundle_manifest=bundle_manifest,
                claims_payload=claims_payload,
            )

            assert (
                json.loads(
                    (bundle_root / "reviewer-summary.json").read_text(encoding="utf-8")
                )
                == reviewer_payload
            )
            assert (bundle_root / "reviewer-summary.md").read_text(
                encoding="utf-8"
            ) == render_bundle_reviewer_summary(reviewer_payload)


def test_repository_docs_evidence_overview_matches_generated_payload() -> None:
    index_root = REPO_ROOT / "evidence-book" / "index"
    index_payload = _load_json(index_root / "evidence-index.json")
    teaching_payload = _load_json(index_root / "teaching-and-migration.json")
    freshness_payload = build_evidence_freshness_report(REPO_ROOT)
    integrity_payload = build_evidence_integrity_report(REPO_ROOT)
    coverage_payload = build_evidence_coverage_gap_report(REPO_ROOT)
    closure_payload = build_closure_criteria(REPO_ROOT)
    scorecard_payload = build_evidence_maturity_scorecard(REPO_ROOT)

    overview = render_docs_evidence_overview(
        index_payload=index_payload,
        teaching_payload=teaching_payload,
        freshness_payload=freshness_payload,
        integrity_payload=integrity_payload,
        coverage_payload=coverage_payload,
        closure_payload=closure_payload,
        scorecard_payload=scorecard_payload,
    )

    assert (REPO_ROOT / DOCS_EVIDENCE_OVERVIEW).read_text(encoding="utf-8") == overview
