from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from .book import evidence_book_root, write_evidence_book_index
from .coverage import (
    COVERAGE_GAPS_JSON,
    COVERAGE_GAPS_MARKDOWN,
    build_evidence_coverage_gap_report,
    encode_evidence_coverage_gap_report,
    render_evidence_coverage_gap_report,
)
from .freshness import (
    FRESHNESS_REPORT_JSON,
    FRESHNESS_REPORT_MARKDOWN,
    build_evidence_freshness_report,
    encode_evidence_freshness_report,
    render_evidence_freshness_report,
)
from .integrity import (
    INTEGRITY_REPORT_JSON,
    INTEGRITY_REPORT_MARKDOWN,
    build_evidence_integrity_report,
    encode_evidence_integrity_report,
    render_evidence_integrity_report,
)
from .reviewer import (
    REVIEWER_SUMMARY_JSON,
    REVIEWER_SUMMARY_MARKDOWN,
    build_bundle_reviewer_summary,
    encode_bundle_reviewer_summary,
    render_bundle_reviewer_summary,
)
from .study_registry import (
    EvidenceStudyBuildReport,
    EvidenceStudyRegistration,
    EvidenceStudyRerunReport,
    build_registered_study,
    rerun_selected_evidence,
    study_registrations,
)


DOCS_EVIDENCE_OVERVIEW = Path("docs") / "02-evidence-book" / "index.md"


@dataclass(slots=True)
class EvidenceBookRefreshReport:
    book_root: Path
    reviewer_summary_count: int
    updated_paths: list[Path]


@dataclass(slots=True)
class EvidenceBookStudyBuildReport:
    study_report: EvidenceStudyBuildReport
    refresh_report: EvidenceBookRefreshReport


@dataclass(slots=True)
class EvidenceBookStudyRerunReport:
    rerun_report: EvidenceStudyRerunReport
    refresh_report: EvidenceBookRefreshReport


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def _write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def write_bundle_reviewer_summaries(repo_root: Path | str) -> list[Path]:
    repo_root = Path(repo_root)
    book_root = evidence_book_root(repo_root)
    written_paths: list[Path] = []
    for bundle_root in sorted(book_root.glob("studies/*/evidence-*")):
        if not bundle_root.is_dir():
            continue
        study_manifest = _load_json(bundle_root.parent / "study.json")
        bundle_manifest = _load_json(bundle_root / "manifest.json")
        claims_path = bundle_root / "claims.json"
        claim_verdicts_path = bundle_root / "claim_verdicts.json"
        if claims_path.exists():
            claims_payload: dict[str, object] | list[object] | None = _load_json(
                claims_path
            )
        elif claim_verdicts_path.exists():
            claims_payload = json.loads(claim_verdicts_path.read_text(encoding="utf-8"))
        else:
            claims_payload = None
        payload = build_bundle_reviewer_summary(
            study_manifest=study_manifest,
            bundle_manifest=bundle_manifest,
            claims_payload=claims_payload,
        )
        reviewer_json_path = bundle_root / REVIEWER_SUMMARY_JSON
        reviewer_md_path = bundle_root / REVIEWER_SUMMARY_MARKDOWN
        _write_text(reviewer_json_path, encode_bundle_reviewer_summary(payload))
        _write_text(reviewer_md_path, render_bundle_reviewer_summary(payload))
        written_paths.extend([reviewer_json_path, reviewer_md_path])
    return written_paths


def render_docs_evidence_overview(
    *,
    index_payload: dict[str, object],
    teaching_payload: dict[str, object],
    freshness_payload: dict[str, object],
    integrity_payload: dict[str, object],
    coverage_payload: dict[str, object],
) -> str:
    lines = [
        "---",
        "title: Evidence Book",
        "audience: public",
        "type: evidence-navigation",
        "status: active",
        "owner: bijux-phylogenetics-dev",
        "last_reviewed: 2026-05-10",
        "---",
        "",
        "# Evidence Book",
        "",
        "The evidence-book is the repository trust surface for governed parity,",
        "teaching, migration, and explicit coverage boundaries.",
        "",
        f"- studies: `{index_payload['study_count']}`",
        f"- evidence bundles: `{index_payload['evidence_count']}`",
        f"- teaching studies: `{teaching_payload['teaching_study_count']}`",
        f"- migration studies: `{teaching_payload['migration_study_count']}`",
        f"- freshness statuses: `{', '.join(f'{key}={value}' for key, value in freshness_payload['freshness_status_counts'].items())}`",
        f"- coverage gaps: `{coverage_payload['coverage_gap_count']}` debt entries",
        "",
        "## Review Surfaces",
        "",
        "- repository evidence root: [`evidence-book/README.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/README.md)",
        "- teaching and migration index: [`evidence-book/index/teaching-and-migration.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/teaching-and-migration.md)",
        "- parity dashboard: [`evidence-book/index/parity-dashboard.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/parity-dashboard.md)",
        "- freshness report: [`evidence-book/index/freshness-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/freshness-report.md)",
        "- integrity report: [`evidence-book/index/integrity-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/integrity-report.md)",
        "- coverage gaps: [`evidence-book/index/coverage-gaps.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/coverage-gaps.md)",
        "",
        "## Studies",
        "",
    ]
    for study in index_payload["studies"]:
        lines.append(f"### {study['study_title']}")
        lines.append("")
        lines.append(f"- study id: `{study['study_id']}`")
        lines.append(f"- categories: `{', '.join(study['study_categories'])}`")
        lines.append(f"- bundles: `{study['bundle_count']}`")
        lines.append(
            f"- study manifest: [GitHub](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/studies/{study['study_id']}/study.json)"
        )
        lines.append("")
    lines.extend(
        [
            "## Integrity",
            "",
            f"- tracked bundles: `{integrity_payload['bundle_count']}`",
            f"- action required: `{integrity_payload['action_required_count']}`",
            "",
            "## Coverage",
            "",
            f"- family gaps: `{coverage_payload['family_gap_count']}`",
            f"- debt gaps: `{coverage_payload['coverage_gap_count']}`",
            "",
        ]
    )
    return "\n".join(lines)


def write_docs_evidence_overview(repo_root: Path | str) -> Path:
    repo_root = Path(repo_root)
    book_root = evidence_book_root(repo_root)
    docs_path = repo_root / DOCS_EVIDENCE_OVERVIEW
    index_payload = _load_json(book_root / "index" / "evidence-index.json")
    teaching_payload = _load_json(book_root / "index" / "teaching-and-migration.json")
    freshness_payload = _load_json(book_root / "index" / FRESHNESS_REPORT_JSON)
    integrity_payload = _load_json(book_root / "index" / INTEGRITY_REPORT_JSON)
    coverage_payload = _load_json(book_root / "index" / COVERAGE_GAPS_JSON)
    return _write_text(
        docs_path,
        render_docs_evidence_overview(
            index_payload=index_payload,
            teaching_payload=teaching_payload,
            freshness_payload=freshness_payload,
            integrity_payload=integrity_payload,
            coverage_payload=coverage_payload,
        ),
    )


def refresh_evidence_book(repo_root: Path | str) -> EvidenceBookRefreshReport:
    repo_root = Path(repo_root)
    book_root = evidence_book_root(repo_root)
    updated_paths: list[Path] = []
    updated_paths.extend(write_evidence_book_index(repo_root))
    reviewer_paths = write_bundle_reviewer_summaries(repo_root)
    updated_paths.extend(reviewer_paths)
    updated_paths.extend(write_evidence_book_index(repo_root))

    index_root = book_root / "index"
    freshness_payload = build_evidence_freshness_report(repo_root)
    integrity_payload = build_evidence_integrity_report(repo_root)
    coverage_payload = build_evidence_coverage_gap_report(repo_root)
    updated_paths.extend(
        [
            _write_text(
                index_root / FRESHNESS_REPORT_JSON,
                encode_evidence_freshness_report(freshness_payload),
            ),
            _write_text(
                index_root / FRESHNESS_REPORT_MARKDOWN,
                render_evidence_freshness_report(freshness_payload),
            ),
            _write_text(
                index_root / INTEGRITY_REPORT_JSON,
                encode_evidence_integrity_report(integrity_payload),
            ),
            _write_text(
                index_root / INTEGRITY_REPORT_MARKDOWN,
                render_evidence_integrity_report(integrity_payload),
            ),
            _write_text(
                index_root / COVERAGE_GAPS_JSON,
                encode_evidence_coverage_gap_report(coverage_payload),
            ),
            _write_text(
                index_root / COVERAGE_GAPS_MARKDOWN,
                render_evidence_coverage_gap_report(coverage_payload),
            ),
        ]
    )
    updated_paths.extend(write_evidence_book_index(repo_root))
    updated_paths.append(write_docs_evidence_overview(repo_root))
    return EvidenceBookRefreshReport(
        book_root=book_root,
        reviewer_summary_count=len(reviewer_paths) // 2,
        updated_paths=sorted(set(updated_paths)),
    )


def build_evidence_book_study(
    repo_root: Path | str, study_id: str
) -> EvidenceBookStudyBuildReport:
    report = build_registered_study(repo_root, study_id)
    refresh_report = refresh_evidence_book(repo_root)
    return EvidenceBookStudyBuildReport(
        study_report=report,
        refresh_report=refresh_report,
    )


def rerun_evidence_book_selection(
    repo_root: Path | str, study_id: str, evidence_ids: list[str]
) -> EvidenceBookStudyRerunReport:
    report = rerun_selected_evidence(repo_root, study_id, evidence_ids)
    refresh_report = refresh_evidence_book(repo_root)
    return EvidenceBookStudyRerunReport(
        rerun_report=report,
        refresh_report=refresh_report,
    )


def list_registered_evidence_studies(
    repo_root: Path | str,
) -> tuple[EvidenceStudyRegistration, ...]:
    return study_registrations(repo_root)
