from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from .book import evidence_book_root, write_evidence_book_index
from .closure import (
    ANALYTICAL_SURFACE_COVERAGE_JSON,
    ANALYTICAL_SURFACE_COVERAGE_MARKDOWN,
    CLAIM_REAUDIT_JSON,
    CLAIM_REAUDIT_MARKDOWN,
    CLOSURE_CRITERIA_JSON,
    CLOSURE_CRITERIA_MARKDOWN,
    COMPLETION_GATES_JSON,
    COMPLETION_GATES_MARKDOWN,
    EVIDENCE_MATURITY_SCORECARD_JSON,
    EVIDENCE_MATURITY_SCORECARD_MARKDOWN,
    EVIDENCE_REVIEW_RITUAL_JSON,
    EVIDENCE_REVIEW_RITUAL_MARKDOWN,
    build_analytical_surface_coverage,
    build_claim_reaudit,
    build_closure_criteria,
    build_completion_gates,
    build_evidence_maturity_scorecard,
    build_evidence_review_ritual,
    render_analytical_surface_coverage,
    render_claim_reaudit,
    render_closure_criteria,
    render_completion_gates,
    render_evidence_maturity_scorecard,
    render_evidence_review_ritual,
)
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
from .study_contracts import load_study_contract
from .study_registry import (
    EvidenceStudyBuildReport,
    EvidenceStudyRegistration,
    EvidenceStudyRerunReport,
    build_registered_study,
    rerun_selected_evidence,
    study_registrations,
)

DOCS_EVIDENCE_OVERVIEW = (
    Path("docs") / "public" / "phylogenetics-evidence-book" / "index.md"
)


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


@dataclass(slots=True)
class EvidenceBookSelectionBuildReport:
    study_id: str
    selected_evidence_ids: list[str]
    updated_paths: list[str]
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
        study_manifest = load_study_contract(bundle_root.parent)
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
    closure_payload: dict[str, object],
    scorecard_payload: dict[str, object],
) -> str:
    foundational_status = next(
        criterion["current_status"]
        for criterion in closure_payload["criteria"]
        if criterion["criterion_id"] == "foundational-numerical-trust"
    )
    reviewer_status = next(
        criterion["current_status"]
        for criterion in closure_payload["criteria"]
        if criterion["criterion_id"] == "reviewer-readiness"
    )
    study_page_by_id = {
        "primate-longevity-signal": "pcm1-primate-longevity-signal",
        "primate-pgls-and-signal": "pcm2-primate-pgls-and-signal",
        "primate-continuous-and-categorical-regression": (
            "pcm3-primate-continuous-and-categorical-regression"
        ),
        "cichlid-bayesian-linear-and-phylogenetic-mixed-models": (
            "pcm4-cichlid-bayesian-linear-and-phylogenetic-mixed-models"
        ),
        "cichlid-phylogenetic-generalized-linear-mixed-models": (
            "pcm5-cichlid-phylogenetic-generalized-linear-mixed-models"
        ),
    }
    study_question_by_id = {
        "primate-longevity-signal": (
            "How do primate traits and tree preparation become a usable signal and"
            " ancestral-state workflow?"
        ),
        "primate-pgls-and-signal": (
            "How do baseline GLS, PGLS, transformed trees, and signal tests behave"
            " on the same primate dataset?"
        ),
        "primate-continuous-and-categorical-regression": (
            "How should continuous predictors, categorical predictors, and"
            " interactions be interpreted together in comparative regression?"
        ),
        "cichlid-bayesian-linear-and-phylogenetic-mixed-models": (
            "Where does the cichlid workflow move from standard regression into"
            " Bayesian mixed models with a phylogenetic effect?"
        ),
        "cichlid-phylogenetic-generalized-linear-mixed-models": (
            "How do count models, binary-style models, repeated measures, and"
            " posterior signal estimates fit together in one study family?"
        ),
    }
    lines = [
        "---",
        "title: Evidence Book",
        "audience: public",
        "type: evidence-navigation",
        "status: active",
        "owner: bijux-phylogenetics-dev",
        "last_reviewed: 2026-07-04",
        "---",
        "",
        "# Evidence Book",
        "",
        "This section explains the comparative study record behind the public",
        "phylogenetics runtime.",
        "",
        "The five study families come from course material that the repository has",
        "turned into checked-in, reviewable evidence. All five now have governed",
        "parity-style bundle IDs. Some studies already carry stronger numerical",
        "parity records, while others keep parity verdicts open because the",
        "scientific workflow still depends on an external engine or an unclosed",
        "lecture-specific boundary.",
        "",
        "```mermaid",
        "flowchart LR",
        "    A[Source lecture] --> B[Public study guide]",
        "    B --> C[Checked-in evidence dossier]",
        "    C --> D[Freshness and integrity reports]",
        "    C --> E[Study-specific evidence bundles]",
        "```",
        "",
        "## Start Here",
        "",
        "- [evidence book reading guide](reading-guide.md)",
        "- checked-in study dossiers:",
        "  [`evidence-book/README.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/README.md)",
        "- repository-wide review reports:",
        "  [`evidence-book/index/`](https://github.com/bijux/bijux-phylogenetics/tree/main/evidence-book/index)",
        "",
        "Treat those checked-in reports as the current source of truth when you need",
        "live repository status rather than the stable narrative on this page.",
        "",
        "## At A Glance",
        "",
        f"- studies: `{index_payload['study_count']}`",
        f"- evidence bundles: `{index_payload['evidence_count']}`",
        f"- teaching studies: `{teaching_payload['teaching_study_count']}`",
        f"- migration studies: `{teaching_payload['migration_study_count']}`",
        f"- freshness statuses: `{', '.join(f'{key}={value}' for key, value in freshness_payload['freshness_status_counts'].items())}`",
        f"- coverage gaps still tracked: `{coverage_payload['coverage_gap_count']}`",
        f"- foundational numerical trust: `{foundational_status}`",
        f"- reviewer readiness: `{reviewer_status}`",
        f"- maturity tier: `{scorecard_payload['maturity_tier']}`",
        "",
        "## Study Atlas",
        "",
        "| Study | Main question | Evidence IDs | Current repository state |",
        "| --- | --- | --- | --- |",
    ]
    for study in index_payload["studies"]:
        study_id = str(study["study_id"])
        study_page = study_page_by_id[study_id]
        state = "Documented study family with open limits"
        evidence_ids = [str(entry["evidence_id"]) for entry in study["evidence"]]
        evidence_span = (
            f"`{evidence_ids[0]}` to `{evidence_ids[-1]}` ({len(evidence_ids)} bundles)"
            if evidence_ids
            else "none"
        )
        if study_id == "primate-longevity-signal":
            state = "Governed parity bundles deployed; plotting boundary remains explicit"
        elif study_id == "primate-pgls-and-signal":
            state = "Governed parity bundles deployed; one coverage boundary remains explicit"
        elif study_id == "primate-continuous-and-categorical-regression":
            state = "Governed parity bundle program deployed; lecture-specific verdicts remain open"
        elif study_id.startswith("cichlid-"):
            state = "Governed parity bundle program deployed; posterior workflow still external"
        lines.append(
            f"| [{study['study_title']}](./{study_page}/index.md) | {study_question_by_id[study_id]} | {evidence_span} | {state} |"
        )
    lines.extend(
        [
            "",
            "## When To Use These Pages",
            "",
            "- Use the public study guides when you want the scientific story in plain language.",
            "- Use the checked-in evidence-book when you need provenance, registries, or bundle-level proof.",
            "- Use the index reports when you need freshness, integrity, or repository-wide gap tracking.",
            "",
            "## Review Reports",
            "",
            "- teaching and migration index: [`evidence-book/index/teaching-and-migration.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/teaching-and-migration.md)",
            "- parity dashboard: [`evidence-book/index/parity-dashboard.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/parity-dashboard.md)",
            "- freshness report: [`evidence-book/index/freshness-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/freshness-report.md)",
            "- integrity report: [`evidence-book/index/integrity-report.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/integrity-report.md)",
            "- coverage gaps: [`evidence-book/index/coverage-gaps.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/coverage-gaps.md)",
            "- claim re-audit: [`evidence-book/index/claim-reaudit.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/claim-reaudit.md)",
            "- closure criteria: [`evidence-book/index/closure-criteria.md`](https://github.com/bijux/bijux-phylogenetics/blob/main/evidence-book/index/closure-criteria.md)",
            "",
            "The study record is useful because it keeps strong results visible without",
            "pretending that every workflow is fully reproduced inside the repository.",
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
    closure_payload = _load_json(book_root / "index" / CLOSURE_CRITERIA_JSON)
    scorecard_payload = _load_json(
        book_root / "index" / EVIDENCE_MATURITY_SCORECARD_JSON
    )
    return _write_text(
        docs_path,
        render_docs_evidence_overview(
            index_payload=index_payload,
            teaching_payload=teaching_payload,
            freshness_payload=freshness_payload,
            integrity_payload=integrity_payload,
            coverage_payload=coverage_payload,
            closure_payload=closure_payload,
            scorecard_payload=scorecard_payload,
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
    claim_reaudit_payload = build_claim_reaudit(repo_root)
    analytical_surface_payload = build_analytical_surface_coverage(repo_root)
    closure_payload = build_closure_criteria(repo_root)
    maturity_scorecard_payload = build_evidence_maturity_scorecard(repo_root)
    review_ritual_payload = build_evidence_review_ritual(repo_root)
    completion_gates_payload = build_completion_gates(repo_root)
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
            _write_text(
                index_root / CLAIM_REAUDIT_JSON,
                json.dumps(claim_reaudit_payload, indent=2, sort_keys=True) + "\n",
            ),
            _write_text(
                index_root / CLAIM_REAUDIT_MARKDOWN,
                render_claim_reaudit(claim_reaudit_payload),
            ),
            _write_text(
                index_root / ANALYTICAL_SURFACE_COVERAGE_JSON,
                json.dumps(analytical_surface_payload, indent=2, sort_keys=True) + "\n",
            ),
            _write_text(
                index_root / ANALYTICAL_SURFACE_COVERAGE_MARKDOWN,
                render_analytical_surface_coverage(analytical_surface_payload),
            ),
            _write_text(
                index_root / CLOSURE_CRITERIA_JSON,
                json.dumps(closure_payload, indent=2, sort_keys=True) + "\n",
            ),
            _write_text(
                index_root / CLOSURE_CRITERIA_MARKDOWN,
                render_closure_criteria(closure_payload),
            ),
            _write_text(
                index_root / EVIDENCE_MATURITY_SCORECARD_JSON,
                json.dumps(maturity_scorecard_payload, indent=2, sort_keys=True) + "\n",
            ),
            _write_text(
                index_root / EVIDENCE_MATURITY_SCORECARD_MARKDOWN,
                render_evidence_maturity_scorecard(maturity_scorecard_payload),
            ),
            _write_text(
                index_root / EVIDENCE_REVIEW_RITUAL_JSON,
                json.dumps(review_ritual_payload, indent=2, sort_keys=True) + "\n",
            ),
            _write_text(
                index_root / EVIDENCE_REVIEW_RITUAL_MARKDOWN,
                render_evidence_review_ritual(review_ritual_payload),
            ),
            _write_text(
                index_root / COMPLETION_GATES_JSON,
                json.dumps(completion_gates_payload, indent=2, sort_keys=True) + "\n",
            ),
            _write_text(
                index_root / COMPLETION_GATES_MARKDOWN,
                render_completion_gates(completion_gates_payload),
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


def build_evidence_book_selection(
    repo_root: Path | str, study_id: str, evidence_ids: list[str]
) -> EvidenceBookSelectionBuildReport:
    rerun_report = rerun_selected_evidence(repo_root, study_id, evidence_ids)
    refresh_report = refresh_evidence_book(repo_root)
    return EvidenceBookSelectionBuildReport(
        study_id=rerun_report.study_id,
        selected_evidence_ids=rerun_report.selected_evidence_ids,
        updated_paths=rerun_report.updated_paths,
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
