from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

from .book import (
    build_evidence_book_index,
    build_evidence_false_confidence_audit,
    build_evidence_mismatch_archive,
)
from .coverage import build_evidence_coverage_gap_report
from .freshness import build_evidence_freshness_report
from .integrity import build_evidence_integrity_report


CLAIM_REAUDIT_JSON = "claim-reaudit.json"
CLAIM_REAUDIT_MARKDOWN = "claim-reaudit.md"
ANALYTICAL_SURFACE_COVERAGE_JSON = "analytical-surface-coverage.json"
ANALYTICAL_SURFACE_COVERAGE_MARKDOWN = "analytical-surface-coverage.md"
CLOSURE_CRITERIA_JSON = "closure-criteria.json"
CLOSURE_CRITERIA_MARKDOWN = "closure-criteria.md"
EVIDENCE_MATURITY_SCORECARD_JSON = "evidence-maturity-scorecard.json"
EVIDENCE_MATURITY_SCORECARD_MARKDOWN = "evidence-maturity-scorecard.md"
EVIDENCE_REVIEW_RITUAL_JSON = "evidence-review-ritual.json"
EVIDENCE_REVIEW_RITUAL_MARKDOWN = "evidence-review-ritual.md"
COMPLETION_GATES_JSON = "completion-gates.json"
COMPLETION_GATES_MARKDOWN = "completion-gates.md"


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object at {path}")
    return payload


def _bundle_entry_map(repo_root: Path) -> dict[tuple[str, str], dict[str, object]]:
    index_payload = build_evidence_book_index(repo_root)
    return {
        (str(entry["study_id"]), str(entry["evidence_id"])): entry
        for entry in index_payload["evidence"]
        if isinstance(entry, dict)
    }


def _surface_definitions() -> tuple[dict[str, object], ...]:
    return (
        {
            "surface_id": "tree-diagnostics",
            "surface_title": "Tree diagnostics and tree-data correspondence",
            "runtime_surfaces": [
                "validate",
                "inspect",
                "topology",
                "dataset crosswalk",
            ],
            "linked_evidence_ids": [
                ("primate-longevity-signal", "evidence-006"),
                ("primate-longevity-signal", "evidence-007"),
                ("primate-longevity-signal", "evidence-008"),
            ],
            "coverage_status": "bounded",
            "coverage_summary": (
                "The evidence-book proves governed lecture-derived tree import, "
                "diagnostics, and tree-data correspondence behavior, but it does not yet "
                "close the broader generic tree validation surface."
            ),
            "known_gaps": [
                "No dedicated governed evidence family yet covers generic malformed-tree fixtures outside the Lund-derived studies.",
                "Topology and branch-length comparison workflows remain runtime-tested but not yet mapped into first-class evidence families.",
            ],
        },
        {
            "surface_id": "alignment-diagnostics",
            "surface_title": "Alignment diagnostics and trimming",
            "runtime_surfaces": [
                "alignment trim",
                "alignment coding",
                "alignment translate",
                "alignment identity-matrix",
            ],
            "linked_evidence_ids": [],
            "coverage_status": "uncovered",
            "coverage_summary": (
                "The runtime exposes a broad alignment surface, but the evidence-book "
                "does not yet contain governed parity or trust bundles for it."
            ),
            "known_gaps": [
                "No evidence family yet covers trimming, coding-quality, translation, or identity-matrix behavior.",
            ],
        },
        {
            "surface_id": "comparative-analysis",
            "surface_title": "Comparative models, signal testing, and regression",
            "runtime_surfaces": [
                "comparative readiness",
                "comparative signal",
                "comparative pgls",
                "ancestral compare",
            ],
            "linked_evidence_ids": [
                ("primate-longevity-signal", "evidence-001"),
                ("primate-longevity-signal", "evidence-002"),
                ("primate-longevity-signal", "evidence-003"),
                ("primate-longevity-signal", "evidence-004"),
                ("primate-longevity-signal", "evidence-005"),
                ("primate-longevity-signal", "evidence-009"),
                ("primate-pgls-and-signal", "evidence-001"),
                ("primate-pgls-and-signal", "evidence-002"),
                ("primate-pgls-and-signal", "evidence-003"),
                ("primate-pgls-and-signal", "evidence-004"),
                ("primate-pgls-and-signal", "evidence-005"),
                ("primate-pgls-and-signal", "evidence-006"),
            ],
            "coverage_status": "bounded",
            "coverage_summary": (
                "Comparative workflows are the strongest evidence-grounded surface in the "
                "repository, but transformed-tree and ancestral parity still retain explicit "
                "coverage boundaries and the broader comparative surface is not yet closed."
            ),
            "known_gaps": [
                "Transformed-tree EB parity is still explicit debt rather than closed evidence.",
                "Ancestral-mode parity from the Lund PCM2 lecture is still intentionally unclaimed.",
            ],
        },
        {
            "surface_id": "distance-analysis",
            "surface_title": "Distance analysis and distance maturity",
            "runtime_surfaces": [
                "distance matrix workflows",
                "distance maturity reports",
            ],
            "linked_evidence_ids": [],
            "coverage_status": "uncovered",
            "coverage_summary": (
                "Distance workflows have runtime and maturity tests, but they are not yet "
                "mapped into governed evidence-book families."
            ),
            "known_gaps": [
                "No checked-in evidence bundle yet proves governed distance outputs or cross-tool parity.",
            ],
        },
        {
            "surface_id": "reporting-surfaces",
            "surface_title": "Reviewer-facing reporting and evidence rendering",
            "runtime_surfaces": [
                "report tree",
                "report dataset",
                "figure packaging",
                "reviewer summaries",
            ],
            "linked_evidence_ids": [
                ("primate-longevity-signal", "evidence-001"),
            ],
            "coverage_status": "bounded",
            "coverage_summary": (
                "Reviewer summaries, catalog indexes, and governed study reports are "
                "checked in and portable, but the broader runtime report families are not "
                "yet all represented by evidence families."
            ),
            "known_gaps": [
                "Generic HTML reporting surfaces still depend more on runtime tests than on evidence-book-backed coverage.",
            ],
        },
    )


def build_analytical_surface_coverage(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    bundle_map = _bundle_entry_map(repo_root)
    surface_entries: list[dict[str, object]] = []
    status_counts: Counter[str] = Counter()
    for definition in _surface_definitions():
        linked_evidence = []
        for study_id, evidence_id in definition["linked_evidence_ids"]:
            bundle_entry = bundle_map.get((study_id, evidence_id))
            if bundle_entry is None:
                continue
            linked_evidence.append(
                {
                    "study_id": study_id,
                    "evidence_id": evidence_id,
                    "evidence_title": bundle_entry["evidence_title"],
                    "relative_path": bundle_entry["relative_path"],
                    "verdict_status": bundle_entry["verdict_status"],
                }
            )
        status_counts.update([str(definition["coverage_status"])])
        surface_entries.append(
            {
                "surface_id": definition["surface_id"],
                "surface_title": definition["surface_title"],
                "runtime_surfaces": definition["runtime_surfaces"],
                "coverage_status": definition["coverage_status"],
                "coverage_summary": definition["coverage_summary"],
                "linked_evidence_count": len(linked_evidence),
                "linked_evidence": linked_evidence,
                "known_gaps": definition["known_gaps"],
            }
        )
    return {
        "schema_version": 1,
        "surface_count": len(surface_entries),
        "coverage_status_counts": dict(sorted(status_counts.items())),
        "surfaces": surface_entries,
    }


def render_analytical_surface_coverage(payload: dict[str, object]) -> str:
    lines = [
        "# Analytical Surface Coverage",
        "",
        "This map records which major runtime analysis surfaces are already backed",
        "by governed evidence-book bundles and which remain bounded or uncovered.",
        "",
        f"- analytical surfaces: `{payload['surface_count']}`",
        "",
        "## Coverage Status Counts",
        "",
    ]
    for status, count in payload["coverage_status_counts"].items():
        lines.append(f"- `{status}`: `{count}`")
    lines.extend(["", "## Surfaces", ""])
    for surface in payload["surfaces"]:
        lines.append(f"### {surface['surface_title']}")
        lines.append("")
        lines.append(f"- surface id: `{surface['surface_id']}`")
        lines.append(f"- coverage: `{surface['coverage_status']}`")
        lines.append(f"- linked evidence bundles: `{surface['linked_evidence_count']}`")
        lines.append(f"- summary: {surface['coverage_summary']}")
        lines.append("")
        lines.append("Runtime surfaces:")
        for runtime_surface in surface["runtime_surfaces"]:
            lines.append(f"- `{runtime_surface}`")
        lines.append("")
        lines.append("Linked evidence:")
        if surface["linked_evidence"]:
            for entry in surface["linked_evidence"]:
                lines.append(
                    f"- `{entry['study_id']}/{entry['evidence_id']}` — `{entry['verdict_status']}`"
                )
        else:
            lines.append("- none yet")
        lines.append("")
        lines.append("Known gaps:")
        for gap in surface["known_gaps"]:
            lines.append(f"- {gap}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_claim_reaudit(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    coverage_payload = build_analytical_surface_coverage(repo_root)
    surface_map = {
        str(surface["surface_id"]): surface
        for surface in coverage_payload["surfaces"]
        if isinstance(surface, dict)
    }
    claims = [
        {
            "claim_id": "repository-overview-runtime-breadth",
            "relative_path": "README.md",
            "claim_text": (
                "The repository provides runtime surfaces for phylogenetic validation, "
                "alignment diagnostics, comparative analysis, and reviewable reporting."
            ),
            "linked_surface_ids": [
                "tree-diagnostics",
                "alignment-diagnostics",
                "comparative-analysis",
                "reporting-surfaces",
            ],
        },
        {
            "claim_id": "repository-overview-evidence-book",
            "relative_path": "README.md",
            "claim_text": (
                "The evidence-book is a first-class checked-in trust surface rather than "
                "a transient artifact dump."
            ),
            "linked_surface_ids": [
                "comparative-analysis",
                "reporting-surfaces",
            ],
        },
        {
            "claim_id": "docs-home-runtime-workbench",
            "relative_path": "docs/index.md",
            "claim_text": (
                "The repository is a reproducible phylogenetics workbench across trees, "
                "alignments, evidence capture, and reports."
            ),
            "linked_surface_ids": [
                "tree-diagnostics",
                "alignment-diagnostics",
                "reporting-surfaces",
            ],
        },
        {
            "claim_id": "repository-overview-comparative-scope",
            "relative_path": "docs/01-bijux-phylogenetics/index.md",
            "claim_text": (
                "Comparative analysis, ancestral-state workflows, and evidence bundles "
                "are real current repository surfaces."
            ),
            "linked_surface_ids": [
                "comparative-analysis",
                "reporting-surfaces",
            ],
        },
        {
            "claim_id": "repository-overview-alignment-scope",
            "relative_path": "docs/01-bijux-phylogenetics/index.md",
            "claim_text": (
                "Alignment trimming, translation, and identity workflows exist in the "
                "runtime surface."
            ),
            "linked_surface_ids": [
                "alignment-diagnostics",
            ],
        },
        {
            "claim_id": "repository-overview-distance-scope",
            "relative_path": "docs/01-bijux-phylogenetics/index.md",
            "claim_text": (
                "Distance workflows and their maturity surfaces exist in the runtime."
            ),
            "linked_surface_ids": [
                "distance-analysis",
            ],
        },
    ]
    status_counts: Counter[str] = Counter()
    downgraded_claim_count = 0
    audited_claims = []
    for claim in claims:
        linked_surfaces = [
            surface_map[surface_id] for surface_id in claim["linked_surface_ids"]
        ]
        statuses = {str(surface["coverage_status"]) for surface in linked_surfaces}
        if statuses == {"evidence_grounded"}:
            audit_status = "supported"
            recommendation = (
                "Current evidence coverage supports this claim without downgrade."
            )
        elif "uncovered" in statuses:
            audit_status = "downgraded_to_bounded"
            downgraded_claim_count += 1
            recommendation = (
                "Keep the runtime existence claim, but describe analytical trust as bounded "
                "until uncovered surfaces gain governed evidence families."
            )
        else:
            audit_status = "supported_with_limits"
            recommendation = (
                "Keep the claim, but pair it with an explicit evidence-book coverage link so "
                "reviewers can see the remaining boundaries."
            )
        status_counts.update([audit_status])
        audited_claims.append(
            {
                "claim_id": claim["claim_id"],
                "relative_path": claim["relative_path"],
                "claim_text": claim["claim_text"],
                "audit_status": audit_status,
                "linked_surface_ids": claim["linked_surface_ids"],
                "linked_surface_statuses": {
                    surface["surface_id"]: surface["coverage_status"]
                    for surface in linked_surfaces
                },
                "recommendation": recommendation,
            }
        )
    return {
        "schema_version": 1,
        "claim_count": len(audited_claims),
        "downgraded_claim_count": downgraded_claim_count,
        "audit_status_counts": dict(sorted(status_counts.items())),
        "claims": audited_claims,
    }


def render_claim_reaudit(payload: dict[str, object]) -> str:
    lines = [
        "# Claim Re-Audit",
        "",
        "This re-audit revisits broad repository claims against the current",
        "evidence-book reality and downgrades them wherever portable closure is still missing.",
        "",
        f"- audited claims: `{payload['claim_count']}`",
        f"- downgraded claims: `{payload['downgraded_claim_count']}`",
        "",
        "## Status Counts",
        "",
    ]
    for status, count in payload["audit_status_counts"].items():
        lines.append(f"- `{status}`: `{count}`")
    lines.extend(["", "## Claims", ""])
    for claim in payload["claims"]:
        lines.append(f"- `{claim['claim_id']}` — `{claim['audit_status']}`")
        lines.append(f"  Path: `{claim['relative_path']}`")
        lines.append(f"  Claim: {claim['claim_text']}")
        lines.append(f"  Recommendation: {claim['recommendation']}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_closure_criteria(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    coverage_payload = build_analytical_surface_coverage(repo_root)
    mismatch_archive = build_evidence_mismatch_archive(repo_root)
    coverage_gaps = build_evidence_coverage_gap_report(repo_root)
    freshness_report = build_evidence_freshness_report(repo_root)
    integrity_report = build_evidence_integrity_report(repo_root)
    false_confidence_audit = build_evidence_false_confidence_audit(repo_root)
    surface_map = {
        str(surface["surface_id"]): surface
        for surface in coverage_payload["surfaces"]
        if isinstance(surface, dict)
    }
    unresolved_mismatch_count = int(
        mismatch_archive["verdict_counts"].get("mismatch_unexplained", 0)
    )
    uncovered_surface_count = int(
        coverage_payload["coverage_status_counts"].get("uncovered", 0)
    )
    bounded_surface_count = int(
        coverage_payload["coverage_status_counts"].get("bounded", 0)
    )
    freshness_current_count = int(
        freshness_report["freshness_status_counts"].get("current", 0)
    )
    criteria = [
        {
            "criterion_id": "foundational-numerical-trust",
            "current_status": "bounded",
            "pass_when": [
                "Major runtime analytical surfaces no longer include uncovered foundational categories.",
                "Unresolved scalar mismatches are closed or explicitly relegated to bounded, non-foundational surfaces.",
                "Coverage-gap debt for core comparative and tree workflows has been reduced to bounded, explained exceptions.",
            ],
            "blocking_factors": [
                f"{uncovered_surface_count} analytical surfaces are still uncovered by governed evidence families.",
                f"{unresolved_mismatch_count} unresolved scalar mismatch rows remain visible in the mismatch archive.",
                f"{coverage_gaps['coverage_gap_count']} coverage-gap debt entries remain reviewer-visible.",
            ],
            "supporting_evidence": [
                "analytical-surface-coverage.json",
                "mismatch-archive.json",
                "coverage-gaps.json",
            ],
        },
        {
            "criterion_id": "visualization-competition",
            "current_status": "blocked_by_foundational_numerical_trust",
            "pass_when": [
                "Foundational numerical trust is already closed for the compared workflow family.",
                "Figure-comparison evidence bundles exist for the specific visual surface being claimed competitive.",
            ],
            "blocking_factors": [
                "Foundational numerical trust is still bounded rather than closed.",
                "Plot-only scientific debt remains visible for the PCM1 visual surfaces.",
            ],
            "supporting_evidence": [
                "scientific-debt-register.json",
                "claim-reaudit.json",
            ],
        },
        {
            "criterion_id": "performance-competition",
            "current_status": "blocked_by_foundational_numerical_trust",
            "pass_when": [
                "Foundational numerical trust is already closed for the workflow family.",
                "Representative benchmark evidence bundles compare Bijux against the reference toolchain under governed inputs.",
            ],
            "blocking_factors": [
                "No governed benchmark evidence families yet exist in the evidence-book.",
                "Foundational numerical trust is still bounded rather than closed.",
            ],
            "supporting_evidence": [
                "analytical-surface-coverage.json",
                "completion-gates.json",
            ],
        },
        {
            "criterion_id": "migration-usefulness",
            "current_status": "bounded",
            "pass_when": [
                "Migration studies cover more than one major analytical family for each public migration claim.",
                "Side-by-side R and Bijux mappings exist for the claimed workflow family and link to governed Evidence IDs.",
            ],
            "blocking_factors": [
                "Migration evidence is currently strongest for Lund-derived comparative workflows only.",
                "Alignment and distance migration claims would currently outrun the evidence-book.",
            ],
            "supporting_evidence": [
                "teaching-and-migration.json",
                "analytical-surface-coverage.json",
            ],
        },
        {
            "criterion_id": "reviewer-readiness",
            "current_status": "bounded",
            "pass_when": [
                "A skeptical reviewer can find linked evidence status for every major analytical surface.",
                "Coverage gaps, unresolved mismatches, and bounded claims remain explicit in one navigation path.",
                "High-level docs link directly to the evidence-book closure surfaces.",
            ],
            "blocking_factors": [
                f"{bounded_surface_count} surfaces still rely on bounded rather than evidence-grounded coverage.",
                f"{false_confidence_audit['action_required_count']} false-confidence surface findings would block reviewer readiness if nonzero.",
            ],
            "supporting_evidence": [
                "claim-reaudit.json",
                "analytical-surface-coverage.json",
                "coverage-gaps.json",
            ],
        },
    ]
    return {
        "schema_version": 1,
        "criteria_count": len(criteria),
        "freshness_current_count": freshness_current_count,
        "integrity_tracked_count": int(integrity_report["bundle_count"]),
        "surface_status_snapshot": {
            surface_id: surface["coverage_status"]
            for surface_id, surface in surface_map.items()
        },
        "criteria": criteria,
    }


def render_closure_criteria(payload: dict[str, object]) -> str:
    lines = [
        "# Closure Criteria",
        "",
        "These criteria define what the repository must prove before it can claim",
        "hard closure for foundational trust, visualization competition, performance competition,",
        "migration usefulness, or reviewer readiness.",
        "",
        f"- criteria: `{payload['criteria_count']}`",
        f"- current freshness-tracked bundles: `{payload['freshness_current_count']}`",
        f"- integrity-tracked bundles: `{payload['integrity_tracked_count']}`",
        "",
    ]
    for criterion in payload["criteria"]:
        lines.append(f"## {criterion['criterion_id']}")
        lines.append("")
        lines.append(f"- current status: `{criterion['current_status']}`")
        lines.append("")
        lines.append("Pass when:")
        for row in criterion["pass_when"]:
            lines.append(f"- {row}")
        lines.append("")
        lines.append("Blocking factors:")
        for row in criterion["blocking_factors"]:
            lines.append(f"- {row}")
        lines.append("")
        lines.append("Supporting evidence:")
        for row in criterion["supporting_evidence"]:
            lines.append(f"- `{row}`")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_maturity_scorecard(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    coverage_payload = build_analytical_surface_coverage(repo_root)
    closure_payload = build_closure_criteria(repo_root)
    portability_payload = _load_json(
        repo_root / "evidence-book" / "index" / "portability-audit.json"
    )
    freshness_payload = build_evidence_freshness_report(repo_root)
    integrity_payload = build_evidence_integrity_report(repo_root)
    mismatch_archive = build_evidence_mismatch_archive(repo_root)
    review_components = [
        {
            "component_id": "portability",
            "score": 3,
            "max_score": 3,
            "status": "strong",
            "summary": "Portable locators are enforced and the current audit has zero action-required entries.",
        },
        {
            "component_id": "parity-depth",
            "score": 2,
            "max_score": 3,
            "status": "bounded",
            "summary": "Comparative parity is deep, but foundational closure is still bounded by uncovered surfaces and one unresolved mismatch.",
        },
        {
            "component_id": "failure-transparency",
            "score": 3,
            "max_score": 3,
            "status": "strong",
            "summary": "Mismatch, verdict workflow, coverage-gap, and scientific-debt outputs keep failures visible instead of hidden.",
        },
        {
            "component_id": "coverage-breadth",
            "score": 1,
            "max_score": 3,
            "status": "weak",
            "summary": "Major runtime surfaces still include uncovered alignment and distance families.",
        },
        {
            "component_id": "reviewer-navigation",
            "score": 3,
            "max_score": 3,
            "status": "strong",
            "summary": "Reviewer summaries, docs navigation, and governed index outputs make the current trust surface inspectable.",
        },
        {
            "component_id": "regeneration-discipline",
            "score": 3,
            "max_score": 3,
            "status": "strong",
            "summary": "Workbench refresh, freshness, integrity, and partial rerun support make the evidence-book reproducible as a maintained system.",
        },
    ]
    total_score = sum(component["score"] for component in review_components)
    max_score = sum(component["max_score"] for component in review_components)
    maturity_tier = "reviewable_but_incomplete"
    return {
        "schema_version": 1,
        "maturity_tier": maturity_tier,
        "total_score": total_score,
        "max_score": max_score,
        "coverage_status_counts": coverage_payload["coverage_status_counts"],
        "open_unresolved_mismatch_count": int(
            mismatch_archive["verdict_counts"].get("mismatch_unexplained", 0)
        ),
        "freshness_status_counts": freshness_payload["freshness_status_counts"],
        "integrity_action_required_count": int(
            integrity_payload["action_required_count"]
        ),
        "portability_action_required_count": int(
            portability_payload["action_required_count"]
        ),
        "foundational_numerical_trust_status": next(
            criterion["current_status"]
            for criterion in closure_payload["criteria"]
            if criterion["criterion_id"] == "foundational-numerical-trust"
        ),
        "components": review_components,
    }


def render_evidence_maturity_scorecard(payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Maturity Scorecard",
        "",
        "This scorecard summarizes how mature the repository evidence program is",
        "across portability, parity depth, failure transparency, coverage, reviewability,",
        "and regeneration discipline.",
        "",
        f"- maturity tier: `{payload['maturity_tier']}`",
        f"- score: `{payload['total_score']}/{payload['max_score']}`",
        f"- foundational numerical trust: `{payload['foundational_numerical_trust_status']}`",
        "",
        "## Components",
        "",
    ]
    for component in payload["components"]:
        lines.append(
            f"- `{component['component_id']}` — `{component['status']}` `{component['score']}/{component['max_score']}`"
        )
        lines.append(f"  {component['summary']}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_review_ritual(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    scorecard = build_evidence_maturity_scorecard(repo_root)
    ritual_steps = [
        {
            "ritual_id": "weekly-evidence-triage",
            "cadence": "weekly",
            "goal": "Keep freshness, integrity, and open gaps visible before they drift into release-time surprises.",
            "commands": [
                "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c \"from pathlib import Path; from bijux_phylogenetics.evidence.workbench import refresh_evidence_book; refresh_evidence_book(Path('.'))\"",
                "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c \"from pathlib import Path; from bijux_phylogenetics.evidence.book import validate_evidence_book; report = validate_evidence_book(Path('.')); print({'valid': report.valid, 'issue_count': len(report.issues)})\"",
            ],
            "required_inputs": [
                "freshness-report.json",
                "integrity-report.json",
                "coverage-gaps.json",
                "mismatch-archive.json",
            ],
        },
        {
            "ritual_id": "pre-release-closure-review",
            "cadence": "before release tags",
            "goal": "Re-check broad confidence claims against the current closure state before public release messaging goes out.",
            "commands": [
                "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 pytest packages/bijux-phylogenetics/tests/test_evidence_book_contract.py packages/bijux-phylogenetics/tests/test_evidence_book_repository.py packages/bijux-phylogenetics/tests/test_evidence_book_outputs.py packages/bijux-phylogenetics/tests/test_evidence_workbench.py",
                "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 python -c \"from pathlib import Path; from bijux_phylogenetics.evidence.closure import build_claim_reaudit, build_evidence_maturity_scorecard; print(build_claim_reaudit(Path('.'))['downgraded_claim_count']); print(build_evidence_maturity_scorecard(Path('.'))['maturity_tier'])\"",
            ],
            "required_inputs": [
                "claim-reaudit.json",
                "closure-criteria.json",
                "evidence-maturity-scorecard.json",
                "completion-gates.json",
            ],
        },
        {
            "ritual_id": "new-study-intake-review",
            "cadence": "for every new evidence study family",
            "goal": "Force new study intake to declare evidence ownership, closure effects, and unresolved boundaries up front.",
            "commands": [
                "bijux-phylogenetics evidence book studies --json",
                "bijux-phylogenetics evidence book build --json",
                "bijux-phylogenetics evidence book validate --json",
            ],
            "required_inputs": [
                "analytical-surface-coverage.json",
                "completion-gates.json",
                "scientific-debt-register.json",
            ],
        },
    ]
    return {
        "schema_version": 1,
        "ritual_count": len(ritual_steps),
        "maturity_tier_snapshot": scorecard["maturity_tier"],
        "rituals": ritual_steps,
    }


def render_evidence_review_ritual(payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Review Ritual",
        "",
        "This ritual keeps parity gaps, stale bundles, and broad confidence claims",
        "under recurring maintainers review instead of treating the evidence-book as one-off work.",
        "",
        f"- rituals: `{payload['ritual_count']}`",
        f"- current maturity tier snapshot: `{payload['maturity_tier_snapshot']}`",
        "",
    ]
    for ritual in payload["rituals"]:
        lines.append(f"## {ritual['ritual_id']}")
        lines.append("")
        lines.append(f"- cadence: `{ritual['cadence']}`")
        lines.append(f"- goal: {ritual['goal']}")
        lines.append("")
        lines.append("Required inputs:")
        for entry in ritual["required_inputs"]:
            lines.append(f"- `{entry}`")
        lines.append("")
        lines.append("Commands:")
        for entry in ritual["commands"]:
            lines.append(f"- `{entry}`")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_completion_gates(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    coverage_payload = build_analytical_surface_coverage(repo_root)
    scorecard = build_evidence_maturity_scorecard(repo_root)
    gate_entries = []
    status_counts: Counter[str] = Counter()
    for surface in coverage_payload["surfaces"]:
        coverage_status = str(surface["coverage_status"])
        if coverage_status == "uncovered":
            completion_state = "not_ready"
            missing_conditions = [
                "No governed evidence family is linked to this analytical surface yet.",
                "Reviewer-readable parity or trust bundles are still absent.",
            ]
        else:
            completion_state = "bounded"
            missing_conditions = [
                "Broader surface coverage still depends on bounded or partial evidence families.",
                "The repository scorecard still classifies the overall program as reviewable but incomplete.",
            ]
        status_counts.update([completion_state])
        gate_entries.append(
            {
                "surface_id": surface["surface_id"],
                "surface_title": surface["surface_title"],
                "completion_state": completion_state,
                "linked_evidence_count": surface["linked_evidence_count"],
                "minimum_requirements": [
                    "indexed evidence bundles",
                    "reviewer-readable summaries",
                    "portable checked-in artifacts",
                    "source-grounded provenance",
                ],
                "missing_conditions": missing_conditions,
            }
        )
    return {
        "schema_version": 1,
        "gate_count": len(gate_entries),
        "scorecard_maturity_tier": scorecard["maturity_tier"],
        "completion_state_counts": dict(sorted(status_counts.items())),
        "gates": gate_entries,
    }


def render_completion_gates(payload: dict[str, object]) -> str:
    lines = [
        "# Completion Gates",
        "",
        "These gates prevent the repository from calling a major analytical surface",
        "complete before the evidence-book actually supports that claim.",
        "",
        f"- gates: `{payload['gate_count']}`",
        f"- maturity tier: `{payload['scorecard_maturity_tier']}`",
        "",
        "## Completion State Counts",
        "",
    ]
    for state, count in payload["completion_state_counts"].items():
        lines.append(f"- `{state}`: `{count}`")
    lines.extend(["", "## Gates", ""])
    for gate in payload["gates"]:
        lines.append(f"### {gate['surface_title']}")
        lines.append("")
        lines.append(f"- surface id: `{gate['surface_id']}`")
        lines.append(f"- completion state: `{gate['completion_state']}`")
        lines.append(f"- linked evidence bundles: `{gate['linked_evidence_count']}`")
        lines.append("")
        lines.append("Minimum requirements:")
        for requirement in gate["minimum_requirements"]:
            lines.append(f"- {requirement}")
        lines.append("")
        lines.append("Missing conditions:")
        for condition in gate["missing_conditions"]:
            lines.append(f"- {condition}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"
