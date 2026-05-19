from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Any

from ..portability import (
    audit_payload_path_values,
    classify_locator_kind,
    render_portability_rules_markdown,
)
from ..study_contracts import load_study_contract
from ..teaching import (
    TEACHING_AND_MIGRATION_INDEX_FILENAME,
    TEACHING_AND_MIGRATION_SUMMARY_FILENAME,
    build_teaching_and_migration_index,
    render_teaching_and_migration_index_markdown,
)
from .layout import (
    EVIDENCE_BOOK_DIRNAME,
    EVIDENCE_BUNDLE_MANIFEST,
    EVIDENCE_CATALOG,
    EVIDENCE_CLAIM_MAP,
    EVIDENCE_FALSE_CONFIDENCE_AUDIT,
    EVIDENCE_FALSE_CONFIDENCE_SUMMARY,
    EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
    EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
    EVIDENCE_ID_PATTERN,
    EVIDENCE_INDEX,
    EVIDENCE_INDEX_DIRNAME,
    EVIDENCE_MISMATCH_ARCHIVE,
    EVIDENCE_MISMATCH_SUMMARY,
    EVIDENCE_PARITY_DASHBOARD,
    EVIDENCE_PARITY_SUMMARY,
    EVIDENCE_PORTABILITY_AUDIT,
    EVIDENCE_PORTABILITY_SUMMARY,
    EVIDENCE_REGENERATION_CONTRACT,
    EVIDENCE_REGENERATION_SUMMARY,
    EVIDENCE_SCIENTIFIC_DEBT_REGISTER,
    EVIDENCE_SCIENTIFIC_DEBT_SUMMARY,
    EVIDENCE_VERDICT_WORKFLOWS,
    EVIDENCE_VERDICT_WORKFLOWS_SUMMARY,
    evidence_book_root,
    load_bundle_claim_rows as _load_bundle_claim_rows,
    load_json as _load_json,
    relative_to as _relative_to,
    study_paths as _study_paths,
)
from .models import (
    EvidenceBookValidationIssue as EvidenceBookValidationIssue,
    EvidenceBookValidationReport as EvidenceBookValidationReport,
)
from .validation import validate_evidence_book


def build_evidence_book_index(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    report = validate_evidence_book(
        repo_root,
        require_index_outputs=False,
        require_generated_bundle_outputs=False,
    )
    if not report.valid:
        messages = "; ".join(
            f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
        )
        raise ValueError(f"evidence-book is invalid: {messages}")

    evidence_entries: list[dict[str, object]] = []
    studies: list[dict[str, object]] = []
    for study_root in _study_paths(root):
        study_manifest = load_study_contract(study_root)
        study_entries: list[dict[str, object]] = []
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = _load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            verdict = manifest["verdict"]
            entry = {
                "study_id": study_manifest["study_id"],
                "study_title": study_manifest["study_title"],
                "evidence_id": manifest["evidence_id"],
                "evidence_title": manifest["evidence_title"],
                "owner_package": manifest["owner_package"],
                "relative_path": _relative_to(root, bundle_root).as_posix(),
                "claim_tags": manifest["claim_tags"],
                "comparison_mode": manifest["comparison_mode"],
                "verdict_status": verdict["status"],
                "verdict_summary": verdict["summary"],
            }
            study_entries.append(entry)
            evidence_entries.append(entry)
        studies.append(
            {
                "study_id": study_manifest["study_id"],
                "study_title": study_manifest["study_title"],
                "summary": study_manifest["summary"],
                "owner_package": study_manifest["owner_package"],
                "study_categories": study_manifest["study_categories"],
                "bundle_count": len(study_entries),
                "evidence": study_entries,
            }
        )

    verdict_counts = Counter(str(entry["verdict_status"]) for entry in evidence_entries)
    return {
        "schema_version": 1,
        "root": EVIDENCE_BOOK_DIRNAME,
        "study_count": len(studies),
        "evidence_count": len(evidence_entries),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "studies": studies,
        "evidence": evidence_entries,
    }


def build_evidence_claim_map(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    report = validate_evidence_book(
        repo_root,
        require_index_outputs=False,
        require_generated_bundle_outputs=False,
    )
    if not report.valid:
        messages = "; ".join(
            f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
        )
        raise ValueError(f"evidence-book is invalid: {messages}")

    claims_by_id: dict[str, dict[str, object]] = {}
    for study_root in _study_paths(root):
        study_manifest = load_study_contract(study_root)
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = _load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            for row in _load_bundle_claim_rows(bundle_root):
                claim_id = row["claim_id"]
                entry = claims_by_id.setdefault(
                    claim_id,
                    {
                        "claim_id": claim_id,
                        "claim_title": row.get("claim_title", claim_id),
                        "evidence": [],
                    },
                )
                entry["evidence"].append(
                    {
                        "study_id": study_manifest["study_id"],
                        "study_title": study_manifest["study_title"],
                        "evidence_id": manifest["evidence_id"],
                        "relative_path": _relative_to(root, bundle_root).as_posix(),
                        "owner_package": manifest["owner_package"],
                        "bundle_verdict_status": manifest["verdict"]["status"],
                        "claim_verdict": row.get("verdict"),
                        "source_fragments": row.get("source_fragments", []),
                    }
                )

    claims = [claims_by_id[key] for key in sorted(claims_by_id)]
    return {
        "schema_version": 1,
        "claim_count": len(claims),
        "claims": claims,
    }


def build_evidence_parity_dashboard(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    report = validate_evidence_book(
        repo_root,
        require_index_outputs=False,
        require_generated_bundle_outputs=False,
    )
    if not report.valid:
        messages = "; ".join(
            f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
        )
        raise ValueError(f"evidence-book is invalid: {messages}")

    studies: list[dict[str, object]] = []
    total_row_count = 0
    verdict_counts: Counter[str] = Counter()
    comparison_kind_counts: Counter[str] = Counter()

    for study_root in _study_paths(root):
        study_manifest = load_study_contract(study_root)
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            scalar_table_path = bundle_root / "results" / "scalar-parity-table.json"
            if not scalar_table_path.exists():
                continue
            manifest = _load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            scalar_table = _load_json(scalar_table_path)
            total_row_count += int(scalar_table.get("row_count", 0))
            verdict_counts.update(
                {
                    verdict: int(count)
                    for verdict, count in scalar_table.get("verdict_counts", {}).items()
                }
            )
            comparison_kind_counts.update(
                str(row["comparison_kind"])
                for row in scalar_table.get("rows", [])
                if isinstance(row, dict) and isinstance(row.get("comparison_kind"), str)
            )
            studies.append(
                {
                    "study_id": study_manifest["study_id"],
                    "study_title": study_manifest["study_title"],
                    "evidence_id": manifest["evidence_id"],
                    "relative_path": _relative_to(root, bundle_root).as_posix(),
                    "bundle_verdict_status": manifest["verdict"]["status"],
                    "scalar_row_count": scalar_table["row_count"],
                    "scalar_verdict_counts": scalar_table["verdict_counts"],
                    "comparison_kind_counts": Counter(
                        str(row["comparison_kind"])
                        for row in scalar_table.get("rows", [])
                        if isinstance(row, dict)
                        and isinstance(row.get("comparison_kind"), str)
                    ),
                    "parity_expectation_counts": {},
                    "family_verdicts": [],
                }
            )

    normalized_studies = []
    for entry in studies:
        expectation_counter = entry["parity_expectation_counts"]
        normalized_studies.append(
            {
                **entry,
                "comparison_kind_counts": dict(
                    sorted(entry["comparison_kind_counts"].items())
                ),
                "parity_expectation_counts": dict(sorted(expectation_counter.items())),
            }
        )

    return {
        "schema_version": 1,
        "study_count": len(normalized_studies),
        "scalar_row_count": total_row_count,
        "scalar_verdict_counts": dict(sorted(verdict_counts.items())),
        "comparison_kind_counts": dict(sorted(comparison_kind_counts.items())),
        "parity_expectation_counts": {},
        "studies": normalized_studies,
    }


def _iter_scalar_parity_tables(
    repo_root: Path,
) -> list[tuple[dict[str, object], dict[str, object], Path, dict[str, object]]]:
    root = evidence_book_root(repo_root)
    report = validate_evidence_book(
        repo_root,
        require_index_outputs=False,
        require_generated_bundle_outputs=False,
    )
    if not report.valid:
        messages = "; ".join(
            f"{issue.path.as_posix()}: {issue.message}" for issue in report.issues
        )
        raise ValueError(f"evidence-book is invalid: {messages}")

    tables: list[
        tuple[dict[str, object], dict[str, object], Path, dict[str, object]]
    ] = []
    for study_root in _study_paths(root):
        study_manifest = load_study_contract(study_root)
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            scalar_table_path = bundle_root / "results" / "scalar-parity-table.json"
            if not scalar_table_path.exists():
                continue
            tables.append(
                (
                    study_manifest,
                    _load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST),
                    bundle_root,
                    _load_json(scalar_table_path),
                )
            )
    return tables


def build_evidence_mismatch_archive(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    mismatches: list[dict[str, object]] = []
    verdict_counts: Counter[str] = Counter()
    for (
        study_manifest,
        manifest,
        bundle_root,
        scalar_table,
    ) in _iter_scalar_parity_tables(repo_root):
        for row in scalar_table.get("rows", []):
            if not isinstance(row, dict):
                continue
            verdict = row.get("verdict")
            if verdict not in {"mismatch_explained", "mismatch_unexplained"}:
                continue
            verdict_counts.update([str(verdict)])
            mismatches.append(
                {
                    "archive_id": (
                        f"{study_manifest['study_id']}-{manifest['evidence_id']}-{row['row_id']}"
                    ),
                    "study_id": study_manifest["study_id"],
                    "study_title": study_manifest["study_title"],
                    "evidence_id": manifest["evidence_id"],
                    "relative_path": _relative_to(root, bundle_root).as_posix(),
                    "row_id": row["row_id"],
                    "method_family": row["method_family"],
                    "metric_name": row["metric_name"],
                    "verdict": verdict,
                    "resolution_state": (
                        "explained" if verdict == "mismatch_explained" else "open"
                    ),
                    "r_value": row.get("r_value"),
                    "bijux_value": row.get("bijux_value"),
                    "observed_abs_diff": row.get("observed_abs_diff"),
                    "tolerance_abs_diff": row.get("tolerance_abs_diff"),
                    "explanation_kind": row.get("explanation_kind"),
                    "verdict_explanation": row.get("verdict_explanation"),
                }
            )
    mismatches.sort(key=lambda entry: str(entry["archive_id"]))
    return {
        "schema_version": 1,
        "mismatch_count": len(mismatches),
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "mismatches": mismatches,
    }


def render_evidence_mismatch_archive(payload: dict[str, object]) -> str:
    lines = [
        "# Mismatch Archive",
        "",
        "This archive tracks every scalar parity row that is still represented as",
        "`mismatch_explained` or `mismatch_unexplained` inside the evidence-book.",
        "",
        f"- mismatch rows: `{payload['mismatch_count']}`",
        "",
    ]
    verdict_counts = payload.get("verdict_counts", {})
    if isinstance(verdict_counts, dict) and verdict_counts:
        lines.append("## Verdict Counts")
        lines.append("")
        for verdict, count in verdict_counts.items():
            lines.append(f"- `{verdict}`: `{count}`")
        lines.append("")

    lines.append("## Entries")
    lines.append("")
    for entry in payload["mismatches"]:
        lines.append(
            f"- `{entry['archive_id']}` — `{entry['verdict']}` in `{entry['relative_path']}`"
        )
        lines.append(
            f"  Metric: `{entry['metric_name']}` (`{entry['method_family']}`), diff=`{entry['observed_abs_diff']}`"
        )
        if entry.get("verdict_explanation"):
            lines.append(f"  Explanation: {entry['verdict_explanation']}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_verdict_workflows(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    mismatch_archive = build_evidence_mismatch_archive(repo_root)
    explained_entries = [
        entry
        for entry in mismatch_archive["mismatches"]
        if entry["verdict"] == "mismatch_explained"
    ]
    unexplained_entries = [
        entry
        for entry in mismatch_archive["mismatches"]
        if entry["verdict"] == "mismatch_unexplained"
    ]

    not_comparable_entries: list[dict[str, object]] = []
    for study_root in _study_paths(root):
        study_manifest = load_study_contract(study_root)
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = _load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            for row in _load_bundle_claim_rows(bundle_root):
                if row.get("verdict") != "not_comparable":
                    continue
                not_comparable_entries.append(
                    {
                        "entry_id": (
                            f"{study_manifest['study_id']}-{manifest['evidence_id']}-{row['claim_id']}"
                        ),
                        "study_id": study_manifest["study_id"],
                        "study_title": study_manifest["study_title"],
                        "evidence_id": manifest["evidence_id"],
                        "relative_path": _relative_to(root, bundle_root).as_posix(),
                        "claim_id": row["claim_id"],
                        "claim_title": row.get("claim_title", row["claim_id"]),
                        "summary": row.get("summary"),
                        "source_fragments": row.get("source_fragments", []),
                    }
                )
    not_comparable_entries.sort(key=lambda entry: str(entry["entry_id"]))

    workflows = [
        {
            "verdict_status": "mismatch_explained",
            "workflow_rule": "Keep the row visible, record the explanation kind explicitly, and do not promote it to a full match unless the reference source becomes more precise.",
            "entry_count": len(explained_entries),
            "entries": explained_entries,
        },
        {
            "verdict_status": "mismatch_unexplained",
            "workflow_rule": "Keep the row visible, treat it as open scientific debt, and require explicit closure rather than silent tolerance inflation.",
            "entry_count": len(unexplained_entries),
            "entries": unexplained_entries,
        },
        {
            "verdict_status": "not_comparable",
            "workflow_rule": "Keep the boundary explicit, point to the governing claim and evidence bundle, and do not restate it as a pass/fail parity surface until the runtime owns the comparison.",
            "entry_count": len(not_comparable_entries),
            "entries": not_comparable_entries,
        },
    ]
    return {
        "schema_version": 1,
        "workflow_count": len(workflows),
        "workflows": workflows,
    }


def render_evidence_verdict_workflows(payload: dict[str, object]) -> str:
    lines = [
        "# Verdict Workflows",
        "",
        "This index explains how `mismatch_explained`, `mismatch_unexplained`, and",
        "`not_comparable` evidence states are supposed to stay visible.",
        "",
    ]
    for workflow in payload["workflows"]:
        lines.append(f"## {workflow['verdict_status']}")
        lines.append("")
        lines.append(workflow["workflow_rule"])
        lines.append("")
        lines.append(f"- entries: `{workflow['entry_count']}`")
        lines.append("")
        for entry in workflow["entries"]:
            identifier = entry.get("archive_id", entry.get("entry_id"))
            lines.append(f"- `{identifier}`")
            if entry.get("relative_path"):
                lines.append(f"  Path: `{entry['relative_path']}`")
            if entry.get("verdict_explanation"):
                lines.append(f"  Explanation: {entry['verdict_explanation']}")
            if entry.get("claim_title"):
                lines.append(f"  Claim: {entry['claim_title']}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_false_confidence_audit(repo_root: Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    audited_surfaces = [
        {
            "surface_id": "repository-readme",
            "relative_path": "README.md",
            "disallowed_phrases": [
                "comparative validation studies",
            ],
        },
        {
            "surface_id": "evidence-book-readme",
            "relative_path": "evidence-book/README.md",
            "disallowed_phrases": [
                "comparative validation against",
            ],
        },
        {
            "surface_id": "primate-longevity-readme",
            "relative_path": "evidence-book/studies/primate-longevity-signal/README.md",
            "disallowed_phrases": [
                "validated against established R",
            ],
        },
        {
            "surface_id": "primate-pgls-readme",
            "relative_path": "evidence-book/studies/primate-pgls-and-signal/README.md",
            "disallowed_phrases": [
                "fully validated",
            ],
        },
    ]
    entries = []
    action_required_count = 0
    for surface in audited_surfaces:
        surface_path = repo_root / surface["relative_path"]
        if not surface_path.exists():
            entries.append(
                {
                    "surface_id": surface["surface_id"],
                    "relative_path": surface["relative_path"],
                    "status": "not_present",
                    "matched_phrases": [],
                    "review_rule": "Surface is outside the current fixture repository shape and is therefore not audited here.",
                }
            )
            continue
        text = surface_path.read_text(encoding="utf-8")
        matched_phrases = [
            phrase for phrase in surface["disallowed_phrases"] if phrase in text
        ]
        status = "action_required" if matched_phrases else "clear"
        if status == "action_required":
            action_required_count += 1
        entries.append(
            {
                "surface_id": surface["surface_id"],
                "relative_path": surface["relative_path"],
                "status": status,
                "matched_phrases": matched_phrases,
                "review_rule": "High-level public evidence surfaces must describe parity and trust boundaries honestly instead of implying broader validation than the evidence-book currently proves.",
            }
        )
    return {
        "schema_version": 1,
        "surface_count": len(entries),
        "action_required_count": action_required_count,
        "surfaces": entries,
    }


def render_evidence_false_confidence_audit(payload: dict[str, object]) -> str:
    lines = [
        "# False Confidence Audit",
        "",
        "This audit checks high-level public evidence surfaces for phrases that would",
        "overstate confidence beyond the current governed parity state.",
        "",
        f"- audited surfaces: `{payload['surface_count']}`",
        f"- action required: `{payload['action_required_count']}`",
        "",
        "## Surfaces",
        "",
    ]
    for surface in payload["surfaces"]:
        lines.append(
            f"- `{surface['surface_id']}` — `{surface['status']}` at `{surface['relative_path']}`"
        )
        if surface["matched_phrases"]:
            lines.append(
                "  Matched phrases: "
                + ", ".join(f"`{phrase}`" for phrase in surface["matched_phrases"])
            )
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_portability_audit(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    excluded_index_files = {
        EVIDENCE_PORTABILITY_AUDIT,
        EVIDENCE_PORTABILITY_SUMMARY,
        EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
        EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
        EVIDENCE_REGENERATION_CONTRACT,
        EVIDENCE_REGENERATION_SUMMARY,
    }
    json_file_count = 0
    locator_kind_counts: Counter[str] = Counter()
    issues: list[dict[str, object]] = []
    report_like_file_count = 0
    for path in sorted(
        candidate for candidate in root.rglob("*") if candidate.is_file()
    ):
        if (
            path.parent.name == EVIDENCE_INDEX_DIRNAME
            and path.name in excluded_index_files
        ):
            continue
        relative_path = _relative_to(root, path).as_posix()
        if path.suffix in {".json", ".md", ".csv", ".tsv", ".nwk"}:
            report_like_file_count += 1
        if path.suffix != ".json":
            continue
        json_file_count += 1
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        if isinstance(payload, (dict, list)):
            for issue in audit_payload_path_values(
                payload,
                relative_file_path=relative_path,
            ):
                issues.append(
                    {
                        "relative_file_path": issue.relative_file_path,
                        "json_pointer": issue.json_pointer,
                        "value": issue.value,
                        "issue_kind": issue.issue_kind,
                        "message": issue.message,
                    }
                )
            stack: list[dict[str, Any] | list[Any]] = [payload]
            while stack:
                current = stack.pop()
                if isinstance(current, dict):
                    for value in current.values():
                        if isinstance(value, str):
                            locator_kind = classify_locator_kind(value)
                            if locator_kind != "non_path_text":
                                locator_kind_counts.update([locator_kind])
                        elif isinstance(value, (dict, list)):
                            stack.append(value)
                elif isinstance(current, list):
                    for value in current:
                        if isinstance(value, str):
                            locator_kind = classify_locator_kind(value)
                            if locator_kind != "non_path_text":
                                locator_kind_counts.update([locator_kind])
                        elif isinstance(value, (dict, list)):
                            stack.append(value)
    return {
        "schema_version": 1,
        "rule_count": 3,
        "rules_markdown": render_portability_rules_markdown(),
        "audited_json_file_count": json_file_count,
        "report_like_file_count": report_like_file_count,
        "action_required_count": len(issues),
        "locator_kind_counts": dict(sorted(locator_kind_counts.items())),
        "issues": issues,
    }


def render_evidence_portability_audit(payload: dict[str, object]) -> str:
    lines = [
        "# Portability Audit",
        "",
        "This audit enforces portable path semantics for checked-in evidence payloads.",
        "",
        f"- audited json files: `{payload['audited_json_file_count']}`",
        f"- report-like files: `{payload['report_like_file_count']}`",
        f"- action required: `{payload['action_required_count']}`",
        "",
        "## Rules",
        "",
        *render_portability_rules_markdown().splitlines()[2:],
        "",
    ]
    locator_kind_counts = payload.get("locator_kind_counts", {})
    if isinstance(locator_kind_counts, dict) and locator_kind_counts:
        lines.append("## Locator Kinds")
        lines.append("")
        for locator_kind, count in locator_kind_counts.items():
            lines.append(f"- `{locator_kind}`: `{count}`")
        lines.append("")
    lines.append("## Issues")
    lines.append("")
    if not payload["issues"]:
        lines.append("- none")
    else:
        for issue in payload["issues"]:
            lines.append(
                f"- `{issue['relative_file_path']}` `{issue['json_pointer']}` `{issue['issue_kind']}`"
            )
            lines.append(f"  Value: `{issue['value']}`")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_fragile_example_audit(repo_root: Path) -> dict[str, object]:
    scientific_debt = build_evidence_scientific_debt_register(repo_root)
    fragile_kinds = {
        "artifact_only",
        "coverage_gap",
        "instability",
        "model-boundary",
        "plot_only",
        "seeded_input_only",
        "workflow_only",
    }
    fragile_entries = [
        debt for debt in scientific_debt["debts"] if debt["debt_kind"] in fragile_kinds
    ]
    counts = Counter(str(entry["debt_kind"]) for entry in fragile_entries)
    return {
        "schema_version": 1,
        "fragile_example_count": len(fragile_entries),
        "fragility_kind_counts": dict(sorted(counts.items())),
        "examples": fragile_entries,
    }


def render_evidence_fragile_example_audit(payload: dict[str, object]) -> str:
    lines = [
        "# Fragile Example Audit",
        "",
        "This audit lists evidence bundles or bundle fragments that still depend on",
        "narrow assumptions, seeded inputs, plot-only interpretations, or explicit",
        "coverage boundaries.",
        "",
        f"- fragile examples: `{payload['fragile_example_count']}`",
        "",
    ]
    counts = payload.get("fragility_kind_counts", {})
    if isinstance(counts, dict) and counts:
        lines.append("## Fragility Kinds")
        lines.append("")
        for kind, count in counts.items():
            lines.append(f"- `{kind}`: `{count}`")
        lines.append("")
    lines.append("## Examples")
    lines.append("")
    for example in payload["examples"]:
        lines.append(f"- `{example['debt_id']}` — `{example['debt_kind']}`")
        if example.get("relative_path"):
            lines.append(f"  Path: `{example['relative_path']}`")
        if example.get("detail"):
            lines.append(f"  Detail: {example['detail']}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_regeneration_contract(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    studies: list[dict[str, object]] = []
    for study_root in _study_paths(root):
        study_manifest = load_study_contract(study_root)
        source_paths: list[str] = []
        generated_paths: list[str] = []
        for path in sorted(
            candidate for candidate in study_root.rglob("*") if candidate.is_file()
        ):
            relative_path = _relative_to(Path(repo_root), path).as_posix()
            if any(
                parent.name in {"reference", "provenance", "datasets"}
                for parent in path.parents
            ):
                source_paths.append(relative_path)
            else:
                generated_paths.append(relative_path)
        studies.append(
            {
                "study_id": study_manifest["study_id"],
                "study_title": study_manifest["study_title"],
                "build_script_path": None,
                "rerun_command": (
                    "UV_PROJECT_ENVIRONMENT=artifacts/root/venv uv run --python 3.11 "
                    "python -m bijux_phylogenetics.command_line evidence book build "
                    f"--study-id {study_manifest['study_id']}"
                ),
                "bundle_ids": [
                    path.name
                    for path in sorted(
                        candidate
                        for candidate in study_root.iterdir()
                        if candidate.is_dir()
                        and EVIDENCE_ID_PATTERN.fullmatch(candidate.name)
                    )
                ],
                "source_paths": source_paths,
                "generated_paths": generated_paths,
            }
        )
    return {
        "schema_version": 1,
        "study_count": len(studies),
        "studies": studies,
    }


def render_evidence_regeneration_contract(payload: dict[str, object]) -> str:
    lines = [
        "# Regeneration Contract",
        "",
        "This contract records which study files are governed sources, which are",
        "generated durable outputs, and how each study is rerun.",
        "",
        f"- studies: `{payload['study_count']}`",
        "",
    ]
    for study in payload["studies"]:
        lines.append(f"## {study['study_title']}")
        lines.append("")
        lines.append(f"- study id: `{study['study_id']}`")
        lines.append(f"- build script: `{study['build_script_path']}`")
        if study["rerun_command"] is not None:
            lines.append(f"- rerun command: `{study['rerun_command']}`")
        lines.append(f"- bundle ids: `{', '.join(study['bundle_ids'])}`")
        lines.append(f"- source path count: `{len(study['source_paths'])}`")
        lines.append(f"- generated path count: `{len(study['generated_paths'])}`")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_evidence_scientific_debt_register(repo_root: Path) -> dict[str, object]:
    root = evidence_book_root(repo_root)
    mismatch_archive = build_evidence_mismatch_archive(repo_root)
    verdict_workflows = build_evidence_verdict_workflows(repo_root)
    debts: list[dict[str, object]] = []

    for mismatch in mismatch_archive["mismatches"]:
        if mismatch["verdict"] != "mismatch_unexplained":
            continue
        debts.append(
            {
                "debt_id": mismatch["archive_id"],
                "debt_kind": "unresolved_mismatch",
                "study_id": mismatch["study_id"],
                "evidence_id": mismatch["evidence_id"],
                "relative_path": mismatch["relative_path"],
                "detail": (
                    f"{mismatch['metric_name']} remains unresolved with observed absolute "
                    f"difference {mismatch['observed_abs_diff']}."
                ),
                "evidence": [
                    f"verdict={mismatch['verdict']}",
                    f"r_value={mismatch['r_value']}",
                    f"bijux_value={mismatch['bijux_value']}",
                ],
            }
        )

    not_comparable_workflow = next(
        workflow
        for workflow in verdict_workflows["workflows"]
        if workflow["verdict_status"] == "not_comparable"
    )
    for entry in not_comparable_workflow["entries"]:
        debts.append(
            {
                "debt_id": entry["entry_id"],
                "debt_kind": "coverage_gap",
                "study_id": entry["study_id"],
                "evidence_id": entry["evidence_id"],
                "relative_path": entry["relative_path"],
                "detail": entry.get("summary")
                or "This analytical surface is still tracked as not comparable.",
                "evidence": [
                    f"claim_id={entry['claim_id']}",
                    *[
                        f"source_fragment={fragment}"
                        for fragment in entry.get("source_fragments", [])
                    ],
                ],
            }
        )

    for study_root in _study_paths(root):
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            debt_register_path = bundle_root / "scientific_debt_register.json"
            if debt_register_path.exists():
                payload = _load_json(debt_register_path)
                for index, debt in enumerate(payload.get("debts", []), start=1):
                    if not isinstance(debt, dict):
                        continue
                    identifier = (
                        debt.get("debt_id")
                        or debt.get("block_id")
                        or f"{study_root.name}-{bundle_root.name}-bundle-debt-{index}"
                    )
                    debts.append(
                        {
                            "debt_id": identifier,
                            "debt_kind": debt.get("debt_kind", "bundle-debt"),
                            "study_id": study_root.name,
                            "evidence_id": bundle_root.name,
                            "relative_path": _relative_to(
                                root, debt_register_path
                            ).as_posix(),
                            "detail": debt.get("detail"),
                            "evidence": debt.get("evidence", []),
                        }
                    )
            for json_path in sorted(bundle_root.glob("*.json")):
                if json_path.name in {
                    "manifest.json",
                    "claims.json",
                    "scientific_debt_register.json",
                }:
                    continue
                try:
                    payload = _load_json(json_path)
                except (ValueError, json.JSONDecodeError):
                    continue
                for debt in payload.get("scientific_debt_entries", []):
                    if not isinstance(debt, dict):
                        continue
                    debts.append(
                        {
                            "debt_id": debt.get("debt_id")
                            or f"{study_root.name}-{bundle_root.name}-{json_path.stem}",
                            "debt_kind": debt.get("debt_kind", "bundle-inline-debt"),
                            "study_id": study_root.name,
                            "evidence_id": bundle_root.name,
                            "relative_path": _relative_to(root, json_path).as_posix(),
                            "detail": debt.get("detail"),
                            "evidence": debt.get("evidence", []),
                        }
                    )

    debts.sort(key=lambda entry: str(entry["debt_id"]))
    debt_kind_counts = Counter(str(entry["debt_kind"]) for entry in debts)
    return {
        "schema_version": 1,
        "debt_count": len(debts),
        "debt_kind_counts": dict(sorted(debt_kind_counts.items())),
        "debts": debts,
    }


def render_evidence_scientific_debt_register(payload: dict[str, object]) -> str:
    lines = [
        "# Scientific Debt Register",
        "",
        "This register centralizes unresolved parity gaps, explicit coverage",
        "boundaries, and reviewer-visible trust weaknesses across the evidence-book.",
        "",
        f"- debt entries: `{payload['debt_count']}`",
        "",
    ]
    debt_kind_counts = payload.get("debt_kind_counts", {})
    if isinstance(debt_kind_counts, dict) and debt_kind_counts:
        lines.append("## Debt Kinds")
        lines.append("")
        for debt_kind, count in debt_kind_counts.items():
            lines.append(f"- `{debt_kind}`: `{count}`")
        lines.append("")

    lines.append("## Entries")
    lines.append("")
    for debt in payload["debts"]:
        lines.append(f"- `{debt['debt_id']}` — `{debt['debt_kind']}`")
        if debt.get("relative_path"):
            lines.append(f"  Path: `{debt['relative_path']}`")
        if debt.get("detail"):
            lines.append(f"  Detail: {debt['detail']}")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_evidence_parity_dashboard(dashboard_payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Parity Dashboard",
        "",
        f"- studies with scalar parity tables: `{dashboard_payload['study_count']}`",
        f"- scalar parity rows: `{dashboard_payload['scalar_row_count']}`",
        "",
    ]
    verdict_counts = dashboard_payload.get("scalar_verdict_counts", {})
    if isinstance(verdict_counts, dict) and verdict_counts:
        lines.append("## Scalar Verdict Counts")
        lines.append("")
        for verdict, count in verdict_counts.items():
            lines.append(f"- `{verdict}`: `{count}`")
        lines.append("")

    expectation_counts = dashboard_payload.get("parity_expectation_counts", {})
    if isinstance(expectation_counts, dict) and expectation_counts:
        lines.append("## Parity Expectations")
        lines.append("")
        for expectation, count in expectation_counts.items():
            lines.append(f"- `{expectation}`: `{count}`")
        lines.append("")

    comparison_kind_counts = dashboard_payload.get("comparison_kind_counts", {})
    if isinstance(comparison_kind_counts, dict) and comparison_kind_counts:
        lines.append("## Comparison Kinds")
        lines.append("")
        for comparison_kind, count in comparison_kind_counts.items():
            lines.append(f"- `{comparison_kind}`: `{count}`")
        lines.append("")

    lines.append("## Study Summary")
    lines.append("")
    for study in dashboard_payload["studies"]:
        lines.append(f"### {study['study_title']}")
        lines.append("")
        lines.append(f"- evidence id: `{study['evidence_id']}`")
        lines.append(f"- path: `{study['relative_path']}`")
        lines.append(f"- bundle verdict: `{study['bundle_verdict_status']}`")
        lines.append(f"- scalar rows: `{study['scalar_row_count']}`")
        if study["scalar_verdict_counts"]:
            counts = ", ".join(
                f"{key}={value}"
                for key, value in study["scalar_verdict_counts"].items()
            )
            lines.append(f"- scalar verdict counts: {counts}")
        if study["comparison_kind_counts"]:
            counts = ", ".join(
                f"{key}={value}"
                for key, value in study["comparison_kind_counts"].items()
            )
            lines.append(f"- comparison kinds: {counts}")
        if study["parity_expectation_counts"]:
            counts = ", ".join(
                f"{key}={value}"
                for key, value in study["parity_expectation_counts"].items()
            )
            lines.append(f"- parity expectations: {counts}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_evidence_catalog(index_payload: dict[str, object]) -> str:
    lines = [
        "# Evidence Catalog",
        "",
        "This catalog is generated from `index/evidence-index.json` and lists every",
        "governed evidence bundle in the repository evidence-book.",
        "",
        f"- studies: `{index_payload['study_count']}`",
        f"- evidence bundles: `{index_payload['evidence_count']}`",
        "",
    ]
    verdict_counts = index_payload.get("verdict_counts", {})
    if isinstance(verdict_counts, dict) and verdict_counts:
        lines.append("## Verdict Counts")
        lines.append("")
        for verdict, count in verdict_counts.items():
            lines.append(f"- `{verdict}`: `{count}`")
        lines.append("")

    lines.append("## Studies")
    lines.append("")
    for study in index_payload["studies"]:
        lines.append(f"### {study['study_title']}")
        lines.append("")
        lines.append(f"- study id: `{study['study_id']}`")
        lines.append(f"- owner package: `{study['owner_package']}`")
        lines.append(f"- categories: `{', '.join(study['study_categories'])}`")
        lines.append(f"- bundle count: `{study['bundle_count']}`")
        lines.append(f"- summary: {study['summary']}")
        lines.append("")
        for entry in study["evidence"]:
            lines.append(
                f"- `{entry['evidence_id']}` — {entry['evidence_title']} "
                f"(`{entry['verdict_status']}`)"
            )
            lines.append(f"  Path: `{entry['relative_path']}`")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_evidence_book_index(repo_root: Path) -> tuple[Path, Path]:
    root = evidence_book_root(repo_root)
    index_root = root / EVIDENCE_INDEX_DIRNAME
    index_root.mkdir(parents=True, exist_ok=True)
    teaching_migration_path = index_root / TEACHING_AND_MIGRATION_INDEX_FILENAME
    teaching_migration_summary_path = (
        index_root / TEACHING_AND_MIGRATION_SUMMARY_FILENAME
    )
    teaching_migration_payload = build_teaching_and_migration_index(
        [],
        [],
        [],
    )
    teaching_migration_path.write_text(
        json.dumps(teaching_migration_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    teaching_migration_summary_path.write_text(
        render_teaching_and_migration_index_markdown(teaching_migration_payload),
        encoding="utf-8",
    )
    payload = build_evidence_book_index(repo_root)
    claim_map_payload = build_evidence_claim_map(repo_root)
    parity_dashboard_payload = build_evidence_parity_dashboard(repo_root)
    mismatch_archive_payload = build_evidence_mismatch_archive(repo_root)
    verdict_workflows_payload = build_evidence_verdict_workflows(repo_root)
    false_confidence_payload = build_evidence_false_confidence_audit(repo_root)
    scientific_debt_payload = build_evidence_scientific_debt_register(repo_root)
    fragile_examples_payload = build_evidence_fragile_example_audit(repo_root)
    regeneration_payload = build_evidence_regeneration_contract(repo_root)
    index_path = index_root / EVIDENCE_INDEX
    catalog_path = index_root / EVIDENCE_CATALOG
    claim_map_path = index_root / EVIDENCE_CLAIM_MAP
    parity_dashboard_path = index_root / EVIDENCE_PARITY_DASHBOARD
    parity_summary_path = index_root / EVIDENCE_PARITY_SUMMARY
    mismatch_archive_path = index_root / EVIDENCE_MISMATCH_ARCHIVE
    mismatch_summary_path = index_root / EVIDENCE_MISMATCH_SUMMARY
    verdict_workflows_path = index_root / EVIDENCE_VERDICT_WORKFLOWS
    verdict_workflows_summary_path = index_root / EVIDENCE_VERDICT_WORKFLOWS_SUMMARY
    false_confidence_audit_path = index_root / EVIDENCE_FALSE_CONFIDENCE_AUDIT
    false_confidence_summary_path = index_root / EVIDENCE_FALSE_CONFIDENCE_SUMMARY
    scientific_debt_path = index_root / EVIDENCE_SCIENTIFIC_DEBT_REGISTER
    scientific_debt_summary_path = index_root / EVIDENCE_SCIENTIFIC_DEBT_SUMMARY
    portability_path = index_root / EVIDENCE_PORTABILITY_AUDIT
    portability_summary_path = index_root / EVIDENCE_PORTABILITY_SUMMARY
    fragile_examples_path = index_root / EVIDENCE_FRAGILE_EXAMPLE_AUDIT
    fragile_examples_summary_path = index_root / EVIDENCE_FRAGILE_EXAMPLE_SUMMARY
    regeneration_path = index_root / EVIDENCE_REGENERATION_CONTRACT
    regeneration_summary_path = index_root / EVIDENCE_REGENERATION_SUMMARY
    index_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    catalog_path.write_text(
        render_evidence_catalog(payload),
        encoding="utf-8",
    )
    claim_map_path.write_text(
        json.dumps(claim_map_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    parity_dashboard_path.write_text(
        json.dumps(parity_dashboard_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    parity_summary_path.write_text(
        render_evidence_parity_dashboard(parity_dashboard_payload),
        encoding="utf-8",
    )
    mismatch_archive_path.write_text(
        json.dumps(mismatch_archive_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    mismatch_summary_path.write_text(
        render_evidence_mismatch_archive(mismatch_archive_payload),
        encoding="utf-8",
    )
    verdict_workflows_path.write_text(
        json.dumps(verdict_workflows_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    verdict_workflows_summary_path.write_text(
        render_evidence_verdict_workflows(verdict_workflows_payload),
        encoding="utf-8",
    )
    false_confidence_audit_path.write_text(
        json.dumps(false_confidence_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    false_confidence_summary_path.write_text(
        render_evidence_false_confidence_audit(false_confidence_payload),
        encoding="utf-8",
    )
    scientific_debt_path.write_text(
        json.dumps(scientific_debt_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    scientific_debt_summary_path.write_text(
        render_evidence_scientific_debt_register(scientific_debt_payload),
        encoding="utf-8",
    )
    portability_payload = build_evidence_portability_audit(repo_root)
    portability_path.write_text(
        json.dumps(portability_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    portability_summary_path.write_text(
        render_evidence_portability_audit(portability_payload),
        encoding="utf-8",
    )
    fragile_examples_path.write_text(
        json.dumps(fragile_examples_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    fragile_examples_summary_path.write_text(
        render_evidence_fragile_example_audit(fragile_examples_payload),
        encoding="utf-8",
    )
    regeneration_path.write_text(
        json.dumps(regeneration_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    regeneration_summary_path.write_text(
        render_evidence_regeneration_contract(regeneration_payload),
        encoding="utf-8",
    )
    return index_path, catalog_path
