from __future__ import annotations

import json
from pathlib import Path

from ..portability import render_portability_rules_markdown
from ..teaching import (
    TEACHING_AND_MIGRATION_INDEX_FILENAME,
    TEACHING_AND_MIGRATION_SUMMARY_FILENAME,
    build_teaching_and_migration_index,
    render_teaching_and_migration_index_markdown,
)
from .aggregation import (
    build_evidence_book_index,
    build_evidence_claim_map,
    build_evidence_parity_dashboard,
)
from .audits import (
    build_evidence_false_confidence_audit as _build_evidence_false_confidence_audit,
    build_evidence_fragile_example_audit as _build_evidence_fragile_example_audit,
    build_evidence_mismatch_archive as _build_evidence_mismatch_archive,
    build_evidence_portability_audit as _build_evidence_portability_audit,
    build_evidence_regeneration_contract as _build_evidence_regeneration_contract,
    build_evidence_scientific_debt_register as _build_evidence_scientific_debt_register,
    build_evidence_verdict_workflows as _build_evidence_verdict_workflows,
)
from .layout import (
    EVIDENCE_CATALOG,
    EVIDENCE_CLAIM_MAP,
    EVIDENCE_FALSE_CONFIDENCE_AUDIT,
    EVIDENCE_FALSE_CONFIDENCE_SUMMARY,
    EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
    EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
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
)
from .models import (
    EvidenceBookValidationIssue as EvidenceBookValidationIssue,
    EvidenceBookValidationReport as EvidenceBookValidationReport,
)
from .validation import validate_evidence_book as validate_evidence_book


def build_evidence_mismatch_archive(repo_root: Path) -> dict[str, object]:
    return _build_evidence_mismatch_archive(repo_root)


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
    return _build_evidence_verdict_workflows(repo_root)


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
    return _build_evidence_false_confidence_audit(repo_root)


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
    return _build_evidence_portability_audit(repo_root)


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
    return _build_evidence_fragile_example_audit(repo_root)


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
    return _build_evidence_regeneration_contract(repo_root)


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
    return _build_evidence_scientific_debt_register(repo_root)


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
