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
from .aggregation import (
    iter_scalar_parity_tables,
    validated_evidence_book_root,
)
from .layout import (
    EVIDENCE_BUNDLE_MANIFEST,
    EVIDENCE_FRAGILE_EXAMPLE_AUDIT,
    EVIDENCE_FRAGILE_EXAMPLE_SUMMARY,
    EVIDENCE_ID_PATTERN,
    EVIDENCE_INDEX_DIRNAME,
    EVIDENCE_PORTABILITY_AUDIT,
    EVIDENCE_PORTABILITY_SUMMARY,
    EVIDENCE_REGENERATION_CONTRACT,
    EVIDENCE_REGENERATION_SUMMARY,
    load_bundle_claim_rows,
    load_json,
    relative_to,
    study_paths,
)


def build_evidence_mismatch_archive(repo_root: Path) -> dict[str, object]:
    root = validated_evidence_book_root(repo_root)
    mismatches: list[dict[str, object]] = []
    verdict_counts: Counter[str] = Counter()
    for (
        study_manifest,
        manifest,
        bundle_root,
        scalar_table,
    ) in iter_scalar_parity_tables(repo_root):
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
                    "relative_path": relative_to(root, bundle_root).as_posix(),
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


def build_evidence_verdict_workflows(repo_root: Path) -> dict[str, object]:
    root = validated_evidence_book_root(repo_root)
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
    for study_root in study_paths(root):
        study_manifest = load_study_contract(study_root)
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            for row in load_bundle_claim_rows(bundle_root):
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
                        "relative_path": relative_to(root, bundle_root).as_posix(),
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


def build_evidence_portability_audit(repo_root: Path) -> dict[str, object]:
    root = validated_evidence_book_root(repo_root)
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
        relative_path = relative_to(root, path).as_posix()
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


def build_evidence_regeneration_contract(repo_root: Path) -> dict[str, object]:
    root = validated_evidence_book_root(repo_root)
    studies: list[dict[str, object]] = []
    for study_root in study_paths(root):
        study_manifest = load_study_contract(study_root)
        source_paths: list[str] = []
        generated_paths: list[str] = []
        for path in sorted(
            candidate for candidate in study_root.rglob("*") if candidate.is_file()
        ):
            relative_path = relative_to(Path(repo_root), path).as_posix()
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


def build_evidence_scientific_debt_register(repo_root: Path) -> dict[str, object]:
    root = validated_evidence_book_root(repo_root)
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

    for study_root in study_paths(root):
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            debt_register_path = bundle_root / "scientific_debt_register.json"
            if debt_register_path.exists():
                payload = load_json(debt_register_path)
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
                            "relative_path": relative_to(
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
                    payload = load_json(json_path)
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
                            "relative_path": relative_to(root, json_path).as_posix(),
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
