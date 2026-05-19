from __future__ import annotations

from collections import Counter
from pathlib import Path

from ..study_contracts import load_study_contract
from .layout import (
    EVIDENCE_BOOK_DIRNAME,
    EVIDENCE_BUNDLE_MANIFEST,
    EVIDENCE_ID_PATTERN,
    evidence_book_root,
    load_bundle_claim_rows,
    load_json,
    relative_to,
    study_paths,
)
from .validation import validate_evidence_book


def validated_evidence_book_root(repo_root: Path) -> Path:
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
        raise ValueError(f"{EVIDENCE_BOOK_DIRNAME} is invalid: {messages}")
    return root


def build_evidence_book_index(repo_root: Path) -> dict[str, object]:
    root = validated_evidence_book_root(repo_root)
    evidence_entries: list[dict[str, object]] = []
    studies: list[dict[str, object]] = []
    for study_root in study_paths(root):
        study_manifest = load_study_contract(study_root)
        study_entries: list[dict[str, object]] = []
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            verdict = manifest["verdict"]
            entry = {
                "study_id": study_manifest["study_id"],
                "study_title": study_manifest["study_title"],
                "evidence_id": manifest["evidence_id"],
                "evidence_title": manifest["evidence_title"],
                "owner_package": manifest["owner_package"],
                "relative_path": relative_to(root, bundle_root).as_posix(),
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
    root = validated_evidence_book_root(repo_root)
    claims_by_id: dict[str, dict[str, object]] = {}
    for study_root in study_paths(root):
        study_manifest = load_study_contract(study_root)
        for bundle_root in sorted(
            path
            for path in study_root.iterdir()
            if path.is_dir() and EVIDENCE_ID_PATTERN.fullmatch(path.name)
        ):
            manifest = load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST)
            for row in load_bundle_claim_rows(bundle_root):
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
                        "relative_path": relative_to(root, bundle_root).as_posix(),
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
    root = validated_evidence_book_root(repo_root)
    studies: list[dict[str, object]] = []
    total_row_count = 0
    verdict_counts: Counter[str] = Counter()
    comparison_kind_counts: Counter[str] = Counter()

    for (
        study_manifest,
        manifest,
        bundle_root,
        scalar_table,
    ) in iter_scalar_parity_tables(repo_root):
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
                "relative_path": relative_to(root, bundle_root).as_posix(),
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


def iter_scalar_parity_tables(
    repo_root: Path,
) -> list[tuple[dict[str, object], dict[str, object], Path, dict[str, object]]]:
    root = validated_evidence_book_root(repo_root)
    tables: list[
        tuple[dict[str, object], dict[str, object], Path, dict[str, object]]
    ] = []
    for study_root in study_paths(root):
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
                    load_json(bundle_root / EVIDENCE_BUNDLE_MANIFEST),
                    bundle_root,
                    load_json(scalar_table_path),
                )
            )
    return tables
