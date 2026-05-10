from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_portability_audit,
    validate_evidence_book,
    write_evidence_book_index,
)


def _write_portable_fixture(root: Path) -> Path:
    book_root = root / "evidence-book"
    study_root = book_root / "studies" / "portable-study"
    bundle_root = study_root / "evidence-001"
    bundle_root.mkdir(parents=True, exist_ok=True)
    (book_root / "README.md").write_text("# Evidence Book\n", encoding="utf-8")
    (study_root / "README.md").write_text("# Portable Study\n", encoding="utf-8")
    (study_root / "study.json").write_text(
        json.dumps(
            {
                "study_id": "portable-study",
                "study_title": "Portable Study",
                "summary": "Fixture-backed portability study.",
                "owner_package": "bijux-phylogenetics",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (bundle_root / "README.md").write_text("# Evidence 001\n", encoding="utf-8")
    (bundle_root / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "portable-study",
                "evidence_id": "evidence-001",
                "evidence_title": "Portable bundle",
                "summary": "Portable path semantics fixture.",
                "owner_package": "bijux-phylogenetics",
                "claim_ids": ["portable-claim"],
                "source_basis": [
                    {
                        "kind": "repository-fixture",
                        "label": "portable source",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/trees/example_tree.nwk",
                    }
                ],
                "freshness": {
                    "last_generated_on": "2026-05-10",
                    "governed_code_paths": [
                        "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/book.py"
                    ],
                    "source_basis_locators": [
                        "packages/bijux-phylogenetics/tests/fixtures/trees/example_tree.nwk"
                    ],
                },
                "ownership": {
                    "owner_package": "bijux-phylogenetics",
                    "analytical_surfaces": ["portability"],
                },
                "claim_tags": ["portability"],
                "verdict": {
                    "status": "matched",
                    "summary": "Portable fixture output matches expectations.",
                },
                "limitations": ["Fixture scope only."],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (bundle_root / "claims.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "portable-study",
                "evidence_id": "evidence-001",
                "claim_count": 1,
                "claims": [
                    {
                        "claim_id": "portable-claim",
                        "claim_title": "Portable claim",
                        "summary": "Portable fixture claim.",
                        "verdict": "matched",
                        "evidence_ids": ["evidence-001"],
                        "source_fragments": [],
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return root


def test_validate_evidence_book_rejects_workstation_absolute_locator(tmp_path: Path) -> None:
    repo_root = _write_portable_fixture(tmp_path)
    write_evidence_book_index(repo_root)

    manifest_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "portable-study"
        / "evidence-001"
        / "manifest.json"
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["source_basis"][0]["locator"] = "/Users/bijan/private/reference.csv"
    payload["freshness"]["source_basis_locators"] = ["/Users/bijan/private/reference.csv"]
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_evidence_book(repo_root)

    assert report.valid is False
    assert any("workstation-local absolute paths" in issue.message for issue in report.issues)


def test_validate_evidence_book_rejects_parent_traversal_locator(tmp_path: Path) -> None:
    repo_root = _write_portable_fixture(tmp_path)
    write_evidence_book_index(repo_root)

    manifest_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "portable-study"
        / "evidence-001"
        / "manifest.json"
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["source_basis"][0]["locator"] = "../bijux-pollenomics/packages/example.csv"
    payload["freshness"]["source_basis_locators"] = [
        "../bijux-pollenomics/packages/example.csv"
    ]
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_evidence_book(repo_root)
    portability_audit = build_evidence_portability_audit(repo_root)

    assert report.valid is False
    assert any("parent-directory traversal" in issue.message for issue in report.issues)
    assert portability_audit["action_required_count"] >= 1
