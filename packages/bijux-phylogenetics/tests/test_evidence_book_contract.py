from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import validate_evidence_book
from bijux_phylogenetics.evidence.bundle_artifacts import build_bundle_governed_artifacts
from bijux_phylogenetics.evidence.bundle_contracts import (
    ARTIFACT_JSON_FILENAMES,
    RESULT_ARTIFACT_JSON_FILENAMES,
)
from bijux_phylogenetics.evidence.workbench import refresh_evidence_book


def _write_governed_bundle_fixture(repo_root: Path) -> Path:
    study_root = repo_root / "evidence-book" / "studies" / "demo-study"
    bundle_root = study_root / "evidence-001"
    (repo_root / "evidence-book" / "index").mkdir(parents=True, exist_ok=True)
    (study_root / "datasets").mkdir(parents=True, exist_ok=True)
    (study_root / "provenance").mkdir(parents=True, exist_ok=True)
    (study_root / "reference").mkdir(parents=True, exist_ok=True)
    bundle_root.mkdir(parents=True, exist_ok=True)
    (repo_root / "evidence-book" / "README.md").write_text("# Evidence Book\n", encoding="utf-8")
    (study_root / "README.md").write_text(
        "# Demo Study\n\nHuman-first evidence fixture.\n", encoding="utf-8"
    )
    (study_root / "datasets" / "registry.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "demo-study",
                "dataset_count": 1,
                "datasets": [
                    {
                        "dataset_id": "dataset-001",
                        "kind": "repository-fixture",
                        "label": "Demo fixture",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
                        "schema_summary": "Demo dataset fixture.",
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (study_root / "provenance" / "sources.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "demo-study",
                "source_count": 1,
                "intake_policy": "repository-owned-source",
                "sources": [
                    {
                        "source_id": "demo-fixture",
                        "kind": "repository-fixture",
                        "label": "Demo fixture",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
                        "read_only": True,
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (bundle_root / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "demo-study",
                "evidence_id": "evidence-001",
                "evidence_title": "Demo bundle",
                "summary": "Demo governed bundle.",
                "owner_package": "bijux-phylogenetics",
                "claim_ids": ["demo-claim"],
                "source_basis": [
                    {
                        "kind": "repository-fixture",
                        "label": "Demo fixture",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
                    }
                ],
                "freshness": {
                    "last_generated_on": "2026-05-10",
                    "governed_code_paths": [
                        "packages/bijux-phylogenetics/src/bijux_phylogenetics"
                    ],
                    "source_basis_locators": [
                        "packages/bijux-phylogenetics/tests/fixtures/demo.tsv"
                    ],
                },
                "ownership": {
                    "owner_package": "bijux-phylogenetics",
                    "analytical_surfaces": ["demo-surface"],
                },
                "claim_tags": ["demo"],
                "comparison_mode": "direct_parity",
                "verdict": {
                    "status": "matched",
                    "summary": "Fixture is aligned.",
                },
                "limitations": [],
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
                "study_id": "demo-study",
                "evidence_id": "evidence-001",
                "claim_count": 1,
                "claims": [
                    {
                        "claim_id": "demo-claim",
                        "claim_title": "Demo claim",
                        "summary": "Demo summary.",
                        "verdict": "matched",
                        "evidence_ids": ["evidence-001"],
                        "source_fragments": ["demo-fragment"],
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (bundle_root / "README.md").write_text(
        "# Demo Evidence\n\nHuman-first evidence fixture.\n",
        encoding="utf-8",
    )
    for relative_path, payload in build_bundle_governed_artifacts(repo_root, bundle_root).items():
        target = bundle_root / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        if relative_path in ARTIFACT_JSON_FILENAMES or relative_path in RESULT_ARTIFACT_JSON_FILENAMES:
            target.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        else:
            assert isinstance(payload, str)
            target.write_text(payload, encoding="utf-8")
    refresh_evidence_book(repo_root)
    return repo_root


def test_validate_evidence_book_accepts_minimal_human_first_layout(tmp_path: Path) -> None:
    repo_root = _write_governed_bundle_fixture(tmp_path)

    report = validate_evidence_book(repo_root)

    assert report.valid is True
    assert report.issues == []


def test_validate_evidence_book_rejects_study_root_machine_files(tmp_path: Path) -> None:
    repo_root = _write_governed_bundle_fixture(tmp_path)
    study_root = repo_root / "evidence-book" / "studies" / "demo-study"
    (study_root / "study.json").write_text("{}\n", encoding="utf-8")

    report = validate_evidence_book(repo_root, require_index_outputs=False)

    assert report.valid is False
    assert any("study root may only contain" in issue.message for issue in report.issues)


def test_validate_evidence_book_rejects_bundle_root_machine_files(tmp_path: Path) -> None:
    repo_root = _write_governed_bundle_fixture(tmp_path)
    bundle_root = (
        repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-001"
    )
    (bundle_root / "comparison.json").write_text("{}\n", encoding="utf-8")

    report = validate_evidence_book(repo_root, require_index_outputs=False)

    assert report.valid is False
    assert any("move machine outputs under results/" in issue.message for issue in report.issues)
