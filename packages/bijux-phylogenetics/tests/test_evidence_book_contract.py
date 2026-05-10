from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_claim_map,
    build_evidence_book_index,
    build_evidence_false_confidence_audit,
    build_evidence_mismatch_archive,
    build_evidence_parity_dashboard,
    build_evidence_verdict_workflows,
    render_evidence_catalog,
    render_evidence_false_confidence_audit,
    render_evidence_mismatch_archive,
    render_evidence_parity_dashboard,
    render_evidence_verdict_workflows,
    validate_evidence_book,
    write_evidence_book_index,
)


def _write_book_fixture(root: Path) -> Path:
    book_root = root / "evidence-book"
    study_root = book_root / "studies" / "taxon-trust"
    bundle_root = study_root / "evidence-001"
    index_root = book_root / "index"
    bundle_root.mkdir(parents=True, exist_ok=True)
    index_root.mkdir(parents=True, exist_ok=True)
    (book_root / "README.md").write_text("# Evidence Book\n", encoding="utf-8")
    (study_root / "README.md").write_text("# Taxon Trust\n", encoding="utf-8")
    (study_root / "study.json").write_text(
        json.dumps(
            {
                "study_id": "taxon-trust",
                "study_title": "Taxon Trust",
                "summary": "Fixture-backed taxon evidence.",
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
                "study_id": "taxon-trust",
                "evidence_id": "evidence-001",
                "evidence_title": "Taxon workflow review bundle",
                "summary": "Validates taxon-workflow trust surfaces.",
                "owner_package": "bijux-phylogenetics",
                "claim_ids": ["taxonomy-review"],
                "source_basis": [
                    {
                        "kind": "repository-fixture",
                        "label": "example taxon workflow tree",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk",
                    }
                ],
                "freshness": {
                    "last_generated_on": "2026-05-10",
                    "governed_code_paths": [
                        "packages/bijux-phylogenetics/src/bijux_phylogenetics"
                    ],
                    "source_basis_locators": [
                        "packages/bijux-phylogenetics/tests/fixtures/trees/example_taxon_workflow_tree.nwk"
                    ],
                },
                "ownership": {
                    "owner_package": "bijux-phylogenetics",
                    "analytical_surfaces": ["taxonomy"],
                },
                "claim_tags": ["taxonomy", "review"],
                "verdict": {
                    "status": "matched",
                    "summary": "Observed output matches the checked-in fixture expectations.",
                },
                "limitations": ["Covers one workflow family only."],
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
                "study_id": "taxon-trust",
                "evidence_id": "evidence-001",
                "claim_count": 1,
                "claims": [
                    {
                        "claim_id": "taxonomy-review",
                        "claim_title": "Taxonomy review",
                        "summary": "Fixture-backed taxonomy review contract.",
                        "verdict": "matched",
                        "evidence_ids": ["evidence-001"],
                        "source_fragments": ["fixture-taxonomy-review"],
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


def test_validate_evidence_book_accepts_governed_layout(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)
    write_evidence_book_index(repo_root)

    report = validate_evidence_book(repo_root)

    assert report.valid is True
    assert report.issues == []
    assert [path.name for path in report.bundle_paths] == ["evidence-001"]


def test_validate_evidence_book_rejects_missing_manifest_fields(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)
    manifest_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "taxon-trust"
        / "evidence-001"
        / "manifest.json"
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    del payload["verdict"]
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report = validate_evidence_book(repo_root, require_index_outputs=False)

    assert report.valid is False
    assert any("missing keys: verdict" in issue.message for issue in report.issues)


def test_write_evidence_book_index_renders_catalog_from_index(tmp_path: Path) -> None:
    repo_root = _write_book_fixture(tmp_path)

    index_path, catalog_path = write_evidence_book_index(repo_root)
    payload = build_evidence_book_index(repo_root)
    catalog = render_evidence_catalog(payload)
    claim_map = build_evidence_claim_map(repo_root)
    parity_dashboard = build_evidence_parity_dashboard(repo_root)
    mismatch_archive = build_evidence_mismatch_archive(repo_root)
    verdict_workflows = build_evidence_verdict_workflows(repo_root)
    false_confidence_audit = build_evidence_false_confidence_audit(repo_root)
    claim_map_path = repo_root / "evidence-book" / "index" / "claim-map.json"
    parity_dashboard_path = repo_root / "evidence-book" / "index" / "parity-dashboard.json"
    parity_summary_path = repo_root / "evidence-book" / "index" / "parity-dashboard.md"
    mismatch_archive_path = repo_root / "evidence-book" / "index" / "mismatch-archive.json"
    mismatch_summary_path = repo_root / "evidence-book" / "index" / "mismatch-archive.md"
    verdict_workflows_path = repo_root / "evidence-book" / "index" / "verdict-workflows.json"
    verdict_workflows_summary_path = repo_root / "evidence-book" / "index" / "verdict-workflows.md"
    false_confidence_audit_path = repo_root / "evidence-book" / "index" / "false-confidence-audit.json"
    false_confidence_summary_path = repo_root / "evidence-book" / "index" / "false-confidence-audit.md"

    assert index_path.exists()
    assert catalog_path.exists()
    assert claim_map_path.exists()
    assert parity_dashboard_path.exists()
    assert parity_summary_path.exists()
    assert mismatch_archive_path.exists()
    assert mismatch_summary_path.exists()
    assert verdict_workflows_path.exists()
    assert verdict_workflows_summary_path.exists()
    assert false_confidence_audit_path.exists()
    assert false_confidence_summary_path.exists()
    assert payload["study_count"] == 1
    assert payload["evidence_count"] == 1
    assert "Taxon Trust" in catalog
    assert catalog_path.read_text(encoding="utf-8") == catalog
    assert json.loads(claim_map_path.read_text(encoding="utf-8")) == claim_map
    assert json.loads(parity_dashboard_path.read_text(encoding="utf-8")) == parity_dashboard
    assert parity_summary_path.read_text(encoding="utf-8") == render_evidence_parity_dashboard(
        parity_dashboard
    )
    assert json.loads(mismatch_archive_path.read_text(encoding="utf-8")) == mismatch_archive
    assert mismatch_summary_path.read_text(encoding="utf-8") == render_evidence_mismatch_archive(
        mismatch_archive
    )
    assert json.loads(verdict_workflows_path.read_text(encoding="utf-8")) == verdict_workflows
    assert verdict_workflows_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_verdict_workflows(verdict_workflows)
    assert json.loads(
        false_confidence_audit_path.read_text(encoding="utf-8")
    ) == false_confidence_audit
    assert false_confidence_summary_path.read_text(
        encoding="utf-8"
    ) == render_evidence_false_confidence_audit(false_confidence_audit)
