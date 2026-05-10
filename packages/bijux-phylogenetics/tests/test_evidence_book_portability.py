from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import validate_evidence_book


def test_validate_evidence_book_rejects_workstation_paths_in_study_provenance(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    study_root = repo_root / "evidence-book" / "studies" / "demo-study"
    bundle_root = study_root / "evidence-001"
    (repo_root / "evidence-book" / "README.md").parent.mkdir(parents=True, exist_ok=True)
    (repo_root / "evidence-book" / "README.md").write_text("# Evidence Book\n", encoding="utf-8")
    (repo_root / "evidence-book" / "index").mkdir(parents=True, exist_ok=True)
    (study_root / "README.md").write_text("# Demo Study\n\nFixture study.\n", encoding="utf-8")
    (study_root / "datasets").mkdir(parents=True, exist_ok=True)
    (study_root / "reference").mkdir(parents=True, exist_ok=True)
    (study_root / "provenance").mkdir(parents=True, exist_ok=True)
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
                        "schema_summary": "Fixture dataset.",
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
                        "locator": "/Users/demo/private.tsv",
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
    bundle_root.mkdir(parents=True, exist_ok=True)
    for filename, text in {
        "README.md": "# Evidence 001\n",
        "analysis.py": "from __future__ import annotations\n",
        "reference.R": "#!/usr/bin/env Rscript\n",
        "report.md": "# Report\n",
    }.items():
        (bundle_root / filename).write_text(text, encoding="utf-8")
    for filename, payload in {
        "manifest.json": {
            "schema_version": 1,
            "study_id": "demo-study",
            "evidence_id": "evidence-001",
            "evidence_title": "Demo bundle",
            "summary": "Fixture bundle.",
            "owner_package": "bijux-phylogenetics",
            "claim_ids": ["demo-claim"],
            "source_basis": [
                {
                    "kind": "repository-fixture",
                    "label": "Demo fixture",
                    "locator": "/Users/demo/private.tsv",
                }
            ],
            "freshness": {
                "last_generated_on": "2026-05-10",
                "governed_code_paths": [
                    "packages/bijux-phylogenetics/src/bijux_phylogenetics"
                ],
                "source_basis_locators": ["/Users/demo/private.tsv"],
            },
            "ownership": {
                "owner_package": "bijux-phylogenetics",
                "analytical_surfaces": ["demo-surface"],
            },
            "claim_tags": ["demo"],
            "comparison_mode": "direct_parity",
            "verdict": {"status": "matched", "summary": "Fixture aligned."},
            "limitations": [],
        },
        "claims.json": {
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
        "checks.json": {"schema_version": 1},
        "provenance.json": {"schema_version": 1},
        "inputs.manifest.json": {"schema_version": 1, "inputs": []},
        "results/manifest.json": {"schema_version": 1},
    }.items():
        target = bundle_root / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    (bundle_root / "results" / "README.md").write_text("# Results\n", encoding="utf-8")
    (bundle_root / "results" / "reviewer-summary.json").write_text("{}\n", encoding="utf-8")
    (bundle_root / "results" / "reviewer-summary.md").write_text("# Reviewer Summary\n", encoding="utf-8")

    report = validate_evidence_book(repo_root, require_index_outputs=False)

    assert report.valid is False
    assert any("/Users/" in issue.message or "/Users/" in issue.path.as_posix() for issue in report.issues)
