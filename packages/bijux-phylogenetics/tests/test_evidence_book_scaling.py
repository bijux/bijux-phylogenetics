from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import validate_evidence_book
from bijux_phylogenetics.evidence.bundle_artifacts import (
    build_bundle_governed_artifacts,
)
from bijux_phylogenetics.evidence.bundle_contracts import (
    ARTIFACT_JSON_FILENAMES,
    RESULT_ARTIFACT_JSON_FILENAMES,
)
from bijux_phylogenetics.evidence.reviewer import (
    build_bundle_reviewer_summary,
    encode_bundle_reviewer_summary,
    render_bundle_reviewer_summary,
)
from bijux_phylogenetics.evidence.study_contracts import load_study_contract


def test_validate_evidence_book_scales_across_many_minimal_bundles(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    (repo_root / "evidence-book" / "README.md").parent.mkdir(
        parents=True, exist_ok=True
    )
    (repo_root / "evidence-book" / "README.md").write_text(
        "# Evidence Book\n", encoding="utf-8"
    )
    (repo_root / "evidence-book" / "index").mkdir(parents=True, exist_ok=True)

    for study_index in range(5):
        study_root = (
            repo_root / "evidence-book" / "studies" / f"demo-study-{study_index}"
        )
        study_root.mkdir(parents=True, exist_ok=True)
        (study_root / "README.md").write_text(
            f"# Demo Study {study_index}\n\nFixture study.\n", encoding="utf-8"
        )
        (study_root / "datasets").mkdir(parents=True, exist_ok=True)
        (study_root / "provenance").mkdir(parents=True, exist_ok=True)
        (study_root / "reference").mkdir(parents=True, exist_ok=True)
        (study_root / "datasets" / "registry.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "study_id": f"demo-study-{study_index}",
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
                    "study_id": f"demo-study-{study_index}",
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
        for evidence_index in range(10):
            bundle_root = study_root / f"evidence-{evidence_index + 1:03d}"
            bundle_root.mkdir(parents=True, exist_ok=True)
            (bundle_root / "README.md").write_text(
                f"# Evidence {evidence_index + 1:03d}\n", encoding="utf-8"
            )
            (bundle_root / "manifest.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "study_id": f"demo-study-{study_index}",
                        "evidence_id": f"evidence-{evidence_index + 1:03d}",
                        "evidence_title": "Demo bundle",
                        "summary": "Fixture bundle.",
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
                        "verdict": {"status": "matched", "summary": "Fixture aligned."},
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
                        "study_id": f"demo-study-{study_index}",
                        "evidence_id": f"evidence-{evidence_index + 1:03d}",
                        "claim_count": 1,
                        "claims": [
                            {
                                "claim_id": "demo-claim",
                                "claim_title": "Demo claim",
                                "summary": "Demo summary.",
                                "verdict": "matched",
                                "evidence_ids": [f"evidence-{evidence_index + 1:03d}"],
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
            reviewer_summary = build_bundle_reviewer_summary(
                study_manifest=load_study_contract(study_root),
                bundle_manifest=json.loads(
                    (bundle_root / "manifest.json").read_text(encoding="utf-8")
                ),
                claims_payload=json.loads(
                    (bundle_root / "claims.json").read_text(encoding="utf-8")
                ),
            )
            for relative_path, payload in build_bundle_governed_artifacts(
                repo_root, bundle_root
            ).items():
                target = bundle_root / relative_path
                target.parent.mkdir(parents=True, exist_ok=True)
                if (
                    relative_path in ARTIFACT_JSON_FILENAMES
                    or relative_path in RESULT_ARTIFACT_JSON_FILENAMES
                ):
                    target.write_text(
                        json.dumps(payload, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                    )
                else:
                    assert isinstance(payload, str)
                    target.write_text(payload, encoding="utf-8")
            (bundle_root / "results" / "reviewer-summary.json").write_text(
                encode_bundle_reviewer_summary(reviewer_summary),
                encoding="utf-8",
            )
            (bundle_root / "results" / "reviewer-summary.md").write_text(
                render_bundle_reviewer_summary(reviewer_summary),
                encoding="utf-8",
            )

    report = validate_evidence_book(repo_root, require_index_outputs=False)
    assert report.valid is True
    assert len(report.bundle_paths) == 50
