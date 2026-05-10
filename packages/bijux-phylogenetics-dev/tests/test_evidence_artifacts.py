from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics_dev.quality.evidence_artifacts import (
    check_evidence_artifacts,
    sync_evidence_artifacts,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _minimal_repo(tmp_path: Path) -> Path:
    repo_root = tmp_path / "repo"
    bundle_root = (
        repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-001"
    )
    _write(
        repo_root / "evidence-book" / "studies" / "demo-study" / "study.json",
        json.dumps(
            {
                "study_id": "demo-study",
                "study_title": "Demo Study",
                "summary": "Small governed study for bundle artifact checks.",
                "owner_package": "bijux-phylogenetics",
                "study_categories": ["scientific-validation"],
                "provenance_descriptor_locator": "evidence-book/studies/demo-study/provenance/sources.json",
                "dataset_registry_locator": "evidence-book/studies/demo-study/datasets/registry.json",
                "source_intake_policy": "repository-owned-source",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    _write(
        repo_root / "evidence-book" / "studies" / "demo-study" / "provenance" / "sources.json",
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
    )
    _write(
        repo_root / "evidence-book" / "studies" / "demo-study" / "datasets" / "registry.json",
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
                        "schema_summary": "Small fixture for artifact rendering.",
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    _write(
        bundle_root / "manifest.json",
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "demo-study",
                "evidence_id": "evidence-001",
                "evidence_title": "Demo bundle",
                "summary": "Demo bundle for local artifact checks.",
                "owner_package": "bijux-phylogenetics",
                "claim_ids": ["demo-claim"],
                "source_basis": [
                    {
                        "kind": "repository-fixture",
                        "label": "Demo fixture",
                        "locator": "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
                    },
                    {
                        "kind": "repository-reference",
                        "label": "Demo parity payload",
                        "locator": "evidence-book/studies/demo-study/evidence-001/parity.json",
                    },
                ],
                "freshness": {
                    "last_generated_on": "2026-05-10",
                    "governed_code_paths": [
                        "packages/bijux-phylogenetics/src/bijux_phylogenetics"
                    ],
                    "source_basis_locators": [
                        "packages/bijux-phylogenetics/tests/fixtures/demo.tsv",
                        "evidence-book/studies/demo-study/evidence-001/parity.json",
                    ],
                },
                "ownership": {
                    "owner_package": "bijux-phylogenetics",
                    "analytical_surfaces": ["demo-analysis"],
                },
                "claim_tags": ["demo"],
                "comparison_mode": "direct_parity",
                "verdict": {"status": "matched", "summary": "Demo verdict."},
                "limitations": ["Fixture-backed demo only."],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    _write(
        bundle_root / "claims.json",
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
                        "summary": "Demo claim summary.",
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
    )
    _write(bundle_root / "README.md", "# Demo Bundle\n")
    _write(bundle_root / "reviewer-summary.json", "{}\n")
    _write(bundle_root / "reviewer-summary.md", "# Reviewer Summary\n")
    _write(bundle_root / "parity.json", '{"status": "ok"}\n')
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "fixtures"
        / "demo.tsv",
        "species\tvalue\nA\t1\n",
    )
    return repo_root


def test_sync_evidence_artifacts_writes_local_bundle_surfaces(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)

    written = sync_evidence_artifacts(repo_root)

    written_paths = [path.relative_to(repo_root).as_posix() for path in written]
    assert "evidence-book/studies/demo-study/evidence-001/reference.R" in written_paths
    assert "evidence-book/studies/demo-study/evidence-001/analysis.py" in written_paths
    assert "evidence-book/studies/demo-study/evidence-001/checks.json" in written_paths
    assert "evidence-book/studies/demo-study/evidence-001/report.md" in written_paths
    assert "evidence-book/studies/demo-study/evidence-001/provenance.json" in written_paths


def test_check_evidence_artifacts_flags_stale_bundle_surface(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)
    sync_evidence_artifacts(repo_root)
    report_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / "report.md"
    )
    report_path.write_text("# stale\n", encoding="utf-8")

    mismatches = check_evidence_artifacts(repo_root)

    assert mismatches == [
        "evidence-book/studies/demo-study/evidence-001/report.md: stale governed local artifact"
    ]


def test_repository_evidence_artifacts_are_synchronized() -> None:
    assert check_evidence_artifacts(REPO_ROOT) == []
