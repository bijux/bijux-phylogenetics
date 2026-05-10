from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.evidence.book import (
    build_evidence_book_index,
    validate_evidence_book,
    write_evidence_book_index,
)
from bijux_phylogenetics.evidence.scaffolding import (
    EvidenceBundleTemplateSpec,
    write_evidence_bundle_template,
)


def _write_scaling_fixture(root: Path, *, study_count: int, bundles_per_study: int) -> Path:
    book_root = root / "evidence-book"
    studies_root = book_root / "studies"
    studies_root.mkdir(parents=True, exist_ok=True)
    (book_root / "README.md").write_text("# Evidence Book\n", encoding="utf-8")
    for study_index in range(1, study_count + 1):
        study_id = f"scaling-study-{study_index:02d}"
        study_root = studies_root / study_id
        study_root.mkdir(parents=True, exist_ok=True)
        (study_root / "README.md").write_text(f"# {study_id}\n", encoding="utf-8")
        (study_root / "study.json").write_text(
            json.dumps(
                {
                    "study_id": study_id,
                    "study_title": f"Scaling Study {study_index:02d}",
                    "summary": "Synthetic scaling fixture study.",
                    "owner_package": "bijux-phylogenetics",
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        for bundle_index in range(1, bundles_per_study + 1):
            evidence_id = f"evidence-{bundle_index:03d}"
            write_evidence_bundle_template(
                study_root / evidence_id,
                EvidenceBundleTemplateSpec(
                    study_id=study_id,
                    evidence_id=evidence_id,
                    evidence_title=f"Synthetic Bundle {bundle_index:03d}",
                    summary="Synthetic scaling fixture bundle.",
                    owner_package="bijux-phylogenetics",
                    claim_id=f"{study_id}-claim-{bundle_index:03d}",
                    claim_title="Synthetic claim",
                    claim_summary="Synthetic scaling fixture claim.",
                    analytical_surfaces=("scaling",),
                    claim_tags=("scaling", "template"),
                    source_basis_locators=(
                        "packages/bijux-phylogenetics/tests/fixtures/trees/example_tree.nwk",
                    ),
                    governed_code_paths=(
                        "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/scaffolding.py",
                    ),
                    limitation="Synthetic scaling fixture only.",
                    freshness_date="2026-05-10",
                ),
            )
    return root


def test_write_evidence_book_index_scales_across_many_bundles(tmp_path: Path) -> None:
    repo_root = _write_scaling_fixture(tmp_path, study_count=4, bundles_per_study=30)

    write_evidence_book_index(repo_root)
    report = validate_evidence_book(repo_root)
    payload = build_evidence_book_index(repo_root)

    assert report.valid is True
    assert payload["study_count"] == 4
    assert payload["evidence_count"] == 120
    assert len(payload["evidence"]) == 120
