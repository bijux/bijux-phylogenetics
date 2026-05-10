from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics_dev.quality.evidence_inputs import (
    INPUT_MANIFEST_FILENAME,
    build_inputs_manifest,
    check_inputs_manifests,
    sync_inputs_manifests,
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
        bundle_root / "manifest.json",
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "demo-study",
                "evidence_id": "evidence-001",
                "evidence_title": "Demo bundle",
                "summary": "Demo bundle for input manifest checks.",
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
                        "label": "Local parity payload",
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
                "comparison_mode": "direct_r_parity",
                "verdict": {"status": "matched", "summary": "Demo verdict."},
                "limitations": [],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    _write(
        repo_root
        / "packages"
        / "bijux-phylogenetics"
        / "tests"
        / "fixtures"
        / "demo.tsv",
        "species\tvalue\nA\t1\n",
    )
    _write(bundle_root / "parity.json", '{"status": "ok"}\n')
    _write(bundle_root / "reference_table.csv", "species,value\nA,1\n")
    _write(bundle_root / "results" / "README.md", "# Results\n")
    _write(
        bundle_root / "results" / "manifest.json",
        json.dumps(
            {
                "schema_version": 1,
                "study_id": "demo-study",
                "evidence_id": "evidence-001",
                "results_directory": "results",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )
    return repo_root


def test_sync_inputs_manifests_writes_bundle_companion_files(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)

    written = sync_inputs_manifests(repo_root)

    assert [path.relative_to(repo_root).as_posix() for path in written] == [
        "evidence-book/studies/demo-study/evidence-001/inputs.manifest.json"
    ]
    payload = json.loads(written[0].read_text(encoding="utf-8"))
    assert payload["source_input_count"] == 2
    assert payload["governed_local_artifact_count"] == 4
    assert payload["local_input_count"] == 1
    assert {entry["locator"] for entry in payload["local_inputs"]} == {
        "evidence-book/studies/demo-study/evidence-001/reference_table.csv",
    }
    assert payload["local_inputs"][0]["input_class"] == "copied-reference-fragment"


def test_sync_inputs_manifests_supports_single_bundle_selection(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)

    written = sync_inputs_manifests(
        repo_root,
        study_id="demo-study",
        evidence_id="evidence-001",
    )

    assert [path.relative_to(repo_root).as_posix() for path in written] == [
        "evidence-book/studies/demo-study/evidence-001/inputs.manifest.json"
    ]


def test_sync_inputs_manifests_ignores_study_registry_files(tmp_path: Path) -> None:
    repo_root = _minimal_repo(tmp_path)
    _write(
        repo_root / "evidence-book" / "studies" / "demo-study" / "evidence-registry.json",
        "{}\n",
    )

    written = sync_inputs_manifests(repo_root)

    assert [path.relative_to(repo_root).as_posix() for path in written] == [
        "evidence-book/studies/demo-study/evidence-001/inputs.manifest.json"
    ]


def test_check_inputs_manifests_flags_stale_bundle_companion_files(
    tmp_path: Path,
) -> None:
    repo_root = _minimal_repo(tmp_path)
    sync_inputs_manifests(repo_root)
    manifest_path = (
        repo_root
        / "evidence-book"
        / "studies"
        / "demo-study"
        / "evidence-001"
        / INPUT_MANIFEST_FILENAME
    )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["local_input_count"] = 99
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    mismatches = check_inputs_manifests(repo_root)

    assert mismatches == [
        "evidence-book/studies/demo-study/evidence-001/inputs.manifest.json: stale input manifest"
    ]


def test_repository_inputs_manifests_are_synchronized() -> None:
    assert check_inputs_manifests(REPO_ROOT) == []


def test_repository_bundle_inputs_manifest_exists_for_every_bundle() -> None:
    bundle_roots = sorted(
        path
        for path in (REPO_ROOT / "evidence-book" / "studies").glob("*/evidence-*")
        if path.is_dir()
    )

    assert bundle_roots
    for bundle_root in bundle_roots:
        payload = build_inputs_manifest(REPO_ROOT, bundle_root)
        manifest_path = bundle_root / INPUT_MANIFEST_FILENAME
        assert manifest_path.is_file()
        actual = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert actual["study_id"] == payload["study_id"]
        assert actual["evidence_id"] == payload["evidence_id"]
