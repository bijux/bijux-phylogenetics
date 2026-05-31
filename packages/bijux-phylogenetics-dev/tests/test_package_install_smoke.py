from __future__ import annotations

import os
from pathlib import Path

from bijux_phylogenetics_dev.quality.package_install_smoke import (
    select_artifact_paths,
    validate_alignment_payload,
    validate_example_input_probe,
    validate_pgls_payload,
    validate_resource_probe,
    validate_tree_package_payload,
)


def test_select_artifact_paths_supports_wheel_sdist_and_both(tmp_path: Path) -> None:
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir()
    wheel = dist_dir / "bijux_phylogenetics-0.1.0-py3-none-any.whl"
    sdist = dist_dir / "bijux_phylogenetics-0.1.0.tar.gz"
    wheel.write_text("wheel", encoding="utf-8")
    sdist.write_text("sdist", encoding="utf-8")

    assert select_artifact_paths(dist_dir, "wheel") == [("wheel", wheel)]
    assert select_artifact_paths(dist_dir, "sdist") == [("sdist", sdist)]
    assert select_artifact_paths(dist_dir, "both") == [
        ("wheel", wheel),
        ("sdist", sdist),
    ]


def test_select_artifact_paths_prefers_the_newest_matching_artifact(
    tmp_path: Path,
) -> None:
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir()
    older_wheel = dist_dir / "bijux_phylogenetics-0.1.0-py3-none-any.whl"
    newer_wheel = dist_dir / "bijux_phylogenetics-0.1.1-py3-none-any.whl"
    older_wheel.write_text("older", encoding="utf-8")
    newer_wheel.write_text("newer", encoding="utf-8")
    older_stat = older_wheel.stat()
    os.utime(older_wheel, ns=(older_stat.st_atime_ns, older_stat.st_mtime_ns))
    os.utime(
        newer_wheel,
        ns=(older_stat.st_atime_ns + 1, older_stat.st_mtime_ns + 1),
    )

    assert select_artifact_paths(dist_dir, "wheel") == [("wheel", newer_wheel)]


def test_validate_resource_probe_rejects_repository_source_paths(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    source_root = (
        repo_root / "packages" / "bijux-phylogenetics" / "src" / "bijux_phylogenetics"
    )
    resource_root = source_root / "resources"
    for relative in (
        ("examples", "alignments", "example_alignment.fasta"),
        ("examples", "trees", "example_tree.nwk"),
        ("datasets", "mammals", "primate_comparative", "tree.nwk"),
        ("datasets", "mammals", "primate_comparative", "traits.csv"),
        (
            "datasets",
            "pathogens",
            "rabies_cross_host_geography_panel",
            "workflow-config.json",
        ),
    ):
        path = resource_root.joinpath(*relative)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("sentinel\n", encoding="utf-8")
    report = {
        "package_root": source_root.as_posix(),
        "resource_paths": {
            "example_alignment": (
                resource_root / "examples" / "alignments" / "example_alignment.fasta"
            ).as_posix(),
            "example_tree": (
                resource_root / "examples" / "trees" / "example_tree.nwk"
            ).as_posix(),
            "primate_tree": (
                resource_root
                / "datasets"
                / "mammals"
                / "primate_comparative"
                / "tree.nwk"
            ).as_posix(),
            "primate_traits": (
                resource_root
                / "datasets"
                / "mammals"
                / "primate_comparative"
                / "traits.csv"
            ).as_posix(),
            "rabies_workflow_config": (
                resource_root
                / "datasets"
                / "pathogens"
                / "rabies_cross_host_geography_panel"
                / "workflow-config.json"
            ).as_posix(),
        },
        "missing_resources": [],
    }

    issues = validate_resource_probe(report, repo_root)

    issue_codes = {issue.code for issue in issues}
    assert "source-tree-package-root" in issue_codes
    assert "source-tree-resource-path" in issue_codes


def test_validate_alignment_and_pgls_payloads_check_expected_shapes() -> None:
    alignment_issues = validate_alignment_payload(
        {
            "status": "ok",
            "metrics": {
                "sequence_count": 4,
                "detected_type": "dna",
                "selected_type": "dna",
            },
        }
    )
    pgls_issues = validate_pgls_payload(
        {
            "status": "ok",
            "metrics": {
                "taxon_count": 75,
                "predictor_count": 1,
                "coefficient_count": 2,
            },
        }
    )

    assert alignment_issues == []
    assert pgls_issues == []


def test_validate_example_input_probe_rejects_source_tree_paths(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    source_root = (
        repo_root / "packages" / "bijux-phylogenetics" / "src" / "bijux_phylogenetics"
    )
    copied_dir = source_root / "resources" / "examples" / "copied"
    copied_dir.mkdir(parents=True, exist_ok=True)
    copied_paths = {
        "alignment": (copied_dir / "example_alignment.fasta").as_posix(),
        "alt_tree": (copied_dir / "example_tree_alt.nwk").as_posix(),
        "metadata": (copied_dir / "example_metadata.tsv").as_posix(),
        "traits": (copied_dir / "example_traits.tsv").as_posix(),
        "tree": (copied_dir / "example_tree.nwk").as_posix(),
    }
    report = {
        "destination": copied_dir.as_posix(),
        "copied_paths": copied_paths,
    }
    for path in copied_paths.values():
        Path(path).write_text("sentinel\n", encoding="utf-8")

    issues = validate_example_input_probe(report, repo_root)

    issue_codes = {issue.code for issue in issues}
    assert issue_codes == {"source-tree-example-input"}


def test_validate_tree_package_payload_requires_full_output_bundle(
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "tree-package"
    out_dir.mkdir()
    for name in (
        "tree-report.html",
        "tree-image.svg",
        "support-table.tsv",
        "clade-table.tsv",
        "branch-stats.tsv",
    ):
        (out_dir / name).write_text(name, encoding="utf-8")

    issues = validate_tree_package_payload(
        {"status": "ok", "metrics": {"tip_count": 4}},
        out_dir,
    )

    assert len(issues) == 1
    assert issues[0].code == "tree-package-output-missing"
    assert issues[0].path.endswith("tree-report.manifest.json")
