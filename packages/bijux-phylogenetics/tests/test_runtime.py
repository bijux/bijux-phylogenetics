from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory
from bijux_phylogenetics.identity import IDENTITY
from bijux_phylogenetics.reports.service import annotate_tree_against_table, render_phylogenetics_report


FIXTURES = Path(__file__).parent / "fixtures"


def test_package_identity_matches_canonical_names() -> None:
    assert bijux_phylogenetics.__name__ == "bijux_phylogenetics"
    assert IDENTITY.package_name == "bijux-phylogenetics"
    assert IDENTITY.import_name == "bijux_phylogenetics"
    assert IDENTITY.cli_name == "bijux-phylogenetics"
    assert "bijux phylogenetics" == IDENTITY.umbrella_command


def test_validate_tree_path_reports_expected_counts() -> None:
    report = validate_tree_path(FIXTURES / "example_tree.nwk")
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.rooted is True
    assert report.ultrametric is True


def test_compare_tree_paths_reports_nonzero_distance() -> None:
    report = compare_tree_paths(FIXTURES / "example_tree.nwk", FIXTURES / "example_tree_alt.nwk")
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.robinson_foulds_distance > 0


def test_annotate_tree_against_table_finds_missing_and_extra_taxa() -> None:
    report = annotate_tree_against_table(FIXTURES / "example_tree.nwk", FIXTURES / "example_traits.tsv")
    assert report.linked_taxa == 3
    assert report.missing_from_table == ["D"]
    assert report.extra_table_entries == ["E"]


def test_render_phylogenetics_report_writes_html(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    result = render_phylogenetics_report(
        tree_path=FIXTURES / "example_tree.nwk",
        alignment_path=FIXTURES / "example_alignment.fasta",
        traits_path=FIXTURES / "example_traits.tsv",
        metadata_path=FIXTURES / "example_traits.tsv",
        out_path=output,
    )
    assert result.output_path == output
    assert output.exists()
    assert "Bijux Phylogenetics Report" in output.read_text(encoding="utf-8")


def test_bundle_directory_copies_files_and_manifest(tmp_path: Path) -> None:
    run_root = tmp_path / "run"
    run_root.mkdir()
    (run_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    report = bundle_directory(run_root, tmp_path / "bundle")
    assert report.file_count == 1
    manifest = (tmp_path / "bundle" / "manifest.json").read_text(encoding="utf-8")
    assert "artifact.txt" in manifest


def test_cli_validate_json_output(capsys) -> None:
    exit_code = main(["validate", str(FIXTURES / "example_tree.nwk"), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["tip_count"] == 4


def test_cli_commands_json_lists_registered_taxonomy(capsys) -> None:
    exit_code = main(["commands", "--format", "json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    command_names = [item["name"] for item in payload]
    assert exit_code == 0
    assert command_names == [
        "inspect",
        "validate",
        "normalize",
        "compare",
        "annotate",
        "diagnose",
        "render",
        "report",
        "evidence",
        "adapter",
    ]
