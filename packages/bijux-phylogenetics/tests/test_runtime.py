from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory
from bijux_phylogenetics.errors import InvalidBranchLengthError, UnsupportedTreeFormatError
from bijux_phylogenetics.identity import IDENTITY
from bijux_phylogenetics.io.nexus import load_nexus
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.io.phyloxml import load_phyloxml
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
    assert report.source_format == "newick"


def test_newick_loader_raises_invalid_branch_length_error() -> None:
    try:
        loads_newick("((A:abc,B:0.2):0.3,C:0.4);")
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


def test_nexus_loader_reads_translation_block_fixture() -> None:
    tree = load_nexus(FIXTURES / "example_tree.nex")
    assert tree.source_format == "nexus"
    assert tree.tip_names == ["A", "B", "C", "D"]
    assert tree.tip_count == 4


def test_phyloxml_loader_reads_annotated_tree_fixture() -> None:
    tree = load_phyloxml(FIXTURES / "example_tree.phyloxml")
    assert tree.source_format == "phyloxml"
    assert tree.tip_names == ["A", "B", "C"]
    assert tree.tip_count == 3


def test_validate_cli_reports_unsupported_format_error(tmp_path: Path, capsys) -> None:
    path = tmp_path / "tree.txt"
    path.write_text("not a tree\n", encoding="utf-8")
    exit_code = main(["validate", str(path), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"] == [
        {
            "code": UnsupportedTreeFormatError.code,
            "message": f"unsupported tree format for {path}",
        }
    ]


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
    assert payload["status"] == "ok"
    assert payload["command"] == "validate"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["data"]["tip_count"] == 4


def test_cli_commands_json_lists_registered_taxonomy(capsys) -> None:
    exit_code = main(["commands", "--format", "json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    command_names = [item["name"] for item in payload["data"]["commands"]]
    assert exit_code == 0
    assert payload["status"] == "ok"
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


def test_cli_report_json_output_uses_result_envelope(tmp_path: Path, capsys) -> None:
    output = tmp_path / "report.html"
    exit_code = main(
        [
            "report",
            "--tree",
            str(FIXTURES / "example_tree.nwk"),
            "--alignment",
            str(FIXTURES / "example_alignment.fasta"),
            "--traits",
            str(FIXTURES / "example_traits.tsv"),
            "--metadata",
            str(FIXTURES / "example_traits.tsv"),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "report"
    assert payload["outputs"] == [str(output)]
