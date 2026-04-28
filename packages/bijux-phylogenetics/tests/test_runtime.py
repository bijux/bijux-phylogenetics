from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.core.metadata import inspect_metadata_table
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path, validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory
from bijux_phylogenetics.errors import (
    DuplicateTaxonError,
    InvalidBranchLengthError,
    MetadataJoinError,
    UnnamedTipError,
    UnsupportedTreeFormatError,
)
from bijux_phylogenetics.identity import IDENTITY
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.io.nexus import load_nexus
from bijux_phylogenetics.io.phyloxml import load_phyloxml
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.reports.service import annotate_tree_against_table, render_phylogenetics_report


FIXTURES = Path(__file__).parent / "fixtures"


def test_package_identity_matches_canonical_names() -> None:
    assert bijux_phylogenetics.__name__ == "bijux_phylogenetics"
    assert IDENTITY.package_name == "bijux-phylogenetics"
    assert IDENTITY.import_name == "bijux_phylogenetics"
    assert IDENTITY.cli_name == "bijux-phylogenetics"
    assert "bijux phylogenetics" == IDENTITY.umbrella_command


def test_taxon_labels_preserve_raw_names_and_normalized_keys() -> None:
    tree = loads_newick("('Homo sapiens':0.1,'NCBI|123/45':0.2,'A.B-1':0.3);")
    assert [(taxon.raw, taxon.key) for taxon in tree.tip_taxa] == [
        ("Homo sapiens", "Homo_sapiens"),
        ("NCBI|123/45", "NCBI_123_45"),
        ("A.B-1", "A.B-1"),
    ]
    assert dumps_newick(tree) == "(A.B-1:0.3,'Homo sapiens':0.1,'NCBI|123/45':0.2);"


def test_metadata_inspect_reports_taxon_contract() -> None:
    report = inspect_metadata_table(FIXTURES / "example_metadata.tsv")
    assert report.format == "tsv"
    assert report.row_count == 4
    assert report.column_count == 3
    assert report.taxon_column == "taxon"
    assert report.taxa == ["A", "B", "C", "D"]


def test_metadata_inspect_rejects_duplicate_taxa() -> None:
    try:
        inspect_metadata_table(FIXTURES / "example_metadata_duplicate.tsv")
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_metadata_inspect_rejects_missing_requested_taxon_column() -> None:
    try:
        inspect_metadata_table(FIXTURES / "example_metadata_missing_taxon.csv", taxon_column="taxon")
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
        assert error.message == "missing taxon column 'taxon'"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_validate_tree_path_reports_expected_counts() -> None:
    report = validate_tree_path(FIXTURES / "example_tree.nwk")
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.rooted is True
    assert report.ultrametric is True
    assert report.source_format == "newick"


def test_validate_tree_path_rejects_duplicate_tip_labels_by_default() -> None:
    try:
        validate_tree_path(FIXTURES / "example_tree_duplicate.nwk")
    except DuplicateTaxonError as error:
        assert error.code == "duplicate_taxon_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected DuplicateTaxonError")


def test_validate_tree_path_warns_for_unnamed_tips_in_non_strict_mode() -> None:
    report = validate_tree_path(FIXTURES / "example_tree_unnamed_tip.nwk")
    assert report.missing_taxa == 1
    assert "tree contains unnamed tips" in report.warnings


def test_validate_tree_path_rejects_unnamed_tips_in_strict_mode() -> None:
    try:
        validate_tree_path(FIXTURES / "example_tree_unnamed_tip.nwk", strict=True)
    except UnnamedTipError as error:
        assert error.code == "unnamed_tip_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected UnnamedTipError")


def test_validate_tree_path_rejects_negative_branch_lengths_by_default() -> None:
    try:
        validate_tree_path(FIXTURES / "example_tree_negative_length.nwk")
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


def test_inspect_tree_path_returns_normalized_json_summary_contract() -> None:
    report = inspect_tree_path(FIXTURES / "example_tree.nwk")
    assert report.tip_count == 4
    assert report.node_count == 7
    assert report.internal_node_count == 3
    assert report.edge_count == 6
    assert report.clade_count == 3
    assert report.has_branch_lengths is True
    assert report.is_binary is True
    assert report.max_depth == 2
    assert report.taxa == ["A", "B", "C", "D"]


def test_inspect_tree_path_distinguishes_rooted_and_unrooted_fixtures() -> None:
    rooted = inspect_tree_path(FIXTURES / "example_tree.nwk")
    unrooted = inspect_tree_path(FIXTURES / "example_tree_unrooted.nwk")
    assert rooted.rooted is True
    assert unrooted.rooted is False


def test_inspect_tree_path_reports_exact_polytomy_nodes() -> None:
    report = inspect_tree_path(FIXTURES / "example_tree_polytomy.nwk")
    assert report.is_binary is False
    assert report.polytomy_count == 1
    assert report.polytomy_nodes == ["A|B|C"]


def test_inspect_tree_path_classifies_branch_length_completeness() -> None:
    complete = inspect_tree_path(FIXTURES / "example_tree.nwk")
    partial = inspect_tree_path(FIXTURES / "example_tree_partial_lengths.nwk")
    absent = inspect_tree_path(FIXTURES / "example_tree_no_lengths.nwk")
    assert complete.branch_length_status == "complete"
    assert partial.branch_length_status == "partial"
    assert absent.branch_length_status == "absent"


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


def test_detect_tree_format_uses_filename_suffixes() -> None:
    assert detect_tree_format(Path("x.nwk")) == "newick"
    assert detect_tree_format(Path("x.nex")) == "nexus"
    assert detect_tree_format(Path("x.phyloxml")) == "phyloxml"


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


def test_validate_cli_can_allow_duplicate_tip_labels(capsys) -> None:
    exit_code = main(["validate", str(FIXTURES / "example_tree_duplicate.nwk"), "--allow-duplicates", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["duplicate_taxa"] == ["A"]


def test_validate_cli_strict_mode_rejects_unnamed_tips(capsys) -> None:
    exit_code = main(["validate", str(FIXTURES / "example_tree_unnamed_tip.nwk"), "--strict", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"] == [
        {
            "code": UnnamedTipError.code,
            "message": "tree contains 1 unnamed tip labels",
        }
    ]


def test_validate_cli_can_allow_negative_branch_lengths(capsys) -> None:
    exit_code = main(
        ["validate", str(FIXTURES / "example_tree_negative_length.nwk"), "--allow-negative-branches", "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["negative_branch_lengths"] == 1


def test_cli_inspect_accepts_explicit_tree_format(capsys) -> None:
    exit_code = main(["inspect", str(FIXTURES / "example_tree.nex"), "--format", "nexus", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["source_format"] == "nexus"
    assert payload["data"]["node_count"] == 7
    assert payload["data"]["edge_count"] == 6
    assert payload["data"]["clade_count"] == 3
    assert payload["data"]["taxa"] == ["A", "B", "C", "D"]
    assert payload["metrics"]["tip_count"] == 4


def test_cli_normalize_writes_canonical_newick(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized.nwk"
    exit_code = main(
        ["normalize", str(FIXTURES / "example_tree.nex"), "--format", "nexus", "--out", str(output), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert output.read_text(encoding="utf-8").strip() == "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"


def test_cli_normalize_taxa_writes_mapping_file(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized-taxa.nwk"
    mapping = tmp_path / "normalized-taxa.tsv"
    exit_code = main(
        [
            "normalize-taxa",
            str(FIXTURES / "example_tree_labels.nwk"),
            "--policy",
            "spaces-to-underscores",
            "--out",
            str(output),
            "--mapping-out",
            str(mapping),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["renamed_taxa"] == 2
    assert output.read_text(encoding="utf-8").strip() == "(A.B-1:0.3,Homo_sapiens:0.1,Mus_musculus:0.2);"
    assert mapping.read_text(encoding="utf-8") == (
        "raw_label\tnormalized_label\n"
        "Homo sapiens\tHomo_sapiens\n"
        "Mus musculus\tMus_musculus\n"
    )


def test_compare_tree_paths_reports_nonzero_distance() -> None:
    report = compare_tree_paths(FIXTURES / "example_tree.nwk", FIXTURES / "example_tree_alt.nwk")
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.robinson_foulds_distance > 0


def test_annotate_tree_against_table_finds_missing_and_extra_taxa() -> None:
    report = annotate_tree_against_table(FIXTURES / "example_tree.nwk", FIXTURES / "example_traits.tsv")
    assert report.linked_taxa == 3
    assert report.missing_from_table == ["D"]
    assert report.extra_table_entries == ["E"]


def test_cli_metadata_inspect_json_output(capsys) -> None:
    exit_code = main(["metadata", "inspect", str(FIXTURES / "example_metadata.tsv"), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "metadata"
    assert payload["data"]["taxon_column"] == "taxon"
    assert payload["metrics"]["taxon_count"] == 4


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
        "metadata",
        "inspect",
        "validate",
        "normalize",
        "normalize-taxa",
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
