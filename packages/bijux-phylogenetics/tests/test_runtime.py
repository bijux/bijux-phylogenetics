from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.compare.topology import compare_branch_lengths, compare_support_values, compare_tree_paths
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.core.alignment import AlignmentSummary
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.core.metadata import inspect_metadata_table
from bijux_phylogenetics.core.pruning import prune_alignment_to_tree, prune_tree_to_alignment, prune_tree_to_taxa
from bijux_phylogenetics.core.traits import detect_unusable_trait_columns, link_tree_to_traits, validate_traits_table
from bijux_phylogenetics.diagnostics.root_to_tip import compute_root_to_tip_distances
from bijux_phylogenetics.diagnostics.root_to_tip import diagnose_ultrametricity
from bijux_phylogenetics.diagnostics.validation import diagnose_tree_path, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory, validate_bundle
from bijux_phylogenetics.errors import (
    AlignmentTaxonMismatchError,
    DuplicateTaxonError,
    EngineUnavailableError,
    InvalidBranchLengthError,
    InvalidAlignmentError,
    MetadataJoinError,
    NonUltrametricTreeError,
    UnnamedTipError,
    UnsupportedTreeFormatError,
    UnrootedTreeError,
)
from bijux_phylogenetics.identity import IDENTITY
from bijux_phylogenetics.io.newick import dumps_newick, loads_newick
from bijux_phylogenetics.io.nexus import load_nexus
from bijux_phylogenetics.io.phyloxml import load_phyloxml
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.io.fasta import link_alignment_to_tree, load_fasta_alignment, summarise_fasta
from bijux_phylogenetics.io.fasta import (
    detect_sequences_with_excessive_missing_data,
    detect_sites_with_excessive_missing_data,
    write_fasta_alignment,
)
from bijux_phylogenetics.render.svg import render_tree_svg
from bijux_phylogenetics.reports.service import (
    annotate_tree_against_table,
    render_dataset_report,
    render_phylo_inputs_report,
    render_phylogenetics_report,
    render_tree_report,
)


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


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
    report = inspect_metadata_table(fixture("example_metadata.tsv"))
    assert report.format == "tsv"
    assert report.row_count == 4
    assert report.column_count == 3
    assert report.taxon_column == "taxon"
    assert report.taxa == ["A", "B", "C", "D"]
    assert [(row.name, row.missing_count, row.completeness_fraction) for row in report.column_completeness] == [
        ("taxon", 0, 1.0),
        ("species", 0, 1.0),
        ("location", 0, 1.0),
    ]


def test_inspect_environment_reports_available_and_optional_dependencies() -> None:
    report = inspect_environment()
    status_by_name = {item.name: item for item in report.dependencies}
    assert report.python_version
    assert report.host_platform
    assert status_by_name["biopython"].available is True
    assert "dendropy" in status_by_name


def test_metadata_inspect_rejects_duplicate_taxa() -> None:
    try:
        inspect_metadata_table(fixture("example_metadata_duplicate.tsv"))
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_metadata_inspect_rejects_missing_requested_taxon_column() -> None:
    try:
        inspect_metadata_table(fixture("example_metadata_missing_taxon.csv"), taxon_column="taxon")
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
        assert error.message == "missing taxon column 'taxon'"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_traits_validate_infers_numeric_and_categorical_schema() -> None:
    report = validate_traits_table(fixture("example_traits_validate.tsv"))
    assert report.taxon_column == "taxon"
    assert [
        (column.name, column.kind, column.missing_count, column.missing_fraction)
        for column in report.trait_columns
    ] == [
        ("height_cm", "numeric", 0, 0.0),
        ("habitat", "categorical", 0, 0.0),
        ("status", "categorical", 1, 0.25),
    ]


def test_traits_detect_unusable_columns_by_missingness() -> None:
    columns = detect_unusable_trait_columns(
        fixture("example_traits_validate.tsv"),
        missingness_threshold=0.2,
    )
    assert [(column.name, column.missing_fraction) for column in columns] == [("status", 0.25)]


def test_traits_link_reports_mismatch_and_usable_taxa() -> None:
    report = link_tree_to_traits(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
    assert report.tree_taxa == 4
    assert report.trait_taxa == 4
    assert report.linked_taxa == 3
    assert report.usable_taxa == ["A", "B", "C"]
    assert report.missing_from_traits == ["D"]
    assert report.extra_trait_taxa == ["E"]


def test_traits_link_strict_mode_rejects_mismatch() -> None:
    try:
        link_tree_to_traits(fixture("example_tree.nwk"), fixture("example_traits.tsv"), strict=True)
    except MetadataJoinError as error:
        assert error.code == "metadata_join_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected MetadataJoinError")


def test_prune_tree_to_taxa_writes_expected_tip_set() -> None:
    tree, report = prune_tree_to_taxa(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
    assert tree.tip_names == ["A", "B", "C"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]


def test_prune_alignment_to_tree_keeps_exact_tree_taxa() -> None:
    records, report = prune_alignment_to_tree(
        fixture("example_alignment_extra_taxon.fasta"),
        fixture("example_tree.nwk"),
    )
    assert [record.identifier for record in records] == ["A", "B", "C", "D"]
    assert report.original_sequence_count == 5
    assert report.kept_ids == ["A", "B", "C", "D"]
    assert report.removed_ids == ["E"]


def test_prune_tree_to_alignment_keeps_exact_alignment_taxa() -> None:
    tree, report = prune_tree_to_alignment(
        fixture("example_tree.nwk"),
        fixture("example_alignment_mismatch.fasta"),
    )
    assert tree.tip_names == ["A", "B", "C"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]
    assert report.taxon_column == "identifier"


def test_alignment_inspect_reports_core_diagnostics() -> None:
    report = summarise_fasta(fixture("example_alignment.fasta"))
    assert isinstance(report, AlignmentSummary)
    assert report.sequence_count == 4
    assert report.alignment_length == 8
    assert report.ids == ["A", "B", "C", "D"]
    assert report.missing_data_fraction == 0.0
    assert report.gap_fraction == 0.0
    assert report.constant_site_count == 6
    assert report.variable_site_count == 2
    assert report.parsimony_informative_site_count == 2
    assert [(row.identifier, row.missing_fraction) for row in report.per_sequence_missingness] == [
        ("A", 0.0),
        ("B", 0.0),
        ("C", 0.0),
        ("D", 0.0),
    ]


def test_alignment_inspect_reports_per_sequence_missingness() -> None:
    report = summarise_fasta(fixture("example_alignment_missingness.fasta"))
    assert report.sequence_count == 3
    assert report.alignment_length == 6
    assert report.constant_site_count == 6
    assert report.variable_site_count == 0
    assert report.parsimony_informative_site_count == 0
    assert report.missing_data_fraction == 4 / 18
    assert [(row.identifier, row.missing_fraction) for row in report.per_sequence_missingness] == [
        ("A", 2 / 6),
        ("B", 2 / 6),
        ("C", 0.0),
    ]
    assert [(row.position, row.missing_fraction) for row in report.per_site_missingness] == [
        (1, 0.0),
        (2, 0.0),
        (3, 0.0),
        (4, 0.0),
        (5, 2 / 3),
        (6, 2 / 3),
    ]
    assert report.all_gap_columns == []
    assert report.all_missing_columns == []


def test_alignment_inspect_reports_site_missingness_and_empty_columns() -> None:
    report = summarise_fasta(fixture("example_alignment_site_missingness.fasta"))
    assert report.alignment_length == 5
    assert [(row.position, row.missing_fraction) for row in report.per_site_missingness] == [
        (1, 0.0),
        (2, 0.0),
        (3, 1.0),
        (4, 1.0),
        (5, 0.5),
    ]
    assert report.all_gap_columns == [2]
    assert report.all_missing_columns == [3, 4]


def test_alignment_inspect_rejects_unequal_lengths() -> None:
    try:
        summarise_fasta(fixture("example_alignment_invalid_lengths.fasta"))
    except InvalidAlignmentError as error:
        assert error.code == "invalid_alignment_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidAlignmentError")


def test_alignment_link_reports_exact_mismatch() -> None:
    report = link_alignment_to_tree(fixture("example_tree.nwk"), fixture("example_alignment.fasta"))
    assert report.tree_taxa == 4
    assert report.alignment_ids == 4
    assert report.linked_taxa == 4
    assert report.missing_from_alignment == []
    assert report.extra_alignment_ids == []


def test_alignment_link_strict_mode_rejects_mismatch() -> None:
    try:
        link_alignment_to_tree(fixture("example_tree.nwk"), fixture("example_alignment_mismatch.fasta"), strict=True)
    except AlignmentTaxonMismatchError as error:
        assert error.code == "alignment_taxon_mismatch_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected AlignmentTaxonMismatchError")


def test_write_fasta_alignment_preserves_record_order_and_sequences(tmp_path: Path) -> None:
    records = load_fasta_alignment(fixture("example_alignment.fasta"))
    output = tmp_path / "alignment.fasta"
    write_fasta_alignment(output, records)
    assert output.read_text(encoding="utf-8") == (
        ">A\nACTGACTG\n"
        ">B\nACTGACTA\n"
        ">C\nACTGACGG\n"
        ">D\nACTGACGA\n"
    )
    assert load_fasta_alignment(output) == records


def test_alignment_detects_sequences_with_excessive_missing_data() -> None:
    rows = detect_sequences_with_excessive_missing_data(
        fixture("example_alignment_missingness.fasta"),
        threshold=0.3,
    )
    assert [(row.identifier, row.missing_fraction) for row in rows] == [
        ("A", 2 / 6),
        ("B", 2 / 6),
    ]


def test_alignment_detects_sites_with_excessive_missing_data() -> None:
    rows = detect_sites_with_excessive_missing_data(
        fixture("example_alignment_site_missingness.fasta"),
        threshold=0.4,
    )
    assert [(row.position, row.missing_fraction) for row in rows] == [
        (3, 1.0),
        (4, 1.0),
        (5, 0.5),
    ]


def test_build_run_manifest_captures_checksums_and_environment(tmp_path: Path) -> None:
    input_path = tmp_path / "input.txt"
    output_path = tmp_path / "output.txt"
    input_path.write_text("input\n", encoding="utf-8")
    output_path.write_text("output\n", encoding="utf-8")
    manifest = build_run_manifest(
        command="inspect",
        arguments=["inspect", str(input_path), "--json"],
        input_paths=[input_path],
        output_paths=[output_path],
    )
    manifest_path = write_run_manifest(tmp_path / "run.manifest.json", manifest)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["command"] == "inspect"
    assert payload["arguments"] == ["inspect", str(input_path), "--json"]
    assert payload["input_paths"] == [str(input_path)]
    assert payload["output_paths"] == [str(output_path)]
    assert payload["input_checksums"][str(input_path)]
    assert payload["output_checksums"][str(output_path)]
    assert payload["python_version"]
    assert payload["host_platform"]


def test_validate_tree_path_reports_expected_counts() -> None:
    report = validate_tree_path(fixture("example_tree.nwk"))
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.rooted is True
    assert report.ultrametric is True
    assert report.source_format == "newick"


def test_validate_tree_path_rejects_duplicate_tip_labels_by_default() -> None:
    try:
        validate_tree_path(fixture("example_tree_duplicate.nwk"))
    except DuplicateTaxonError as error:
        assert error.code == "duplicate_taxon_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected DuplicateTaxonError")


def test_validate_tree_path_warns_for_unnamed_tips_in_non_strict_mode() -> None:
    report = validate_tree_path(fixture("example_tree_unnamed_tip.nwk"))
    assert report.missing_taxa == 1
    assert "tree contains unnamed tips" in report.warnings


def test_validate_tree_path_rejects_unnamed_tips_in_strict_mode() -> None:
    try:
        validate_tree_path(fixture("example_tree_unnamed_tip.nwk"), strict=True)
    except UnnamedTipError as error:
        assert error.code == "unnamed_tip_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected UnnamedTipError")


def test_validate_tree_path_rejects_negative_branch_lengths_by_default() -> None:
    try:
        validate_tree_path(fixture("example_tree_negative_length.nwk"))
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


def test_validate_tree_path_can_require_rooted_tree() -> None:
    try:
        validate_tree_path(fixture("example_tree_unrooted.nwk"), require_rooted=True)
    except UnrootedTreeError as error:
        assert error.code == "unrooted_tree_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected UnrootedTreeError")


def test_validate_tree_path_can_require_ultrametric_tree() -> None:
    try:
        validate_tree_path(fixture("example_tree_ladderized.nwk"), require_ultrametric=True)
    except NonUltrametricTreeError as error:
        assert error.code == "non_ultrametric_tree_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected NonUltrametricTreeError")


def test_validate_tree_path_warns_for_zero_length_branches() -> None:
    report = validate_tree_path(fixture("example_tree_zero_lengths.nwk"))
    assert report.zero_length_branches == 3
    assert "tree contains zero-length branches" in report.warnings


def test_inspect_tree_path_returns_normalized_json_summary_contract() -> None:
    report = inspect_tree_path(fixture("example_tree.nwk"))
    assert report.tip_count == 4
    assert report.node_count == 7
    assert report.internal_node_count == 3
    assert report.edge_count == 6
    assert report.clade_count == 3
    assert report.has_branch_lengths is True
    assert report.is_binary is True
    assert report.is_ultrametric is True
    assert report.zero_length_branch_count == 0
    assert report.max_depth == 2
    assert report.mean_depth == 2.0
    assert report.imbalance_summary == "balanced"
    assert report.cherry_count == 2
    assert report.warnings == []
    assert report.taxa == ["A", "B", "C", "D"]


def test_inspect_tree_path_distinguishes_ladderized_shape() -> None:
    report = inspect_tree_path(fixture("example_tree_ladderized.nwk"))
    assert report.max_depth == 3
    assert report.mean_depth == 2.25
    assert report.imbalance_summary == "ladderized"
    assert report.cherry_count == 1


def test_inspect_tree_path_distinguishes_rooted_and_unrooted_fixtures() -> None:
    rooted = inspect_tree_path(fixture("example_tree.nwk"))
    unrooted = inspect_tree_path(fixture("example_tree_unrooted.nwk"))
    assert rooted.rooted is True
    assert unrooted.rooted is False


def test_inspect_tree_path_reports_exact_polytomy_nodes() -> None:
    report = inspect_tree_path(fixture("example_tree_polytomy.nwk"))
    assert report.is_binary is False
    assert report.polytomy_count == 1
    assert report.polytomy_nodes == ["A|B|C"]


def test_inspect_tree_path_classifies_branch_length_completeness() -> None:
    complete = inspect_tree_path(fixture("example_tree.nwk"))
    partial = inspect_tree_path(fixture("example_tree_partial_lengths.nwk"))
    absent = inspect_tree_path(fixture("example_tree_no_lengths.nwk"))
    assert complete.branch_length_status == "complete"
    assert partial.branch_length_status == "partial"
    assert absent.branch_length_status == "absent"
    assert partial.has_branch_lengths is True
    assert partial.warnings == ["tree contains partial branch lengths"]
    assert absent.has_branch_lengths is False
    assert absent.warnings == ["tree contains no branch lengths"]


def test_newick_loader_raises_invalid_branch_length_error() -> None:
    try:
        loads_newick("((A:abc,B:0.2):0.3,C:0.4);")
    except InvalidBranchLengthError as error:
        assert error.code == "invalid_branch_length_error"
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidBranchLengthError")


def test_nexus_loader_reads_translation_block_fixture() -> None:
    tree = load_nexus(fixture("example_tree.nex"))
    assert tree.source_format == "nexus"
    assert tree.tip_names == ["A", "B", "C", "D"]
    assert tree.tip_count == 4


def test_phyloxml_loader_reads_annotated_tree_fixture() -> None:
    tree = load_phyloxml(fixture("example_tree.phyloxml"))
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
    exit_code = main(["validate", str(fixture("example_tree_duplicate.nwk")), "--allow-duplicates", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["duplicate_taxa"] == ["A"]


def test_validate_cli_strict_mode_rejects_unnamed_tips(capsys) -> None:
    exit_code = main(["validate", str(fixture("example_tree_unnamed_tip.nwk")), "--strict", "--json"])
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
        ["validate", str(fixture("example_tree_negative_length.nwk")), "--allow-negative-branches", "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["negative_branch_lengths"] == 1


def test_validate_cli_can_require_rooted_and_ultrametric_typed_errors(capsys) -> None:
    exit_code = main(["validate", str(fixture("example_tree_unrooted.nwk")), "--require-rooted", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == UnrootedTreeError.code

    exit_code = main(["validate", str(fixture("example_tree_ladderized.nwk")), "--require-ultrametric", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == NonUltrametricTreeError.code


def test_cli_inspect_accepts_explicit_tree_format(capsys) -> None:
    exit_code = main(["inspect", str(fixture("example_tree.nex")), "--format", "nexus", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["source_format"] == "nexus"
    assert payload["data"]["node_count"] == 7
    assert payload["data"]["edge_count"] == 6
    assert payload["data"]["clade_count"] == 3
    assert payload["data"]["mean_depth"] == 2.0
    assert payload["data"]["imbalance_summary"] == "balanced"
    assert payload["metrics"]["cherry_count"] == 2
    assert payload["data"]["taxa"] == ["A", "B", "C", "D"]
    assert payload["metrics"]["tip_count"] == 4


def test_cli_normalize_writes_canonical_newick(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized.nwk"
    exit_code = main(
        ["normalize", str(fixture("example_tree.nex")), "--format", "nexus", "--out", str(output), "--json"]
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
            str(fixture("example_tree_labels.nwk")),
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
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.robinson_foulds_distance > 0


def test_compare_support_values_pairs_shared_clades() -> None:
    report = compare_support_values(fixture("example_tree_support_left.nwk"), fixture("example_tree_support_right.nwk"))
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert [(row.split_id, row.left_support, row.right_support) for row in report.shared_clades] == [
        ("A|B", 95.0, 90.0),
        ("C|D", 88.0, 85.0),
    ]


def test_compare_branch_lengths_reports_delta_ratio_and_missing_lengths() -> None:
    scaled = compare_branch_lengths(fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_right.nwk"))
    missing = compare_branch_lengths(fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_missing.nwk"))
    assert [(row.split_id, row.delta, row.ratio) for row in scaled.shared_splits] == [
        ("A|B", 0.2, 2.0),
        ("C|D", 0.1, 2.0),
    ]
    assert [(row.split_id, row.left_length, row.right_length, row.delta, row.ratio) for row in missing.shared_splits] == [
        ("A|B", 0.2, None, None, None),
        ("C|D", 0.1, 0.2, 0.1, 2.0),
    ]


def test_build_tree_comparison_report_writes_html_with_checksums(tmp_path: Path) -> None:
    output = tmp_path / "compare.html"
    result = build_tree_comparison_report(
        fixture("example_tree_support_left.nwk"),
        fixture("example_tree_support_right.nwk"),
        out_path=output,
    )
    html = output.read_text(encoding="utf-8")
    assert result.output_path == output
    assert "Bijux Tree Comparison Report" in html
    assert "input-checksums" in html
    assert "support-comparison" in html


def test_render_tree_svg_writes_static_tree_image(tmp_path: Path) -> None:
    output = tmp_path / "tree.svg"
    result = render_tree_svg(fixture("example_tree.nwk"), out_path=output)
    svg = output.read_text(encoding="utf-8")
    assert result.output_path == output
    assert result.format == "svg"
    assert "<svg" in svg
    assert "A" in svg and "D" in svg


def test_render_tree_svg_can_use_metadata_labels(tmp_path: Path) -> None:
    output = tmp_path / "annotated.svg"
    result = render_tree_svg(
        fixture("example_tree.nwk"),
        out_path=output,
        labels={"A": "Alpha species", "B": "Beta species", "C": "Gamma species"},
    )
    svg = output.read_text(encoding="utf-8")
    assert "Alpha species" in svg
    assert "Beta species" in svg
    assert result.missing_metadata_labels == ["D"]


def test_diagnose_tree_path_combines_inspection_and_validation() -> None:
    report = diagnose_tree_path(fixture("example_tree.nwk"))
    assert report.inspection.tip_count == 4
    assert report.validation.tip_count == 4


def test_compute_root_to_tip_distances_reports_one_row_per_tip() -> None:
    report = compute_root_to_tip_distances(fixture("example_tree.nwk"))
    assert [(row.tip, row.distance) for row in report.distances] == [
        ("A", 0.30000000000000004),
        ("B", 0.30000000000000004),
        ("C", 0.30000000000000004),
        ("D", 0.30000000000000004),
    ]


def test_diagnose_ultrametricity_reports_max_deviation() -> None:
    ultrametric = diagnose_ultrametricity(fixture("example_tree.nwk"), tolerance=1e-6)
    non_ultrametric = diagnose_ultrametricity(fixture("example_tree_ladderized.nwk"), tolerance=1e-6)
    assert ultrametric.ultrametric is True
    assert ultrametric.max_deviation == 0.0
    assert non_ultrametric.ultrametric is False
    assert non_ultrametric.max_deviation == 0.2


def test_annotate_tree_against_table_finds_missing_and_extra_taxa() -> None:
    report = annotate_tree_against_table(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
    assert report.linked_taxa == 3
    assert report.annotated_taxa == ["A", "B", "C"]
    assert report.missing_from_table == ["D"]
    assert report.extra_table_entries == ["E"]


def test_cli_metadata_inspect_json_output(capsys) -> None:
    exit_code = main(["metadata", "inspect", str(fixture("example_metadata.tsv")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "metadata"
    assert payload["data"]["taxon_column"] == "taxon"
    assert payload["metrics"]["taxon_count"] == 4


def test_cli_env_inspect_json_output(capsys) -> None:
    exit_code = main(["env", "inspect", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "env"
    assert payload["metrics"]["dependency_count"] >= 1


def test_cli_annotate_writes_annotation_json(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree.annotated.json"
    exit_code = main(
        [
            "annotate",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_traits.tsv")),
            "--taxon-column",
            "taxon",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    artifact = json.loads(output.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert artifact["annotated_taxa"] == ["A", "B", "C"]


def test_cli_traits_validate_json_output(capsys) -> None:
    exit_code = main(["traits", "validate", str(fixture("example_traits_validate.tsv")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "traits"
    assert payload["metrics"]["trait_column_count"] == 3
    assert payload["data"]["trait_columns"][0]["kind"] == "numeric"


def test_cli_traits_link_json_output(capsys) -> None:
    exit_code = main(
        ["traits", "link", str(fixture("example_tree.nwk")), str(fixture("example_traits.tsv")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["linked_taxa"] == 3
    assert payload["data"]["missing_from_traits"] == ["D"]


def test_cli_traits_link_strict_mode_returns_typed_error(capsys) -> None:
    exit_code = main(
        [
            "traits",
            "link",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--strict",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == MetadataJoinError.code


def test_cli_prune_writes_tree_and_pruned_taxa_report(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree.pruned.nwk"
    exit_code = main(
        [
            "prune",
            str(fixture("example_tree.nwk")),
            "--keep-from",
            str(fixture("example_traits.tsv")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert output.read_text(encoding="utf-8").strip() == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert payload["data"]["removed_taxa"] == ["D"]
    assert (tmp_path / "pruned_taxa.tsv").read_text(encoding="utf-8") == "taxon\nD\n"


def test_cli_alignment_inspect_json_output(capsys) -> None:
    exit_code = main(["alignment", "inspect", str(fixture("example_alignment.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "alignment"
    assert payload["metrics"]["alignment_length"] == 8
    assert payload["data"]["variable_site_count"] == 2


def test_cli_alignment_link_json_output(capsys) -> None:
    exit_code = main(
        ["alignment", "link", str(fixture("example_tree.nwk")), str(fixture("example_alignment.fasta")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_alignment_link_strict_mode_returns_typed_error(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "link",
            str(fixture("example_tree.nwk")),
            str(fixture("example_alignment_mismatch.fasta")),
            "--strict",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == AlignmentTaxonMismatchError.code


def test_cli_compare_support_json_output(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "support",
            str(fixture("example_tree_support_left.nwk")),
            str(fixture("example_tree_support_right.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["shared_clades"] == 2
    assert payload["data"]["shared_clades"][0]["split_id"] == "A|B"


def test_cli_compare_branch_lengths_json_output(capsys) -> None:
    exit_code = main(
        [
            "compare",
            "branch-lengths",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_branch_lengths_right.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["shared_splits"] == 2
    assert payload["data"]["shared_splits"][0]["ratio"] == 2.0


def test_cli_compare_report_json_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "compare.html"
    exit_code = main(
        [
            "compare",
            "report",
            str(fixture("example_tree_support_left.nwk")),
            str(fixture("example_tree_support_right.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert "Bijux Tree Comparison Report" in output.read_text(encoding="utf-8")


def test_cli_render_writes_svg_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree.svg"
    exit_code = main(["render", str(fixture("example_tree.nwk")), "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert "<svg" in output.read_text(encoding="utf-8")


def test_cli_render_with_metadata_labels_reports_missing_taxa(tmp_path: Path, capsys) -> None:
    output = tmp_path / "annotated.svg"
    exit_code = main(
        [
            "render",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--label-column",
            "species",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["missing_metadata_labels"] == []
    assert "Alpha species" in output.read_text(encoding="utf-8")


def test_cli_render_with_partial_metadata_warns_for_missing_labels(tmp_path: Path, capsys) -> None:
    output = tmp_path / "partial.svg"
    exit_code = main(
        [
            "render",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_traits.tsv")),
            "--label-column",
            "value",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["warnings"] == ["D"]


def test_render_phylogenetics_report_writes_html(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    result = render_phylogenetics_report(
        tree_path=fixture("example_tree.nwk"),
        alignment_path=fixture("example_alignment.fasta"),
        traits_path=fixture("example_traits.tsv"),
        metadata_path=fixture("example_traits.tsv"),
        out_path=output,
    )
    assert result.output_path == output
    assert output.exists()
    assert "Bijux Phylogenetics Report" in output.read_text(encoding="utf-8")


def test_render_tree_report_writes_embedded_manifest(tmp_path: Path) -> None:
    output = tmp_path / "tree-report.html"
    result = render_tree_report(tree_path=fixture("example_tree.nwk"), out_path=output)
    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "tree"
    assert result.machine_manifest["report_kind"] == "tree"
    assert result.machine_manifest["input_paths"] == [str(fixture("example_tree.nwk"))]
    assert 'id="bijux-report-manifest"' in text
    assert "Bijux Tree Report" in text


def test_render_dataset_report_writes_metadata_sections(tmp_path: Path) -> None:
    output = tmp_path / "dataset-report.html"
    result = render_dataset_report(
        tree_path=fixture("example_tree.nwk"),
        metadata_path=fixture("example_metadata.tsv"),
        traits_path=fixture("example_traits.tsv"),
        out_path=output,
    )
    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "dataset"
    assert result.metadata_linkage is not None
    assert result.machine_manifest["sections"] == [
        "tree-validation",
        "tree-inspection",
        "metadata-linkage",
        "traits-linkage",
    ]
    assert "Bijux Dataset Report" in text


def test_render_phylo_inputs_report_writes_alignment_sections(tmp_path: Path) -> None:
    output = tmp_path / "phylo-inputs-report.html"
    result = render_phylo_inputs_report(
        tree_path=fixture("example_tree.nwk"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=output,
    )
    text = output.read_text(encoding="utf-8")
    assert result.report_kind == "phylo-inputs"
    assert result.alignment is not None
    assert result.alignment_linkage is not None
    assert result.machine_manifest["sections"] == [
        "tree-validation",
        "tree-inspection",
        "alignment-summary",
        "alignment-linkage",
    ]
    assert "Bijux Phylo Inputs Report" in text


def test_bundle_directory_copies_files_and_manifest(tmp_path: Path) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "tree.nwk").write_text("(A:0.1,B:0.1);\n", encoding="utf-8")
    (outputs_root / "report.html").write_text("<html></html>\n", encoding="utf-8")
    report = bundle_directory([inputs_root], [outputs_root], tmp_path / "bundle")
    assert report.file_count == 2
    assert report.input_file_count == 1
    assert report.output_file_count == 1
    manifest = (tmp_path / "bundle" / "manifest.json").read_text(encoding="utf-8")
    assert "tree.nwk" in manifest
    assert (tmp_path / "bundle" / "checksums.tsv").exists()
    assert (tmp_path / "bundle" / "environment.json").exists()
    assert (tmp_path / "bundle" / "README.md").exists()


def test_validate_bundle_accepts_matching_checksums(tmp_path: Path) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    (outputs_root / "summary.txt").write_text("result\n", encoding="utf-8")
    bundle_root = tmp_path / "bundle"
    bundle_directory([inputs_root], [outputs_root], bundle_root)
    report = validate_bundle(bundle_root)
    assert report.valid is True
    assert report.mismatches == []


def test_validate_bundle_detects_checksum_drift(tmp_path: Path) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    (outputs_root / "summary.txt").write_text("result\n", encoding="utf-8")
    bundle_root = tmp_path / "bundle"
    bundle_directory([inputs_root], [outputs_root], bundle_root)
    (bundle_root / "outputs" / "outputs" / "summary.txt").write_text("drift\n", encoding="utf-8")
    report = validate_bundle(bundle_root)
    assert report.valid is False
    assert report.mismatches[0].reason in {"checksum_mismatch", "size_mismatch"}


def test_cli_validate_json_output(capsys) -> None:
    exit_code = main(["validate", str(fixture("example_tree.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "validate"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["data"]["tip_count"] == 4


def test_cli_inspect_reports_zero_length_branch_count(capsys) -> None:
    exit_code = main(["inspect", str(fixture("example_tree_zero_lengths.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["zero_length_branch_count"] == 3
    assert payload["data"]["zero_length_branch_count"] == 3
    assert payload["warnings"] == ["tree contains zero-length branches"]


def test_cli_diagnose_json_output(capsys) -> None:
    exit_code = main(["diagnose", str(fixture("example_tree.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "diagnose"
    assert payload["metrics"]["cherry_count"] == 2
    assert payload["data"]["inspection"]["imbalance_summary"] == "balanced"


def test_cli_diagnose_distances_writes_tsv(tmp_path: Path, capsys) -> None:
    output = tmp_path / "distances.tsv"
    exit_code = main(["diagnose", "distances", str(fixture("example_tree.nwk")), "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["outputs"] == [str(output)]
    assert output.read_text(encoding="utf-8") == (
        "tip\tdistance\n"
        "A\t0.3\n"
        "B\t0.3\n"
        "C\t0.3\n"
        "D\t0.3\n"
    )


def test_cli_diagnose_ultrametric_reports_tolerance_and_deviation(capsys) -> None:
    exit_code = main(["diagnose", "ultrametric", str(fixture("example_tree_ladderized.nwk")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["data"]["ultrametric"] is False
    assert payload["metrics"]["max_deviation"] == 0.2


def test_cli_validate_writes_run_manifest(tmp_path: Path, capsys) -> None:
    manifest = tmp_path / "validate.manifest.json"
    exit_code = main(["validate", str(fixture("example_tree.nwk")), "--json", "--manifest", str(manifest)])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert payload["outputs"] == [str(manifest)]
    assert manifest_payload["command"] == "validate"
    assert manifest_payload["arguments"] == [
        "validate",
        str(fixture("example_tree.nwk")),
        "--json",
        "--manifest",
        str(manifest),
    ]
    assert manifest_payload["input_checksums"][str(fixture("example_tree.nwk"))]


def test_cli_normalize_includes_manifest_in_output_list(tmp_path: Path, capsys) -> None:
    output = tmp_path / "normalized.nwk"
    manifest = tmp_path / "normalize.manifest.json"
    exit_code = main(
        [
            "normalize",
            str(fixture("example_tree.nex")),
            "--format",
            "nexus",
            "--out",
            str(output),
            "--json",
            "--manifest",
            str(manifest),
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["outputs"] == [str(output), str(manifest)]


def test_cli_commands_json_lists_registered_taxonomy(capsys) -> None:
    exit_code = main(["commands", "--format", "json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    command_names = [item["name"] for item in payload["data"]["commands"]]
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert command_names == [
        "env",
        "metadata",
        "traits",
        "prune",
        "alignment",
        "inspect",
        "validate",
        "normalize",
        "normalize-taxa",
        "compare",
        "annotate",
        "diagnose",
        "render",
        "report",
        "demo",
        "evidence",
        "adapter",
    ]


def test_run_capability_demo_creates_expected_artifacts(tmp_path: Path) -> None:
    result = run_capability_demo(tmp_path / "demo")
    assert result.tree_report.exists()
    assert result.dataset_report.exists()
    assert result.phylo_inputs_report.exists()
    assert result.comparison_report.exists()
    assert result.evidence_bundle.exists()
    assert result.capability_summary.exists()


def test_cli_evidence_bundle_and_validate_json_output(tmp_path: Path, capsys) -> None:
    inputs_root = tmp_path / "inputs"
    outputs_root = tmp_path / "outputs"
    bundle_root = tmp_path / "bundle"
    inputs_root.mkdir()
    outputs_root.mkdir()
    (inputs_root / "artifact.txt").write_text("evidence\n", encoding="utf-8")
    (outputs_root / "summary.txt").write_text("result\n", encoding="utf-8")

    exit_code = main(
        [
            "evidence",
            "bundle",
            "--inputs",
            str(inputs_root),
            "--outputs",
            str(outputs_root),
            "--out",
            str(bundle_root),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    bundle_payload = json.loads(captured.out)
    assert exit_code == 0
    assert bundle_payload["status"] == "ok"
    assert bundle_payload["metrics"]["input_file_count"] == 1
    assert bundle_payload["metrics"]["output_file_count"] == 1

    exit_code = main(["evidence", "validate", str(bundle_root), "--json"])
    captured = capsys.readouterr()
    validate_payload = json.loads(captured.out)
    assert exit_code == 0
    assert validate_payload["status"] == "ok"
    assert validate_payload["data"]["valid"] is True


def test_cli_demo_run_json_output_reports_generated_artifacts(tmp_path: Path, capsys) -> None:
    output = tmp_path / "demo"
    exit_code = main(["demo", "run", "--out", str(output), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 6
    assert payload["data"]["tree_report"] == str(output / "reports" / "tree-report.html")
    assert payload["data"]["evidence_bundle"] == str(output / "evidence-pack")


def test_cli_report_json_output_uses_result_envelope(tmp_path: Path, capsys) -> None:
    output = tmp_path / "report.html"
    exit_code = main(
        [
            "report",
            "tree",
            str(fixture("example_tree.nwk")),
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
    assert payload["data"]["report_kind"] == "tree"


def test_cli_report_dataset_json_output_uses_dataset_contract(tmp_path: Path, capsys) -> None:
    output = tmp_path / "dataset-report.html"
    exit_code = main(
        [
            "report",
            "dataset",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--traits",
            str(fixture("example_traits.tsv")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["data"]["report_kind"] == "dataset"
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_report_phylo_inputs_json_output_uses_alignment_contract(tmp_path: Path, capsys) -> None:
    output = tmp_path / "phylo-inputs-report.html"
    exit_code = main(
        [
            "report",
            "phylo-inputs",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--alignment",
            str(fixture("example_alignment.fasta")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["command"] == "report"
    assert payload["data"]["report_kind"] == "phylo-inputs"
    assert payload["metrics"]["alignment_length"] == 8
    assert payload["metrics"]["linked_taxa"] == 4


def test_cli_adapter_returns_typed_engine_error(capsys) -> None:
    exit_code = main(["adapter", "iqtree", "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == EngineUnavailableError.code
