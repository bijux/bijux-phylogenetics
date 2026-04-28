from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics
from bijux_phylogenetics.cli import main
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_clade_sets,
    compare_support_values,
    compare_tree_paths,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
    write_tree_comparison_table,
)
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.core.alignment import AlignmentSummary
from bijux_phylogenetics.core.dataset import summarize_dataset_readiness
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.core.metadata import inspect_metadata_table, join_table_to_taxa
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_alignment_to_tree,
    prune_tree_to_alignment,
    prune_tree_to_requested_taxa,
    prune_tree_to_taxa,
)
from bijux_phylogenetics.core.topology import (
    collapse_branches_below_length,
    extract_named_clade,
    ladderize_tree,
    sort_tree_tips_alphabetically,
)
from bijux_phylogenetics.core.taxonomy import inspect_tree_taxa_safety, normalize_tree_taxa, write_taxon_mapping
from bijux_phylogenetics.core.traits import (
    detect_missing_trait_values,
    detect_unusable_trait_columns,
    link_tree_to_traits,
    prune_traits_to_tree,
    validate_traits_table,
)
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
    build_alignment_quality_report,
    compute_pairwise_sequence_identity_matrix,
    compute_amino_acid_composition,
    compute_nucleotide_composition,
    detect_composition_outlier_sequences,
    detect_identical_duplicate_sequences,
    detect_invalid_alignment_characters,
    detect_near_duplicate_sequences,
    infer_alignment_alphabet,
    detect_sequences_with_excessive_missing_data,
    detect_sites_with_excessive_missing_data,
    inspect_coding_alignment,
    remove_all_gap_columns,
    remove_all_missing_columns,
    remove_sequences_above_missingness_threshold,
    translate_coding_alignment,
    trim_alignment,
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
    assert tree.branch_lengths() == [0.1, 0.2, 0.3]
    assert tree.terminal_branch_lengths() == [
        ("Homo sapiens", 0.1),
        ("NCBI|123/45", 0.2),
        ("A.B-1", 0.3),
    ]
    assert dumps_newick(tree) == "(A.B-1:0.3,'Homo sapiens':0.1,'NCBI|123/45':0.2);"


def test_normalize_tree_taxa_reports_rename_mapping() -> None:
    tree = loads_newick("('Homo sapiens':0.1,'Mus musculus':0.2,A:0.3);")
    normalized_tree, report = normalize_tree_taxa(tree, policy="spaces-to-underscores")
    assert normalized_tree.tip_names == ["Homo_sapiens", "Mus_musculus", "A"]
    assert [(rename.raw_label, rename.normalized_label) for rename in report.renamed_taxa] == [
        ("Homo sapiens", "Homo_sapiens"),
        ("Mus musculus", "Mus_musculus"),
    ]


def test_taxon_safety_reports_unsafe_labels_and_normalization_collisions(tmp_path: Path) -> None:
    tree = loads_newick(
        "('Homo sapiens':0.1,Homo_sapiens:0.2,'NCBI/123':0.3,'Quoted''Name':0.4,A:0.5);"
    )
    report = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    assert [(entry.raw_label, entry.normalized_label, entry.reasons) for entry in report.unsafe_taxa] == [
        ("Homo sapiens", "Homo_sapiens", ["contains whitespace", "collides with another label after normalization"]),
        ("Homo_sapiens", "Homo_sapiens", ["collides with another label after normalization"]),
        ("NCBI/123", "NCBI/123", ["contains slash characters"]),
        ("Quoted'Name", "Quoted'Name", ["contains quote characters"]),
    ]
    assert [(entry.normalized_label, entry.raw_labels) for entry in report.collisions] == [
        ("Homo_sapiens", ["Homo sapiens", "Homo_sapiens"])
    ]

    mapping_path = tmp_path / "taxon-mapping.tsv"
    write_taxon_mapping(mapping_path, normalize_tree_taxa(tree, policy="spaces-to-underscores")[1].renamed_taxa)
    assert mapping_path.read_text(encoding="utf-8") == (
        "raw_label\tnormalized_label\n"
        "Homo sapiens\tHomo_sapiens\n"
    )


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


def test_join_table_to_taxa_returns_tip_by_tip_metadata_rows() -> None:
    report = join_table_to_taxa(["A", "B", "Z"], fixture("example_metadata.tsv"))
    assert [(row.taxon, row.matched, row.values.get("species", "")) for row in report.joined_rows] == [
        ("A", True, "Alpha species"),
        ("B", True, "Beta species"),
        ("Z", False, ""),
    ]
    assert report.missing_from_metadata == ["Z"]
    assert report.extra_metadata_taxa == ["C", "D"]


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


def test_detect_missing_trait_values_reports_taxon_and_column() -> None:
    report = detect_missing_trait_values(fixture("example_traits_validate.tsv"))
    assert [(item.taxon, item.trait) for item in report.missing_values] == [("C", "status")]


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


def test_prune_traits_to_tree_keeps_tree_order_for_overlapping_taxa() -> None:
    rows, report = prune_traits_to_tree(
        fixture("example_tree.nwk"),
        fixture("example_traits.tsv"),
    )
    assert [row["taxon"] for row in rows] == ["A", "B", "C"]
    assert report.original_row_count == 4
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["E"]


def test_traits_prune_cli_writes_pruned_table_in_tree_order(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "pruned-traits.tsv"
    exit_code = main(
        [
            "traits",
            "prune",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "taxon\tvalue\n"
        "A\t1.2\n"
        "B\t1.4\n"
        "C\t1.8\n"
    )
    assert payload["data"]["kept_taxa"] == ["A", "B", "C"]
    assert payload["data"]["removed_taxa"] == ["E"]


def test_traits_missing_cli_reports_taxon_and_column(capsys) -> None:
    exit_code = main(["traits", "missing", str(fixture("example_traits_validate.tsv")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["missing_value_count"] == 1
    assert payload["data"]["missing_values"] == [{"taxon": "C", "trait": "status"}]


def test_dataset_readiness_reports_ready_comparative_inputs() -> None:
    report = summarize_dataset_readiness(
        fixture("example_tree.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits_validate.tsv"),
    )
    assert report.ready_for_comparative_analysis is True
    assert report.analysis_taxa == ["A", "B", "C", "D"]
    assert report.missing_metadata_taxa == []
    assert report.missing_trait_taxa == []
    assert report.metadata_only_taxa == []
    assert report.trait_only_taxa == []
    assert report.unusable_trait_columns == []
    assert report.blockers == []
    assert report.warnings == []


def test_dataset_readiness_reports_linkage_blockers() -> None:
    report = summarize_dataset_readiness(
        fixture("example_tree.nwk"),
        fixture("example_metadata.tsv"),
        fixture("example_traits.tsv"),
    )
    assert report.ready_for_comparative_analysis is False
    assert report.analysis_taxa == ["A", "B", "C"]
    assert report.missing_trait_taxa == ["D"]
    assert report.trait_only_taxa == ["E"]
    assert report.blockers == ["trait table is missing one or more tree taxa"]
    assert report.warnings == ["trait table contains taxa absent from the tree"]


def test_prune_tree_to_taxa_writes_expected_tip_set() -> None:
    tree, report = prune_tree_to_taxa(fixture("example_tree.nwk"), fixture("example_traits.tsv"))
    assert tree.tip_names == ["A", "B", "C"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1):0.2,C:0.3);"
    assert report.kept_taxa == ["A", "B", "C"]
    assert report.removed_taxa == ["D"]


def test_prune_tree_to_requested_taxa_reports_absent_requests() -> None:
    tree, report = prune_tree_to_requested_taxa(
        fixture("example_tree.nwk"),
        ["A", "C", "Z"],
    )
    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
    assert report.requested_taxa == ["A", "C", "Z"]
    assert report.kept_taxa == ["A", "C"]
    assert report.removed_taxa == ["B", "D"]
    assert report.absent_requested_taxa == ["Z"]


def test_drop_tree_taxa_excludes_exact_requested_tips() -> None:
    tree, report = drop_tree_taxa(
        fixture("example_tree.nwk"),
        ["B", "D", "Z"],
    )
    assert tree.tip_names == ["A", "C"]
    assert dumps_newick(tree) == "(A:0.3,C:0.3);"
    assert report.requested_taxa == ["B", "D", "Z"]
    assert report.kept_taxa == ["A", "C"]
    assert report.removed_taxa == ["B", "D"]
    assert report.absent_requested_taxa == ["Z"]


def test_prune_cli_accepts_explicit_taxon_keep_lists(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "pruned-tree.nwk"
    pruned_taxa_path = tmp_path / "removed.tsv"
    exit_code = main(
        [
            "prune",
            str(fixture("example_tree.nwk")),
            "--taxa",
            "A",
            "C",
            "Z",
            "--out",
            str(output_path),
            "--pruned-taxa-out",
            str(pruned_taxa_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(A:0.3,C:0.3);\n"
    assert pruned_taxa_path.read_text(encoding="utf-8") == "taxon\nB\nD\n"
    assert payload["data"]["absent_requested_taxa"] == ["Z"]
    assert payload["data"]["kept_taxa"] == ["A", "C"]


def test_prune_cli_accepts_explicit_taxon_exclusion_lists(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "pruned-tree.nwk"
    pruned_taxa_path = tmp_path / "removed.tsv"
    exit_code = main(
        [
            "prune",
            str(fixture("example_tree.nwk")),
            "--exclude-taxa",
            "B",
            "D",
            "Z",
            "--out",
            str(output_path),
            "--pruned-taxa-out",
            str(pruned_taxa_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == "(A:0.3,C:0.3);\n"
    assert pruned_taxa_path.read_text(encoding="utf-8") == "taxon\nB\nD\n"
    assert payload["data"]["absent_requested_taxa"] == ["Z"]
    assert payload["data"]["removed_taxa"] == ["B", "D"]


def test_extract_named_clade_returns_exact_descendant_subtree() -> None:
    tree, report = extract_named_clade(
        fixture("example_tree_named_clades.nwk"),
        clade_name="Mammals",
    )
    assert tree.tip_names == ["A", "B"]
    assert dumps_newick(tree) == "(A:0.1,B:0.1)Mammals;"
    assert report.clade_name == "Mammals"
    assert report.tip_count == 2
    assert report.taxa == ["A", "B"]


def test_collapse_branches_below_length_turns_short_internal_edges_into_polytomies() -> None:
    tree, report = collapse_branches_below_length(
        fixture("example_tree_collapse_threshold.nwk"),
        threshold=0.05,
    )
    assert tree.tip_names == ["A", "B", "C", "D"]
    assert dumps_newick(tree) == "((A:0.1,B:0.1,C:0.2):0.3,D:0.4);"
    assert report.threshold == 0.05
    assert report.collapsed_clades == ["A|B"]


def test_ladderize_tree_orders_larger_subtrees_first() -> None:
    tree, report = ladderize_tree(fixture("example_tree_ordering.nwk"))
    assert tree.tip_names == ["X", "Y", "Z", "A", "B"]
    assert [len(child.children) if child.children else 1 for child in tree.root.children] == [3, 2]
    assert report.strategy == "ladderize"
    assert report.tip_order == ["X", "Y", "Z", "A", "B"]


def test_sort_tree_tips_alphabetically_preserves_topology_with_stable_tip_order() -> None:
    tree, report = sort_tree_tips_alphabetically(fixture("example_tree_ordering.nwk"))
    assert tree.tip_names == ["A", "B", "X", "Y", "Z"]
    assert [len(child.children) if child.children else 1 for child in tree.root.children] == [2, 3]
    assert report.strategy == "alphabetical"
    assert report.tip_order == ["A", "B", "X", "Y", "Z"]


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
    assert report.inferred_alphabet == "dna"
    assert report.nucleotide_composition == {"A": 0.3125, "C": 0.25, "G": 0.25, "T": 0.1875}
    assert report.whole_alignment_gc_content == 0.5
    assert [(row.identifier, row.missing_fraction) for row in report.per_sequence_missingness] == [
        ("A", 0.0),
        ("B", 0.0),
        ("C", 0.0),
        ("D", 0.0),
    ]
    assert [(row.identifier, row.gc_fraction) for row in report.per_sequence_gc_content] == [
        ("A", 0.5),
        ("B", 0.375),
        ("C", 0.625),
        ("D", 0.5),
    ]


def test_alignment_detects_sequence_alphabet_types() -> None:
    dna_records = load_fasta_alignment(fixture("example_alignment.fasta"))
    protein_records = load_fasta_alignment(fixture("example_alignment_protein.fasta"))
    assert infer_alignment_alphabet(dna_records) == "dna"
    assert infer_alignment_alphabet(protein_records) == "protein"


def test_alignment_detects_invalid_characters_for_declared_alphabet() -> None:
    invalid = detect_invalid_alignment_characters(
        fixture("example_alignment_invalid_dna.fasta"),
        alphabet="dna",
    )
    assert [(row.identifier, row.position, row.character) for row in invalid] == [("A", 5, "Z")]


def test_alignment_reports_nucleotide_and_amino_acid_composition() -> None:
    dna_records = load_fasta_alignment(fixture("example_alignment.fasta"))
    protein_records = load_fasta_alignment(fixture("example_alignment_protein.fasta"))
    assert compute_nucleotide_composition(dna_records, alphabet="dna") == {
        "A": 0.3125,
        "C": 0.25,
        "G": 0.25,
        "T": 0.1875,
    }
    assert compute_amino_acid_composition(protein_records, alphabet="protein") == {
        "F": 0.083333333333333,
        "I": 0.083333333333333,
        "K": 0.083333333333333,
        "L": 0.125,
        "M": 0.25,
        "R": 0.041666666666667,
        "T": 0.125,
        "V": 0.041666666666667,
        "W": 0.125,
        "Y": 0.041666666666667,
    }


def test_alignment_detects_gc_composition_outliers() -> None:
    outliers = detect_composition_outlier_sequences(
        fixture("example_alignment_gc_outlier.fasta"),
        deviation_threshold=0.2,
    )
    assert [(row.identifier, row.deviation) for row in outliers] == [("C", 1.0)]


def test_alignment_detects_identical_and_near_duplicate_sequences() -> None:
    duplicates = detect_identical_duplicate_sequences(fixture("example_alignment_duplicates.fasta"))
    near_duplicates = detect_near_duplicate_sequences(
        fixture("example_alignment_duplicates.fasta"),
        identity_threshold=0.875,
    )
    assert [(group.identifiers, group.sequence) for group in duplicates] == [
        (["A", "B"], "ACTGACTG")
    ]
    assert [(pair.left_identifier, pair.right_identifier, pair.identity, pair.comparable_sites) for pair in near_duplicates] == [
        ("A", "C", 0.875, 8),
        ("A", "D", 0.875, 8),
        ("B", "C", 0.875, 8),
        ("B", "D", 0.875, 8),
        ("C", "D", 0.875, 8),
    ]


def test_alignment_quality_report_collects_composition_duplicates_and_warnings() -> None:
    report = build_alignment_quality_report(fixture("example_alignment_duplicates.fasta"))
    assert report.sequence_count == 4
    assert report.alignment_length == 8
    assert report.variable_site_count == 1
    assert report.inferred_alphabet == "dna"
    assert report.invalid_characters == []
    assert report.duplicate_sequence_groups
    assert report.near_duplicate_pairs == []
    assert report.warnings == ["alignment contains identical duplicate sequences"]


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


def test_alignment_removes_all_gap_columns() -> None:
    records, report = remove_all_gap_columns(fixture("example_alignment_site_missingness.fasta"))
    assert [record.sequence for record in records] == ["AN?T", "CN?N", "GN?A", "TN?N"]
    assert report.original_alignment_length == 5
    assert report.trimmed_alignment_length == 4
    assert [(column.position, column.reason) for column in report.removed_columns] == [(2, "all-gap")]


def test_alignment_removes_all_missing_columns() -> None:
    records, report = remove_all_missing_columns(fixture("example_alignment_site_missingness.fasta"))
    assert [record.sequence for record in records] == ["A-T", "C-N", "G-A", "T-N"]
    assert report.original_alignment_length == 5
    assert report.trimmed_alignment_length == 3
    assert [(column.position, column.reason) for column in report.removed_columns] == [
        (3, "all-missing"),
        (4, "all-missing"),
    ]


def test_alignment_removes_sequences_above_missingness_threshold() -> None:
    records, report = remove_sequences_above_missingness_threshold(
        fixture("example_alignment_missingness.fasta"),
        threshold=0.3,
    )
    assert [record.identifier for record in records] == ["C"]
    assert report.original_sequence_count == 3
    assert report.trimmed_sequence_count == 1
    assert [(row.identifier, row.missing_fraction, row.reason) for row in report.removed_sequences] == [
        ("A", 2 / 6, "missingness-threshold"),
        ("B", 2 / 6, "missingness-threshold"),
    ]


def test_alignment_trimming_report_combines_sequence_and_column_transforms() -> None:
    records, report = trim_alignment(
        fixture("example_alignment_trim.fasta"),
        sequence_missingness_threshold=0.3,
    )
    assert [(record.identifier, record.sequence) for record in records] == [("B", "ACTG")]
    assert report.original_sequence_count == 3
    assert report.trimmed_sequence_count == 1
    assert report.original_alignment_length == 6
    assert report.trimmed_alignment_length == 4
    assert [(column.position, column.reason) for column in report.removed_columns] == [
        (3, "all-gap"),
        (4, "all-missing"),
    ]
    assert [(row.identifier, row.reason) for row in report.removed_sequences] == [
        ("A", "missingness-threshold"),
        ("C", "missingness-threshold"),
    ]


def test_alignment_identity_matrix_reports_pairs_and_comparable_sites() -> None:
    report = compute_pairwise_sequence_identity_matrix(fixture("example_alignment_duplicates.fasta"))
    assert report.identifiers == ["A", "B", "C", "D"]
    assert [(pair.left_identifier, pair.right_identifier, pair.identity, pair.comparable_sites) for pair in report.pairs] == [
        ("A", "A", 1.0, 8),
        ("A", "B", 1.0, 8),
        ("A", "C", 0.875, 8),
        ("A", "D", 0.875, 8),
        ("B", "B", 1.0, 8),
        ("B", "C", 0.875, 8),
        ("B", "D", 0.875, 8),
        ("C", "C", 1.0, 8),
        ("C", "D", 0.875, 8),
        ("D", "D", 1.0, 8),
    ]


def test_coding_alignment_reports_frameshift_like_sequences_and_stop_codons() -> None:
    diagnostics = inspect_coding_alignment(fixture("example_alignment_coding.fasta"))
    assert diagnostics.sequence_count == 4
    assert [(row.identifier, row.comparable_length, row.remainder) for row in diagnostics.frameshift_like_sequences] == [
        ("C", 8, 2)
    ]
    assert [(row.identifier, row.codon_index, row.nucleotide_start, row.codon, row.terminal) for row in diagnostics.stop_codons] == [
        ("A", 3, 7, "TAA", True),
        ("D", 2, 4, "TAG", False),
    ]


def test_translate_coding_alignment_emits_amino_acid_records() -> None:
    records, report = translate_coding_alignment(fixture("example_alignment_coding.fasta"))
    assert [(record.identifier, record.sequence) for record in records] == [
        ("A", "ME*"),
        ("B", "MEW"),
        ("C", "MXW"),
        ("D", "M*W"),
    ]
    assert report.translated_sequence_count == 4
    assert report.source_alignment_length == 9
    assert report.translated_alignment_length == 3
    assert report.stop_codon_count == 2
    assert report.frameshift_like_sequence_count == 1


def test_cli_alignment_trim_writes_trimmed_fasta_and_report(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "trimmed.fasta"
    exit_code = main(
        [
            "alignment",
            "trim",
            str(fixture("example_alignment_trim.fasta")),
            "--out",
            str(output_path),
            "--sequence-missingness-threshold",
            "0.3",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == ">B\nACTG\n"
    assert payload["metrics"]["removed_column_count"] == 2
    assert payload["metrics"]["removed_sequence_count"] == 2


def test_cli_alignment_identity_matrix_writes_tsv(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "identity.tsv"
    exit_code = main(
        [
            "alignment",
            "identity-matrix",
            str(fixture("example_alignment_duplicates.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tidentity\tcomparable_sites",
        "A\tA\t1\t8",
        "A\tB\t1\t8",
        "A\tC\t0.875\t8",
    ]
    assert payload["metrics"]["pair_count"] == 10


def test_cli_alignment_coding_reports_frameshifts_and_stops(capsys) -> None:
    exit_code = main(["alignment", "coding", str(fixture("example_alignment_coding.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["frameshift_like_sequence_count"] == 1
    assert payload["metrics"]["stop_codon_count"] == 2


def test_cli_alignment_translate_writes_amino_acid_alignment(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "translated.fasta"
    exit_code = main(
        [
            "alignment",
            "translate",
            str(fixture("example_alignment_coding.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        ">A\nME*\n"
        ">B\nMEW\n"
        ">C\nMXW\n"
        ">D\nM*W\n"
    )
    assert payload["metrics"]["translated_sequence_count"] == 4
    assert payload["metrics"]["stop_codon_count"] == 2


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
    assert report.branch_length_summary is not None
    assert (
        report.branch_length_summary.count,
        report.branch_length_summary.minimum,
        report.branch_length_summary.maximum,
        report.branch_length_summary.mean,
        report.branch_length_summary.median,
        report.branch_length_summary.first_quartile,
        report.branch_length_summary.third_quartile,
    ) == (6, 0.1, 0.2, 0.15, 0.15, 0.1, 0.2)
    assert report.tree_diameter == 0.6
    assert report.zero_length_branch_count == 0
    assert report.max_depth == 2
    assert report.mean_depth == 2.0
    assert report.colless_imbalance_index == 0.0
    assert report.normalized_colless_imbalance == 0.0
    assert report.sackin_imbalance_index == 8
    assert report.unusually_imbalanced is False
    assert report.long_branch_taxa == []
    assert report.star_like is False
    assert report.comb_like is False
    assert report.tree_quality_score == 100.0
    assert report.tree_quality_warnings == []
    assert report.imbalance_summary == "balanced"
    assert report.cherry_count == 2
    assert report.warnings == []
    assert report.taxa == ["A", "B", "C", "D"]


def test_inspect_tree_path_distinguishes_ladderized_shape() -> None:
    report = inspect_tree_path(fixture("example_tree_ladderized.nwk"))
    assert report.tree_diameter == 0.4
    assert report.max_depth == 3
    assert report.mean_depth == 2.25
    assert report.colless_imbalance_index == 3.0
    assert report.normalized_colless_imbalance == 1.0
    assert report.sackin_imbalance_index == 9
    assert report.unusually_imbalanced is True
    assert report.long_branch_taxa == []
    assert report.star_like is False
    assert report.comb_like is True
    assert report.tree_quality_score == 75.0
    assert [warning.code for warning in report.tree_quality_warnings] == [
        "unusually_imbalanced",
        "comb_like",
    ]
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
    assert report.colless_imbalance_index is None
    assert report.normalized_colless_imbalance is None
    assert report.unusually_imbalanced is None
    assert report.long_branch_taxa == []
    assert report.star_like is False
    assert report.comb_like is False
    assert report.tree_quality_score == 90.0
    assert [warning.code for warning in report.tree_quality_warnings] == ["polytomies"]
    assert report.polytomy_count == 1
    assert report.polytomy_nodes == ["A|B|C"]


def test_inspect_tree_path_detects_long_branch_taxa() -> None:
    report = inspect_tree_path(fixture("example_tree_long_branch.nwk"))
    assert report.long_branch_taxa == ["A"]
    assert report.star_like is False
    assert report.tree_quality_score == 85.0
    assert [warning.code for warning in report.tree_quality_warnings] == ["long_branches"]


def test_inspect_tree_path_detects_star_like_tree() -> None:
    report = inspect_tree_path(fixture("example_tree_star.nwk"))
    assert report.star_like is True
    assert report.long_branch_taxa == []
    assert report.tree_quality_score == 80.0
    assert [warning.code for warning in report.tree_quality_warnings] == ["polytomies", "star_like"]


def test_inspect_tree_path_classifies_branch_length_completeness() -> None:
    complete = inspect_tree_path(fixture("example_tree.nwk"))
    partial = inspect_tree_path(fixture("example_tree_partial_lengths.nwk"))
    absent = inspect_tree_path(fixture("example_tree_no_lengths.nwk"))
    assert complete.branch_length_status == "complete"
    assert partial.branch_length_status == "partial"
    assert absent.branch_length_status == "absent"
    assert complete.branch_length_summary is not None
    assert partial.branch_length_summary is not None
    assert absent.branch_length_summary is None
    assert complete.tree_diameter == 0.6
    assert partial.tree_diameter is None
    assert absent.tree_diameter is None
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
    assert payload["data"]["tree_diameter"] == 1.6
    assert payload["data"]["mean_depth"] == 2.0
    assert payload["data"]["colless_imbalance_index"] == 0.0
    assert payload["data"]["sackin_imbalance_index"] == 8
    assert payload["data"]["imbalance_summary"] == "balanced"
    assert payload["data"]["tree_quality_score"] == 100.0
    assert payload["metrics"]["cherry_count"] == 2
    assert payload["metrics"]["tree_diameter"] == 1.6
    assert payload["metrics"]["tree_quality_score"] == 100.0
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


def test_prune_trees_to_shared_taxa_keeps_identical_tip_sets() -> None:
    left, right, report = prune_trees_to_shared_taxa(
        fixture("example_tree.nwk"),
        fixture("example_tree_overlap.nwk"),
    )
    assert left.tip_names == ["A", "B", "C"]
    assert right.tip_names == ["A", "B", "C"]
    assert report.shared_taxa == ["A", "B", "C"]
    assert report.left_only_taxa == ["D"]
    assert report.right_only_taxa == ["E"]


def test_compare_tree_paths_reports_identical_topology_boolean() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree.nwk"))
    assert report.topology_equal is True
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is False
    assert report.same_topology_different_branch_lengths is False
    assert report.robinson_foulds_distance == 0


def test_compare_tree_paths_reports_different_topology_boolean() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_topology_diff.nwk"))
    assert report.topology_equal is False
    assert report.same_unrooted_topology is False
    assert report.same_taxa_different_rooting is False


def test_compare_clade_sets_reports_shared_and_unique_clades() -> None:
    report = compare_clade_sets(fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert report.shared_clades == ["A|B"]
    assert report.left_only_clades == ["C|D"]
    assert report.right_only_clades == ["A|B|C"]


def test_compare_tree_paths_detects_same_topology_with_different_branch_lengths() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_branch_lengths_right.nwk"))
    assert report.topology_equal is True
    assert report.same_topology_different_branch_lengths is True


def test_compare_tree_paths_detects_same_taxa_with_different_rooting() -> None:
    report = compare_tree_paths(fixture("example_tree.nwk"), fixture("example_tree_rooting_diff.nwk"))
    assert report.topology_equal is False
    assert report.same_unrooted_topology is True
    assert report.same_taxa_different_rooting is True


def test_detect_clade_changes_reports_lost_and_gained_sets() -> None:
    report = detect_clade_changes(fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert report.lost_clades == ["C|D"]
    assert report.gained_clades == ["A|B|C"]


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


def test_write_tree_comparison_table_writes_one_row_per_compared_split(tmp_path: Path) -> None:
    output = tmp_path / "comparison.tsv"
    write_tree_comparison_table(output, fixture("example_tree.nwk"), fixture("example_tree_alt.nwk"))
    assert output.read_text(encoding="utf-8") == (
        "split_id\tcomparison_status\tshared_clade\tleft_support\tright_support\tleft_length\tright_length\tlength_delta\tlength_ratio\n"
        "A|B\tshared\ttrue\t\t\t0.2\t0.1\t-0.1\t0.5\n"
        "A|B|C\tright_only\tfalse\t\t\t\t\t\t\n"
        "C|D\tleft_only\tfalse\t\t\t\t\t\t\n"
    )


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
    assert "clade-comparison" in html
    assert "clade-changes" in html
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
    assert [(row.taxon, row.matched) for row in report.joined_rows] == [
        ("A", True),
        ("B", True),
        ("C", True),
        ("D", False),
    ]


def test_cli_metadata_inspect_json_output(capsys) -> None:
    exit_code = main(["metadata", "inspect", str(fixture("example_metadata.tsv")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["command"] == "metadata"
    assert payload["data"]["taxon_column"] == "taxon"
    assert payload["metrics"]["taxon_count"] == 4


def test_cli_annotate_can_write_joined_tip_rows(tmp_path: Path, capsys) -> None:
    report_path = tmp_path / "annotation.json"
    joined_path = tmp_path / "annotation.tsv"
    exit_code = main(
        [
            "annotate",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_metadata.tsv")),
            "--out",
            str(report_path),
            "--joined-out",
            str(joined_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["joined_rows"][3]["taxon"] == "D"
    assert joined_path.read_text(encoding="utf-8") == (
        "taxon\tmatched\tspecies\tlocation\n"
        "A\ttrue\tAlpha species\tSweden\n"
        "B\ttrue\tBeta species\tNorway\n"
        "C\ttrue\tGamma species\tDenmark\n"
        "D\ttrue\tDelta species\tFinland\n"
    )


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
    assert payload["metrics"]["alphabet"] == "dna"
    assert payload["metrics"]["duplicate_group_count"] == 0
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


def test_cli_alignment_quality_json_output(capsys) -> None:
    exit_code = main(["alignment", "quality", str(fixture("example_alignment_duplicates.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["duplicate_group_count"] == 1
    assert payload["metrics"]["near_duplicate_count"] == 0
    assert payload["warnings"] == ["alignment contains identical duplicate sequences"]


def test_cli_alignment_composition_json_output(capsys) -> None:
    exit_code = main(["alignment", "composition", str(fixture("example_alignment.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["alphabet"] == "dna"
    assert payload["data"]["whole_alignment_gc_content"] == 0.5
    assert payload["data"]["nucleotide_composition"] == {
        "A": 0.3125,
        "C": 0.25,
        "G": 0.25,
        "T": 0.1875,
    }


def test_cli_alignment_alphabet_json_output(capsys) -> None:
    exit_code = main(["alignment", "alphabet", str(fixture("example_alignment_protein.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["alphabet"] == "protein"
    assert payload["data"]["inferred_alphabet"] == "protein"


def test_cli_alignment_gc_json_output(capsys) -> None:
    exit_code = main(["alignment", "gc", str(fixture("example_alignment.fasta")), "--json"])
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["whole_alignment_gc_content"] == 0.5
    assert payload["data"]["per_sequence_gc_content"][1] == {"gc_fraction": 0.375, "identifier": "B"}


def test_cli_alignment_invalid_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "invalid",
            str(fixture("example_alignment_invalid_dna.fasta")),
            "--alphabet",
            "dna",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["invalid_character_count"] == 1
    assert payload["data"] == [{"character": "Z", "identifier": "A", "position": 5}]


def test_cli_alignment_duplicates_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "duplicates",
            str(fixture("example_alignment_duplicates.fasta")),
            "--identity-threshold",
            "0.875",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["duplicate_group_count"] == 1
    assert payload["metrics"]["near_duplicate_count"] == 5
    assert payload["data"]["duplicate_sequence_groups"][0]["identifiers"] == ["A", "B"]


def test_cli_alignment_outliers_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "outliers",
            str(fixture("example_alignment_gc_outlier.fasta")),
            "--deviation-threshold",
            "0.2",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["composition_outlier_count"] == 1
    assert payload["data"] == [{"deviation": 1.0, "identifier": "C"}]


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


def test_cli_compare_clades_json_output(capsys) -> None:
    exit_code = main(
        ["compare", "clades", str(fixture("example_tree.nwk")), str(fixture("example_tree_alt.nwk")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["shared_clades"] == ["A|B"]
    assert payload["data"]["left_only_clades"] == ["C|D"]
    assert payload["data"]["right_only_clades"] == ["A|B|C"]


def test_cli_compare_changes_json_output(capsys) -> None:
    exit_code = main(
        ["compare", "changes", str(fixture("example_tree.nwk")), str(fixture("example_tree_alt.nwk")), "--json"]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["lost_clades"] == ["C|D"]
    assert payload["data"]["gained_clades"] == ["A|B|C"]


def test_cli_compare_prune_writes_shared_taxon_trees(tmp_path: Path, capsys) -> None:
    output_dir = tmp_path / "shared"
    exit_code = main(
        [
            "compare",
            "prune",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_overlap.nwk")),
            "--out",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["data"]["shared_taxa"] == ["A", "B", "C"]
    assert (output_dir / "left-shared.nwk").read_text(encoding="utf-8") == "((A:0.1,B:0.1):0.2,C:0.3);\n"
    assert (output_dir / "right-shared.nwk").read_text(encoding="utf-8") == "((A:0.1,B:0.1):0.2,C:0.3);\n"


def test_cli_compare_table_writes_tsv_output(tmp_path: Path, capsys) -> None:
    output = tmp_path / "comparison.tsv"
    exit_code = main(
        [
            "compare",
            "table",
            str(fixture("example_tree.nwk")),
            str(fixture("example_tree_alt.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["table_rows"] == 3
    assert output.read_text(encoding="utf-8").startswith("split_id\tcomparison_status\tshared_clade\t")
    assert payload["data"]["table_path"] == str(output)


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
    assert result.dataset_readiness is not None
    assert result.machine_manifest["sections"] == [
        "tree-validation",
        "tree-inspection",
        "metadata-linkage",
        "traits-linkage",
        "trait-missing-values",
        "dataset-readiness",
    ]
    assert result.trait_missing_values is not None
    assert result.trait_missing_values.missing_values == []
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
    assert result.alignment_quality is not None
    assert result.alignment_linkage is not None
    assert result.machine_manifest["sections"] == [
        "tree-validation",
        "tree-inspection",
        "alignment-summary",
        "alignment-quality",
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
    assert payload["metrics"]["tree_diameter"] == 0.6
    assert payload["metrics"]["tree_quality_score"] == 100.0
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
