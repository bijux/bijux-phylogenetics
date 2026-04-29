from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.pruning import prune_tree_to_taxa
from bijux_phylogenetics.core.taxon_workflows import (
    build_taxon_stability_report,
    build_taxon_workflow_loss_report,
    load_taxon_run_source,
)
from bijux_phylogenetics.core.taxonomy import (
    audit_tree_taxon_synonyms,
    inspect_tree_taxon_namespaces,
    resolve_tree_taxon_synonyms,
    write_synonym_resolution_mapping,
)
from bijux_phylogenetics.io.trees import load_tree


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


def test_audit_tree_taxon_synonyms_detects_candidates_and_ambiguity() -> None:
    tree = load_tree(fixture("example_taxonomy_tree.nwk"))
    report = audit_tree_taxon_synonyms(tree, fixture("example_taxon_synonyms.tsv"))

    assert [(row.raw_label, row.accepted_label) for row in report.candidates] == [
        ("Felis_concolor", "Puma_concolor"),
    ]
    assert [(row.raw_label, row.accepted_labels) for row in report.ambiguous_mappings] == [
        ("Jaguar", ["Panthera_onca", "Panthera_pardus"]),
    ]


def test_resolve_tree_taxon_synonyms_is_reversible_and_rejects_ambiguous_rows(tmp_path: Path) -> None:
    tree = load_tree(fixture("example_taxonomy_tree.nwk"))
    resolved_tree, report = resolve_tree_taxon_synonyms(
        tree,
        synonym_table_path=fixture("example_taxon_synonyms.tsv"),
    )

    assert "Puma_concolor" in resolved_tree.tip_names
    assert "Felis_concolor" not in resolved_tree.tip_names
    assert "Jaguar" in resolved_tree.tip_names
    assert [row.raw_label for row in report.renamed_taxa] == ["Felis_concolor"]
    assert len(report.ambiguous_mappings) == 1

    mapping_path = write_synonym_resolution_mapping(tmp_path / "synonyms.tsv", report)
    assert "accepted_label" in mapping_path.read_text(encoding="utf-8")


def test_inspect_tree_taxon_namespaces_reports_mixed_identifier_styles() -> None:
    tree = load_tree(fixture("example_taxonomy_tree.nwk"))
    report = inspect_tree_taxon_namespaces(tree)

    assert report.mixed_namespaces is True
    assert report.namespace_counts["accession_id"] == 1
    assert report.namespace_counts["species_name"] == 2
    assert report.namespace_counts["sample_id"] == 1
    assert report.namespace_counts["isolate_id"] == 1
    assert any("multiple explicit taxon namespaces" in warning for warning in report.warnings)


def test_prune_tree_to_taxa_reports_information_loss() -> None:
    _, report = prune_tree_to_taxa(
        fixture("example_taxon_workflow_tree.nwk"),
        fixture("example_taxon_workflow_metadata.csv"),
    )

    assert report.information_loss.lost_taxa_count == 1
    assert report.information_loss.lost_taxa_fraction == 0.25
    assert report.information_loss.lost_metadata_count == 1
    assert report.information_loss.lost_branch_length > 0.0


def test_build_taxon_workflow_loss_report_tracks_first_loss_stage() -> None:
    report = build_taxon_workflow_loss_report(
        fixture("example_taxon_workflow_tree.nwk"),
        fixture("example_taxon_workflow_metadata.csv"),
        fixture("example_taxon_workflow_traits.csv"),
        alignment_path=fixture("example_taxon_workflow_alignment.fasta"),
        filtered_alignment_path=fixture("example_taxon_workflow_filtered_alignment.fasta"),
        inference_tree_path=fixture("example_taxon_workflow_inference.nwk"),
        reported_taxa_path=fixture("example_taxon_workflow_reported.csv"),
    )

    rows = {row.taxon: row for row in report.rows}
    assert rows["A"].first_loss_stage is None
    assert rows["B"].first_loss_stage == "alignment_filtering"
    assert rows["C"].first_loss_stage == "trait_missingness"
    assert rows["D"].first_loss_stage == "alignment"
    assert report.loss_stage_counts == {
        "alignment": 1,
        "alignment_filtering": 1,
        "trait_missingness": 1,
    }


def test_build_taxon_stability_report_compares_named_sources() -> None:
    report = build_taxon_stability_report(
        [
            load_taxon_run_source(label="tree", path=fixture("example_taxon_workflow_tree.nwk")),
            load_taxon_run_source(label="alignment", path=fixture("example_taxon_workflow_alignment.fasta")),
            load_taxon_run_source(label="filtered", path=fixture("example_taxon_workflow_filtered_alignment.fasta")),
        ]
    )

    assert report.shared_taxa == ["A", "C"]
    assert sorted(report.unstable_taxa) == ["B", "D"]
    rows = {row.taxon: row for row in report.rows}
    assert rows["A"].retention_fraction == 1.0
    assert rows["B"].retention_fraction == 2 / 3
    assert rows["D"].retention_fraction == 1 / 3
