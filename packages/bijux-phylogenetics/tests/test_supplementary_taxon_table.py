from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.reports import write_supplementary_taxon_table

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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_write_supplementary_taxon_table_tracks_analysis_and_reporting_loss(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-taxa.tsv"

    result = write_supplementary_taxon_table(
        output_path,
        tree_path=fixture("example_taxon_workflow_tree.nwk"),
        metadata_path=fixture("example_taxon_workflow_metadata.csv"),
        traits_path=fixture("example_taxon_workflow_traits.csv"),
        alignment_path=fixture("example_taxon_workflow_alignment.fasta"),
        filtered_alignment_path=fixture(
            "example_taxon_workflow_filtered_alignment.fasta"
        ),
        inference_tree_path=fixture("example_taxon_workflow_inference.nwk"),
        reported_taxa_path=fixture("example_taxon_workflow_reported.csv"),
    )

    assert result.output_path == output_path
    assert result.row_count == 4
    assert result.analysis_included_count == 2
    assert result.analysis_excluded_count == 2
    assert result.reporting_retained_count == 1
    assert result.reporting_dropped_count == 3
    assert result.metadata_column_count == 1
    assert result.trait_column_count == 2

    rows = {row["taxon"]: row for row in read_tsv(output_path)}
    assert rows["A"]["analysis_status"] == "included"
    assert rows["A"]["reporting_status"] == "retained"
    assert rows["A"]["metadata_group"] == "g1"
    assert rows["A"]["trait_trait_alpha"] == "1"
    assert rows["B"]["analysis_status"] == "excluded"
    assert rows["B"]["analysis_exclusion_reason"] == "absent_from_metadata"
    assert rows["B"]["analysis_first_failed_surface"] == "metadata"
    assert rows["B"]["reporting_status"] == "dropped"
    assert rows["B"]["workflow_first_loss_stage"] == "alignment_filtering"
    assert "removed_during_alignment_filtering" in rows["B"]["workflow_loss_reasons"]
    assert rows["C"]["analysis_status"] == "included"
    assert rows["C"]["reporting_status"] == "dropped"
    assert rows["C"]["reporting_loss_reason"] == (
        "trait_missingness:one_or_more_trait_values_missing"
    )
    assert rows["C"]["trait_trait_beta"] == ""
    assert rows["D"]["analysis_status"] == "excluded"
    assert rows["D"]["analysis_exclusion_reason"] == "absent_from_alignment"
    assert rows["D"]["reporting_status"] == "dropped"
    assert rows["D"]["workflow_first_loss_stage"] == "alignment"


def test_write_supplementary_taxon_table_prefixes_metadata_and_trait_columns(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "tree.nwk"
    metadata_path = tmp_path / "metadata.tsv"
    traits_path = tmp_path / "traits.tsv"
    output_path = tmp_path / "supplementary-taxa.tsv"

    tree_path.write_text("(Alpha:1,Beta:1);", encoding="utf-8")
    metadata_path.write_text(
        "taxon\tstatus\ttaxonomy_id\nAlpha\tretained\t123\nBeta\texcluded\t456\n",
        encoding="utf-8",
    )
    traits_path.write_text(
        "taxon\tstatus\tmass_g\nAlpha\tlarge\t10.5\nBeta\tsmall\t9.1\n",
        encoding="utf-8",
    )

    result = write_supplementary_taxon_table(
        output_path,
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
    )

    assert "metadata_status" in result.columns
    assert "trait_status" in result.columns
    rows = {row["taxon"]: row for row in read_tsv(output_path)}
    assert rows["Alpha"]["metadata_status"] == "retained"
    assert rows["Alpha"]["trait_status"] == "large"
    assert rows["Beta"]["external_taxonomy_ids"] == "metadata.taxonomy_id=456"
