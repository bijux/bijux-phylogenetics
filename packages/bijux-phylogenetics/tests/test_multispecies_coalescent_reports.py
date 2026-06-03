from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.simulation import (
    simulate_multispecies_coalescent_gene_tree,
    write_multispecies_coalescent_branch_table,
    write_multispecies_coalescent_event_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_multispecies_coalescent_ledgers_emit_reviewable_rows(
    tmp_path: Path,
) -> None:
    _gene_tree, report = simulate_multispecies_coalescent_gene_tree(
        fixture("trees", "multispecies_coalescent_species_tree_3_taxa.nwk"),
        sample_count_table_path=fixture(
            "metadata", "multispecies_coalescent_sample_counts_3_taxa.tsv"
        ),
        population_size_table_path=fixture(
            "metadata", "multispecies_coalescent_population_sizes_3_taxa.tsv"
        ),
        seed=7,
    )

    event_path = write_multispecies_coalescent_event_table(
        tmp_path / "multispecies-coalescent-events.tsv",
        report,
    )
    branch_path = write_multispecies_coalescent_branch_table(
        tmp_path / "multispecies-coalescent-branches.tsv",
        report,
    )

    event_text = event_path.read_text(encoding="utf-8")
    branch_text = branch_path.read_text(encoding="utf-8")

    assert event_path.is_file()
    assert branch_path.is_file()
    assert event_text.startswith(
        "event_index\tspecies_branch\tbranch_role\tdescendant_species\tpopulation_size\t"
    )
    assert "\n1\tA\ttip-branch\tA\t0.05\t" in event_text
    assert "A__1|A__2" in event_text
    assert branch_text.startswith(
        "species_branch\tbranch_role\tdescendant_species\tbranch_duration\tpopulation_size\t"
    )
    assert "\nA|B\tinternal-branch\tA|B\t1\t1000000\t2\t0\t2\t1\ttrue\n" in branch_text
    assert branch_text.rstrip().endswith(
        "A|B|C\troot-population\tA|B|C\t\t0.5\t3\t2\t1\t0\tfalse"
    )
