from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.influence import write_taxon_influence_table


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_write_taxon_influence_table_writes_ranked_leave_one_out_rows(
    tmp_path: Path,
) -> None:
    output = tmp_path / "taxon-influence.tsv"

    write_taxon_influence_table(
        output,
        fixture("example_tree_taxon_influence_left.nwk"),
        fixture("example_tree_taxon_influence_right.nwk"),
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith(
        "influence_rank\ttaxon\tretained_taxa\tretained_taxon_count"
    )
    columns = lines[1].split("\t")
    assert columns[0:4] == ["1", "C", "A|B|D|E", "4"]
    assert columns[4] == "false"
    assert columns[5] == "true"
    assert columns[10:13] == ["2", "0", "-2"]
    assert columns[-3:] == ["true", "true", "3.333333"]
