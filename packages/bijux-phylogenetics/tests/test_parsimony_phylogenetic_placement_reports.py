from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.parsimony import (
    place_parsimony_queries,
    write_parsimony_placement_alternative_table,
    write_parsimony_placement_artifacts,
    write_parsimony_placement_run_json,
    write_parsimony_placement_summary_table,
    write_parsimony_placement_tree_set,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def _placement_report():
    return place_parsimony_queries(
        fixture("placement_reference_tree_4_taxa.nwk"),
        fixture("placement_reference_matrix.tsv"),
        fixture("placement_query_matrix.tsv"),
    )


def test_write_parsimony_placement_summary_table_writes_best_rows(
    tmp_path: Path,
) -> None:
    report = _placement_report()
    output_path = tmp_path / "summary.tsv"

    write_parsimony_placement_summary_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "query_id\tcharacter_count\tbest_edge_id\tbest_child_name\tbest_descendant_taxa\tbest_total_steps\tbest_additional_steps\tbest_total_weighted_score\tbest_additional_weighted_score\tcandidate_placement_count\tequally_best_placement_count"
    )
    assert lines[1].startswith("Q_C\t2\t")
    assert "\tC\tC\t2\t0\t2.0\t0.0\t6\t1" in lines[1]
    assert lines[2].startswith("Q_TIE\t2\t")
    assert "\t\tA|B\t2\t0\t2.0\t0.0\t6\t4" in lines[2]


def test_write_parsimony_placement_alternative_table_writes_ranked_rows(
    tmp_path: Path,
) -> None:
    report = _placement_report()
    output_path = tmp_path / "alternative_placements.tsv"

    write_parsimony_placement_alternative_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "query_id\tplacement_rank\tedge_id\tchild_name\tdescendant_taxa\ttotal_steps\tadditional_steps\ttotal_weighted_score\tadditional_weighted_score\tis_equally_best"
    )
    assert len(lines) == 13
    assert lines[1].startswith("Q_C\t1\t")
    assert lines[7].startswith("Q_TIE\t1\t")


def test_write_parsimony_placement_run_json_and_artifacts_materialize_outputs(
    tmp_path: Path,
) -> None:
    report = _placement_report()
    run_json_path = tmp_path / "run.json"

    write_parsimony_placement_run_json(run_json_path, report)
    payload = json.loads(run_json_path.read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-placement"
    assert payload["query_count"] == 2
    assert payload["edge_count"] == 6

    tree_set_path = write_parsimony_placement_tree_set(
        tmp_path / "equally_best_placements.nwk",
        report,
    )
    assert len(load_newick_tree_set(tree_set_path)) == 5

    outputs = write_parsimony_placement_artifacts(tmp_path / "placement-run", report)
    assert set(outputs) == {
        "equally_best_tree_path",
        "summary_path",
        "alternative_path",
        "run_json_path",
    }
    assert len(load_newick_tree_set(outputs["equally_best_tree_path"])) == 5
    assert outputs["summary_path"].is_file()
    assert outputs["alternative_path"].is_file()
    assert outputs["run_json_path"].is_file()
