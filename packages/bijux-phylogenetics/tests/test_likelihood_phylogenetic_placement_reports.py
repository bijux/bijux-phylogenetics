from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.phylo.likelihood import (
    place_queries_by_likelihood_from_alignment,
    write_likelihood_placement_alternative_table,
    write_likelihood_placement_artifacts,
    write_likelihood_placement_run_json,
    write_likelihood_placement_summary_table,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def _placement_report():
    return place_queries_by_likelihood_from_alignment(
        fixture("trees", "likelihood_placement_reference_tree_4_taxa.nwk"),
        fixture("alignments", "likelihood_placement_reference_alignment_4_taxa.fasta"),
        fixture("alignments", "likelihood_placement_query_alignment_2_taxa.fasta"),
    )


def test_write_likelihood_placement_summary_table_writes_best_edge_rows(
    tmp_path: Path,
) -> None:
    report = _placement_report()
    output_path = tmp_path / "summary.tsv"

    write_likelihood_placement_summary_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "query_id\tsite_count\tpattern_count\tbest_edge_id\tbest_child_name\tbest_descendant_taxa\tbest_original_branch_length\tbest_proximal_length\tbest_distal_length\tbest_pendant_length\tbest_log_likelihood\tbest_likelihood_weight_ratio\tcandidate_placement_count\tequally_best_placement_count"
    )
    assert lines[1].startswith("Q_A\t12\t")
    assert "\tA\tA\t0.1\t" in lines[1]
    assert lines[2].startswith("Q_C\t12\t")
    assert "\tC\tC\t0.1\t" in lines[2]


def test_write_likelihood_placement_alternative_table_writes_ranked_rows(
    tmp_path: Path,
) -> None:
    report = _placement_report()
    output_path = tmp_path / "alternative_placements.tsv"

    write_likelihood_placement_alternative_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "query_id\tplacement_rank\tedge_id\tchild_name\tdescendant_taxa\toriginal_branch_length\toptimized_proximal_length\toptimized_distal_length\toptimized_pendant_length\tlog_likelihood\tlikelihood_weight_ratio\tfunction_evaluation_count\toptimization_pass_count\tconverged"
    )
    assert len(lines) == 13
    assert lines[1].startswith("Q_A\t1\t")
    assert lines[7].startswith("Q_C\t1\t")


def test_write_likelihood_placement_run_json_and_artifacts_materialize_governed_outputs(
    tmp_path: Path,
) -> None:
    report = _placement_report()
    run_json_path = tmp_path / "run.json"

    write_likelihood_placement_run_json(run_json_path, report)
    payload = json.loads(run_json_path.read_text(encoding="utf-8"))
    assert payload["model_name"] == "JC69"
    assert payload["query_count"] == 2
    assert payload["edge_count"] == 6

    outputs = write_likelihood_placement_artifacts(tmp_path / "placement-run", report)
    assert set(outputs) == {
        "alternative_path",
        "best_tree_path",
        "run_json_path",
        "summary_path",
    }
    assert len(load_newick_tree_set(outputs["best_tree_path"])) == 2
    assert outputs["summary_path"].is_file()
    assert outputs["alternative_path"].is_file()
    assert outputs["run_json_path"].is_file()
