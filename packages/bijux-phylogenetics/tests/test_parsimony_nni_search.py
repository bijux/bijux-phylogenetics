from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    search_parsimony_nni,
    write_parsimony_nni_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_nni_search_surface() -> None:
    assert parsimony_api.search_parsimony_nni is search_parsimony_nni
    assert parsimony_api.write_parsimony_nni_artifacts is write_parsimony_nni_artifacts


def test_parsimony_nni_search_accepts_one_score_improving_move_and_stops() -> None:
    report = search_parsimony_nni(
        fixture("nni_search_start_tree_4_taxa.nwk"),
        fixture("nni_search_matrix.tsv"),
        method="fitch",
    )

    assert report.algorithm == "parsimony-nni-search"
    assert report.method == "fitch"
    assert report.start_tree_newick == "(((A,C),B),D);"
    assert report.start_score == 4.0
    assert report.final_tree_newick == "(((A,B),C),D);"
    assert report.final_score == 2.0
    assert report.accepted_move_count == 1
    assert report.evaluated_neighbor_count == 8
    assert report.stopping_reason == "no-improving-neighbor"
    assert [
        (
            row.event_index,
            row.event_kind,
            row.iteration,
            row.score_before,
            row.score_after,
            row.score_delta,
            row.tree_before_newick,
            row.tree_after_newick,
            row.pivot_branch_id,
            row.sibling_clade_id,
            row.exchanged_clade_id,
            row.stopping_reason,
        )
        for row in report.trace_rows
    ] == [
        (
            1,
            "start",
            0,
            None,
            4.0,
            None,
            None,
            "(((A,C),B),D);",
            None,
            None,
            None,
            None,
        ),
        (
            2,
            "accepted-move",
            1,
            4.0,
            2.0,
            -2.0,
            "(((A,C),B),D);",
            "(((A,B),C),D);",
            "A|C",
            "B",
            "C",
            None,
        ),
        (
            3,
            "final",
            1,
            None,
            2.0,
            None,
            None,
            "(((A,B),C),D);",
            None,
            None,
            None,
            "no-improving-neighbor",
        ),
    ]


def test_write_parsimony_nni_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_parsimony_nni(
        fixture("nni_search_start_tree_4_taxa.nwk"),
        fixture("nni_search_matrix.tsv"),
        method="fitch",
    )

    outputs = write_parsimony_nni_artifacts(tmp_path / "nni-search-run", report)

    assert set(outputs) == {
        "start_tree_path",
        "final_tree_path",
        "trace_path",
        "run_json_path",
    }
    assert (
        outputs["trace_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "event_index\tevent_kind\titeration\tscore_before\tscore_after\tscore_delta\ttree_before_newick\ttree_after_newick\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tstopping_reason\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-nni-search"
    assert payload["method"] == "fitch"
    assert payload["start_tree_newick"] == "(((A,C),B),D);"
    assert payload["start_score"] == 4.0
    assert payload["final_tree_newick"] == "(((A,B),C),D);"
    assert payload["final_score"] == 2.0
    assert payload["accepted_move_count"] == 1
    assert payload["evaluated_neighbor_count"] == 8
    assert payload["stopping_reason"] == "no-improving-neighbor"


def test_parsimony_nni_search_records_clean_stop_when_start_tree_is_already_optimal() -> (
    None
):
    report = search_parsimony_nni(
        loads_newick("(((A,B),C),D);"),
        fixture("nni_search_matrix.tsv"),
        method="fitch",
    )

    assert report.start_tree_newick == "(((A,B),C),D);"
    assert report.final_tree_newick == "(((A,B),C),D);"
    assert report.start_score == 2.0
    assert report.final_score == 2.0
    assert report.accepted_move_count == 0
    assert report.evaluated_neighbor_count == 4
    assert report.stopping_reason == "no-improving-neighbor"
    assert [row.event_kind for row in report.trace_rows] == ["start", "final"]
