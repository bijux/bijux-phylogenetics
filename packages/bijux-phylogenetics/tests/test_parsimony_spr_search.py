from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    search_parsimony_nni,
    search_parsimony_spr,
    write_parsimony_spr_artifacts,
)
from bijux_phylogenetics.parsimony.spr import _iter_spr_move_candidates

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_spr_search_surface() -> None:
    assert parsimony_api.search_parsimony_spr is search_parsimony_spr
    assert parsimony_api.write_parsimony_spr_artifacts is write_parsimony_spr_artifacts


def test_parsimony_spr_candidates_exclude_self_regrafts() -> None:
    candidates = list(
        _iter_spr_move_candidates(
            loads_newick(
                fixture("spr_search_start_tree_5_taxa.nwk").read_text(encoding="utf-8")
            )
        )
    )

    assert candidates
    for candidate in candidates:
        assert candidate.regraft_target_branch_id != candidate.pruned_clade_id
        if candidate.regraft_target_descendant_taxa is not None:
            assert not set(candidate.regraft_target_descendant_taxa).issubset(
                set(candidate.pruned_descendant_taxa)
            )


def test_parsimony_spr_search_improves_one_case_where_nni_stays_stuck() -> None:
    nni_report = search_parsimony_nni(
        fixture("spr_search_start_tree_5_taxa.nwk"),
        fixture("spr_search_matrix.tsv"),
        method="fitch",
    )
    spr_report = search_parsimony_spr(
        fixture("spr_search_start_tree_5_taxa.nwk"),
        fixture("spr_search_matrix.tsv"),
        method="fitch",
    )

    assert nni_report.start_score == 3.0
    assert nni_report.final_score == 3.0
    assert nni_report.accepted_move_count == 0
    assert spr_report.algorithm == "parsimony-spr-search"
    assert spr_report.method == "fitch"
    assert spr_report.start_tree_newick == "((((A,D),B),C),E);"
    assert spr_report.start_score == 3.0
    assert spr_report.final_tree_newick == "((((A,B),C),D),E);"
    assert spr_report.final_score == 2.0
    assert spr_report.accepted_move_count == 1
    assert spr_report.evaluated_neighbor_count == 50
    assert spr_report.stopping_reason == "no-improving-neighbor"
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
            row.pruned_clade_id,
            row.regraft_target_branch_id,
            row.stopping_reason,
        )
        for row in spr_report.trace_rows
    ] == [
        (
            1,
            "start",
            0,
            None,
            3.0,
            None,
            None,
            "((((A,D),B),C),E);",
            None,
            None,
            None,
        ),
        (
            2,
            "accepted-move",
            1,
            3.0,
            2.0,
            -1.0,
            "((((A,D),B),C),E);",
            "((((A,B),C),D),E);",
            "D",
            "A|B|C",
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
            "((((A,B),C),D),E);",
            None,
            None,
            "no-improving-neighbor",
        ),
    ]


def test_write_parsimony_spr_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_parsimony_spr(
        fixture("spr_search_start_tree_5_taxa.nwk"),
        fixture("spr_search_matrix.tsv"),
        method="fitch",
    )

    outputs = write_parsimony_spr_artifacts(tmp_path / "spr-search-run", report)

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
            "event_index\tevent_kind\titeration\tscore_before\tscore_after\tscore_delta\ttree_before_newick\ttree_after_newick\tpruned_clade_id\tregraft_target_branch_id\tstopping_reason\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-spr-search"
    assert payload["method"] == "fitch"
    assert payload["start_tree_newick"] == "((((A,D),B),C),E);"
    assert payload["start_score"] == 3.0
    assert payload["final_tree_newick"] == "((((A,B),C),D),E);"
    assert payload["final_score"] == 2.0
    assert payload["accepted_move_count"] == 1
    assert payload["evaluated_neighbor_count"] == 50
    assert payload["stopping_reason"] == "no-improving-neighbor"
