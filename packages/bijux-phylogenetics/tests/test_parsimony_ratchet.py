from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    run_parsimony_ratchet,
    write_parsimony_ratchet_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_ratchet_surface() -> None:
    assert parsimony_api.run_parsimony_ratchet is run_parsimony_ratchet
    assert (
        parsimony_api.write_parsimony_ratchet_artifacts
        is write_parsimony_ratchet_artifacts
    )


def test_parsimony_ratchet_repeats_same_seed_identically() -> None:
    left_report = run_parsimony_ratchet(
        fixture("ratchet_search_start_tree_5_taxa.nwk"),
        fixture("ratchet_search_matrix.tsv"),
        method="fitch",
        cycle_count=3,
        random_seed=1,
        perturbed_character_count=1,
        perturbation_factor=2.0,
    )
    right_report = run_parsimony_ratchet(
        fixture("ratchet_search_start_tree_5_taxa.nwk"),
        fixture("ratchet_search_matrix.tsv"),
        method="fitch",
        cycle_count=3,
        random_seed=1,
        perturbed_character_count=1,
        perturbation_factor=2.0,
    )

    assert left_report == right_report
    assert left_report.algorithm == "parsimony-ratchet"
    assert left_report.start_tree_newick == "(A,((B,D),(C,E)));"
    assert left_report.start_score == 5.0
    assert left_report.final_tree_newick == "(((A,(B,E)),C),D);"
    assert left_report.final_score == 4.0
    assert left_report.best_tree_newick == "(((A,(B,E)),C),D);"
    assert left_report.best_score == 4.0
    assert [
        (
            row.cycle_index,
            row.start_score,
            row.start_tree_newick,
            row.perturbed_character_ids,
            row.perturbation_factor,
            row.perturbed_score,
            row.perturbed_tree_newick,
            row.perturbed_accepted_move_count,
            row.restored_score,
            row.restored_tree_newick,
            row.restored_accepted_move_count,
            row.best_score_after_cycle,
            row.best_tree_after_cycle,
            row.best_tree_improved,
        )
        for row in left_report.cycle_rows
    ] == [
        (
            1,
            5.0,
            "(A,((B,D),(C,E)));",
            ["c1_terminal_be"],
            2.0,
            5.0,
            "(((A,(B,E)),C),D);",
            2,
            4.0,
            "(((A,(B,E)),C),D);",
            0,
            4.0,
            "(((A,(B,E)),C),D);",
            True,
        ),
        (
            2,
            4.0,
            "(((A,(B,E)),C),D);",
            ["c3_split_cd"],
            2.0,
            5.0,
            "(((A,(B,E)),C),D);",
            0,
            4.0,
            "(((A,(B,E)),C),D);",
            0,
            4.0,
            "(((A,(B,E)),C),D);",
            False,
        ),
        (
            3,
            4.0,
            "(((A,(B,E)),C),D);",
            ["c1_terminal_be"],
            2.0,
            5.0,
            "(((A,(B,E)),C),D);",
            0,
            4.0,
            "(((A,(B,E)),C),D);",
            0,
            4.0,
            "(((A,(B,E)),C),D);",
            False,
        ),
    ]
    assert [
        (
            row.history_index,
            row.cycle_index,
            row.best_score,
            row.best_tree_newick,
        )
        for row in left_report.best_tree_history_rows
    ] == [
        (1, 0, 5.0, "(A,((B,D),(C,E)));"),
        (2, 1, 4.0, "(((A,(B,E)),C),D);"),
    ]


def test_parsimony_ratchet_requires_real_temporary_reweighting_to_escape_local_optimum() -> (
    None
):
    report = run_parsimony_ratchet(
        fixture("ratchet_search_start_tree_5_taxa.nwk"),
        fixture("ratchet_search_matrix.tsv"),
        method="fitch",
        cycle_count=3,
        random_seed=1,
        perturbed_character_count=1,
        perturbation_factor=1.0,
    )

    assert report.start_score == 5.0
    assert report.final_score == 5.0
    assert report.best_score == 5.0
    assert report.final_tree_newick == "(A,((B,D),(C,E)));"
    assert report.best_tree_newick == "(A,((B,D),(C,E)));"
    assert [row.perturbed_character_ids for row in report.cycle_rows] == [
        ["c1_terminal_be"],
        ["c3_split_cd"],
        ["c1_terminal_be"],
    ]
    assert all(row.perturbed_accepted_move_count == 0 for row in report.cycle_rows)
    assert all(row.restored_accepted_move_count == 0 for row in report.cycle_rows)
    assert len(report.best_tree_history_rows) == 1


def test_write_parsimony_ratchet_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = run_parsimony_ratchet(
        fixture("ratchet_search_start_tree_5_taxa.nwk"),
        fixture("ratchet_search_matrix.tsv"),
        method="fitch",
        cycle_count=3,
        random_seed=1,
        perturbed_character_count=1,
        perturbation_factor=2.0,
    )

    outputs = write_parsimony_ratchet_artifacts(tmp_path / "ratchet-run", report)

    assert set(outputs) == {
        "start_tree_path",
        "final_tree_path",
        "best_tree_path",
        "cycle_history_path",
        "best_tree_history_path",
        "run_json_path",
    }
    assert (
        outputs["cycle_history_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "cycle_index\tstart_score\tstart_tree_newick\tperturbed_character_ids\tperturbation_factor\tperturbed_score\tperturbed_tree_newick\tperturbed_accepted_move_count\trestored_score\trestored_tree_newick\trestored_accepted_move_count\tbest_score_after_cycle\tbest_tree_after_cycle\tbest_tree_improved\n"
        )
    )
    assert (
        outputs["best_tree_history_path"]
        .read_text(encoding="utf-8")
        .startswith("history_index\tcycle_index\tbest_score\tbest_tree_newick\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-ratchet"
    assert payload["method"] == "fitch"
    assert payload["cycle_count"] == 3
    assert payload["random_seed"] == 1
    assert payload["perturbed_character_count"] == 1
    assert payload["perturbation_factor"] == 2.0
    assert payload["start_score"] == 5.0
    assert payload["final_score"] == 4.0
    assert payload["best_score"] == 4.0
    assert payload["best_tree_newick"] == "(((A,(B,E)),C),D);"
