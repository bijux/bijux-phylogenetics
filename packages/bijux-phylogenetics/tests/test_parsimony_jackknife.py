from __future__ import annotations

import json
from pathlib import Path

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    bootstrap_parsimony,
    jackknife_parsimony,
    write_parsimony_jackknife_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_jackknife_surface() -> None:
    assert parsimony_api.jackknife_parsimony is jackknife_parsimony
    assert (
        parsimony_api.write_parsimony_jackknife_artifacts
        is write_parsimony_jackknife_artifacts
    )


def test_jackknife_parsimony_repeats_same_seed_identically() -> None:
    left_report = jackknife_parsimony(
        fixture("jackknife_matrix.tsv"),
        method="fitch",
        replicate_count=12,
        random_seed=5,
        retain_probability=0.6,
    )
    right_report = jackknife_parsimony(
        fixture("jackknife_matrix.tsv"),
        method="fitch",
        replicate_count=12,
        random_seed=5,
        retain_probability=0.6,
    )

    assert left_report.algorithm == "parsimony-jackknife"
    assert left_report.method == "fitch"
    assert left_report.candidate_tree_count == 15
    assert left_report.reference_tree_newick == "(((A,B),C),D);"
    assert left_report.reference_score == 7.0
    assert [row.retained_character_count for row in left_report.replicate_rows] == [
        1,
        2,
        5,
        3,
        3,
        3,
        2,
        1,
        5,
        3,
        4,
        2,
    ]
    assert [
        (
            row.replicate_index,
            row.retained_character_ids,
            row.best_score,
            row.optimal_tree_count,
            row.tree_newick,
        )
        for row in left_report.replicate_rows
    ] == [
        (
            row.replicate_index,
            row.retained_character_ids,
            row.best_score,
            row.optimal_tree_count,
            row.tree_newick,
        )
        for row in right_report.replicate_rows
    ]
    assert [
        (
            row.branch_id,
            row.supporting_tree_count,
            row.clade_frequency,
            row.support_percent,
        )
        for row in left_report.clade_support_rows
    ] == [
        ("A|B", 11, 0.916666666666667, 91.6666666666667),
        ("A|B|C", 11, 0.916666666666667, 91.6666666666667),
    ]


def test_jackknife_uses_without_replacement_subsampling_not_bootstrap_draws() -> None:
    jackknife_report = jackknife_parsimony(
        fixture("jackknife_matrix.tsv"),
        method="fitch",
        replicate_count=3,
        random_seed=5,
        retain_probability=0.6,
    )
    bootstrap_report = bootstrap_parsimony(
        fixture("jackknife_matrix.tsv"),
        method="fitch",
        replicate_count=3,
        random_seed=5,
    )

    assert all(
        len(row.retained_character_ids) == len(set(row.retained_character_ids))
        for row in jackknife_report.replicate_rows
    )
    assert any(
        len(row.sampled_character_ids) != len(set(row.sampled_character_ids))
        for row in bootstrap_report.replicate_rows
    )


def test_write_parsimony_jackknife_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = jackknife_parsimony(
        fixture("jackknife_matrix.tsv"),
        method="fitch",
        replicate_count=12,
        random_seed=5,
        retain_probability=0.6,
    )

    outputs = write_parsimony_jackknife_artifacts(tmp_path / "jackknife-run", report)

    assert set(outputs) == {
        "reference_tree_path",
        "replicate_trees_path",
        "replicate_scores_path",
        "retained_characters_path",
        "clade_support_path",
        "consensus_tree_path",
        "clade_frequencies_path",
        "run_json_path",
    }
    assert (
        outputs["replicate_scores_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "replicate_index\tretained_character_count\tbest_score\toptimal_tree_count\ttree_newick\n"
        )
    )
    assert (
        outputs["retained_characters_path"]
        .read_text(encoding="utf-8")
        .startswith("replicate_index\tretained_index\tsource_character_id\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-jackknife"
    assert payload["method"] == "fitch"
    assert payload["replicate_count"] == 12
    assert payload["random_seed"] == 5
    assert payload["retain_probability"] == 0.6
