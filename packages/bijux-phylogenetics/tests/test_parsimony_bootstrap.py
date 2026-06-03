from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    bootstrap_parsimony,
    write_parsimony_bootstrap_artifacts,
)
from bijux_phylogenetics.parsimony.bootstrap import _build_clade_support_rows

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_bootstrap_surface() -> None:
    assert parsimony_api.bootstrap_parsimony is bootstrap_parsimony
    assert (
        parsimony_api.write_parsimony_bootstrap_artifacts
        is write_parsimony_bootstrap_artifacts
    )


def test_bootstrap_parsimony_repeats_same_seed_identically() -> None:
    left_report = bootstrap_parsimony(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
        replicate_count=20,
        random_seed=1,
    )
    right_report = bootstrap_parsimony(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
        replicate_count=20,
        random_seed=1,
    )

    assert left_report.algorithm == "parsimony-bootstrap"
    assert left_report.method == "fitch"
    assert left_report.candidate_tree_count == 15
    assert left_report.reference_score == 5.0
    assert left_report.reference_optimal_tree_count == 5
    assert left_report.reference_tree_newick == "(((A,B),C),D);"
    assert [
        (
            row.replicate_index,
            row.sampled_character_ids,
            row.best_score,
            row.optimal_tree_count,
            row.tree_newick,
        )
        for row in left_report.replicate_rows
    ] == [
        (
            row.replicate_index,
            row.sampled_character_ids,
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
        (
            "A|B",
            17,
            0.85,
            85.0,
        ),
        (
            "A|B|C",
            20,
            1.0,
            100.0,
        ),
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
        (
            row.branch_id,
            row.supporting_tree_count,
            row.clade_frequency,
            row.support_percent,
        )
        for row in right_report.clade_support_rows
    ]


def test_bootstrap_clade_support_rows_use_clade_identity_not_node_index() -> None:
    reference_tree = loads_newick("(((A,B),C),(D,E));")
    replicate_trees = [
        loads_newick("((D,E),((A,B),C));"),
        loads_newick("((D,E),(A,(B,C)));"),
    ]

    rows = _build_clade_support_rows(reference_tree, replicate_trees)

    assert [
        (
            row.branch_id,
            row.descendant_taxa,
            row.supporting_tree_count,
            row.clade_frequency,
        )
        for row in rows
    ] == [
        ("A|B", ["A", "B"], 1, 0.5),
        ("D|E", ["D", "E"], 2, 1.0),
        ("A|B|C", ["A", "B", "C"], 2, 1.0),
    ]


def test_write_parsimony_bootstrap_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = bootstrap_parsimony(
        fixture("bootstrap_matrix.tsv"),
        method="fitch",
        replicate_count=20,
        random_seed=1,
    )

    outputs = write_parsimony_bootstrap_artifacts(tmp_path / "bootstrap-run", report)

    assert set(outputs) == {
        "reference_tree_path",
        "replicate_trees_path",
        "replicate_scores_path",
        "replicate_draws_path",
        "clade_support_path",
        "consensus_tree_path",
        "clade_frequencies_path",
        "run_json_path",
    }
    assert (
        outputs["replicate_scores_path"]
        .read_text(encoding="utf-8")
        .startswith("replicate_index\tbest_score\toptimal_tree_count\ttree_newick\n")
    )
    assert (
        outputs["replicate_draws_path"]
        .read_text(encoding="utf-8")
        .startswith("replicate_index\tdraw_index\tsource_character_id\n")
    )
    assert (
        outputs["clade_support_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tnode_name\tdescendant_taxa\tsupporting_tree_count\tclade_frequency\tsupport_percent\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-bootstrap"
    assert payload["method"] == "fitch"
    assert payload["replicate_count"] == 20
    assert payload["random_seed"] == 1
    assert payload["candidate_tree_count"] == 15
    assert payload["reference_tree_newick"] == "(((A,B),C),D);"
