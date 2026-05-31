from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import loads_newick
import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment,
    resolve_candidate_tree_records,
    write_candidate_tree_site_likelihood_matrix_artifacts,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.runtime.errors import TreeParseError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def write_candidate_tree_set(path: Path, trees: list[str]) -> Path:
    path.write_text("".join(f"{tree}\n" for tree in trees), encoding="utf-8")
    return path


def test_package_likelihood_gateway_exports_candidate_tree_matrix_surface() -> None:
    assert (
        likelihood_api.evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment
        is evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment
    )
    assert (
        likelihood_api.write_candidate_tree_site_likelihood_matrix_artifacts
        is write_candidate_tree_site_likelihood_matrix_artifacts
    )
    assert (
        likelihood_api.resolve_candidate_tree_records is resolve_candidate_tree_records
    )


def test_candidate_tree_site_likelihood_matrix_reports_tree_by_site_rows_and_totals(
    tmp_path: Path,
) -> None:
    tree_set_path = write_candidate_tree_set(
        tmp_path / "candidate-trees.nwk",
        [
            "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);",
        ],
    )
    alignment_path = fixture("alignments", "jc69_site_pattern_alignment.fasta")

    report = evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment(
        tree_set_path,
        alignment_path,
        model_name="jc69",
    )

    assert report.model_name == "JC69"
    assert report.tree_count == 2
    assert report.site_count == 10
    assert report.pattern_count == 6
    assert report.compression_used is True
    assert report.expansion_policy == "candidate-tree-expanded-site-rows"
    assert report.comparison_caution_label == (
        "all candidate trees are rescored under one shared alignment/model surface; prior per-tree fitted-model differences are not preserved in this comparison"
    )
    assert report.parameter_values == {}
    assert [row.candidate_tree_id for row in report.candidate_trees] == [
        "candidate-tree-1",
        "candidate-tree-2",
    ]
    assert report.candidate_trees[
        0
    ].topology_fingerprint == rooted_topology_fingerprint(
        loads_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
    )
    assert report.candidate_trees[
        1
    ].topology_fingerprint == rooted_topology_fingerprint(
        loads_newick("(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);")
    )
    assert len(report.matrix_rows) == 20

    totals_by_tree_id: dict[str, float] = {}
    for row in report.matrix_rows:
        totals_by_tree_id[row.candidate_tree_id] = (
            totals_by_tree_id.get(
                row.candidate_tree_id,
                0.0,
            )
            + row.log_likelihood
        )

    for summary in report.candidate_trees:
        assert math.isclose(
            totals_by_tree_id[summary.candidate_tree_id],
            summary.log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )

    best_log_likelihood = max(row.log_likelihood for row in report.candidate_trees)
    for summary in report.candidate_trees:
        assert math.isclose(
            summary.observed_delta_log_likelihood,
            best_log_likelihood - summary.log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    assert math.isclose(
        min(row.observed_delta_log_likelihood for row in report.candidate_trees),
        0.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    first_tree_total = evaluate_jc69_tree_likelihood(
        loads_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"),
        load_fasta_alignment(alignment_path),
    )
    assert math.isclose(
        report.candidate_trees[0].log_likelihood,
        first_tree_total.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_candidate_tree_site_likelihood_matrix_rejects_incompatible_taxa(
    tmp_path: Path,
) -> None:
    tree_set_path = write_candidate_tree_set(
        tmp_path / "candidate-trees.nwk",
        [
            "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
            "((A:0.1,B:0.1):0.1,(C:0.1,E:0.1):0.1);",
        ],
    )

    try:
        evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment(
            tree_set_path,
            fixture("alignments", "jc69_site_pattern_alignment.fasta"),
            model_name="jc69",
        )
    except ValueError as error:
        assert (
            "candidate-tree-2 is incompatible with the shared alignment/model surface"
            in str(error)
        )
    else:
        raise AssertionError("candidate trees with incompatible taxa must fail")


def test_candidate_tree_site_likelihood_matrix_requires_multiple_candidates() -> None:
    try:
        resolve_candidate_tree_records(
            [loads_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")]
        )
    except TreeParseError as error:
        assert str(error) == (
            "candidate tree site likelihood matrix requires at least two candidate trees"
        )
    else:
        raise AssertionError("candidate tree matrix must require multiple candidates")


def test_write_candidate_tree_site_likelihood_matrix_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    tree_set_path = write_candidate_tree_set(
        tmp_path / "candidate-trees.nwk",
        [
            "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);",
        ],
    )
    report = evaluate_nucleotide_candidate_tree_site_likelihood_matrix_from_alignment(
        tree_set_path,
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )

    outputs = write_candidate_tree_site_likelihood_matrix_artifacts(
        tmp_path / "candidate-tree-matrix-run",
        report,
    )

    assert set(outputs) == {
        "candidate_tree_path",
        "summary_path",
        "matrix_path",
        "run_json_path",
    }
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "candidate_tree_id\tcandidate_tree_label\ttopology_fingerprint\ttree_newick\tlog_likelihood\tobserved_delta_log_likelihood\n"
        )
    )
    assert (
        outputs["matrix_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "model_name\tcandidate_tree_id\tcandidate_tree_label\ttaxon_order\tpattern_id\tpattern_weight\tsite_position\tsite_states\tlog_likelihood\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["model_name"] == "JC69"
    assert payload["tree_count"] == 2
    assert payload["site_count"] == 10
    assert payload["comparison_caution_label"] == (
        "all candidate trees are rescored under one shared alignment/model surface; prior per-tree fitted-model differences are not preserved in this comparison"
    )
    assert len(payload["candidate_trees"]) == 2
    assert set(payload["candidate_trees"][0]) == {
        "candidate_tree_id",
        "candidate_tree_label",
        "topology_fingerprint",
        "tree_newick",
        "log_likelihood",
        "observed_delta_log_likelihood",
    }
    assert len(payload["matrix_rows"]) == 20
