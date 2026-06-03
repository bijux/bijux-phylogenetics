from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_multi_start_from_alignment,
    write_nucleotide_likelihood_multi_start_artifacts,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_nucleotide_likelihood_multi_start_artifacts_records_rank_and_iterations(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_multi_start_artifacts(
        tmp_path / "likelihood-multi-start-run",
        report,
    )

    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "start_tree_source_kind\tstart_tree_source_label\tstart_tree_generation_seed\tsearch_algorithm\tstart_log_likelihood\tfinal_log_likelihood\tfinal_likelihood_rank\tfinal_topology_fingerprint\tsearch_iteration_count\taccepted_move_count\tevaluated_neighbor_count\tbranch_reoptimization_policy\tsubstitution_parameter_policy\tstopping_reason\tbest_run\tstart_tree_newick\tfinal_tree_newick\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert [row["final_likelihood_rank"] for row in payload["run_summaries"]] == [
        1,
        2,
        3,
    ]
    assert [row["search_iteration_count"] for row in payload["run_summaries"]] == [
        3,
        3,
        3,
    ]
