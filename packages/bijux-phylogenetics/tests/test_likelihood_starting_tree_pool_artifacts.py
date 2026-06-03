from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_pool import (
    write_nucleotide_likelihood_starting_tree_pool_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_nucleotide_likelihood_starting_tree_pool_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = build_nucleotide_likelihood_starting_tree_pool_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        random_start_tree_count=2,
        random_start_tree_seed=17,
    )

    outputs = write_nucleotide_likelihood_starting_tree_pool_artifacts(
        tmp_path / "starting-tree-pool",
        report,
    )

    assert set(outputs) == {
        "start_tree_path",
        "score_table_path",
        "run_json_path",
    }
    assert outputs["start_tree_path"].name == "starting_trees.nwk"
    assert outputs["score_table_path"].name == "starting_trees.tsv"
    assert outputs["run_json_path"].name == "run.json"
    assert (
        outputs["score_table_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "model_name\ttree_id\tsource_strategy\tgeneration_seed\tstarting_log_likelihood\ttopology_hash\tsubstitution_parameter_policy\tsubstitution_parameter_values\tsubstitution_parameter_warnings\ttree_newick\n"
        )
    )

    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))

    assert payload["algorithm"] == "nucleotide-likelihood-starting-tree-pool"
    assert payload["model_name"] == "JC69"
    assert payload["random_start_tree_count"] == 2
    assert payload["random_start_tree_seed"] == 17
    assert [row["tree_id"] for row in payload["starting_tree_summaries"]] == [
        "input-tree",
        "likelihood-stepwise-addition-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
