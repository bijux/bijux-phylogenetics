from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick_tree_set
from bijux_phylogenetics.phylo.likelihood import (
    bootstrap_nucleotide_likelihood_tree_inference_from_alignment,
    write_nucleotide_likelihood_bootstrap_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_native_likelihood_bootstrap_is_seed_reproducible() -> None:
    left_report = bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=2,
        start_tree_seed=5,
        replicate_count=4,
        bootstrap_seed=9,
        upper_branch_length_bound=1.0,
    )
    right_report = bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=2,
        start_tree_seed=5,
        replicate_count=4,
        bootstrap_seed=9,
        upper_branch_length_bound=1.0,
    )

    assert left_report == right_report
    assert left_report.algorithm == "nucleotide-likelihood-bootstrap-tree-inference"
    assert left_report.requested_model_name == "jc69"
    assert left_report.selected_reference_model_name == "JC69"
    assert left_report.search_method == "nni"
    assert left_report.taxon_count == 4
    assert left_report.site_count == 12
    assert left_report.pattern_count == 2
    assert left_report.replicate_count == 4
    assert len(left_report.replicate_rows) == 4
    assert all(
        len(row.sampled_site_indices) == left_report.site_count
        for row in left_report.replicate_rows
    )
    assert all(
        0.0 <= row.support_percent <= 100.0
        for row in left_report.clade_support_rows
    )
    assert all(
        row.supporting_tree_count <= left_report.replicate_count
        for row in left_report.clade_support_rows
    )


def test_native_likelihood_bootstrap_writes_governed_outputs(
    tmp_path: Path,
) -> None:
    report = bootstrap_nucleotide_likelihood_tree_inference_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=2,
        start_tree_seed=5,
        replicate_count=4,
        bootstrap_seed=9,
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_bootstrap_artifacts(
        tmp_path / "likelihood-bootstrap-tree",
        report,
    )

    assert set(outputs) == {
        "reference_tree_path",
        "replicate_trees_path",
        "replicate_draws_path",
        "clade_support_path",
        "consensus_tree_path",
        "clade_frequencies_path",
        "run_json_path",
    }
    assert len(load_newick_tree_set(outputs["replicate_trees_path"])) == report.replicate_count
    assert outputs["replicate_draws_path"].read_text(encoding="utf-8").startswith(
        "replicate_index\tsampled_site_indices\treplicate_start_tree_seed\tselected_model_name\tbest_run_source_label\tfinal_log_likelihood\tfinal_topology_fingerprint\taccepted_move_count\tsearch_iteration_count\tfinal_tree_newick\n"
    )
    assert outputs["clade_support_path"].read_text(encoding="utf-8").startswith(
        "branch_id\tnode_label\tdescendant_taxa\tsupporting_tree_count\tclade_frequency\tsupport_percent\n"
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-bootstrap-tree-inference"
    assert payload["selected_reference_model_name"] == "JC69"
    assert len(payload["replicate_rows"]) == 4
    assert len(payload["clade_support_rows"]) == len(report.clade_support_rows)
