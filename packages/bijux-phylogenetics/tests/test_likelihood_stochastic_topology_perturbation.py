from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment,
    write_nucleotide_likelihood_stochastic_topology_perturbation_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_stochastic_topology_perturbation_records_reproducible_nni_path() -> (
    None
):
    left_report = (
        search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment(
            fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
            perturbation_move_family="nni",
            local_search_method="nni",
            perturbation_seed=7,
            perturbation_move_count=2,
            upper_branch_length_bound=1.0,
        )
    )
    right_report = (
        search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment(
            fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
            perturbation_move_family="nni",
            local_search_method="nni",
            perturbation_seed=7,
            perturbation_move_count=2,
            upper_branch_length_bound=1.0,
        )
    )

    assert (
        left_report.algorithm
        == "nucleotide-likelihood-stochastic-topology-perturbation-search"
    )
    assert left_report.perturbation_move_family == "nni"
    assert left_report.local_search_method == "nni"
    assert left_report.perturbation_move_count_requested == 2
    assert left_report.perturbation_move_count_applied == 2
    assert left_report.local_search_stopping_reason == "no-improving-neighbor"
    assert left_report.local_search_accepted_move_count == 3
    assert left_report.local_search_evaluated_neighbor_count == 16
    assert left_report.perturbed_tree_newick == "((A:0.1,D:0.1):0.1,(B:0.1,C:0.1):0.1);"
    assert (
        left_report.local_search_start_tree_newick
        == "((A:0.999999999756461,D:2.43539016857288e-10):2.43539016857288e-10,(B:0.999999999756461,C:2.43539016857288e-10):2.43539016857288e-10);"
    )
    assert left_report.final_tree_newick != left_report.local_search_start_tree_newick
    assert math.isclose(
        left_report.local_search_start_log_likelihood,
        -57.25006765431177,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        left_report.final_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [
        (
            step.step_index,
            step.move_family,
            step.pivot_branch_id,
            step.sibling_clade_id,
            step.exchanged_clade_id,
            step.topology_fingerprint_before,
            step.topology_fingerprint_after,
        )
        for step in left_report.perturbation_steps
    ] == [
        (
            1,
            "nni",
            "A|C",
            "B",
            "A",
            "55d4921cbc851957d00b78341edf0c50186067abfc687a9d667708c07f2d640a",
            "6135ea2cf08b5016b0541f0a80074b31e0b9f1bc565b8411f1801575244545c1",
        ),
        (
            2,
            "nni",
            "A|B|C",
            "D",
            "B|C",
            "6135ea2cf08b5016b0541f0a80074b31e0b9f1bc565b8411f1801575244545c1",
            "18592aa3a645c91f5a28f9857ef6d4c8178481af10a48b3d04f98ec930ae3fd6",
        ),
    ]
    assert [step.tree_after_newick for step in left_report.perturbation_steps] == [
        "((A:0.1,(B:0.1,C:0.1):0.1):0.1,D:0.1);",
        "((A:0.1,D:0.1):0.1,(B:0.1,C:0.1):0.1);",
    ]
    assert [
        (
            step.tree_after_newick,
            step.topology_fingerprint_after,
        )
        for step in left_report.perturbation_steps
    ] == [
        (
            step.tree_after_newick,
            step.topology_fingerprint_after,
        )
        for step in right_report.perturbation_steps
    ]
    assert left_report.final_tree_newick == right_report.final_tree_newick
    assert left_report.final_log_likelihood == right_report.final_log_likelihood


def test_likelihood_stochastic_topology_perturbation_records_topology_changing_spr_path() -> (
    None
):
    report = (
        search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment(
            fixture("trees", "jc69_likelihood_spr_start_tree_5_taxa.nwk"),
            fixture("alignments", "jc69_likelihood_spr_alignment_5_taxa.fasta"),
            model_name="jc69",
            perturbation_move_family="spr",
            local_search_method="nni",
            perturbation_seed=11,
            perturbation_move_count=2,
            upper_branch_length_bound=1.0,
        )
    )

    assert report.perturbation_move_family == "spr"
    assert report.local_search_method == "nni"
    assert report.perturbation_move_count_applied == 2
    assert [step.pruned_clade_id for step in report.perturbation_steps] == [
        "E",
        "A|E",
    ]
    assert [step.regraft_target_branch_id for step in report.perturbation_steps] == [
        "A",
        "B",
    ]
    assert all(
        step.topology_fingerprint_before != step.topology_fingerprint_after
        for step in report.perturbation_steps
    )
    assert (
        report.perturbed_tree_newick
        == "((((A:0.05,E:0.1):0.05,B:0.05):0.05,D:0.1):0.1,C:0.1);"
    )
    assert report.local_search_accepted_move_count == 1
    assert report.local_search_stopping_reason == "no-improving-neighbor"
    assert math.isclose(
        report.local_search_start_log_likelihood,
        -8.64203998404885,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.final_log_likelihood,
        -8.642039983561773,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.final_log_likelihood > report.local_search_start_log_likelihood


def test_write_likelihood_stochastic_topology_perturbation_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = (
        search_nucleotide_likelihood_stochastic_topology_perturbation_from_alignment(
            fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
            perturbation_move_family="nni",
            local_search_method="nni",
            perturbation_seed=7,
            perturbation_move_count=2,
            upper_branch_length_bound=1.0,
        )
    )

    outputs = write_nucleotide_likelihood_stochastic_topology_perturbation_artifacts(
        tmp_path / "likelihood-stochastic-perturbation-run",
        report,
    )

    assert set(outputs) == {
        "input_tree_path",
        "perturbed_tree_path",
        "final_tree_path",
        "perturbation_trace_path",
        "run_json_path",
    }
    assert (
        outputs["perturbation_trace_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "step_index\tmove_family\ttree_before_newick\ttree_after_newick\ttopology_fingerprint_before\ttopology_fingerprint_after\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tpruned_clade_id\tregraft_target_branch_id\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert (
        payload["algorithm"]
        == "nucleotide-likelihood-stochastic-topology-perturbation-search"
    )
    assert payload["perturbation_move_family"] == "nni"
    assert payload["local_search_method"] == "nni"
    assert payload["perturbation_seed"] == 7
    assert payload["perturbation_move_count_applied"] == 2
    assert payload["local_search_stopping_reason"] == "no-improving-neighbor"
    assert payload["local_search_accepted_move_count"] == 3
    assert payload["final_tree_newick"] == report.final_tree_newick
    assert payload["perturbation_steps"][0]["pivot_branch_id"] == "A|C"
    assert payload["perturbation_steps"][1]["exchanged_clade_id"] == "B|C"
