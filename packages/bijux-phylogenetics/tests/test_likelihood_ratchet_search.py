from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_ratchet_from_alignment,
    write_nucleotide_likelihood_ratchet_artifacts,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_ratchet_search_surface() -> None:
    assert (
        likelihood_api.search_nucleotide_likelihood_ratchet_from_alignment
        is search_nucleotide_likelihood_ratchet_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_ratchet_artifacts
        is write_nucleotide_likelihood_ratchet_artifacts
    )


def test_likelihood_ratchet_repeats_same_seed_identically() -> None:
    left_report = search_nucleotide_likelihood_ratchet_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        local_search_method="nni",
        cycle_count=2,
        perturbation_seed=5,
        perturbed_site_count=2,
        perturbation_factor=3,
        upper_branch_length_bound=1.0,
    )
    right_report = search_nucleotide_likelihood_ratchet_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        local_search_method="nni",
        cycle_count=2,
        perturbation_seed=5,
        perturbed_site_count=2,
        perturbation_factor=3,
        upper_branch_length_bound=1.0,
    )

    assert left_report.algorithm == "nucleotide-likelihood-ratchet-search"
    assert left_report.local_search_method == "nni"
    assert left_report.cycle_count == 2
    assert left_report.perturbation_seed == 5
    assert left_report.perturbed_site_count == 2
    assert left_report.perturbation_factor == 3
    assert math.isclose(
        left_report.start_log_likelihood,
        -54.442517349645854,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        left_report.final_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert left_report.best_log_likelihood == left_report.final_log_likelihood
    assert left_report.final_tree_newick == right_report.final_tree_newick
    assert left_report.best_tree_newick == right_report.best_tree_newick
    assert [
        (
            row.cycle_index,
            row.reweighted_site_positions,
            row.temporary_site_weights,
            row.perturbed_score,
            row.restored_score,
            row.best_tree_improved,
        )
        for row in left_report.cycle_rows
    ] == [
        (
            row.cycle_index,
            row.reweighted_site_positions,
            row.temporary_site_weights,
            row.perturbed_score,
            row.restored_score,
            row.best_tree_improved,
        )
        for row in right_report.cycle_rows
    ]


def test_likelihood_ratchet_records_temporary_weights_best_history_and_restored_scores() -> (
    None
):
    report = search_nucleotide_likelihood_ratchet_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        local_search_method="nni",
        cycle_count=2,
        perturbation_seed=5,
        perturbed_site_count=2,
        perturbation_factor=3,
        upper_branch_length_bound=1.0,
    )

    assert (
        report.start_tree_newick
        == "(((A:2.43539016857288e-10,C:0.999999999756461):2.43539016857288e-10,B:2.43539016857288e-10):0.999999999756461,D:0.999999999756461);"
    )
    assert (
        report.final_tree_newick
        == "((A:2.43539016857288e-10,B:2.43539016857288e-10):0.999999999756461,(C:2.43539016857288e-10,D:2.43539016857288e-10):0.999999999756461);"
    )
    assert len(report.cycle_rows) == 2
    assert len(report.best_tree_history_rows) == 2

    first_cycle, second_cycle = report.cycle_rows

    assert first_cycle.reweighted_site_positions == [5, 10]
    assert first_cycle.temporary_site_weights == {5: 3, 10: 3}
    assert first_cycle.perturbed_alignment_length == 16
    assert first_cycle.perturbed_pattern_count == 2
    assert first_cycle.perturbed_search_algorithm == "nucleotide-likelihood-nni-search"
    assert first_cycle.perturbed_accepted_move_count == 2
    assert first_cycle.perturbed_evaluated_neighbor_count == 12
    assert first_cycle.perturbed_stopping_reason == "no-improving-neighbor"
    assert math.isclose(
        first_cycle.perturbed_score,
        -45.51366626396895,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert first_cycle.restored_search_algorithm == "nucleotide-likelihood-nni-search"
    assert first_cycle.restored_accepted_move_count == 0
    assert first_cycle.restored_evaluated_neighbor_count == 4
    assert first_cycle.restored_stopping_reason == "no-improving-neighbor"
    assert math.isclose(
        first_cycle.restored_score,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert first_cycle.best_tree_improved is True

    assert second_cycle.reweighted_site_positions == [6, 12]
    assert second_cycle.temporary_site_weights == {6: 3, 12: 3}
    assert second_cycle.perturbed_accepted_move_count == 0
    assert second_cycle.restored_accepted_move_count == 0
    assert second_cycle.best_tree_improved is False

    assert [row.cycle_index for row in report.best_tree_history_rows] == [0, 1]
    assert report.best_tree_history_rows[1].best_tree_newick == report.best_tree_newick


def test_write_nucleotide_likelihood_ratchet_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = search_nucleotide_likelihood_ratchet_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        local_search_method="nni",
        cycle_count=2,
        perturbation_seed=5,
        perturbed_site_count=2,
        perturbation_factor=3,
        upper_branch_length_bound=1.0,
    )

    outputs = write_nucleotide_likelihood_ratchet_artifacts(
        tmp_path / "likelihood-ratchet-run",
        report,
    )

    assert set(outputs) == {
        "input_tree_path",
        "start_tree_path",
        "final_tree_path",
        "cycle_table_path",
        "best_tree_history_path",
        "run_json_path",
    }
    assert (
        outputs["cycle_table_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "cycle_index\tstart_log_likelihood\tstart_tree_newick\treweighted_site_positions\ttemporary_site_weights\tperturbation_factor\tperturbed_alignment_length\tperturbed_pattern_count\tperturbed_search_algorithm\tperturbed_score\tperturbed_tree_newick\tperturbed_accepted_move_count\tperturbed_evaluated_neighbor_count\tperturbed_stopping_reason\trestored_search_algorithm\trestored_score\trestored_tree_newick\trestored_accepted_move_count\trestored_evaluated_neighbor_count\trestored_stopping_reason\tbest_score_after_cycle\tbest_tree_after_cycle\tbest_tree_improved\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-likelihood-ratchet-search"
    assert payload["local_search_method"] == "nni"
    assert payload["cycle_count"] == 2
    assert payload["perturbation_seed"] == 5
    assert payload["perturbed_site_count"] == 2
    assert payload["perturbation_factor"] == 3
    assert payload["cycle_rows"][0]["reweighted_site_positions"] == [5, 10]
    assert payload["cycle_rows"][0]["temporary_site_weights"] == {"5": 3, "10": 3}
    assert payload["cycle_rows"][0]["best_tree_improved"] is True
    assert payload["best_tree_history_rows"][1]["cycle_index"] == 1
