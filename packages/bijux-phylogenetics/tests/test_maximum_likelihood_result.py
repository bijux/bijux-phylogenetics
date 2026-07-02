from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.likelihood.maximum_likelihood_result import (
    NucleotideMaximumLikelihoodResult,
    build_nucleotide_maximum_likelihood_result,
    load_nucleotide_maximum_likelihood_result_json,
    write_nucleotide_maximum_likelihood_result_json,
)
from bijux_phylogenetics.phylo.likelihood.tree_inference import (
    infer_nucleotide_likelihood_tree_from_alignment,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_maximum_likelihood_result_round_trips_inference_payload(
    tmp_path: Path,
) -> None:
    inference_report = infer_nucleotide_likelihood_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        search_method="nni",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    result = build_nucleotide_maximum_likelihood_result(inference_report)
    written_path = write_nucleotide_maximum_likelihood_result_json(
        tmp_path / "maximum-likelihood-result.json",
        result,
    )
    reloaded_result = load_nucleotide_maximum_likelihood_result_json(written_path)

    assert isinstance(reloaded_result, NucleotideMaximumLikelihoodResult)
    assert reloaded_result.to_dict() == result.to_dict()
    assert reloaded_result.model_name == "JC69"
    assert reloaded_result.final_tree_newick == inference_report.best_final_tree_newick
    assert reloaded_result.final_log_likelihood == inference_report.best_final_log_likelihood
    assert reloaded_result.parameter_values == (
        inference_report.model_selection_report.rows[0].parameter_values
    )
    assert reloaded_result.search_trace_rows
    assert reloaded_result.support_reports.bootstrap is None
    assert reloaded_result.support_reports.sh_like is None
