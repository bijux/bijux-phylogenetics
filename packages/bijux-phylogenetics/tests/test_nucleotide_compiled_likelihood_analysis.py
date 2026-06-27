from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
)
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    evaluate_compiled_finite_state_likelihood_analysis,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    compile_selected_nucleotide_likelihood_analysis,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("model_name", "kwargs", "report_evaluator"),
    [
        ("jc69", {}, evaluate_jc69_tree_likelihood),
        (
            "hky85",
            {
                "kappa": 4.0,
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
            },
            evaluate_hky85_tree_likelihood,
        ),
        (
            "gtr",
            {
                "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 4.5,
                    "AT": 0.8,
                    "CG": 1.6,
                    "CT": 2.4,
                    "GT": 3.1,
                },
            },
            evaluate_gtr_tree_likelihood,
        ),
    ],
)
def test_compiled_nucleotide_analysis_matches_existing_tree_reports(
    model_name: str,
    kwargs: dict[str, object],
    report_evaluator,
) -> None:
    tree = load_tree(fixture("trees", f"{model_name}_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", f"{model_name}_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_selected_nucleotide_likelihood_analysis(
        tree,
        records,
        model_name=model_name,
        owner_name="compiled nucleotide likelihood test",
        **kwargs,
    )

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = report_evaluator(tree, records, **kwargs)

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_nucleotide_analysis_keeps_partition_weights_and_policy() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_selected_nucleotide_likelihood_analysis(
        tree,
        records,
        model_name="jc69",
        owner_name="compiled nucleotide likelihood partition test",
        observation_policy="reject",
        partition_map=(
            LocusPartition(
                name="leading-sites",
                segments=(LocusSegment(start=1, end=2),),
                total_sites=2,
                data_type="DNA",
            ),
            LocusPartition(
                name="trailing-sites",
                segments=(LocusSegment(start=3, end=4),),
                total_sites=2,
                data_type="DNA",
            ),
        ),
        site_weight_by_position={2: 2.0, 4: 0.5},
    )

    assert analysis.partition_count == 2
    assert analysis.missing_data_policy == "reject"
    assert math.isclose(analysis.site_weight_total, 4.5, abs_tol=1e-12)
    assert [partition.partition_name for partition in analysis.partition_contexts] == [
        "leading-sites",
        "trailing-sites",
    ]
