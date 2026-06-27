from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("model_name", "kwargs", "report_evaluator"),
    [
        ("jc69", {}, evaluate_jc69_tree_likelihood),
        ("k80", {"kappa": 4.0}, evaluate_k80_tree_likelihood),
        (
            "f81",
            {"base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}},
            evaluate_f81_tree_likelihood,
        ),
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
def test_selected_nucleotide_surface_matches_existing_tree_reports(
    model_name: str,
    kwargs: dict[str, object],
    report_evaluator,
) -> None:
    tree = load_tree(fixture("trees", f"{model_name}_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", f"{model_name}_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name="finite-state nucleotide protocol test",
        **kwargs,
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = report_evaluator(tree, records, **kwargs)

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
