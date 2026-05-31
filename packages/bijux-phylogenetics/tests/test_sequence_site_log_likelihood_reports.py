from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment,
    evaluate_protein_poisson_tree_likelihood_from_alignment,
    fit_local_clock_likelihood_from_alignment,
    fit_strict_clock_likelihood_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
    write_local_clock_likelihood_artifacts,
    write_site_log_likelihood_table,
    write_strict_clock_likelihood_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_site_log_likelihood_table_supports_fixed_protein_models(
    tmp_path: Path,
) -> None:
    poisson_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"),
    )
    empirical_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
    )

    poisson_path = write_site_log_likelihood_table(
        tmp_path / "poisson-site-log-likelihoods.tsv",
        poisson_report,
    )
    empirical_path = write_site_log_likelihood_table(
        tmp_path / "empirical-site-log-likelihoods.tsv",
        empirical_report,
    )

    assert len(poisson_report.site_log_likelihoods) == poisson_report.site_count
    assert math.isclose(
        _expanded_site_log_likelihood_sum(poisson_path),
        poisson_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert len(empirical_report.site_log_likelihoods) == empirical_report.site_count
    assert math.isclose(
        _expanded_site_log_likelihood_sum(empirical_path),
        empirical_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_write_site_log_likelihood_table_supports_protein_mixture_models(
    tmp_path: Path,
) -> None:
    gamma_report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments",
                "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=0.8,
            category_count=4,
            matrix_label="compact-polar",
        )
    )
    invariant_report = optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        initial_invariant_proportion=0.1,
    )
    gamma_invariant_report = optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        alpha=0.8,
        category_count=4,
        matrix_label="compact-polar",
        initial_invariant_proportion=0.1,
    )

    for output_name, report in (
        ("gamma-site-log-likelihoods.tsv", gamma_report),
        ("invariant-site-log-likelihoods.tsv", invariant_report),
        ("gamma-invariant-site-log-likelihoods.tsv", gamma_invariant_report),
    ):
        output_path = write_site_log_likelihood_table(tmp_path / output_name, report)
        assert math.isclose(
            _expanded_site_log_likelihood_sum(output_path),
            report.log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )


def test_strict_clock_artifacts_include_site_log_likelihood_table(
    tmp_path: Path,
) -> None:
    report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )

    outputs = write_strict_clock_likelihood_artifacts(tmp_path, report)

    assert len(report.site_log_likelihoods) == report.site_count
    assert outputs["site_log_likelihood_path"].name == "site_log_likelihoods.tsv"
    assert math.isclose(
        _expanded_site_log_likelihood_sum(outputs["site_log_likelihood_path"]),
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        _weighted_pattern_log_likelihood_sum(outputs["site_log_likelihood_path"]),
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_local_clock_artifacts_include_site_log_likelihood_table(
    tmp_path: Path,
) -> None:
    report = fit_local_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "local_clock_likelihood_alignment_4_taxa.fasta"),
        fixture("metadata", "local_clock_regimes_4_taxa.tsv"),
    )

    outputs = write_local_clock_likelihood_artifacts(tmp_path, report)

    assert len(report.site_log_likelihoods) == report.site_count
    assert outputs["site_log_likelihood_path"].name == "site_log_likelihoods.tsv"
    assert math.isclose(
        _expanded_site_log_likelihood_sum(outputs["site_log_likelihood_path"]),
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        _weighted_pattern_log_likelihood_sum(outputs["site_log_likelihood_path"]),
        report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _expanded_site_log_likelihood_sum(path: Path) -> float:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return sum(float(row["log_likelihood"]) for row in rows)


def _weighted_pattern_log_likelihood_sum(path: Path) -> float:
    pattern_rows: dict[str, tuple[int, float]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pattern_rows.setdefault(
                row["pattern_id"],
                (int(row["pattern_weight"]), float(row["log_likelihood"])),
            )
    return sum(
        weight * log_likelihood for weight, log_likelihood in pattern_rows.values()
    )


def _compact_polar_rate_matrix() -> numpy.ndarray:
    return _build_empirical_rate_matrix(
        boosted_pairs={
            ("A", "C"): 0.45,
            ("C", "D"): 0.35,
            ("D", "E"): 0.55,
            ("A", "E"): 0.20,
        }
    )


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = _protein_state_index()
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _build_empirical_rate_matrix(
    *,
    boosted_pairs: dict[tuple[str, str], float],
) -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = _protein_state_index()
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _protein_state_order() -> tuple[str, ...]:
    return (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )


def _protein_state_index() -> dict[str, int]:
    return {state: index for index, state in enumerate(_protein_state_order())}
