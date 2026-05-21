from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_signal_test,
    summarize_phylogenetic_signal,
    write_phylogenetic_signal_permutation_table,
    write_phylogenetic_signal_summary_table,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_phylogenetic_signal_test_preserves_permutation_rows() -> None:
    report = compute_phylogenetic_signal_test(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        permutations=7,
        seed=11,
    )
    assert report.permutations == 7
    assert len(report.permutation_rows) == 7
    assert report.permutation_rows[0].permutation_index == 1
    assert all(1 <= row.permutation_index <= 7 for row in report.permutation_rows)
    assert (
        sum(1 for row in report.permutation_rows if row.at_or_above_observed)
        == report.permuted_k_at_or_above_observed
    )
    assert math.isclose(
        report.null_distribution_minimum,
        min(row.permuted_k for row in report.permutation_rows),
    )
    assert math.isclose(
        report.null_distribution_mean,
        sum(row.permuted_k for row in report.permutation_rows)
        / len(report.permutation_rows),
    )
    assert math.isclose(
        report.null_distribution_maximum,
        max(row.permuted_k for row in report.permutation_rows),
    )


def test_summarize_phylogenetic_signal_matches_reference_fixture_cases() -> None:
    fixture_payload = json.loads(
        fixture("comparative_reference_validation.json").read_text(encoding="utf-8")
    )
    expected_by_case = {
        observation["case"]: observation
        for observation in fixture_payload["observations"]
        if observation["case"]
        in {"blomberg-k-example-tree", "pagel-lambda-example-tree"}
    }
    report = summarize_phylogenetic_signal(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        permutations=19,
        seed=7,
    )
    blomberg_observation = expected_by_case["blomberg-k-example-tree"]
    lambda_observation = expected_by_case["pagel-lambda-example-tree"]
    blomberg_tolerance = float(blomberg_observation["tolerance"])
    lambda_tolerance = float(lambda_observation["tolerance"])
    assert math.isclose(
        report.blombergs_k.k,
        blomberg_observation["expected_parameters"]["k"],
        rel_tol=blomberg_tolerance,
        abs_tol=blomberg_tolerance,
    )
    assert math.isclose(
        report.pagels_lambda.lambda_value,
        lambda_observation["expected_parameters"]["lambda_value"],
        rel_tol=lambda_tolerance,
        abs_tol=lambda_tolerance,
    )
    assert math.isclose(
        report.pagels_lambda.log_likelihood,
        lambda_observation["expected_parameters"]["log_likelihood"],
        rel_tol=lambda_tolerance,
        abs_tol=lambda_tolerance,
    )
    assert report.lambda_likelihood_ratio_statistic >= 0.0
    assert 0.0 <= report.lambda_likelihood_ratio_p_value <= 1.0


def test_write_phylogenetic_signal_tables_write_summary_and_permutations(
    tmp_path: Path,
) -> None:
    report = summarize_phylogenetic_signal(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        permutations=5,
        seed=3,
    )
    summary_path = tmp_path / "phylogenetic-signal-summary.tsv"
    permutation_path = tmp_path / "phylogenetic-signal-permutations.tsv"
    write_phylogenetic_signal_summary_table(summary_path, report)
    write_phylogenetic_signal_permutation_table(permutation_path, report)
    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    permutation_rows = permutation_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_count\tblombergs_k")
    assert len(summary_rows) == 2
    assert permutation_rows[0].startswith(
        "trait\tobserved_k\testimated_lambda\tpermutations"
    )
    assert len(permutation_rows) == 6
