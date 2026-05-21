from __future__ import annotations

from bijux_phylogenetics.simulation import validate_geiger_sim_char_reference_examples
from bijux_phylogenetics.simulation.reference.geiger_sim_char import (
    GEIGER_SIM_CHAR_REFERENCE_PAYLOADS,
)


def test_validate_geiger_sim_char_reference_examples_passes() -> None:
    report = validate_geiger_sim_char_reference_examples()

    assert report.all_passed is True
    assert report.case_count == 3
    assert [row.case_id for row in report.observations] == [
        "brownian-internal-long-branch-envelope",
        "speciational-internal-long-branch-envelope",
        "discrete-rate-matrix-internal-long-branch-envelope",
    ]


def test_validate_geiger_sim_char_reference_examples_tracks_expected_envelopes() -> (
    None
):
    report = validate_geiger_sim_char_reference_examples()
    observed_by_case = {row.case_id: row for row in report.observations}

    brownian = observed_by_case["brownian-internal-long-branch-envelope"]
    assert (
        brownian.expected_metrics
        == GEIGER_SIM_CHAR_REFERENCE_PAYLOADS["brownian-internal-long-branch-envelope"][
            "expected_metrics"
        ]
    )
    assert (
        brownian.observed_metrics["tip_distribution"]["A"]["standard_deviation"]
        > brownian.observed_metrics["tip_distribution"]["C"]["standard_deviation"]
    )

    speciational = observed_by_case["speciational-internal-long-branch-envelope"]
    assert (
        abs(
            speciational.observed_metrics["tip_distribution"]["A"]["standard_deviation"]
            - speciational.observed_metrics["tip_distribution"]["C"][
                "standard_deviation"
            ]
        )
        < 0.15
    )

    discrete = observed_by_case["discrete-rate-matrix-internal-long-branch-envelope"]
    assert discrete.observed_metrics["tip_state_frequency"]["A:1"] > 0.3
    assert discrete.observed_metrics["tip_state_frequency"]["C:1"] < 0.1
