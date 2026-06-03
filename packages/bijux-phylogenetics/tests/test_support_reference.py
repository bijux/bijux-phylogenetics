from __future__ import annotations

from bijux_phylogenetics.compare.reference import (
    validate_support_reference_examples,
)


def test_validate_support_reference_examples_reports_passing_cases() -> None:
    report = validate_support_reference_examples()

    assert report.case_count == 6
    assert report.reference_case_count == 4
    assert report.policy_case_count == 2
    assert report.all_passed is True
    assert {observation.category for observation in report.observations} == {
        "reference-pair",
        "clade-mapping-policy",
        "topology-mismatch-policy",
    }


def test_validate_support_reference_examples_keeps_clade_mapping_branch_accurate() -> (
    None
):
    report = validate_support_reference_examples()
    observations = {row.case_id: row for row in report.observations}

    rotation = observations["clade-mapped-support-rotation"]
    bootstrap = observations["iqtree-ufboot-reference"]
    fasttree = observations["fasttree-local-support-reference"]

    assert rotation.passed is True
    assert rotation.observed_metrics["conflicting_clade_count"] == 0
    assert (
        rotation.observed_metrics["shared_support_pairs"]["A|B"]["left_support"] == 97.0
    )
    assert bootstrap.observed_metrics["support_by_clade"] == {
        "A|B": 95.0,
        "C|D": 88.0,
    }
    assert fasttree.observed_metrics["support_by_clade"] == {
        "A|B": 0.91,
        "C|D": 0.73,
    }


def test_validate_support_reference_examples_flags_topology_mismatch_explicitly() -> (
    None
):
    report = validate_support_reference_examples()
    observations = {row.case_id: row for row in report.observations}

    mismatch = observations["bootstrap-posterior-topology-mismatch"]
    posterior = observations["posterior-clade-frequency-reference"]

    assert mismatch.passed is True
    assert mismatch.observed_metrics["topology_mismatch_detected"] is True
    assert mismatch.observed_metrics["topology_mismatch_clade_count"] == 2
    assert mismatch.observed_metrics["bootstrap_support_by_clade"] == {
        "A|B": 0.97,
        "C|D": 0.96,
    }
    assert posterior.observed_metrics["frequency_by_clade"]["A|B"]["frequency"] == (
        0.666666666666667
    )
