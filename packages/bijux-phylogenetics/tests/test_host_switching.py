from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.ecology import (
    summarize_host_switching,
    write_host_state_node_table,
    write_host_switch_branch_table,
    write_host_switch_count_table,
    write_host_switch_exclusion_table,
    write_host_switch_fit_table,
    write_host_switch_summary_table,
    write_unsupported_host_switch_claim_table,
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


def test_summarize_host_switching_reconstructs_hosts_and_counts_switches() -> None:
    report = summarize_host_switching(
        fixture("example_tree.nwk"),
        fixture("example_traits_host.tsv"),
        trait="host",
        model="er",
    )

    assert report.summary.analysis_constraint_mode == "unconstrained"
    assert report.summary.observed_host_count == 3
    assert report.summary.host_switch_count == 2
    assert len(report.node_rows) == 3
    assert [row.transition for row in report.branch_rows if row.changed] == [
        "human->bat",
        "bat->rodent",
    ]
    assert [row.transition for row in report.count_rows] == [
        "bat->rodent",
        "human->bat",
    ]


@pytest.mark.slow
def test_summarize_host_switching_supports_constrained_transition_models() -> None:
    report = summarize_host_switching(
        fixture("example_tree.nwk"),
        fixture("example_traits_host.tsv"),
        trait="host",
        model="ard",
        constraint_path=fixture("example_host_transition_constraints.tsv"),
    )

    assert report.summary.analysis_constraint_mode == "constrained"
    assert len(report.fit_rows) == 2
    assert report.summary.allowed_transition_count == 4
    assert report.summary.forbidden_transition_count == 2
    assert report.summary.preferred_constraint in {"constrained", "unconstrained"}


@pytest.mark.slow
def test_summarize_host_switching_reports_unsupported_unconstrained_claims(
    tmp_path: Path,
) -> None:
    constraint_path = tmp_path / "constraints.tsv"
    constraint_path.write_text(
        (
            "source_host\ttarget_host\ttransition_allowed\n"
            "human\tbat\tallowed\n"
            "bat\thuman\tallowed\n"
            "human\trodent\tallowed\n"
            "rodent\thuman\tallowed\n"
        ),
        encoding="utf-8",
    )

    report = summarize_host_switching(
        fixture("example_tree.nwk"),
        fixture("example_traits_host.tsv"),
        trait="host",
        model="ard",
        constraint_path=constraint_path,
    )

    assert report.summary.unsupported_switch_claim_count >= 1
    assert any(not row.claim_resolved for row in report.unsupported_claim_rows)


@pytest.mark.slow
def test_summarize_host_switching_rejects_unknown_constraint_host(
    tmp_path: Path,
) -> None:
    constraint_path = tmp_path / "ghost.tsv"
    constraint_path.write_text(
        "source_host\ttarget_host\ttransition_allowed\nghost\thuman\tallowed\n",
        encoding="utf-8",
    )
    with pytest.raises(
        ValueError,
        match="source host is not present in the analyzed host vocabulary",
    ):
        summarize_host_switching(
            fixture("example_tree.nwk"),
            fixture("example_traits_host.tsv"),
            trait="host",
            model="ard",
            constraint_path=constraint_path,
        )


@pytest.mark.slow
def test_write_host_switching_tables_emit_expected_ledgers(tmp_path: Path) -> None:
    report = summarize_host_switching(
        fixture("example_tree.nwk"),
        fixture("example_traits_host.tsv"),
        trait="host",
        model="ard",
        constraint_path=fixture("example_host_transition_constraints.tsv"),
    )

    summary_path = write_host_switch_summary_table(tmp_path / "summary.tsv", report)
    nodes_path = write_host_state_node_table(tmp_path / "nodes.tsv", report)
    branches_path = write_host_switch_branch_table(tmp_path / "branches.tsv", report)
    counts_path = write_host_switch_count_table(tmp_path / "counts.tsv", report)
    fits_path = write_host_switch_fit_table(tmp_path / "fits.tsv", report)
    unsupported_path = write_unsupported_host_switch_claim_table(
        tmp_path / "unsupported.tsv",
        report,
    )
    exclusions_path = write_host_switch_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "preferred_constraint" in summary_path.read_text(encoding="utf-8")
    assert "host_probabilities" in nodes_path.read_text(encoding="utf-8")
    assert "certainty_class" in branches_path.read_text(encoding="utf-8")
    assert "certain_switch_count" in counts_path.read_text(encoding="utf-8")
    assert "constraint_mode" in fits_path.read_text(encoding="utf-8")
    assert "claim_resolved" in unsupported_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
