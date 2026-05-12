from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.biogeography import (
    summarize_constrained_geographic_model,
    summarize_constrained_geographic_report,
    write_constrained_geographic_exclusion_table,
    write_constrained_geographic_fit_table,
    write_constrained_geographic_summary_table,
    write_constrained_geographic_transition_table,
    write_unsupported_geographic_transition_claim_table,
)
from bijux_phylogenetics.errors import AncestralReconstructionError

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


def test_summarize_constrained_geographic_model_supports_aliases() -> None:
    for model in ("er", "sym", "ard"):
        report = summarize_constrained_geographic_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            fixture("example_geographic_adjacency.tsv"),
            trait="region",
            model=model,
        )

        assert report.model == model
        assert report.fit_rows
        assert report.transition_rows


def test_summarize_constrained_geographic_model_reports_unsupported_transition_claims(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "simulated-geography.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tnorth\nB\tnorth\nC\tisland\nD\tsouth\n",
        encoding="utf-8",
    )

    report = summarize_constrained_geographic_model(
        fixture("example_tree.nwk"),
        traits_path,
        fixture("example_geographic_adjacency.tsv"),
        trait="region",
        model="ard",
    )
    summary = summarize_constrained_geographic_report(report)

    assert summary.unsupported_transition_claim_count == len(
        report.unsupported_claim_rows
    )
    assert report.unsupported_claim_rows
    assert any(
        {
            row.unconstrained_source_region,
            row.unconstrained_target_region,
        }
        == {"north", "island"}
        for row in report.unsupported_claim_rows
    )


def test_summarize_constrained_geographic_model_rejects_invalid_adjacency_matrix(
    tmp_path: Path,
) -> None:
    adjacency_path = tmp_path / "bad-adjacency.tsv"
    adjacency_path.write_text(
        "region\tnorth\tsouth\nnorth\t0\t1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing rows"):
        summarize_constrained_geographic_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            adjacency_path,
            trait="region",
            model="ard",
        )


def test_write_constrained_geographic_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_constrained_geographic_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        fixture("example_geographic_adjacency.tsv"),
        trait="region",
        model="ard",
    )

    summary_path = write_constrained_geographic_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    fits_path = write_constrained_geographic_fit_table(
        tmp_path / "fits.tsv",
        report,
    )
    transitions_path = write_constrained_geographic_transition_table(
        tmp_path / "transitions.tsv",
        report,
    )
    unsupported_path = write_unsupported_geographic_transition_claim_table(
        tmp_path / "unsupported.tsv",
        report,
    )
    exclusions_path = write_constrained_geographic_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "preferred_constraint" in summary_path.read_text(encoding="utf-8")
    assert "constraint_mode" in fits_path.read_text(encoding="utf-8")
    assert "transition_allowed" in transitions_path.read_text(encoding="utf-8")
    assert "claim_resolved" in unsupported_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_summarize_constrained_geographic_model_rejects_in_tree_unknown_regions(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "bad-geography.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tnorth\nB\tlagoon\nC\tsouth\nD\tisland\n",
        encoding="utf-8",
    )

    with pytest.raises(
        AncestralReconstructionError,
        match="does not define one or more analyzed region labels",
    ):
        summarize_constrained_geographic_model(
            fixture("example_tree.nwk"),
            traits_path,
            fixture("example_geographic_adjacency.tsv"),
            trait="region",
            model="ard",
        )
