from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.beast import (
    detect_impossible_calibration_constraints,
    prepare_beast_time_tree_analysis,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
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


def test_validate_fossil_calibration_table_accepts_named_and_taxon_targets() -> None:
    report = validate_fossil_calibration_table(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_calibrations.tsv"),
    )

    assert report.calibration_count == 2
    assert report.valid_calibration_count == 2
    assert report.invalid_calibration_count == 0
    assert [calibration.target_kind for calibration in report.calibrations] == ["named-clade", "taxa"]
    assert report.calibrations[0].taxa == ["A", "B"]
    assert report.calibrations[1].taxa == ["C", "D"]


def test_detect_impossible_calibration_constraints_reports_unknown_and_invalid_targets() -> None:
    report = detect_impossible_calibration_constraints(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_calibrations_invalid.tsv"),
    )

    assert report.impossible_calibration_ids == ["bad-clade", "bad-empty", "bad-order"]
    assert {issue.code for issue in report.issues} >= {
        "unknown-clade-name",
        "non-monophyletic-target",
        "minimum-exceeds-maximum",
        "missing-target",
        "missing-age-bounds",
    }


def test_validate_tip_dating_metadata_checks_tree_and_alignment_membership() -> None:
    report = validate_tip_dating_metadata(
        fixture("example_tree_named_clades.nwk"),
        fixture("example_tip_dates_invalid.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
    )

    assert report.valid_tip_count == 1
    assert report.invalid_tip_count == 2
    assert report.extra_tip_taxa == ["E"]
    assert report.extra_alignment_taxa == ["E"]
    assert {issue.code for issue in report.issues} >= {"taxon-missing-from-tree", "taxon-missing-from-alignment", "invalid-date"}


def test_prepare_beast_time_tree_analysis_writes_clock_prior_calibrations_and_tip_dates(tmp_path: Path) -> None:
    output_path = tmp_path / "analysis.xml"

    report = prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        output_path,
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        clock_model="relaxed-lognormal",
        tree_prior="birth-death",
        chain_length=200000,
        log_every=500,
    )

    text = output_path.read_text(encoding="utf-8")
    assert report.taxon_count == 4
    assert report.character_count == 8
    assert report.calibration_count == 2
    assert report.tip_date_count == 4
    assert '<clockModel name="relaxed-lognormal" />' in text
    assert '<treePrior name="birth-death" />' in text
    assert '<calibration id="cal-mammals"' in text
    assert '<date taxon="A" value="2012.5" />' in text
