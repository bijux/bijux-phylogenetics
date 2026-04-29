from __future__ import annotations

from pathlib import Path
import json

from bijux_phylogenetics.bayesian.beast import (
    assess_beast_burnin_sensitivity,
    assess_beast_chain_mixing,
    assess_beast_convergence,
    detect_impossible_calibration_constraints,
    parse_beast_log,
    prepare_beast_time_tree_analysis,
    validate_beast_posterior_log,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.bayesian.evidence import build_bayesian_evidence_package
from bijux_phylogenetics.bayesian.reports import render_calibration_audit_report


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


def test_parse_beast_log_and_assess_convergence_return_parameter_summaries() -> None:
    log_report = parse_beast_log(fixture("example_beast.log"))
    convergence = assess_beast_convergence(
        fixture("example_beast.log"),
        ess_threshold=5.0,
        mean_shift_threshold=0.1,
    )

    assert log_report.row_count == 4
    assert log_report.columns == ["posterior", "likelihood", "clockRate", "treeHeight"]
    assert log_report.rows[1].state == 1000
    assert convergence.converged is False
    assert {warning["code"] for warning in convergence.warnings} == {"low-ess", "mean-drift"}


def test_validate_beast_posterior_log_reports_missing_columns_and_nonmonotonic_states(tmp_path: Path) -> None:
    broken_log = tmp_path / "broken-beast.log"
    broken_log.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tclockRate\n"
        "0\t-510.0\t0.0010\n"
        "0\t-509.0\tbad\n",
        encoding="utf-8",
    )

    report = validate_beast_posterior_log(broken_log)

    assert report.valid is False
    assert report.missing_columns == ["likelihood"]
    assert {issue.code for issue in report.issues} >= {
        "missing-required-column",
        "nonmonotonic-state",
        "invalid-parameter-value",
    }


def test_assess_beast_burnin_sensitivity_reports_tree_and_log_shifts() -> None:
    report = assess_beast_burnin_sensitivity(
        fixture("example_tree_set_left.nwk"),
        log_path=fixture("example_beast.log"),
        burnin_fractions=(0.0, 0.25, 0.5),
    )

    assert [row.burnin_fraction for row in report.slices] == [0.0, 0.25, 0.5]
    assert report.slices[0].posterior_mean == -503.0
    assert report.slices[1].posterior_mean == -500.666667
    assert report.slices[2].tree_height_mean == 13.1


def test_assess_beast_chain_mixing_flags_stuck_and_inconsistent_chains(tmp_path: Path) -> None:
    stable = tmp_path / "stable.log"
    stable.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-510.0\t-490.0\t0.0010\t12.0\n"
        "1000\t-509.0\t-489.0\t0.0011\t12.2\n"
        "2000\t-508.0\t-488.0\t0.0012\t12.3\n"
        "3000\t-507.0\t-487.0\t0.0013\t12.4\n",
        encoding="utf-8",
    )
    stuck = tmp_path / "stuck.log"
    stuck.write_text(
        "# BEAST fixture log\n"
        "state\tposterior\tlikelihood\tclockRate\ttreeHeight\n"
        "0\t-500.0\t-480.0\t0.0010\t15.0\n"
        "1000\t-500.0\t-480.0\t0.0010\t15.0\n"
        "2000\t-500.0\t-480.0\t0.0010\t15.0\n"
        "3000\t-500.0\t-480.0\t0.0010\t15.0\n",
        encoding="utf-8",
    )

    report = assess_beast_chain_mixing(
        [stable, stuck],
        ess_threshold=2.0,
        mean_shift_threshold=0.5,
        cross_chain_mean_shift_threshold=1.0,
    )

    assert report.converged is False
    assert {issue.code for issue in report.issues} >= {"stuck-parameter", "inconsistent-chains"}


def test_render_calibration_audit_report_includes_calibration_and_tip_date_sections(tmp_path: Path) -> None:
    output_path = tmp_path / "calibration-audit.html"

    report = render_calibration_audit_report(
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=output_path,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.output_path == output_path
    assert report.invalid_calibration_count == 0
    assert "fossil-calibrations" in html
    assert "impossible-constraints" in html
    assert "tip-dates" in html


def test_build_bayesian_evidence_package_bundles_inputs_outputs_and_reports(tmp_path: Path) -> None:
    config_path = tmp_path / "analysis.xml"
    report_path = tmp_path / "calibration-audit.html"
    diagnostic_path = tmp_path / "diagnostics.json"
    prepare_beast_time_tree_analysis(
        fixture("example_alignment.fasta"),
        config_path,
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
    )
    render_calibration_audit_report(
        tree_path=fixture("example_tree_named_clades.nwk"),
        calibration_path=fixture("example_calibrations.tsv"),
        tip_dates_path=fixture("example_tip_dates.tsv"),
        alignment_path=fixture("example_alignment.fasta"),
        out_path=report_path,
    )
    diagnostic_path.write_text(json.dumps({"warning_count": 0}, indent=2) + "\n", encoding="utf-8")

    bundle = build_bayesian_evidence_package(
        bundle_root=tmp_path / "bayesian-bundle",
        input_paths=[fixture("example_alignment.fasta"), fixture("example_calibrations.tsv"), fixture("example_tip_dates.tsv")],
        config_paths=[config_path],
        tree_paths=[fixture("example_tree_named_clades.nwk")],
        log_paths=[fixture("example_beast.log")],
        diagnostic_paths=[diagnostic_path],
        report_paths=[report_path],
    )

    assert bundle.valid is True
    assert bundle.file_count == 8
    assert bundle.config_count == 1
    assert bundle.tree_count == 1
    assert bundle.log_count == 1
    assert bundle.report_count == 1
