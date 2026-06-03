from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.macroevolution_recovery_suite as macroevolution_recovery_suite_api
from bijux_phylogenetics.datasets.macroevolution_recovery_suite import (
    export_macroevolution_recovery_suite_dataset,
    load_macroevolution_recovery_suite_dataset,
    run_macroevolution_recovery_suite_demo,
    run_macroevolution_recovery_suite_workflow,
    write_macroevolution_recovery_suite_workflow_bundle,
)

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_macroevolution_recovery_suite_dataset_exposes_packaged_surface() -> None:
    dataset = load_macroevolution_recovery_suite_dataset()

    assert dataset.dataset_id == "macroevolution_recovery_suite"
    assert dataset.label == "Macroevolution recovery suite"
    assert dataset.component_count == 3
    assert dataset.geiger_component_count == 2
    assert dataset.max_taxon_count == 24
    assert dataset.total_recovery_case_count == 22
    assert dataset.geiger_recovery_case_count == 11
    assert dataset.truth_threshold_row_count == 11
    assert dataset.reference_output_root.is_dir()
    assert dataset.continuous_panel.dataset_id == "continuous_mode_recovery_panel"
    assert dataset.discrete_panel.dataset_id == "discrete_mode_recovery_panel"
    assert dataset.known_answer_panel.dataset_id == "known_answer_reference_panel"


def test_export_macroevolution_recovery_suite_dataset_copies_component_exports_and_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_macroevolution_recovery_suite_dataset(tmp_path / "dataset")

    assert result.readme_path.is_file()
    assert result.component_root.is_dir()
    assert result.continuous_panel_export.readme_path.is_file()
    assert result.discrete_panel_export.readme_path.is_file()
    assert result.known_answer_panel_export.readme_path.is_file()
    assert result.expected_output_root.is_dir()
    expected_files = {
        path.relative_to(result.expected_output_root)
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    }
    assert Path("workflow-summary.tsv") in expected_files
    assert Path("component-summary.tsv") in expected_files
    assert Path("requirement-summary.tsv") in expected_files
    assert Path("sim-char-summary.tsv") in expected_files
    assert (
        Path("components/continuous-mode-recovery-panel/workflow-summary.tsv")
        in expected_files
    )
    assert (
        Path("components/discrete-mode-recovery-panel/workflow-summary.tsv")
        in expected_files
    )
    assert (
        Path("components/known-answer-reference-panel/workflow-summary.tsv")
        in expected_files
    )


def test_public_runtime_exports_include_macroevolution_recovery_suite() -> None:
    assert (
        macroevolution_recovery_suite_api.load_macroevolution_recovery_suite_dataset
        is load_macroevolution_recovery_suite_dataset
    )
    assert (
        macroevolution_recovery_suite_api.export_macroevolution_recovery_suite_dataset
        is export_macroevolution_recovery_suite_dataset
    )
    assert (
        macroevolution_recovery_suite_api.run_macroevolution_recovery_suite_workflow
        is run_macroevolution_recovery_suite_workflow
    )
    assert (
        macroevolution_recovery_suite_api.write_macroevolution_recovery_suite_workflow_bundle
        is write_macroevolution_recovery_suite_workflow_bundle
    )
    assert (
        macroevolution_recovery_suite_api.run_macroevolution_recovery_suite_demo
        is run_macroevolution_recovery_suite_demo
    )


@pytest.mark.slow
def test_write_macroevolution_recovery_suite_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_macroevolution_recovery_suite_workflow()
    bundle = write_macroevolution_recovery_suite_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated_paths = sorted(
        path for path in bundle.output_root.rglob("*") if path.is_file()
    )
    expected_paths = sorted(path for path in expected_root.rglob("*") if path.is_file())
    assert [path.relative_to(bundle.output_root) for path in generated_paths] == [
        path.relative_to(expected_root) for path in expected_paths
    ]
    generated = {path.relative_to(bundle.output_root): path for path in generated_paths}
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_macroevolution_recovery_suite_workflow_reruns_component_workflows() -> (
    None
):
    report = run_macroevolution_recovery_suite_workflow()

    assert len(report.continuous_panel_workflow.recovery_report.case_reports) == 7
    assert len(report.discrete_panel_workflow.recovery_report.case_reports) == 4
    assert len(report.known_answer_panel_workflow.threshold_evaluation_rows) == 11
    assert report.continuous_component.selection_match_count == 4
    assert report.discrete_component.selection_match_count == 2
    assert report.known_answer_component.truth_threshold_pass_count == 10


@pytest.mark.slow
def test_write_macroevolution_recovery_suite_workflow_bundle_records_explicit_goal_300_coverage(
    tmp_path: Path,
) -> None:
    report = run_macroevolution_recovery_suite_workflow()
    bundle = write_macroevolution_recovery_suite_workflow_bundle(
        tmp_path / "workflow",
        report,
    )

    requirement_lines = bundle.requirement_summary_path.read_text(
        encoding="utf-8"
    ).splitlines()
    assert any(
        "brownian, early-burst, ornstein-uhlenbeck" in line
        for line in requirement_lines
    )
    assert any(
        "lambda-transformed-branch-review, kappa-transformed-branch-review, delta-transformed-branch-review"
        in line
        for line in requirement_lines
    )
    assert any(
        "all-rates-different, equal-rates, symmetric" in line
        for line in requirement_lines
    )
    assert any(
        "true-parameters.tsv" in line
        and "continuous_root_state" in line
        and "ou_sigma_squared" in line
        for line in requirement_lines
    )


@pytest.mark.slow
def test_run_macroevolution_recovery_suite_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_macroevolution_recovery_suite_demo(tmp_path / "demo")

    assert result.dataset.component_count == 3
    assert result.dataset_export.readme_path.is_file()
    assert result.dataset_export.continuous_panel_export.readme_path.is_file()
    assert result.dataset_export.discrete_panel_export.readme_path.is_file()
    assert result.dataset_export.known_answer_panel_export.readme_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.component_summary_path.is_file()
    assert result.workflow_bundle.requirement_summary_path.is_file()
    assert result.workflow_bundle.sim_char_summary_path.is_file()
    assert result.workflow_bundle.component_root.is_dir()
    assert result.overview_path.is_file()
    overview_text = result.overview_path.read_text(encoding="utf-8")
    assert "geiger-backed recovery panels" in overview_text
    assert "known-answer truth thresholds passed" in overview_text


@pytest.mark.slow
def test_cli_demo_macroevolution_recovery_suite_json_output_reports_metrics(
    tmp_path: Path,
    capsys,
) -> None:
    output = tmp_path / "macroevolution-suite-demo"
    exit_code = main(
        [
            "demo",
            "macroevolution-recovery-suite",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 6
    assert payload["metrics"]["component_count"] == 3
    assert payload["metrics"]["geiger_component_count"] == 2
    assert payload["metrics"]["case_count"] == 22
    assert payload["metrics"]["geiger_case_count"] == 11
    assert payload["metrics"]["max_taxon_count"] == 24
    assert payload["metrics"]["selection_review_case_count"] == 6
    assert payload["metrics"]["selection_match_count"] == 6
    assert payload["metrics"]["geiger_selection_match_count"] == 5
    assert payload["metrics"]["governed_value_pass_count"] == 42
    assert payload["metrics"]["governed_value_row_count"] == 46
    assert payload["metrics"]["governed_comparison_row_count"] == 23
    assert payload["metrics"]["truth_threshold_pass_count"] == 10
    assert payload["metrics"]["truth_threshold_row_count"] == 11
    assert payload["metrics"]["sim_char_case_count"] == 3
    assert payload["metrics"]["sim_char_all_passed"] is True
    assert payload["metrics"]["requirement_pass_count"] == 11
    assert payload["metrics"]["requirement_row_count"] == 11
    assert payload["data"]["dataset"]["dataset_id"] == "macroevolution_recovery_suite"
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
