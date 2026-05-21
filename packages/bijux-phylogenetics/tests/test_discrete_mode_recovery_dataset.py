from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.discrete_mode_recovery as discrete_mode_recovery_api
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    export_discrete_mode_recovery_panel_dataset,
    load_discrete_mode_recovery_panel_dataset,
    run_discrete_mode_recovery_panel_demo,
    run_discrete_mode_recovery_panel_workflow,
    write_discrete_mode_recovery_panel_workflow_bundle,
)

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def test_load_discrete_mode_recovery_panel_dataset_exposes_packaged_surface() -> None:
    dataset = load_discrete_mode_recovery_panel_dataset()
    assert dataset.dataset_id == "discrete_mode_recovery_panel"
    assert dataset.label == "Discrete trait-model recovery panel"
    assert dataset.taxon_count == 24
    assert dataset.tree_count == 2
    assert dataset.case_count == 4
    assert dataset.default_tree_path.is_file()
    assert len(dataset.reference_tree_paths) == 2
    assert dataset.simulation_cases_path.is_file()
    assert dataset.reference_output_root.is_dir()
    header = dataset.simulation_cases_path.read_text(encoding="utf-8").splitlines()[0]
    assert "transform" in header.split("\t")
    assert "transform_parameter_value" in header.split("\t")
    assert "parameter_tolerances" in header.split("\t")


@pytest.mark.slow
def test_write_discrete_mode_recovery_panel_workflow_bundle_matches_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_discrete_mode_recovery_panel_workflow()
    bundle = write_discrete_mode_recovery_panel_workflow_bundle(
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
def test_run_discrete_mode_recovery_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_discrete_mode_recovery_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 24
    assert result.dataset.tree_count == 2
    assert result.dataset_export.default_tree_path.is_file()
    assert result.dataset_export.reference_tree_root.is_dir()
    assert result.dataset_export.simulation_cases_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.recovery_summary_path.is_file()
    assert result.workflow_bundle.parameter_recovery_path.is_file()
    assert result.workflow_bundle.parameter_comparison_path.is_file()
    assert result.workflow_bundle.rate_recovery_path.is_file()
    assert result.workflow_bundle.rate_comparison_path.is_file()
    assert result.workflow_bundle.model_choice_path.is_file()
    assert result.workflow_bundle.execution_review_path.is_file()
    assert result.workflow_bundle.warning_review_path.is_file()
    assert result.workflow_bundle.geiger_reference_path.is_file()
    assert result.workflow_bundle.simulated_traits_root.is_dir()
    assert result.overview_path.is_file()
    overview_text = result.overview_path.read_text(encoding="utf-8")
    assert "geiger model-selection matches expectation" in overview_text
    assert "paired transform-parameter comparisons" in overview_text
    assert "all paired rate comparisons" in overview_text


def test_export_discrete_mode_recovery_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_discrete_mode_recovery_panel_dataset(tmp_path / "dataset")
    expected_files = sorted(
        path.relative_to(result.expected_output_root)
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    )
    assert result.readme_path.is_file()
    assert result.default_tree_path.is_file()
    assert result.reference_tree_root.is_dir()
    assert result.simulation_cases_path.is_file()
    assert Path("workflow-summary.tsv") in expected_files
    assert Path("parameter-recovery.tsv") in expected_files
    assert Path("parameter-comparison.tsv") in expected_files
    assert Path("rate-comparison.tsv") in expected_files
    assert Path("execution-review.tsv") in expected_files
    assert (
        Path("simulated-traits/ard-three-state-weak-identification-review.tsv")
        in expected_files
    )
    assert Path("geiger-reference.tsv") in expected_files


def test_public_runtime_exports_include_discrete_mode_recovery_panel() -> None:
    assert (
        discrete_mode_recovery_api.load_discrete_mode_recovery_panel_dataset
        is load_discrete_mode_recovery_panel_dataset
    )
    assert (
        discrete_mode_recovery_api.export_discrete_mode_recovery_panel_dataset
        is export_discrete_mode_recovery_panel_dataset
    )
    assert (
        discrete_mode_recovery_api.run_discrete_mode_recovery_panel_workflow
        is run_discrete_mode_recovery_panel_workflow
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_panel_workflow_bundle
        is write_discrete_mode_recovery_panel_workflow_bundle
    )
    assert (
        discrete_mode_recovery_api.run_discrete_mode_recovery_panel_demo
        is run_discrete_mode_recovery_panel_demo
    )


@pytest.mark.slow
def test_cli_demo_discrete_mode_recovery_panel_json_output_reports_metrics(
    tmp_path: Path,
    capsys,
) -> None:
    output = tmp_path / "discrete-mode-demo"
    exit_code = main(
        [
            "demo",
            "discrete-mode-recovery-panel",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 14
    assert payload["metrics"]["taxon_count"] == 24
    assert payload["metrics"]["tree_count"] == 2
    assert payload["metrics"]["case_count"] == 4
    assert payload["metrics"]["parameter_row_count"] == 0
    assert payload["metrics"]["parameter_comparison_row_count"] == 0
    assert payload["data"]["dataset"]["dataset_id"] == "discrete_mode_recovery_panel"
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
