from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    export_continuous_mode_recovery_panel_dataset,
    load_continuous_mode_recovery_panel_dataset,
    run_continuous_mode_recovery_panel_demo,
    run_continuous_mode_recovery_panel_workflow,
    write_continuous_mode_recovery_panel_workflow_bundle,
)


def test_load_continuous_mode_recovery_panel_dataset_exposes_packaged_surface() -> None:
    dataset = load_continuous_mode_recovery_panel_dataset()
    assert dataset.dataset_id == "continuous_mode_recovery_panel"
    assert dataset.label == "Continuous trait-model recovery panel"
    assert dataset.taxon_count == 12
    assert dataset.case_count == 4
    assert dataset.reference_tree_path.is_file()
    assert dataset.simulation_cases_path.is_file()
    assert dataset.reference_output_root.is_dir()


@pytest.mark.slow
def test_write_continuous_mode_recovery_panel_workflow_bundle_matches_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_continuous_mode_recovery_panel_workflow()
    bundle = write_continuous_mode_recovery_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated_paths = sorted(path for path in bundle.output_root.rglob("*") if path.is_file())
    expected_paths = sorted(path for path in expected_root.rglob("*") if path.is_file())
    assert [path.relative_to(bundle.output_root) for path in generated_paths] == [
        path.relative_to(expected_root) for path in expected_paths
    ]
    for generated_path in generated_paths:
        relative = generated_path.relative_to(bundle.output_root)
        assert generated_path.read_text(encoding="utf-8") == (
            expected_root / relative
        ).read_text(encoding="utf-8")


@pytest.mark.slow
def test_run_continuous_mode_recovery_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_continuous_mode_recovery_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 12
    assert result.dataset_export.reference_tree_path.is_file()
    assert result.dataset_export.simulation_cases_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.recovery_summary_path.is_file()
    assert result.workflow_bundle.parameter_recovery_path.is_file()
    assert result.workflow_bundle.model_choice_path.is_file()
    assert result.workflow_bundle.warning_review_path.is_file()
    assert result.workflow_bundle.simulated_traits_root.is_dir()
    assert result.overview_path.is_file()
    assert "parameter recoveries within tolerance" in result.overview_path.read_text(
        encoding="utf-8"
    )


def test_export_continuous_mode_recovery_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_continuous_mode_recovery_panel_dataset(tmp_path / "dataset")
    expected_files = sorted(
        path.relative_to(result.expected_output_root)
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    )
    assert result.readme_path.is_file()
    assert result.reference_tree_path.is_file()
    assert result.simulation_cases_path.is_file()
    assert len(expected_files) == 9
    assert Path("workflow-summary.tsv") in expected_files
    assert Path("simulated-traits/brownian-sigma-recovery.tsv") in expected_files
