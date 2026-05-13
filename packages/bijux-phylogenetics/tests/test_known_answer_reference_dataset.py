from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.known_answer_reference import (
    export_known_answer_reference_dataset,
    load_known_answer_reference_dataset,
    run_known_answer_reference_demo,
    run_known_answer_reference_workflow,
    write_known_answer_reference_workflow_bundle,
)


def test_load_known_answer_reference_dataset_exposes_packaged_surface() -> None:
    dataset = load_known_answer_reference_dataset()
    assert dataset.dataset_id == "known_answer_reference_panel"
    assert dataset.label == "Known-answer simulation reference panel"
    assert dataset.taxon_count == 8
    assert dataset.sequence_length == 5000
    assert dataset.sequence_type == "dna"
    assert dataset.distance_method == "neighbor-joining"
    assert dataset.distance_model == "p-distance"
    assert dataset.true_tree_path.is_file()
    assert dataset.alignment_path.is_file()
    assert dataset.continuous_traits_path.is_file()
    assert dataset.discrete_traits_path.is_file()
    assert dataset.true_parameters_path.is_file()
    assert dataset.true_continuous_nodes_path.is_file()
    assert dataset.true_discrete_nodes_path.is_file()
    assert dataset.reference_output_root.is_dir()


def test_write_known_answer_reference_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_known_answer_reference_workflow()
    bundle = write_known_answer_reference_workflow_bundle(tmp_path / "workflow", report)
    expected_root = report.dataset.reference_output_root
    generated = {
        bundle.workflow_summary_path.name: bundle.workflow_summary_path,
        bundle.distance_tree_path.name: bundle.distance_tree_path,
        bundle.tree_recovery_path.name: bundle.tree_recovery_path,
        bundle.parameter_recovery_path.name: bundle.parameter_recovery_path,
        bundle.brownian_fit_summary_path.name: bundle.brownian_fit_summary_path,
        bundle.continuous_ancestral_summary_path.name: (
            bundle.continuous_ancestral_summary_path
        ),
        bundle.continuous_ancestral_uncertainty_path.name: (
            bundle.continuous_ancestral_uncertainty_path
        ),
        bundle.continuous_node_recovery_path.name: bundle.continuous_node_recovery_path,
        bundle.discrete_ancestral_summary_path.name: (
            bundle.discrete_ancestral_summary_path
        ),
        bundle.discrete_ancestral_probability_path.name: (
            bundle.discrete_ancestral_probability_path
        ),
        bundle.discrete_node_recovery_path.name: bundle.discrete_node_recovery_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    for name, generated_path in generated.items():
        assert generated_path.read_text(encoding="utf-8") == (
            expected_root / name
        ).read_text(encoding="utf-8")


def test_run_known_answer_reference_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_known_answer_reference_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 8
    assert result.dataset_export.true_tree_path.is_file()
    assert result.dataset_export.alignment_path.is_file()
    assert result.dataset_export.true_parameters_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.distance_tree_path.is_file()
    assert result.workflow_bundle.parameter_recovery_path.is_file()
    assert result.workflow_bundle.continuous_node_recovery_path.is_file()
    assert result.workflow_bundle.discrete_node_recovery_path.is_file()
    assert result.overview_path.is_file()
    assert "unrooted topology" in result.overview_path.read_text(encoding="utf-8")


def test_export_known_answer_reference_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_known_answer_reference_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.true_tree_path.is_file()
    assert result.alignment_path.is_file()
    assert result.continuous_traits_path.is_file()
    assert result.discrete_traits_path.is_file()
    assert result.true_parameters_path.is_file()
    assert result.true_continuous_nodes_path.is_file()
    assert result.true_discrete_nodes_path.is_file()
    assert len(expected_files) == 11
    assert "workflow-summary.tsv" in expected_files
