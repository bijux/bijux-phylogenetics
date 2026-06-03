from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.known_answer_reference as known_answer_reference_api
from bijux_phylogenetics.datasets.known_answer_reference import (
    export_known_answer_reference_dataset,
    load_known_answer_reference_dataset,
    run_known_answer_reference_demo,
    run_known_answer_reference_workflow,
    write_known_answer_reference_workflow_bundle,
)

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
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
    assert dataset.ou_traits_path.is_file()
    assert dataset.discrete_traits_path.is_file()
    assert dataset.host_traits_path.is_file()
    assert dataset.geographic_traits_path.is_file()
    assert dataset.true_parameters_path.is_file()
    assert dataset.true_continuous_nodes_path.is_file()
    assert dataset.true_ou_nodes_path.is_file()
    assert dataset.true_discrete_nodes_path.is_file()
    assert dataset.true_host_nodes_path.is_file()
    assert dataset.true_geographic_nodes_path.is_file()
    assert dataset.true_host_switch_events_path.is_file()
    assert dataset.true_geographic_transition_events_path.is_file()
    assert dataset.recovery_thresholds_path.is_file()
    assert dataset.reference_output_root.is_dir()


@pytest.mark.slow
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
        bundle.ou_fit_summary_path.name: bundle.ou_fit_summary_path,
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
        bundle.host_switch_summary_path.name: bundle.host_switch_summary_path,
        bundle.host_state_nodes_path.name: bundle.host_state_nodes_path,
        bundle.host_switch_branches_path.name: bundle.host_switch_branches_path,
        bundle.host_node_recovery_path.name: bundle.host_node_recovery_path,
        bundle.host_event_recovery_path.name: bundle.host_event_recovery_path,
        bundle.geographic_ancestral_summary_path.name: (
            bundle.geographic_ancestral_summary_path
        ),
        bundle.geographic_state_probability_path.name: (
            bundle.geographic_state_probability_path
        ),
        bundle.geographic_transition_summary_path.name: (
            bundle.geographic_transition_summary_path
        ),
        bundle.geographic_node_recovery_path.name: (
            bundle.geographic_node_recovery_path
        ),
        bundle.geographic_event_recovery_path.name: (
            bundle.geographic_event_recovery_path
        ),
        bundle.threshold_evaluation_path.name: bundle.threshold_evaluation_path,
    }
    assert {path.name for path in expected_root.glob("*")} == set(generated)
    assert_selected_scientific_outputs_equivalent(expected_root, generated)


@pytest.mark.slow
def test_run_known_answer_reference_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_known_answer_reference_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 8
    assert result.dataset_export.true_tree_path.is_file()
    assert result.dataset_export.alignment_path.is_file()
    assert result.dataset_export.ou_traits_path.is_file()
    assert result.dataset_export.host_traits_path.is_file()
    assert result.dataset_export.geographic_traits_path.is_file()
    assert result.dataset_export.true_parameters_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.distance_tree_path.is_file()
    assert result.workflow_bundle.parameter_recovery_path.is_file()
    assert result.workflow_bundle.continuous_node_recovery_path.is_file()
    assert result.workflow_bundle.discrete_node_recovery_path.is_file()
    assert result.workflow_bundle.ou_fit_summary_path.is_file()
    assert result.workflow_bundle.host_event_recovery_path.is_file()
    assert result.workflow_bundle.geographic_event_recovery_path.is_file()
    assert result.workflow_bundle.threshold_evaluation_path.is_file()
    assert result.overview_path.is_file()
    overview = result.overview_path.read_text(encoding="utf-8")
    assert "unrooted topology" in overview
    assert "host event accuracy" in overview
    assert "threshold passes" in overview


def test_export_known_answer_reference_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_known_answer_reference_dataset(tmp_path / "dataset")
    expected_files = {path.name for path in result.expected_output_root.glob("*")}
    assert result.readme_path.is_file()
    assert result.true_tree_path.is_file()
    assert result.alignment_path.is_file()
    assert result.continuous_traits_path.is_file()
    assert result.ou_traits_path.is_file()
    assert result.discrete_traits_path.is_file()
    assert result.host_traits_path.is_file()
    assert result.geographic_traits_path.is_file()
    assert result.true_parameters_path.is_file()
    assert result.true_continuous_nodes_path.is_file()
    assert result.true_ou_nodes_path.is_file()
    assert result.true_discrete_nodes_path.is_file()
    assert result.true_host_nodes_path.is_file()
    assert result.true_geographic_nodes_path.is_file()
    assert result.true_host_switch_events_path.is_file()
    assert result.true_geographic_transition_events_path.is_file()
    assert result.recovery_thresholds_path.is_file()
    assert len(expected_files) == 23
    assert "workflow-summary.tsv" in expected_files
    assert "host-event-recovery.tsv" in expected_files
    assert "geographic-event-recovery.tsv" in expected_files
    assert "recovery-threshold-evaluation.tsv" in expected_files


@pytest.mark.slow
def test_public_runtime_exports_include_known_answer_reference_surface() -> None:
    assert (
        known_answer_reference_api.load_known_answer_reference_dataset
        is load_known_answer_reference_dataset
    )
    assert (
        known_answer_reference_api.export_known_answer_reference_dataset
        is export_known_answer_reference_dataset
    )
    assert (
        known_answer_reference_api.run_known_answer_reference_workflow
        is run_known_answer_reference_workflow
    )
    assert (
        known_answer_reference_api.write_known_answer_reference_workflow_bundle
        is write_known_answer_reference_workflow_bundle
    )
    assert (
        known_answer_reference_api.run_known_answer_reference_demo
        is run_known_answer_reference_demo
    )


@pytest.mark.slow
def test_cli_demo_known_answer_reference_panel_json_output_reports_recovery_metrics(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "known-answer-demo"
    exit_code = main(
        [
            "demo",
            "known-answer-reference-panel",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 41
    assert payload["metrics"]["taxon_count"] == 8
    assert payload["metrics"]["sequence_length"] == 5000
    assert payload["metrics"]["distance_method"] == "neighbor-joining"
    assert payload["metrics"]["distance_model"] == "p-distance"
    assert payload["metrics"]["rooted_topology_equal"] is False
    assert payload["metrics"]["same_unrooted_topology"] is True
    assert payload["metrics"]["same_taxa_different_rooting"] is True
    assert payload["metrics"]["robinson_foulds_distance"] == 3
    assert payload["metrics"]["parameter_row_count"] == 5
    assert payload["metrics"]["threshold_pass_count"] == 10
    assert payload["metrics"]["threshold_row_count"] == 11
    assert payload["metrics"]["discrete_internal_node_accuracy"] == 1.0
    assert payload["metrics"]["host_internal_node_accuracy"] == 1.0
    assert payload["metrics"]["host_event_accuracy"] == 1.0
    assert payload["metrics"]["geographic_internal_node_accuracy"] == 1.0
    assert payload["metrics"]["geographic_event_accuracy"] == 1.0
    assert payload["metrics"]["reference_output_count"] == 23
    assert payload["data"]["dataset"]["dataset_id"] == "known_answer_reference_panel"
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
    assert payload["data"]["dataset_export"]["host_traits_path"] == str(
        output / "dataset" / "host-traits.tsv"
    )
    assert payload["data"]["workflow_bundle"]["host_event_recovery_path"] == str(
        output / "workflow" / "host-event-recovery.tsv"
    )
    assert payload["data"]["workflow_bundle"]["threshold_evaluation_path"] == str(
        output / "workflow" / "recovery-threshold-evaluation.tsv"
    )
