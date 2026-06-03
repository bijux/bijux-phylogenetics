from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.rabies_method_sensitivity as rabies_method_sensitivity
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    export_rabies_method_sensitivity_panel_dataset,
    load_rabies_method_sensitivity_panel_dataset,
    run_rabies_method_sensitivity_panel_demo,
    run_rabies_method_sensitivity_panel_workflow,
    write_rabies_method_sensitivity_panel_workflow_bundle,
)
import bijux_phylogenetics.datasets.rabies_method_sensitivity.workflow as rabies_method_sensitivity_workflow

from .support.scientific_output_assertions import (
    assert_selected_scientific_outputs_equivalent,
)


def _build_stub_dataset(
    *, variant_count: int, parallel_workers: int
) -> rabies_method_sensitivity.RabiesMethodSensitivityPanelDataset:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    return replace(
        dataset,
        parallel_workers=parallel_workers,
        variants=dataset.variants[:variant_count],
    )


def test_load_rabies_method_sensitivity_panel_dataset_exposes_packaged_surface() -> (
    None
):
    dataset = load_rabies_method_sensitivity_panel_dataset()
    assert dataset.dataset_id == "rabies_method_sensitivity_panel"
    assert dataset.label == "Rabies method-sensitivity panel"
    assert dataset.taxon_count == 9
    assert dataset.sequence_type == "dna"
    assert dataset.outgroup_taxa == ("bat_chile_rv108",)
    assert dataset.parallel_workers == 2
    assert len(dataset.variants) == 4
    assert dataset.variants[0].variant_id == "auto-gap-threshold"
    assert dataset.sequences_path.is_file()
    assert dataset.metadata_path.is_file()
    assert dataset.reference_output_root.is_dir()
    assert "MG458305" in dataset.source_accessions


def test_export_rabies_method_sensitivity_panel_dataset_copies_expected_outputs(
    tmp_path: Path,
) -> None:
    result = export_rabies_method_sensitivity_panel_dataset(tmp_path / "dataset")
    expected_files = {
        path.relative_to(result.expected_output_root)
        for path in result.expected_output_root.rglob("*")
        if path.is_file()
    }
    assert result.readme_path.is_file()
    assert result.config_path.is_file()
    assert result.sequences_path.is_file()
    assert result.metadata_path.is_file()
    assert len(expected_files) == 147
    assert Path("parallel-execution-summary.tsv") in expected_files
    assert Path("rabies-method-sensitivity-panel.run.json") in expected_files
    assert Path("rabies-method-sensitivity.manifest.json") in expected_files
    assert Path("slurm-job-plan.tsv") in expected_files
    assert Path("slurm-estimation-assumptions.tsv") in expected_files
    assert Path("slurm-planning-summary.json") in expected_files
    assert Path("slurm-array-partitions.tsv") in expected_files
    assert Path("slurm-array-members.tsv") in expected_files
    assert Path("slurm-array-strategy.json") in expected_files
    assert Path("slurm-job-evidence.tsv") in expected_files
    assert Path("slurm-job-evidence-summary.json") in expected_files
    assert (
        Path("slurm-job-evidence/auto-gap-threshold/job-evidence.json")
        in expected_files
    )
    assert Path("slurm-storage-categories.tsv") in expected_files
    assert Path("slurm-storage-variants.tsv") in expected_files
    assert Path("slurm-storage-report.json") in expected_files
    assert Path("slurm-storage-report.html") in expected_files
    assert Path("slurm-output-explosion-checks.tsv") in expected_files
    assert Path("slurm-output-explosion-variants.tsv") in expected_files
    assert Path("slurm-output-explosion-report.json") in expected_files
    assert Path("slurm-output-explosion-report.html") in expected_files
    assert Path("slurm-tree-retention-checks.tsv") in expected_files
    assert Path("slurm-tree-retention-files.tsv") in expected_files
    assert Path("slurm-tree-retention-policy.json") in expected_files
    assert Path("slurm-tree-retention-policy.html") in expected_files
    assert Path("slurm-merge-checks.tsv") in expected_files
    assert Path("slurm-merge-variants.tsv") in expected_files
    assert Path("slurm-merge-report.json") in expected_files
    assert Path("slurm-merge-report.html") in expected_files
    assert Path("slurm-output-freshness.tsv") in expected_files
    assert Path("slurm-output-freshness-checks.tsv") in expected_files
    assert Path("slurm-output-freshness.json") in expected_files
    assert Path("slurm-job-status.tsv") in expected_files
    assert Path("slurm-partition-status.tsv") in expected_files
    assert Path("slurm-workflow-status.json") in expected_files
    assert Path("slurm-failure-recovery-jobs.tsv") in expected_files
    assert Path("slurm-failure-recovery-partitions.tsv") in expected_files
    assert Path("slurm-failure-recovery-report.json") in expected_files
    assert Path("slurm-failure-recovery-report.html") in expected_files
    assert Path("slurm-arrays/compact-mafft-auto-standard.sbatch") in expected_files
    assert Path("reproducibility-checks.tsv") in expected_files
    assert Path("reproducibility-variants.tsv") in expected_files
    assert Path("reproducibility-audit.json") in expected_files
    assert (
        Path("report-artifacts/rabies-method-sensitivity-report.manifest.json")
        in expected_files
    )
    assert Path("parallel-logs/auto-gap-threshold.log") in expected_files
    assert Path("workflow-summary.tsv") in expected_files
    assert (
        Path("variants/auto-gap-threshold/unrooted-conclusions.tsv") in expected_files
    )
    report_html = (
        result.expected_output_root / "rabies-method-sensitivity-report.html"
    ).read_text(encoding="utf-8")
    assert 'href="workflow-summary.tsv"' in report_html
    assert 'href="slurm-job-plan.tsv"' in report_html
    assert 'href="slurm-array-partitions.tsv"' in report_html
    assert 'href="slurm-job-evidence-summary.json"' in report_html
    assert 'href="slurm-storage-report.json"' in report_html
    assert 'href="slurm-output-explosion-report.json"' in report_html
    assert 'href="slurm-tree-retention-policy.json"' in report_html
    assert 'href="slurm-merge-report.json"' in report_html
    assert 'href="slurm-output-freshness.json"' in report_html
    assert 'href="slurm-workflow-status.json"' in report_html
    assert 'href="slurm-failure-recovery-report.json"' in report_html
    assert 'href="reproducibility-audit.json"' in report_html
    assert (
        "report-artifacts/rabies-method-sensitivity-report.manifest.json" in report_html
    )


@pytest.mark.slow
def test_write_rabies_method_sensitivity_panel_workflow_bundle_matches_packaged_expected_outputs(
    tmp_path: Path,
) -> None:
    report = run_rabies_method_sensitivity_panel_workflow(tmp_path / "run")
    bundle = write_rabies_method_sensitivity_panel_workflow_bundle(
        tmp_path / "workflow",
        report,
    )
    expected_root = report.dataset.reference_output_root
    generated = {
        path.relative_to(bundle.output_root): path
        for path in bundle.output_root.rglob("*")
        if path.is_file()
    }
    expected = {
        path.relative_to(expected_root): path
        for path in expected_root.rglob("*")
        if path.is_file()
    }
    assert set(generated) == set(expected)
    assert_selected_scientific_outputs_equivalent(expected_root, expected)


@pytest.mark.slow
def test_run_rabies_method_sensitivity_panel_demo_materializes_dataset_and_workflow(
    tmp_path: Path,
) -> None:
    result = run_rabies_method_sensitivity_panel_demo(tmp_path / "demo")
    assert result.dataset.taxon_count == 9
    assert result.dataset_export.config_path.is_file()
    assert result.workflow_bundle.workflow_summary_path.is_file()
    assert result.workflow_bundle.variant_summary_path.is_file()
    assert result.workflow_bundle.parallel_summary_path.is_file()
    assert result.workflow_bundle.execution_record_path.is_file()
    assert result.workflow_bundle.manifest_path.is_file()
    assert result.workflow_bundle.report_manifest_path.is_file()
    assert result.workflow_bundle.slurm_job_plan_path.is_file()
    assert result.workflow_bundle.slurm_assumptions_path.is_file()
    assert result.workflow_bundle.slurm_summary_path.is_file()
    assert result.workflow_bundle.slurm_array_partitions_path.is_file()
    assert result.workflow_bundle.slurm_array_members_path.is_file()
    assert result.workflow_bundle.slurm_array_strategy_path.is_file()
    assert result.workflow_bundle.slurm_job_evidence_root.is_dir()
    assert result.workflow_bundle.slurm_job_evidence_index_path.is_file()
    assert result.workflow_bundle.slurm_job_evidence_summary_path.is_file()
    assert result.workflow_bundle.slurm_storage_categories_path.is_file()
    assert result.workflow_bundle.slurm_storage_variants_path.is_file()
    assert result.workflow_bundle.slurm_storage_summary_path.is_file()
    assert result.workflow_bundle.slurm_storage_report_path.is_file()
    assert result.workflow_bundle.slurm_output_explosion_checks_path.is_file()
    assert result.workflow_bundle.slurm_output_explosion_variants_path.is_file()
    assert result.workflow_bundle.slurm_output_explosion_summary_path.is_file()
    assert result.workflow_bundle.slurm_output_explosion_report_path.is_file()
    assert result.workflow_bundle.slurm_tree_retention_checks_path.is_file()
    assert result.workflow_bundle.slurm_tree_retention_files_path.is_file()
    assert result.workflow_bundle.slurm_tree_retention_summary_path.is_file()
    assert result.workflow_bundle.slurm_tree_retention_report_path.is_file()
    assert result.workflow_bundle.slurm_merge_checks_path.is_file()
    assert result.workflow_bundle.slurm_merge_variants_path.is_file()
    assert result.workflow_bundle.slurm_merge_summary_path.is_file()
    assert result.workflow_bundle.slurm_merge_report_path.is_file()
    assert result.workflow_bundle.slurm_output_freshness_path.is_file()
    assert result.workflow_bundle.slurm_output_freshness_checks_path.is_file()
    assert result.workflow_bundle.slurm_output_freshness_summary_path.is_file()
    assert result.workflow_bundle.slurm_job_status_path.is_file()
    assert result.workflow_bundle.slurm_partition_status_path.is_file()
    assert result.workflow_bundle.slurm_workflow_status_path.is_file()
    assert result.workflow_bundle.slurm_failure_recovery_jobs_path.is_file()
    assert result.workflow_bundle.slurm_failure_recovery_partitions_path.is_file()
    assert result.workflow_bundle.slurm_failure_recovery_summary_path.is_file()
    assert result.workflow_bundle.slurm_failure_recovery_report_path.is_file()
    assert result.workflow_bundle.slurm_array_scripts_root.is_dir()
    assert result.workflow_bundle.slurm_job_count == 4
    assert result.workflow_bundle.slurm_total_estimated_core_hours > 0
    assert result.workflow_bundle.slurm_maximum_estimated_memory_mib >= 1024
    assert result.workflow_bundle.slurm_maximum_estimated_wallclock_minutes >= 20
    assert result.workflow_bundle.slurm_total_estimated_scratch_mib > 0
    assert result.workflow_bundle.slurm_total_estimated_output_mib > 0
    assert result.workflow_bundle.slurm_array_partition_count == 2
    assert result.workflow_bundle.slurm_array_script_count == 2
    assert result.workflow_bundle.slurm_array_largest_partition_size == 2
    assert result.workflow_bundle.slurm_job_evidence_file_count > 0
    assert result.workflow_bundle.slurm_job_evidence_total_runtime_seconds > 0
    assert result.workflow_bundle.slurm_job_evidence_total_output_byte_count > 0
    assert result.workflow_bundle.slurm_storage_total_estimated_mib > 0
    assert result.workflow_bundle.slurm_storage_output_byte_count > 0
    assert result.workflow_bundle.slurm_storage_log_byte_count > 0
    assert result.workflow_bundle.slurm_storage_tree_byte_count > 0
    assert result.workflow_bundle.slurm_storage_posterior_sample_byte_count == 0
    assert result.workflow_bundle.slurm_storage_report_byte_count > 0
    assert result.workflow_bundle.slurm_storage_largest_variant_id in {
        "auto-gap-threshold",
        "auto-gappyout",
        "ginsi-gap-threshold",
        "ginsi-gappyout",
    }
    assert result.workflow_bundle.slurm_output_explosion_status == "low"
    assert result.workflow_bundle.slurm_output_explosion_global_issue_count == 0
    assert result.workflow_bundle.slurm_output_explosion_warning_variant_count == 0
    assert result.workflow_bundle.slurm_output_explosion_high_risk_variant_count == 0
    assert result.workflow_bundle.slurm_tree_retention_status == "no_action"
    assert result.workflow_bundle.slurm_tree_set_file_count == 0
    assert result.workflow_bundle.slurm_tree_posterior_sample_file_count == 0
    assert result.workflow_bundle.slurm_tree_thinning_recommended_file_count == 0
    assert result.workflow_bundle.slurm_tree_thinning_required_file_count == 0
    assert result.workflow_bundle.slurm_tree_compression_recommended_file_count == 0
    assert result.workflow_bundle.slurm_tree_compression_required_file_count == 0
    assert result.workflow_bundle.slurm_merge_status == "merge-ready"
    assert result.workflow_bundle.slurm_merge_ready is True
    assert result.workflow_bundle.slurm_mergeable_variant_count == 4
    assert result.workflow_bundle.slurm_merge_failed_check_count == 0
    assert result.workflow_bundle.slurm_output_freshness_check_count > 0
    assert result.workflow_bundle.slurm_output_freshness_failed_check_count == 0
    assert result.workflow_bundle.slurm_fresh_output_job_count == 4
    assert result.workflow_bundle.slurm_stale_output_job_count == 0
    assert result.workflow_bundle.slurm_completed_job_count == 4
    assert result.workflow_bundle.slurm_failed_job_count == 0
    assert result.workflow_bundle.slurm_pending_job_count == 0
    assert result.workflow_bundle.slurm_stale_job_count == 0
    assert result.workflow_bundle.slurm_failure_recovery_status == "clean"
    assert result.workflow_bundle.slurm_failure_recovery_rerunnable_job_count == 0
    assert result.workflow_bundle.slurm_failure_recovery_blocked_job_count == 0
    assert result.workflow_bundle.slurm_failure_recovery_partition_count == 0
    assert result.workflow_bundle.reproducibility_checks_path.is_file()
    assert result.workflow_bundle.reproducibility_variant_audit_path.is_file()
    assert result.workflow_bundle.reproducibility_audit_path.is_file()
    assert result.workflow_bundle.reproducibility_passed is True
    assert result.workflow_bundle.task_logs_root.is_dir()
    assert result.workflow_bundle.report_path.is_file()
    assert result.overview_path.is_file()
    overview = result.overview_path.read_text(encoding="utf-8")
    assert "variants" in overview
    assert "reproducibility audit" in overview


def test_public_runtime_exports_include_rabies_method_sensitivity_surface() -> None:
    assert (
        rabies_method_sensitivity.load_rabies_method_sensitivity_panel_dataset
        is load_rabies_method_sensitivity_panel_dataset
    )
    assert (
        rabies_method_sensitivity.export_rabies_method_sensitivity_panel_dataset
        is export_rabies_method_sensitivity_panel_dataset
    )
    assert (
        rabies_method_sensitivity.run_rabies_method_sensitivity_panel_workflow
        is run_rabies_method_sensitivity_panel_workflow
    )
    assert (
        rabies_method_sensitivity.write_rabies_method_sensitivity_panel_workflow_bundle
        is write_rabies_method_sensitivity_panel_workflow_bundle
    )
    assert (
        rabies_method_sensitivity.run_rabies_method_sensitivity_panel_demo
        is run_rabies_method_sensitivity_panel_demo
    )


def test_run_rabies_method_sensitivity_panel_workflow_writes_execution_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset = _build_stub_dataset(variant_count=1, parallel_workers=1)
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "load_rabies_method_sensitivity_panel_dataset",
        lambda: dataset,
    )
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_run_variant_workflow",
        lambda **_: object(),
    )
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_build_preprocessing_comparison_rows",
        lambda _: [],
    )
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_aggregate_clades",
        lambda *_args, **_kwargs: [],
    )
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_build_conclusion_rows",
        lambda **_: [],
    )

    output_root = tmp_path / "workflow"
    report = run_rabies_method_sensitivity_panel_workflow(output_root)

    execution_record_path = output_root / "rabies-method-sensitivity-panel.run.json"
    payload = json.loads(execution_record_path.read_text(encoding="utf-8"))
    assert report.execution_mode == "serial"
    assert payload["dataset_id"] == dataset.dataset_id
    assert payload["status"] == "succeeded"
    assert payload["parallel_workers"] == 1
    assert payload["execution_mode"] == "serial"
    assert payload["selected_variant_ids"] == [dataset.variants[0].variant_id]
    assert payload["successful_variants"] == [dataset.variants[0].variant_id]
    assert payload["failed_variants"] == []
    assert (
        payload["task_records"][0]["log_path"] == "parallel-logs/auto-gap-threshold.log"
    )


@pytest.mark.slow
def test_cli_demo_rabies_method_sensitivity_panel_json_output_reports_method_review(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "rabies-demo"
    exit_code = main(
        [
            "demo",
            "rabies-method-sensitivity-panel",
            "--out",
            str(output),
            "--mafft-executable",
            "mafft",
            "--trimal-executable",
            "trimal",
            "--iqtree-executable",
            "iqtree2",
            "--fasttree-executable",
            "FastTree",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["artifact_count"] == 54
    assert payload["metrics"]["taxon_count"] == 9
    assert payload["metrics"]["variant_count"] == 4
    assert payload["metrics"]["parallel_workers"] == 2
    assert payload["metrics"]["execution_mode"] == "parallel"
    assert payload["metrics"]["stable_clade_count"] == 2
    assert payload["metrics"]["changed_clade_count"] == 8
    assert payload["metrics"]["preprocessing_change_pair_count"] == 0
    assert payload["metrics"]["rooted_engine_change_variant_count"] == 0
    assert payload["metrics"]["serious_conflict_variant_count"] == 4
    assert payload["metrics"]["report_linked_artifact_count"] == 48
    assert payload["metrics"]["report_html_size_bytes"] > 0
    assert payload["metrics"]["report_linked_artifact_bytes"] > 0
    assert (
        payload["metrics"]["report_total_output_bytes"]
        >= payload["metrics"]["report_html_size_bytes"]
    )
    assert payload["metrics"]["slurm_job_count"] == 4
    assert payload["metrics"]["slurm_total_estimated_core_hours"] > 0
    assert payload["metrics"]["slurm_maximum_estimated_memory_mib"] >= 1024
    assert payload["metrics"]["slurm_maximum_estimated_wallclock_minutes"] >= 20
    assert payload["metrics"]["slurm_total_estimated_scratch_mib"] > 0
    assert payload["metrics"]["slurm_total_estimated_output_mib"] > 0
    assert payload["metrics"]["slurm_array_partition_count"] == 2
    assert payload["metrics"]["slurm_array_script_count"] == 2
    assert payload["metrics"]["slurm_array_largest_partition_size"] == 2
    assert payload["metrics"]["slurm_job_evidence_file_count"] > 0
    assert payload["metrics"]["slurm_job_evidence_total_runtime_seconds"] > 0
    assert payload["metrics"]["slurm_job_evidence_total_output_byte_count"] > 0
    assert payload["metrics"]["slurm_storage_total_estimated_mib"] > 0
    assert payload["metrics"]["slurm_storage_output_byte_count"] > 0
    assert payload["metrics"]["slurm_storage_log_byte_count"] > 0
    assert payload["metrics"]["slurm_storage_tree_byte_count"] > 0
    assert payload["metrics"]["slurm_storage_posterior_sample_byte_count"] == 0
    assert payload["metrics"]["slurm_storage_report_byte_count"] > 0
    assert payload["metrics"]["slurm_storage_largest_variant_id"] in {
        "auto-gap-threshold",
        "auto-gappyout",
        "ginsi-gap-threshold",
        "ginsi-gappyout",
    }
    assert payload["metrics"]["slurm_output_explosion_status"] == "low"
    assert payload["metrics"]["slurm_output_explosion_global_issue_count"] == 0
    assert payload["metrics"]["slurm_output_explosion_warning_variant_count"] == 0
    assert payload["metrics"]["slurm_output_explosion_high_risk_variant_count"] == 0
    assert payload["metrics"]["slurm_tree_retention_status"] == "no_action"
    assert payload["metrics"]["slurm_tree_set_file_count"] == 0
    assert payload["metrics"]["slurm_tree_posterior_sample_file_count"] == 0
    assert payload["metrics"]["slurm_tree_thinning_recommended_file_count"] == 0
    assert payload["metrics"]["slurm_tree_thinning_required_file_count"] == 0
    assert payload["metrics"]["slurm_tree_compression_recommended_file_count"] == 0
    assert payload["metrics"]["slurm_tree_compression_required_file_count"] == 0
    assert payload["metrics"]["slurm_merge_status"] == "merge-ready"
    assert payload["metrics"]["slurm_merge_ready"] is True
    assert payload["metrics"]["slurm_mergeable_variant_count"] == 4
    assert payload["metrics"]["slurm_merge_failed_check_count"] == 0
    assert payload["metrics"]["slurm_output_freshness_check_count"] > 0
    assert payload["metrics"]["slurm_output_freshness_failed_check_count"] == 0
    assert payload["metrics"]["slurm_fresh_output_job_count"] == 4
    assert payload["metrics"]["slurm_stale_output_job_count"] == 0
    assert payload["metrics"]["slurm_completed_job_count"] == 4
    assert payload["metrics"]["slurm_failed_job_count"] == 0
    assert payload["metrics"]["slurm_pending_job_count"] == 0
    assert payload["metrics"]["slurm_stale_job_count"] == 0
    assert payload["metrics"]["slurm_failure_recovery_status"] == "clean"
    assert payload["metrics"]["slurm_failure_recovery_rerunnable_job_count"] == 0
    assert payload["metrics"]["slurm_failure_recovery_blocked_job_count"] == 0
    assert payload["metrics"]["slurm_failure_recovery_partition_count"] == 0
    assert payload["metrics"]["reproducibility_passed"] is True
    assert payload["metrics"]["reproducibility_check_count"] > 0
    assert payload["metrics"]["reproducibility_failed_check_count"] == 0
    assert payload["metrics"]["reproducibility_failed_variant_count"] == 0
    assert payload["metrics"]["reference_output_count"] == 147
    assert payload["data"]["dataset"]["dataset_id"] == "rabies_method_sensitivity_panel"
    assert payload["data"]["workflow_bundle"]["workflow_summary_path"] == str(
        output / "workflow" / "workflow-summary.tsv"
    )
    assert payload["data"]["workflow_bundle"]["parallel_summary_path"] == str(
        output / "workflow" / "parallel-execution-summary.tsv"
    )
    assert payload["data"]["workflow_bundle"]["execution_record_path"] == str(
        output / "workflow" / "rabies-method-sensitivity-panel.run.json"
    )
    assert payload["data"]["workflow_bundle"]["report_manifest_path"] == str(
        output
        / "workflow"
        / "report-artifacts"
        / "rabies-method-sensitivity-report.manifest.json"
    )
    assert payload["data"]["workflow_bundle"]["slurm_job_plan_path"] == str(
        output / "workflow" / "slurm-job-plan.tsv"
    )
    assert payload["data"]["workflow_bundle"]["slurm_array_partitions_path"] == str(
        output / "workflow" / "slurm-array-partitions.tsv"
    )
    assert payload["data"]["workflow_bundle"]["slurm_job_evidence_index_path"] == str(
        output / "workflow" / "slurm-job-evidence.tsv"
    )
    assert payload["data"]["workflow_bundle"]["slurm_job_evidence_summary_path"] == str(
        output / "workflow" / "slurm-job-evidence-summary.json"
    )
    assert payload["data"]["workflow_bundle"]["slurm_storage_summary_path"] == str(
        output / "workflow" / "slurm-storage-report.json"
    )
    assert payload["data"]["workflow_bundle"]["slurm_storage_report_path"] == str(
        output / "workflow" / "slurm-storage-report.html"
    )
    assert payload["data"]["workflow_bundle"][
        "slurm_output_explosion_summary_path"
    ] == str(output / "workflow" / "slurm-output-explosion-report.json")
    assert payload["data"]["workflow_bundle"][
        "slurm_output_explosion_report_path"
    ] == str(output / "workflow" / "slurm-output-explosion-report.html")
    assert payload["data"]["workflow_bundle"][
        "slurm_tree_retention_summary_path"
    ] == str(output / "workflow" / "slurm-tree-retention-policy.json")
    assert payload["data"]["workflow_bundle"][
        "slurm_tree_retention_report_path"
    ] == str(output / "workflow" / "slurm-tree-retention-policy.html")
    assert payload["data"]["workflow_bundle"]["slurm_merge_summary_path"] == str(
        output / "workflow" / "slurm-merge-report.json"
    )
    assert payload["data"]["workflow_bundle"]["slurm_merge_report_path"] == str(
        output / "workflow" / "slurm-merge-report.html"
    )
    assert payload["data"]["workflow_bundle"]["slurm_output_freshness_path"] == str(
        output / "workflow" / "slurm-output-freshness.tsv"
    )
    assert payload["data"]["workflow_bundle"][
        "slurm_output_freshness_summary_path"
    ] == str(output / "workflow" / "slurm-output-freshness.json")
    assert payload["data"]["workflow_bundle"]["slurm_job_status_path"] == str(
        output / "workflow" / "slurm-job-status.tsv"
    )
    assert payload["data"]["workflow_bundle"]["slurm_workflow_status_path"] == str(
        output / "workflow" / "slurm-workflow-status.json"
    )
    assert payload["data"]["workflow_bundle"][
        "slurm_failure_recovery_summary_path"
    ] == str(output / "workflow" / "slurm-failure-recovery-report.json")
    assert payload["data"]["workflow_bundle"][
        "slurm_failure_recovery_report_path"
    ] == str(output / "workflow" / "slurm-failure-recovery-report.html")
    assert payload["data"]["workflow_bundle"]["reproducibility_checks_path"] == str(
        output / "workflow" / "reproducibility-checks.tsv"
    )
    assert payload["data"]["workflow_bundle"]["reproducibility_audit_path"] == str(
        output / "workflow" / "reproducibility-audit.json"
    )
