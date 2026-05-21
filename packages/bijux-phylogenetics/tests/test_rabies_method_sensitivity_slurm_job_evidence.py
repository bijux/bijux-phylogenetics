from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    run_rabies_method_sensitivity_panel_workflow,
)
import bijux_phylogenetics.datasets.rabies_method_sensitivity.bundle as rabies_method_sensitivity_bundle
from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm import (
    build_rabies_method_sensitivity_slurm_planning_report,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.arrays import (
    build_rabies_method_sensitivity_slurm_array_strategy_report,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.job_evidence import (
    write_rabies_method_sensitivity_slurm_job_evidence_bundle,
)


@pytest.mark.slow
def test_write_rabies_method_sensitivity_panel_workflow_bundle_writes_job_evidence_packages(
    tmp_path: Path,
) -> None:
    report = run_rabies_method_sensitivity_panel_workflow(
        tmp_path / "run",
        variant_ids=("auto-gap-threshold",),
    )
    bundle_root = tmp_path / "workflow"
    bundle_root.mkdir(parents=True, exist_ok=True)
    rabies_method_sensitivity_bundle._copy_output(
        report.execution_record_path,
        bundle_root / report.execution_record_path.name,
    )
    rabies_method_sensitivity_bundle._copy_task_logs(
        bundle_root / "parallel-logs",
        report.task_records,
    )
    rabies_method_sensitivity_bundle._write_variant_outputs(
        bundle_root / "variants",
        report.variant_runs,
    )
    workflow_manifest_path = bundle_root / "rabies-method-sensitivity.manifest.json"
    workflow_manifest_path.write_text("{}\n", encoding="utf-8")
    planning_report = build_rabies_method_sensitivity_slurm_planning_report(report)
    array_strategy_report = build_rabies_method_sensitivity_slurm_array_strategy_report(
        planning_report
    )
    evidence_report = write_rabies_method_sensitivity_slurm_job_evidence_bundle(
        bundle_root / "slurm-job-evidence",
        bundle_root=bundle_root,
        dataset_id=report.dataset.dataset_id,
        workflow_prefix=report.dataset.workflow_prefix,
        execution_mode=report.execution_mode,
        parallel_workers=report.parallel_workers,
        task_records=report.task_records,
        variant_runs=report.variant_runs,
        array_strategy_report=array_strategy_report,
        execution_record_path=bundle_root / report.execution_record_path.name,
        workflow_manifest_path=workflow_manifest_path,
    )

    evidence_root = bundle_root / "slurm-job-evidence" / "auto-gap-threshold"
    evidence_json_path = evidence_root / "job-evidence.json"
    evidence_html_path = evidence_root / "job-evidence.html"
    index_path = bundle_root / "slurm-job-evidence.tsv"
    summary_path = bundle_root / "slurm-job-evidence-summary.json"
    assert evidence_json_path.is_file()
    assert evidence_html_path.is_file()
    assert (evidence_root / "task.log").is_file()
    assert (evidence_root / "alignment.manifest.json").is_file()
    assert (evidence_root / "trimming.manifest.json").is_file()
    assert (evidence_root / "inference-comparison.manifest.json").is_file()
    assert (evidence_root / "model-selection.manifest.json").is_file()
    assert (evidence_root / "iqtree-support.manifest.json").is_file()
    assert (evidence_root / "fasttree.manifest.json").is_file()
    assert index_path.is_file()
    assert summary_path.is_file()
    assert evidence_report.job_count == 1
    assert evidence_report.completed_job_count == 1

    payload = json.loads(evidence_json_path.read_text(encoding="utf-8"))
    assert payload["dataset_id"] == "rabies_method_sensitivity_panel"
    assert payload["workflow_prefix"] == "rabies-method-sensitivity-panel"
    assert payload["variant_id"] == "auto-gap-threshold"
    assert payload["status"] == "succeeded"
    assert payload["variant"]["selected_model"]
    assert payload["commands"]["alignment"]
    assert payload["engine_versions"]["mafft"]
    assert payload["output_inventory"]

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["job_count"] == 1
    assert summary["completed_job_count"] == 1
    assert summary["jobs"][0]["evidence_json_path"] == (
        "slurm-job-evidence/auto-gap-threshold/job-evidence.json"
    )

    table_text = index_path.read_text(encoding="utf-8")
    assert table_text.startswith(
        "partition_id\tarray_index\tvariant_id\tstatus\tselected_model\t"
    )
