from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.datasets.rabies_method_sensitivity as rabies_method_sensitivity
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    load_rabies_method_sensitivity_panel_dataset,
    run_rabies_method_sensitivity_panel_workflow,
)
import bijux_phylogenetics.datasets.rabies_method_sensitivity.workflow as rabies_method_sensitivity_workflow


def _build_stub_dataset(
    *, variant_count: int, parallel_workers: int
) -> rabies_method_sensitivity.RabiesMethodSensitivityPanelDataset:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    return replace(
        dataset,
        parallel_workers=parallel_workers,
        variants=dataset.variants[:variant_count],
    )


def test_run_rabies_method_sensitivity_panel_workflow_restricts_selected_variants(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset = _build_stub_dataset(variant_count=4, parallel_workers=2)
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "load_rabies_method_sensitivity_panel_dataset",
        lambda: dataset,
    )
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_run_variant_workflow",
        lambda *, variant, **_: SimpleNamespace(config=variant),
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

    report = run_rabies_method_sensitivity_panel_workflow(
        tmp_path / "workflow",
        parallel_workers=1,
        variant_ids=("ginsi-gap-threshold", "auto-gappyout"),
    )

    assert [variant.variant_id for variant in report.dataset.variants] == [
        "ginsi-gap-threshold",
        "auto-gappyout",
    ]
    assert [record.variant_id for record in report.task_records] == [
        "ginsi-gap-threshold",
        "auto-gappyout",
    ]
    assert report.execution_mode == "serial"


@pytest.mark.slow
def test_cli_demo_rabies_method_sensitivity_panel_json_output_respects_variant_selection(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "rabies-demo-selected"
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
            "--variant-id",
            "ginsi-gap-threshold",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["variant_count"] == 1
    assert payload["metrics"]["execution_mode"] == "serial"
    assert payload["metrics"]["slurm_job_count"] == 1
    assert payload["metrics"]["slurm_array_partition_count"] == 1
    assert payload["metrics"]["slurm_array_largest_partition_size"] == 1
    assert payload["data"]["workflow_bundle"]["slurm_array_strategy_path"] == str(
        output / "workflow" / "slurm-array-strategy.json"
    )
