from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import threading
import time

import pytest

import bijux_phylogenetics.datasets.rabies_method_sensitivity as rabies_method_sensitivity
from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    load_rabies_method_sensitivity_panel_dataset,
    run_rabies_method_sensitivity_panel_workflow,
)
import bijux_phylogenetics.datasets.rabies_method_sensitivity.workflow as rabies_method_sensitivity_workflow
from bijux_phylogenetics.runtime.errors import EngineWorkflowError


def _build_stub_dataset(
    *, variant_count: int, parallel_workers: int
) -> rabies_method_sensitivity.RabiesMethodSensitivityPanelDataset:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    return replace(
        dataset,
        parallel_workers=parallel_workers,
        variants=dataset.variants[:variant_count],
    )


def test_run_rabies_method_sensitivity_panel_workflow_rejects_concurrent_reuse_of_same_output_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset = _build_stub_dataset(variant_count=1, parallel_workers=1)
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "load_rabies_method_sensitivity_panel_dataset",
        lambda: dataset,
    )

    def slow_failure(**_: object) -> object:
        time.sleep(0.3)
        raise EngineWorkflowError(
            "fixture variant failed", code="fixture_variant_failed"
        )

    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_run_variant_workflow",
        slow_failure,
    )

    shared_root = tmp_path / "shared-workflow"
    first_errors: list[EngineWorkflowError] = []

    def run_first() -> None:
        try:
            run_rabies_method_sensitivity_panel_workflow(shared_root)
        except EngineWorkflowError as error:
            first_errors.append(error)

    thread = threading.Thread(target=run_first)
    thread.start()
    marker_path = shared_root / "rabies-method-sensitivity-panel.run.running.json"
    deadline = time.time() + 5.0
    while not marker_path.exists():
        if time.time() >= deadline:
            raise AssertionError(
                "expected rabies method-sensitivity running marker to appear"
            )
        time.sleep(0.01)

    with pytest.raises(EngineWorkflowError) as second_error:
        run_rabies_method_sensitivity_panel_workflow(shared_root)

    thread.join()

    assert second_error.value.code == (
        "rabies_method_sensitivity_workflow_already_running"
    )
    assert first_errors[0].code == "workflow_parallel_task_failed"
    assert marker_path.exists() is False


def test_run_rabies_method_sensitivity_panel_workflow_allows_parallel_distinct_output_roots(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset = _build_stub_dataset(variant_count=1, parallel_workers=1)
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "load_rabies_method_sensitivity_panel_dataset",
        lambda: dataset,
    )

    def slow_success(**_: object) -> object:
        time.sleep(0.2)
        return object()

    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_run_variant_workflow",
        slow_success,
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

    reports: list[object] = []
    errors: list[BaseException] = []
    output_roots = [tmp_path / "left-workflow", tmp_path / "right-workflow"]

    def run_one(output_root: Path) -> None:
        try:
            reports.append(run_rabies_method_sensitivity_panel_workflow(output_root))
        except BaseException as error:  # pragma: no cover - failure is asserted below
            errors.append(error)

    threads = [
        threading.Thread(target=run_one, args=(output_root,))
        for output_root in output_roots
    ]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert errors == []
    assert len(reports) == 2
    for output_root in output_roots:
        assert (output_root / "rabies-method-sensitivity-panel.run.json").is_file()
        assert (
            output_root / "rabies-method-sensitivity-panel.run.running.json"
        ).exists() is False


def test_run_rabies_method_sensitivity_panel_workflow_preserves_successful_outputs_when_one_parallel_variant_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset = _build_stub_dataset(variant_count=2, parallel_workers=2)
    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "load_rabies_method_sensitivity_panel_dataset",
        lambda: dataset,
    )

    def mixed_outcomes(**kwargs: object) -> object:
        variant = kwargs["variant"]
        variant_root = kwargs["variant_root"]
        variant_root.mkdir(parents=True, exist_ok=True)
        if variant.variant_id == dataset.variants[0].variant_id:
            kept_output = variant_root / "kept-output.txt"
            kept_output.write_text("kept\n", encoding="utf-8")
            return object()
        raise EngineWorkflowError(
            "fixture variant failed", code="fixture_variant_failed"
        )

    monkeypatch.setattr(
        rabies_method_sensitivity_workflow,
        "_run_variant_workflow",
        mixed_outcomes,
    )

    output_root = tmp_path / "workflow"
    with pytest.raises(EngineWorkflowError) as error:
        run_rabies_method_sensitivity_panel_workflow(output_root)

    assert error.value.code == "workflow_parallel_task_failed"
    assert error.value.details["successful_variants"] == [
        dataset.variants[0].variant_id
    ]
    assert error.value.details["failed_variants"] == [dataset.variants[1].variant_id]
    assert (
        output_root / "variants" / dataset.variants[0].variant_id / "kept-output.txt"
    ).read_text(encoding="utf-8") == "kept\n"
    assert (
        output_root / "parallel-logs" / f"{dataset.variants[0].variant_id}.log"
    ).is_file()
    assert (
        output_root / "parallel-logs" / f"{dataset.variants[1].variant_id}.log"
    ).is_file()

    payload = json.loads(
        (output_root / "rabies-method-sensitivity-panel.run.json").read_text(
            encoding="utf-8"
        )
    )
    assert payload["status"] == "failed"
    assert payload["successful_variants"] == [dataset.variants[0].variant_id]
    assert payload["failed_variants"] == [dataset.variants[1].variant_id]
