from __future__ import annotations

import pytest

from bijux_phylogenetics.core import copy_example_inputs, run_capability_demo
from bijux_phylogenetics.reports.service import render_tree_report

pytestmark = pytest.mark.smoke


def test_render_tree_report_writes_html_artifact(tmp_path) -> None:
    inputs = copy_example_inputs(tmp_path / "example-inputs")
    output_path = tmp_path / "tree-report.html"

    result = render_tree_report(tree_path=inputs["tree"], out_path=output_path)

    assert result.output_path == output_path
    assert output_path.exists()
    assert "<html" in output_path.read_text(encoding="utf-8").lower()


def test_run_capability_demo_materializes_expected_report_bundle(tmp_path) -> None:
    result = run_capability_demo(tmp_path / "capability-demo")

    assert result.input_root.exists()
    assert result.report_root.exists()
    assert result.tree_report.exists()
    assert result.dataset_report.exists()
    assert result.phylo_inputs_report.exists()
    assert result.comparison_report.exists()
    assert result.capability_summary.exists()
