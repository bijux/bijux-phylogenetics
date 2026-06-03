from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.presentation.posterior_uncertainty import (
    build_posterior_uncertainty_figure_package,
)
from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.reports.service import render_tree_uncertainty_report
from bijux_phylogenetics.runtime.errors import WorkflowBudgetError
from bijux_phylogenetics.trees import write_bootstrap_tree_set_artifacts

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_build_posterior_uncertainty_figure_package_records_processing_budget(
    tmp_path: Path,
) -> None:
    result = build_posterior_uncertainty_figure_package(
        fixture("example_tree_set_left.nwk"),
        out_dir=tmp_path / "uncertainty-package",
        memory_warning_threshold_bytes=1,
    )

    assert result.tree_count == 3
    assert result.processing.runtime_seconds >= 0.0
    assert result.processing.peak_memory_bytes >= 0
    assert "peak memory exceeded" in " ".join(result.budget_report.warning_messages)


def test_write_bootstrap_tree_set_artifacts_rejects_tree_count_budget(
    tmp_path: Path,
) -> None:
    with pytest.raises(WorkflowBudgetError) as error:
        write_bootstrap_tree_set_artifacts(
            fixture("example_tree_set_left.nwk"),
            out_dir=tmp_path / "bootstrap-review",
            max_tree_count=2,
        )

    assert error.value.code == "tree_set_tree_budget_exceeded"


def test_render_tree_uncertainty_report_truncates_budgeted_sections(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-set-report.html"
    report = render_tree_uncertainty_report(
        tree_set_path=fixture("example_tree_set_left.nwk"),
        out_path=output_path,
        max_report_table_rows=1,
        memory_warning_threshold_bytes=1,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.processing.runtime_seconds >= 0.0
    assert report.artifact_root.is_dir()
    assert report.artifact_manifest_path.is_file()
    assert report.methods_summary_path.is_file()
    assert report.methods_summary_warning_count >= 0
    assert "clade-frequencies" in report.budget_report.truncated_section_names
    assert 'href="tree-set-report.artifacts/clade-frequencies.tsv"' in html
    assert "methods-summary-text" in html
    assert "limitations" in html
    assert "preview_rows" in html
    assert "&quot;rows&quot;: [" not in html
    assert report.machine_manifest["budget"]["truncated_section_names"]
    assert report.machine_manifest["limitations"]
    assert (
        report.machine_manifest["linked_artifact_count"] == report.linked_artifact_count
    )
    assert report.total_output_bytes >= report.html_size_bytes


@pytest.mark.slow
def test_render_tree_uncertainty_report_scales_to_large_tree_sets(
    tmp_path: Path,
) -> None:
    large_tree_set = tmp_path / "large-tree-set.nwk"
    source_lines = (
        fixture("example_tree_set_left.nwk").read_text(encoding="utf-8").splitlines()
    )
    large_tree_set.write_text("\n".join(source_lines * 400) + "\n", encoding="utf-8")
    output_path = tmp_path / "large-tree-set-report.html"

    report = render_tree_uncertainty_report(
        tree_set_path=large_tree_set,
        out_path=output_path,
        max_report_table_rows=3,
    )

    html = output_path.read_text(encoding="utf-8")
    assert report.tree_count == 1200
    assert report.linked_artifact_count >= 10
    assert report.processing.peak_memory_bytes >= 0
    assert report.methods_summary_path.is_file()
    assert report.html_size_bytes > 0
    assert report.linked_artifact_bytes > 0
    assert report.total_output_bytes >= report.html_size_bytes
    assert (report.artifact_root / "clade-frequencies.tsv").is_file()
    assert (report.artifact_root / "tree-uncertainty.manifest.json").is_file()
    assert 'href="large-tree-set-report.artifacts/clade-frequencies.tsv"' in html
    assert "limitations" in html
    assert report.machine_manifest["html_size_bytes"] == report.html_size_bytes
    assert report.machine_manifest["total_output_bytes"] == report.total_output_bytes


def test_cli_tree_set_report_returns_structured_tree_budget_error(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "tree-set-report.html"
    exit_code = main(
        [
            "tree-set",
            "report",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--max-tree-count",
            "2",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 2
    assert payload["errors"][0]["code"] == "tree_set_tree_budget_exceeded"


def test_cli_tree_set_package_reports_budget_warning_metrics(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "uncertainty-package"
    exit_code = main(
        [
            "tree-set",
            "package",
            str(fixture("example_tree_set_left.nwk")),
            "--out-dir",
            str(output_dir),
            "--memory-warning-threshold-bytes",
            "1",
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["budget_warning_count"] >= 1
    assert payload["metrics"]["methods_summary_warning_count"] >= 0
    assert payload["metrics"]["runtime_seconds"] >= 0.0
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["artifact_count"] == 15
    assert (output_dir / "uncertainty-review.html").is_file()
    assert (output_dir / "tree-set-uncertainty-methods-summary.md").is_file()


def test_cli_tree_set_report_reports_output_size_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "tree-set-report.html"
    exit_code = main(
        [
            "tree-set",
            "report",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["linked_artifact_count"] >= 10
    assert payload["metrics"]["methods_summary_warning_count"] >= 0
    assert payload["data"]["methods_summary_path"].endswith(
        "tree-set-uncertainty-methods-summary.md"
    )
    assert payload["data"]["methods_summary_warning_count"] >= 0
    assert payload["data"]["limitations"]
    assert payload["metrics"]["html_size_bytes"] > 0
    assert payload["metrics"]["linked_artifact_bytes"] > 0
    assert (
        payload["metrics"]["total_output_bytes"]
        >= payload["metrics"]["html_size_bytes"]
    )
    assert payload["data"]["artifact_manifest_path"] == str(
        output_path.parent
        / "tree-set-report.artifacts"
        / "tree-uncertainty.manifest.json"
    )
