from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.benchmark import (
    NativeMaximumLikelihoodSpeedBenchmarkObservation,
    NativeMaximumLikelihoodSpeedBenchmarkReport,
    NativeMaximumLikelihoodSpeedTrendRow,
)
import bijux_phylogenetics.command_line as command_line_api
from bijux_phylogenetics.command_line import main


def test_cli_benchmark_native_maximum_likelihood_speed_reports_metrics_and_overrides(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    captured_call: dict[str, object] = {}
    report = NativeMaximumLikelihoodSpeedBenchmarkReport(
        replicates=1,
        observations=[
            NativeMaximumLikelihoodSpeedBenchmarkObservation(
                label="tiny-jc69-nni-4-taxa-12-sites",
                size_class="tiny",
                model_name="jc69",
                search_method="nni",
                taxon_count=4,
                site_count=12,
                runtime_seconds=1.25,
                peak_memory_bytes=4096,
                selected_model_name="JC69",
                best_run_source_label="stepwise-addition-tree",
                search_algorithm="nucleotide-likelihood-nni-search",
                final_log_likelihood=-27.213282538844922,
            ),
            NativeMaximumLikelihoodSpeedBenchmarkObservation(
                label="medium-jc69-nni-4-taxa-16-sites",
                size_class="medium",
                model_name="jc69",
                search_method="nni",
                taxon_count=4,
                site_count=16,
                runtime_seconds=2.5,
                peak_memory_bytes=8192,
                selected_model_name="JC69",
                best_run_source_label="stepwise-addition-tree",
                search_algorithm="nucleotide-likelihood-nni-search",
                final_log_likelihood=-31.5,
            ),
        ],
        trend_rows=[
            NativeMaximumLikelihoodSpeedTrendRow(
                model_name="jc69",
                search_method="nni",
                baseline_size_class="tiny",
                comparison_size_class="medium",
                baseline_taxon_count=4,
                comparison_taxon_count=4,
                baseline_site_count=12,
                comparison_site_count=16,
                baseline_runtime_seconds=1.25,
                comparison_runtime_seconds=2.5,
                runtime_delta_seconds=1.25,
                runtime_ratio=2.0,
            )
        ],
        size_class_labels=["tiny", "medium"],
        model_names=["jc69"],
        search_methods=["nni"],
        limitations=["governed test fixture"],
    )

    def fake_benchmark_native_maximum_likelihood_speed(**kwargs):
        captured_call.update(kwargs)
        return report

    monkeypatch.setattr(
        command_line_api,
        "benchmark_native_maximum_likelihood_speed",
        fake_benchmark_native_maximum_likelihood_speed,
    )

    exit_code = main(
        [
            "benchmark",
            "native-maximum-likelihood-speed",
            "--replicates",
            "1",
            "--model",
            "jc69",
            "--search-method",
            "nni",
            "--size-class",
            "tiny",
            "--taxon-count",
            "4",
            "--site-count",
            "12",
            "--size-class",
            "medium",
            "--taxon-count",
            "4",
            "--site-count",
            "16",
            "--start-tree-count",
            "1",
            "--start-tree-seed",
            "17",
            "--upper-branch-length-bound",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert captured_call == {
        "replicates": 1,
        "size_classes": [("tiny", 4, 12), ("medium", 4, 16)],
        "model_names": ["jc69"],
        "search_methods": ["nni"],
        "start_tree_count": 1,
        "start_tree_seed": 17,
        "upper_branch_length_bound": 1.0,
    }
    assert payload["command"] == "benchmark"
    assert payload["metrics"]["replicates"] == 1
    assert payload["metrics"]["observation_count"] == 2
    assert payload["metrics"]["trend_row_count"] == 1
    assert payload["metrics"]["size_class_count"] == 2
    assert payload["metrics"]["model_count"] == 1
    assert payload["metrics"]["search_method_count"] == 1
    assert payload["data"]["size_class_labels"] == ["tiny", "medium"]
    assert payload["data"]["model_names"] == ["jc69"]
    assert payload["data"]["search_methods"] == ["nni"]


def test_cli_benchmark_native_maximum_likelihood_speed_requires_aligned_size_fields(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit, match="2"):
        main(
            [
                "benchmark",
                "native-maximum-likelihood-speed",
                "--size-class",
                "tiny",
                "--taxon-count",
                "4",
                "--json",
            ]
        )
    captured = capsys.readouterr()
    assert (
        "requires the same number of --size-class, --taxon-count, and --site-count values"
        in captured.err
    )
