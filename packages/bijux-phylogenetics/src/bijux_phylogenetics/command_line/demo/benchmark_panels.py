from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.demo.shared import (
    count_expected_output_files,
    emit_demo_result,
    resolve_demo_runner,
)


def add_benchmark_demo_commands(demo_subparsers: Any) -> None:
    demo_real_dataset_benchmark = demo_subparsers.add_parser(
        "real-dataset-macroevolution",
        help="Materialize the packaged Central European plant dataset together with the governed real-dataset macroevolution benchmark bundle.",
    )
    demo_real_dataset_benchmark.add_argument("--out", required=True, type=Path)
    demo_real_dataset_benchmark.add_argument(
        "--json", action="store_true", help="Emit the demo result as JSON."
    )
    _add_manifest_argument(demo_real_dataset_benchmark)


def run_benchmark_demo_command(args: Any) -> int | None:
    if args.demo_command != "real-dataset-macroevolution":
        return None

    result = resolve_demo_runner("run_real_dataset_macroevolution_benchmark_demo")(
        args.out
    )
    outputs = [
        result.dataset_export.readme_path,
        result.dataset_export.tree_path,
        result.dataset_export.traits_path,
        result.benchmark_bundle.review_traits_path,
        result.benchmark_bundle.summary_path,
        result.benchmark_bundle.model_table_path,
        result.benchmark_bundle.alignment_review_path,
        result.benchmark_bundle.parity_table_path,
        result.benchmark_bundle.geiger_reference_path,
        result.overview_path,
    ]
    return emit_demo_result(
        args,
        outputs=outputs,
        metrics={
            "artifact_count": len(outputs),
            "dataset_taxon_count": result.dataset.taxon_count,
            "summary_row_count": 4,
            "native_model_row_count": 8,
            "alignment_review_row_count": 2,
            "parity_row_count": 10,
            "reference_output_count": count_expected_output_files(
                result.dataset_export.expected_output_root
            ),
        },
        data=result,
        output_root=result.output_root,
    )
