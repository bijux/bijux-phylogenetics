from __future__ import annotations

from typing import Any

from bijux_phylogenetics.benchmark import (
    run_real_dataset_macroevolution_benchmark_demo as run_real_dataset_macroevolution_benchmark_demo,
)
from bijux_phylogenetics.command_line.demo.benchmark_panels import (
    add_benchmark_demo_commands,
    run_benchmark_demo_command,
)
from bijux_phylogenetics.command_line.demo.introductory_panels import (
    add_introductory_demo_commands,
    run_introductory_demo_command,
)
from bijux_phylogenetics.command_line.demo.quality_panels import (
    add_quality_demo_commands,
    run_quality_demo_command,
)
from bijux_phylogenetics.command_line.demo.rabies_panels import (
    add_rabies_demo_commands,
    run_rabies_demo_command,
)
from bijux_phylogenetics.command_line.demo.recovery_panels import (
    add_recovery_demo_commands,
    run_recovery_demo_command,
)
from bijux_phylogenetics.command_line.demo.sequence_panels import (
    add_sequence_demo_commands,
    run_sequence_demo_command,
)
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.datasets import (
    run_rabies_cross_host_geography_panel_demo as run_rabies_cross_host_geography_panel_demo,
)


def add_demo_command(subparsers: Any) -> None:
    demo = subparsers.add_parser(
        get_command_spec("demo").name, help=get_command_spec("demo").summary
    )
    demo_subparsers = demo.add_subparsers(dest="demo_command", required=True)
    add_introductory_demo_commands(demo_subparsers)
    add_benchmark_demo_commands(demo_subparsers)
    add_sequence_demo_commands(demo_subparsers)
    add_rabies_demo_commands(demo_subparsers)
    add_quality_demo_commands(demo_subparsers)
    add_recovery_demo_commands(demo_subparsers)


def run_demo_command(args: Any) -> int:
    introductory_exit_code = run_introductory_demo_command(args)
    if introductory_exit_code is not None:
        return introductory_exit_code
    benchmark_exit_code = run_benchmark_demo_command(args)
    if benchmark_exit_code is not None:
        return benchmark_exit_code
    sequence_exit_code = run_sequence_demo_command(args)
    if sequence_exit_code is not None:
        return sequence_exit_code
    rabies_exit_code = run_rabies_demo_command(args)
    if rabies_exit_code is not None:
        return rabies_exit_code
    quality_exit_code = run_quality_demo_command(args)
    if quality_exit_code is not None:
        return quality_exit_code
    recovery_exit_code = run_recovery_demo_command(args)
    if recovery_exit_code is not None:
        return recovery_exit_code
    raise NotImplementedError(f"unsupported demo command: {args.demo_command}")
