from __future__ import annotations

from typing import Any

from .composition import (
    add_alignment_composition_commands,
    run_alignment_composition_command,
)
from .input_validation import (
    add_alignment_input_validation_commands,
    run_alignment_input_validation_command,
)
from .quality import (
    add_alignment_quality_commands,
    run_alignment_quality_command,
)
from .sequence_anomalies import (
    add_alignment_sequence_anomaly_commands,
    run_alignment_sequence_anomaly_command,
)
from .summary import (
    add_alignment_summary_commands,
    run_alignment_summary_command,
)


def add_alignment_review_commands(alignment_subparsers: Any) -> None:
    add_alignment_summary_commands(alignment_subparsers)
    add_alignment_input_validation_commands(alignment_subparsers)
    add_alignment_composition_commands(alignment_subparsers)
    add_alignment_sequence_anomaly_commands(alignment_subparsers)
    add_alignment_quality_commands(alignment_subparsers)


def run_alignment_review_command(args: Any) -> int | None:
    summary_exit_code = run_alignment_summary_command(args)
    if summary_exit_code is not None:
        return summary_exit_code

    input_validation_exit_code = run_alignment_input_validation_command(args)
    if input_validation_exit_code is not None:
        return input_validation_exit_code

    composition_exit_code = run_alignment_composition_command(args)
    if composition_exit_code is not None:
        return composition_exit_code

    sequence_anomaly_exit_code = run_alignment_sequence_anomaly_command(args)
    if sequence_anomaly_exit_code is not None:
        return sequence_anomaly_exit_code

    return run_alignment_quality_command(args)


__all__ = [
    "add_alignment_review_commands",
    "run_alignment_review_command",
]
