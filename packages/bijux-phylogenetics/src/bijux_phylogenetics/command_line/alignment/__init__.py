from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec

from .coding import (
    add_alignment_coding_commands,
    run_alignment_coding_command,
)
from .distance import (
    add_alignment_distance_commands,
    run_alignment_distance_command,
)
from .linkage import (
    add_alignment_linkage_commands,
    run_alignment_linkage_command,
)
from .matrix import (
    add_alignment_matrix_commands,
    run_alignment_matrix_command,
)
from .review import (
    add_alignment_review_commands,
    run_alignment_review_command,
)


def add_alignment_commands(subparsers: Any) -> None:
    alignment = subparsers.add_parser(
        get_command_spec("alignment").name,
        help=get_command_spec("alignment").summary,
    )
    alignment_subparsers = alignment.add_subparsers(
        dest="alignment_command",
        required=True,
    )
    add_alignment_review_commands(alignment_subparsers)
    add_alignment_matrix_commands(alignment_subparsers)
    add_alignment_distance_commands(alignment_subparsers)
    add_alignment_coding_commands(alignment_subparsers)
    add_alignment_linkage_commands(alignment_subparsers)


def run_alignment_command(args: Any) -> int | None:
    review_exit_code = run_alignment_review_command(args)
    if review_exit_code is not None:
        return review_exit_code
    matrix_exit_code = run_alignment_matrix_command(args)
    if matrix_exit_code is not None:
        return matrix_exit_code
    distance_exit_code = run_alignment_distance_command(args)
    if distance_exit_code is not None:
        return distance_exit_code
    coding_exit_code = run_alignment_coding_command(args)
    if coding_exit_code is not None:
        return coding_exit_code
    return run_alignment_linkage_command(args)
