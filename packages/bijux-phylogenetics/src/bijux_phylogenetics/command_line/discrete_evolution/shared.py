from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _split_csv_values

COMMAND_NAME = "discrete-evolution"


def ordered_states(args: Any) -> list[str] | None:
    return _split_csv_values(args.ordered_states) or None


def allowed_states(args: Any) -> list[str] | None:
    if not hasattr(args, "allowed_states"):
        return None
    return _split_csv_values(args.allowed_states) or None


def model_inputs(args: Any) -> list[Path]:
    return [args.tree, args.table]


def render_density_outputs(outputs: list[Path | str], density_result: Any) -> None:
    outputs.extend(
        [density_result.output_path, density_result.svg_path]
        if density_result.format == "html"
        else [density_result.output_path]
    )
