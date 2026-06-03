from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def count_expected_output_entries(root: Path) -> int:
    return sum(1 for _ in root.glob("*"))


def count_expected_output_files(root: Path) -> int:
    return sum(1 for path in root.rglob("*") if path.is_file())


def emit_demo_result(
    args: Any,
    *,
    outputs: list[Path | str],
    metrics: dict[str, Any],
    data: Any,
    output_root: Path,
) -> int:
    finalized_outputs = _finalize_outputs(
        args,
        command="demo",
        inputs=[],
        outputs=outputs,
    )
    if args.json:
        _print_result(
            build_command_result(
                command="demo",
                inputs=[],
                outputs=finalized_outputs,
                metrics=metrics,
                data=data,
            ),
            json_output=True,
        )
        return 0
    print(output_root)
    return 0


def resolve_demo_runner(name: str) -> Any:
    import bijux_phylogenetics.command_line.demo as demo_command_module

    return getattr(demo_command_module, name)
