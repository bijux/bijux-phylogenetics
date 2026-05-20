from __future__ import annotations

from pathlib import Path
from typing import Any


def tree_and_metadata_inputs(args: Any) -> list[Path]:
    inputs = [args.tree]
    if getattr(args, "metadata", None) is not None:
        inputs.append(args.metadata)
    return inputs


def tree_metadata_traits_inputs(args: Any) -> list[Path]:
    inputs = tree_and_metadata_inputs(args)
    if getattr(args, "traits", None) is not None:
        inputs.append(args.traits)
    return inputs
