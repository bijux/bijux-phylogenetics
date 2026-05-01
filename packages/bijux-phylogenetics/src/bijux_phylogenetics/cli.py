from __future__ import annotations

import argparse
from collections.abc import Sequence
import sys
from typing import Any

from .command_line import build_parser as build_command_parser
from .command_line import run_command as dispatch_command

__all__ = ["build_parser", "main", "run_command"]


def build_parser() -> argparse.ArgumentParser:
    """Build the canonical command-line argument parser."""
    return build_command_parser()


def run_command(
    args: Any,
    *,
    parser: argparse.ArgumentParser | None = None,
) -> int:
    """Dispatch parsed command-line arguments through the CLI handler registry."""
    if parser is None:
        parser = build_parser()
    return dispatch_command(args, parser=parser)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = list(argv) if argv is not None else list(sys.argv[1:])
    return run_command(args, parser=parser)


if __name__ == "__main__":
    raise SystemExit(main())
