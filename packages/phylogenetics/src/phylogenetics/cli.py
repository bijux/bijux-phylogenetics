"""Compatibility CLI entrypoint for the alias distribution."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from typing import Any

from bijux_phylogenetics.cli import build_parser as build_runtime_parser
from bijux_phylogenetics.cli import run_command as dispatch_runtime_command

__all__ = ["build_parser", "main", "run_command"]


def build_parser() -> argparse.ArgumentParser:
    """Build the canonical command-line parser for the alias package."""
    parser = build_runtime_parser()
    parser.prog = "phylogenetics"
    return parser


def run_command(
    args: Any,
    *,
    parser: argparse.ArgumentParser | None = None,
) -> int:
    """Dispatch parsed command-line arguments through runtime handlers."""
    if parser is None:
        parser = build_parser()
    return dispatch_runtime_command(args, parser=parser)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_command(args, parser=parser)


if __name__ == "__main__":
    raise SystemExit(main())

