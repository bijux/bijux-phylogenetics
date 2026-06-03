"""Compatibility CLI facade for the canonical command-line bootstrap."""

from __future__ import annotations

from .command_line.bootstrap import build_parser, main, run_command

__all__ = ["build_parser", "main", "run_command"]
