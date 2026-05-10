from __future__ import annotations

import argparse
from importlib import import_module
from pathlib import Path
import sys
from typing import Any, cast
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[2] / "bijux-phylogenetics" / "src"),
)

from bijux_phylogenetics.command_line.registry import (
    get_command_spec as get_runtime_command_spec,
)
from bijux_phylogenetics.comparative.models import (
    ComparativeModelComparisonReport as RuntimeComparativeModelComparisonReport,
)
from phylogenetic import ComparativeDataset, __version__
from phylogenetic.cli import build_parser

get_alias_command_spec = cast(
    Any,
    import_module("phylogenetic.command_line.registry").get_command_spec,
)
AliasComparativeModelComparisonReport = cast(
    Any,
    import_module("phylogenetic.comparative.models").ComparativeModelComparisonReport,
)


class PhylogeneticCompatibilityTests(unittest.TestCase):
    def test_import_surface_re_exports_runtime_api(self) -> None:
        self.assertTrue(callable(ComparativeDataset))
        self.assertIsInstance(__version__, str)

    def test_cli_parser_supports_runtime_commands(self) -> None:
        parser = build_parser()

        self.assertEqual(parser.prog, "phylogenetic")
        alignment_args = parser.parse_args(["alignment", "profiles", "--json"])

        self.assertEqual(alignment_args.command, "alignment")
        self.assertEqual(alignment_args.alignment_command, "profiles")

    def test_cli_parser_command_set_matches_runtime_parser(self) -> None:
        alias_parser = build_parser()
        runtime_parser = build_runtime_parser()

        alias_commands = _extract_command_choices(alias_parser)
        runtime_commands = _extract_command_choices(runtime_parser)

        self.assertEqual(alias_commands, runtime_commands)

    def test_command_line_registry_module_aliases_runtime_module_identity(self) -> None:
        self.assertIs(get_alias_command_spec, get_runtime_command_spec)

    def test_nested_runtime_types_keep_identity_under_alias_imports(self) -> None:
        self.assertIs(
            AliasComparativeModelComparisonReport,
            RuntimeComparativeModelComparisonReport,
        )


def _extract_command_choices(parser: argparse.ArgumentParser) -> set[str]:
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):  # noqa: SLF001
            return {str(name) for name in action.choices}
    return set()


def build_runtime_parser() -> argparse.ArgumentParser:
    from bijux_phylogenetics.cli import build_parser as build_runtime_parser_impl

    return build_runtime_parser_impl()


if __name__ == "__main__":
    unittest.main()
