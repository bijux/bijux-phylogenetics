from __future__ import annotations

from phylogenetics.cli import build_parser


def test_alias_parser_uses_alias_prog() -> None:
    parser = build_parser()
    assert parser.prog == "phylogenetics"

