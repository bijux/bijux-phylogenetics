from __future__ import annotations

import re
from pathlib import Path

from bijux_phylogenetics.errors import InvalidBranchLengthError
from bijux_phylogenetics.io.biopython import load_biophylo, loads_biophylo

_BRANCH_LENGTH_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")


def _validate_branch_lengths(text: str) -> None:
    index = 0
    while True:
        index = text.find(":", index)
        if index == -1:
            return
        cursor = index + 1
        while cursor < len(text) and text[cursor].isspace():
            cursor += 1
        start = cursor
        while cursor < len(text) and text[cursor] not in ",();":
            cursor += 1
        raw_value = text[start:cursor].strip()
        if not raw_value:
            raise InvalidBranchLengthError(f"missing branch length in Newick string at position {start}")
        if not _BRANCH_LENGTH_PATTERN.fullmatch(raw_value):
            raise InvalidBranchLengthError(f"invalid branch length '{raw_value}' in Newick string")
        index = cursor


def loads_newick(text: str):
    """Parse a Newick string into a minimal tree object."""
    _validate_branch_lengths(text)
    return loads_biophylo(text, source_format="newick")


def load_newick(path: Path):
    """Load a Newick tree from disk."""
    return load_biophylo(path, source_format="newick")
