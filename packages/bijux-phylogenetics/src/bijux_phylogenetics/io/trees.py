from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick
from bijux_phylogenetics.io.nexus import load_nexus
from bijux_phylogenetics.io.phyloxml import load_phyloxml
from bijux_phylogenetics.runtime.errors import UnsupportedTreeFormatError

TreeFormat = str

_SUFFIX_FORMATS: dict[str, TreeFormat] = {
    ".nwk": "newick",
    ".newick": "newick",
    ".tree": "newick",
    ".tre": "newick",
    ".nex": "nexus",
    ".nexus": "nexus",
    ".xml": "phyloxml",
    ".phyloxml": "phyloxml",
}


def detect_tree_format(path: Path) -> TreeFormat:
    """Detect the tree format from filename and file content."""
    if not path.exists():
        suffix = path.suffix.lower()
        if suffix in _SUFFIX_FORMATS:
            return _SUFFIX_FORMATS[suffix]
        raise FileNotFoundError(f"tree file not found: {path}")

    prefix = path.read_text(encoding="utf-8")[:256].lstrip().lower()
    if prefix.startswith("#nexus"):
        return "nexus"
    if prefix.startswith(("<?xml", "<phyloxml")):
        return "phyloxml"
    if prefix.startswith("("):
        return "newick"
    suffix = path.suffix.lower()
    if suffix in _SUFFIX_FORMATS:
        return _SUFFIX_FORMATS[suffix]
    raise UnsupportedTreeFormatError(f"unsupported tree format for {path}")


def load_tree(path: Path, *, source_format: TreeFormat | None = None):
    """Load a tree using either explicit or auto-detected format."""
    tree_format = source_format or detect_tree_format(path)
    if tree_format == "newick":
        return load_newick(path)
    if tree_format == "nexus":
        return load_nexus(path)
    if tree_format == "phyloxml":
        return load_phyloxml(path)
    raise UnsupportedTreeFormatError(f"unsupported tree format for {path}")
