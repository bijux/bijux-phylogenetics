from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.biopython import load_biophylo


def load_nexus(path: Path):
    """Load a NEXUS tree from disk."""
    return load_biophylo(path, source_format="nexus")
