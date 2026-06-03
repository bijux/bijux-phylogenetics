from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick, write_newick_tree_set
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodEqualBestTreeReport,
)


def write_nucleotide_likelihood_best_tree_set(
    path: Path,
    report: NucleotideLikelihoodEqualBestTreeReport,
) -> Path:
    """Write one ordered retained-best tree set for a native likelihood search."""
    return write_newick_tree_set(
        path,
        [loads_newick(row.tree_newick) for row in report.rows],
    )
