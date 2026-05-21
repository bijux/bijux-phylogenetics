from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .budgets import (
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    enforce_tree_set_tree_budget,
)
from .consensus import (
    _build_consensus_tree_with_threshold,
    compute_consensus_tree,
    compute_consensus_tree_with_threshold,
    compute_strict_consensus_tree,
    write_consensus_tree,
)
from .distances import (
    _build_tree_distance_matrix_report,
    compute_tree_distance_matrix,
    write_tree_distance_matrix,
)
from .clade_support import (
    _build_clade_frequency_report,
    _support_classification,
    compute_clade_frequency_table,
    compute_reference_tree_clade_support,
    write_clade_frequency_table,
    write_reference_tree_clade_support_table,
)
from .contracts import ConsensusTreeReport, TreeDistanceMatrixReport, TreeDistancePair
from .contracts import TreeSetReport as TreeSetReport
from .inventory import (
    _TreeSetAnalysis,
    _analyze_tree_set,
    _require_exact_taxa,
    _require_tree_set,
    _validate_same_taxa,
    load_tree_set,
)
from .topology import (
    _clade_counts,
    _clade_signature,
    _clades_conflict,
    _format_clade,
    _rooted_topology_id,
    _tree_distance,
)



