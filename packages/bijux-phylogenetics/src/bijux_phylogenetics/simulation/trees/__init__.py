"""Tree simulation workflows."""

from .generation import (
    simulate_birth_death_trees,
    simulate_coalescent_tree,
    simulate_coalescent_trees,
    simulate_random_tree,
    simulate_random_trees,
    write_coalescent_skyline_table,
    write_coalescent_waiting_time_table,
    write_simulated_tree,
    write_tree_set,
    write_tree_simulation_envelope_table,
    write_tree_simulation_record_table,
)
from .multispecies_coalescent import simulate_multispecies_coalescent_gene_tree
from .multispecies_coalescent_reports import (
    write_multispecies_coalescent_branch_table,
    write_multispecies_coalescent_event_table,
)

__all__ = [
    "simulate_birth_death_trees",
    "simulate_coalescent_tree",
    "simulate_coalescent_trees",
    "simulate_multispecies_coalescent_gene_tree",
    "simulate_random_tree",
    "simulate_random_trees",
    "write_coalescent_skyline_table",
    "write_coalescent_waiting_time_table",
    "write_multispecies_coalescent_branch_table",
    "write_multispecies_coalescent_event_table",
    "write_simulated_tree",
    "write_tree_set",
    "write_tree_simulation_envelope_table",
    "write_tree_simulation_record_table",
]
