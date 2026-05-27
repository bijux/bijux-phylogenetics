from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Jc69BranchLengthOptimizationStep:
    """One bounded branch-length update inside a fixed-topology JC69 search."""

    optimization_pass: int
    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    starting_branch_length: float
    optimized_branch_length: float
    starting_log_likelihood: float
    optimized_log_likelihood: float
    accepted: bool


@dataclass(slots=True)
class Jc69TreeLikelihoodReport:
    """Native JC69 likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    log_likelihood: float


@dataclass(slots=True)
class Jc69BranchLengthOptimizationReport:
    """Fixed-topology JC69 branch-length optimization summary."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    initial_tree_newick: str
    optimized_tree_newick: str
    initial_log_likelihood: float
    optimized_log_likelihood: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool
    lower_branch_length_bound: float
    upper_branch_length_bound: float
    steps: list[Jc69BranchLengthOptimizationStep]
