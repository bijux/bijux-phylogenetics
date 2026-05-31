from __future__ import annotations

import math
import random

from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodStartingTreePoolReport,
    NucleotideLikelihoodStartingTreeSummary,
)

_SUPPORTED_STARTING_TREE_SELECTION_POLICIES = frozenset(
    {"all", "best", "random-k", "strategy-priority"}
)
_SUPPORTED_STARTING_TREE_STRATEGIES = frozenset(
    {"input-tree", "likelihood-stepwise-addition-tree", "random-tree"}
)
_DEFAULT_STARTING_TREE_STRATEGY_PRIORITY = (
    "likelihood-stepwise-addition-tree",
    "input-tree",
    "random-tree",
)


def validate_nucleotide_likelihood_starting_tree_selection_policy(policy: str) -> str:
    """Validate the supported selection policies over one scored starting-tree pool."""
    normalized_policy = policy.strip().lower()
    if normalized_policy not in _SUPPORTED_STARTING_TREE_SELECTION_POLICIES:
        raise ValueError(
            "starting_tree_selection_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_STARTING_TREE_SELECTION_POLICIES))
        )
    return normalized_policy


def validate_nucleotide_likelihood_starting_tree_selection_count(
    policy: str,
    selected_start_tree_count: int | None,
) -> int | None:
    """Validate optional selection counts against the declared starting-tree policy."""
    normalized_policy = validate_nucleotide_likelihood_starting_tree_selection_policy(
        policy
    )
    if normalized_policy in {"all", "best"}:
        if selected_start_tree_count is not None:
            raise ValueError(
                "selected_start_tree_count is supported only for 'random-k' and 'strategy-priority'"
            )
        return None
    if selected_start_tree_count is None:
        if normalized_policy == "random-k":
            raise ValueError(
                "selected_start_tree_count is required when starting_tree_selection_policy is 'random-k'"
            )
        return None
    if selected_start_tree_count < 1:
        raise ValueError("selected_start_tree_count must be at least one when provided")
    return selected_start_tree_count


def validate_nucleotide_likelihood_starting_tree_strategy_priority(
    strategy_priority: tuple[str, ...] | list[str] | None,
) -> tuple[str, ...]:
    """Normalize one declared strategy-priority order over starting-tree sources."""
    if strategy_priority is None:
        return _DEFAULT_STARTING_TREE_STRATEGY_PRIORITY
    normalized_priority = tuple(strategy.strip().lower() for strategy in strategy_priority)
    if not normalized_priority:
        raise ValueError("strategy_priority must not be empty when provided")
    if len(set(normalized_priority)) != len(normalized_priority):
        raise ValueError("strategy_priority must not repeat strategies")
    unsupported = sorted(
        strategy
        for strategy in normalized_priority
        if strategy not in _SUPPORTED_STARTING_TREE_STRATEGIES
    )
    if unsupported:
        raise ValueError(
            "strategy_priority contains unsupported strategies: "
            + ", ".join(unsupported)
        )
    return normalized_priority


def select_nucleotide_likelihood_starting_tree_pool(
    report: NucleotideLikelihoodStartingTreePoolReport,
    *,
    starting_tree_selection_policy: str,
    selected_start_tree_count: int | None = None,
    selection_seed: int = 1,
    strategy_priority: tuple[str, ...] | list[str] | None = None,
) -> list[NucleotideLikelihoodStartingTreeSummary]:
    """Select one explicit subset of scored starting trees from one deterministic pool."""
    normalized_policy = validate_nucleotide_likelihood_starting_tree_selection_policy(
        starting_tree_selection_policy
    )
    validated_selection_count = (
        validate_nucleotide_likelihood_starting_tree_selection_count(
            normalized_policy,
            selected_start_tree_count,
        )
    )
    pool_rows = list(report.starting_tree_summaries)
    if not pool_rows:
        raise ValueError("starting_tree_summaries must not be empty")
    if normalized_policy == "all":
        return pool_rows
    if normalized_policy == "best":
        return [_best_starting_tree_summary(pool_rows)]
    if normalized_policy == "random-k":
        return _select_random_starting_tree_subset(
            pool_rows,
            selected_start_tree_count=validated_selection_count,
            selection_seed=selection_seed,
        )
    validated_priority = validate_nucleotide_likelihood_starting_tree_strategy_priority(
        strategy_priority
    )
    return _select_strategy_priority_starting_tree_subset(
        pool_rows,
        strategy_priority=validated_priority,
        selected_start_tree_count=validated_selection_count,
    )


def _best_starting_tree_summary(
    pool_rows: list[NucleotideLikelihoodStartingTreeSummary],
) -> NucleotideLikelihoodStartingTreeSummary:
    best_row = pool_rows[0]
    for candidate in pool_rows[1:]:
        if _prefer_starting_tree_summary(candidate, best_row):
            best_row = candidate
    return best_row


def _select_random_starting_tree_subset(
    pool_rows: list[NucleotideLikelihoodStartingTreeSummary],
    *,
    selected_start_tree_count: int | None,
    selection_seed: int,
) -> list[NucleotideLikelihoodStartingTreeSummary]:
    if selected_start_tree_count is None:
        raise ValueError(
            "selected_start_tree_count is required when starting_tree_selection_policy is 'random-k'"
        )
    if selected_start_tree_count > len(pool_rows):
        raise ValueError(
            "selected_start_tree_count must not exceed the scored starting-tree pool size"
        )
    selected_indexes = sorted(
        random.Random(selection_seed).sample(
            range(len(pool_rows)),
            selected_start_tree_count,
        )
    )
    return [pool_rows[index] for index in selected_indexes]


def _select_strategy_priority_starting_tree_subset(
    pool_rows: list[NucleotideLikelihoodStartingTreeSummary],
    *,
    strategy_priority: tuple[str, ...],
    selected_start_tree_count: int | None,
) -> list[NucleotideLikelihoodStartingTreeSummary]:
    selected_rows: list[NucleotideLikelihoodStartingTreeSummary] = []
    for strategy in strategy_priority:
        strategy_rows = [
            row for row in pool_rows if row.source_strategy == strategy
        ]
        if not strategy_rows:
            continue
        selected_rows.append(_best_starting_tree_summary(strategy_rows))
        if (
            selected_start_tree_count is not None
            and len(selected_rows) >= selected_start_tree_count
        ):
            break
    if not selected_rows:
        raise ValueError(
            "strategy_priority did not match any starting-tree strategies in the scored pool"
        )
    return selected_rows


def _prefer_starting_tree_summary(
    left: NucleotideLikelihoodStartingTreeSummary,
    right: NucleotideLikelihoodStartingTreeSummary,
) -> bool:
    if (
        left.starting_log_likelihood > right.starting_log_likelihood
        and not math.isclose(
            left.starting_log_likelihood,
            right.starting_log_likelihood,
        )
    ):
        return True
    if (
        right.starting_log_likelihood > left.starting_log_likelihood
        and not math.isclose(
            left.starting_log_likelihood,
            right.starting_log_likelihood,
        )
    ):
        return False
    if left.topology_hash != right.topology_hash:
        return left.topology_hash < right.topology_hash
    if left.tree_id != right.tree_id:
        return left.tree_id < right.tree_id
    return left.tree_newick < right.tree_newick
