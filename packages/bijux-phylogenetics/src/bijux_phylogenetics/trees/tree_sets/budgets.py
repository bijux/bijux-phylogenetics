from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import WorkflowBudgetError

from .contracts import TreeSetWorkflowBudget, TreeSetWorkflowBudgetReport


def _validate_budget_limit(
    value: int | None,
    *,
    name: str,
) -> int | None:
    if value is None:
        return None
    if value < 1:
        raise ValueError(f"{name} must be at least 1, got {value}")
    return value


def build_tree_set_workflow_budget(
    *,
    max_tree_count: int | None = None,
    max_report_table_rows: int | None = None,
    memory_warning_threshold_bytes: int | None = None,
) -> TreeSetWorkflowBudget:
    """Normalize one reviewer-facing resource budget for tree-set workflows."""
    validated_threshold = (
        None
        if memory_warning_threshold_bytes is None
        else _validate_budget_limit(
            memory_warning_threshold_bytes,
            name="memory_warning_threshold_bytes",
        )
    )
    return TreeSetWorkflowBudget(
        max_tree_count=_validate_budget_limit(
            max_tree_count,
            name="max_tree_count",
        ),
        max_report_table_rows=_validate_budget_limit(
            max_report_table_rows,
            name="max_report_table_rows",
        ),
        memory_warning_threshold_bytes=validated_threshold,
    )


def enforce_tree_set_tree_budget(
    *,
    tree_count: int,
    budget: TreeSetWorkflowBudget,
    workflow_name: str,
    source_path: Path,
) -> None:
    """Reject tree-set workflows that exceed an explicit input-size budget."""
    if budget.max_tree_count is None or tree_count <= budget.max_tree_count:
        return
    raise WorkflowBudgetError(
        (
            f"{workflow_name} budget allows at most {budget.max_tree_count} trees, "
            f"but {source_path} contains {tree_count}"
        ),
        code="tree_set_tree_budget_exceeded",
        details={
            "workflow_name": workflow_name,
            "source_path": str(source_path),
            "tree_count": tree_count,
            "max_tree_count": budget.max_tree_count,
        },
    )


def build_tree_set_budget_report(
    *,
    budget: TreeSetWorkflowBudget,
    peak_memory_bytes: int,
    truncated_section_names: list[str] | None = None,
) -> TreeSetWorkflowBudgetReport:
    """Summarize how one tree-set workflow budget was applied."""
    warning_messages: list[str] = []
    if (
        budget.memory_warning_threshold_bytes is not None
        and peak_memory_bytes > budget.memory_warning_threshold_bytes
    ):
        warning_messages.append(
            "peak memory exceeded the configured workflow warning threshold"
        )
    truncated_names = sorted(set(truncated_section_names or []))
    if truncated_names:
        warning_messages.append(
            "reviewer-facing sections were truncated to the configured row limit"
        )
    return TreeSetWorkflowBudgetReport(
        max_tree_count=budget.max_tree_count,
        max_report_table_rows=budget.max_report_table_rows,
        memory_warning_threshold_bytes=budget.memory_warning_threshold_bytes,
        truncated_section_names=truncated_names,
        warning_messages=warning_messages,
    )
