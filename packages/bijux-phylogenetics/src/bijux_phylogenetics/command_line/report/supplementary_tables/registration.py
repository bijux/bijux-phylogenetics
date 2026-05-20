from __future__ import annotations

from typing import Any

from .ancestral_states import (
    add_ancestral_state_supplementary_table_commands,
    run_ancestral_state_supplementary_table_command,
)
from .batch_summary import (
    add_batch_summary_supplementary_table_commands,
    run_batch_summary_supplementary_table_command,
)
from .comparative_models import (
    add_comparative_model_supplementary_table_commands,
    run_comparative_model_supplementary_table_command,
)
from .diversification import (
    add_diversification_supplementary_table_commands,
    run_diversification_supplementary_table_command,
)
from .model_selection import (
    add_model_selection_supplementary_table_commands,
    run_model_selection_supplementary_table_command,
)
from .study_inputs import (
    add_study_input_supplementary_table_commands,
    run_study_input_supplementary_table_command,
)
from .tree_review import (
    add_tree_review_supplementary_table_commands,
    run_tree_review_supplementary_table_command,
)


def add_supplementary_table_report_commands(report_subparsers: Any) -> None:
    add_study_input_supplementary_table_commands(report_subparsers)
    add_tree_review_supplementary_table_commands(report_subparsers)
    add_model_selection_supplementary_table_commands(report_subparsers)
    add_comparative_model_supplementary_table_commands(report_subparsers)
    add_ancestral_state_supplementary_table_commands(report_subparsers)
    add_diversification_supplementary_table_commands(report_subparsers)
    add_batch_summary_supplementary_table_commands(report_subparsers)


def run_supplementary_table_report_command(args: Any) -> int | None:
    study_input_result = run_study_input_supplementary_table_command(args)
    if study_input_result is not None:
        return study_input_result

    tree_review_result = run_tree_review_supplementary_table_command(args)
    if tree_review_result is not None:
        return tree_review_result

    model_selection_result = run_model_selection_supplementary_table_command(args)
    if model_selection_result is not None:
        return model_selection_result

    comparative_model_result = run_comparative_model_supplementary_table_command(args)
    if comparative_model_result is not None:
        return comparative_model_result

    ancestral_state_result = run_ancestral_state_supplementary_table_command(args)
    if ancestral_state_result is not None:
        return ancestral_state_result

    diversification_result = run_diversification_supplementary_table_command(args)
    if diversification_result is not None:
        return diversification_result

    return run_batch_summary_supplementary_table_command(args)
