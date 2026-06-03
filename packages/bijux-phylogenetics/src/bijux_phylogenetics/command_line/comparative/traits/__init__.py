from __future__ import annotations

from typing import Any

from .anova import (
    add_phylogenetic_anova_review_command,
    run_phylogenetic_anova_review_command,
)
from .clade_traits import (
    add_clade_trait_review_command,
    run_clade_trait_review_command,
)
from .imputation import (
    add_trait_imputation_review_command,
    run_trait_imputation_review_command,
)
from .outliers import (
    add_trait_outlier_review_command,
    run_trait_outlier_review_command,
)
from .residuals import (
    add_phylogenetic_residual_review_command,
    run_phylogenetic_residual_review_command,
)


def add_comparative_trait_commands(comparative_subparsers: Any) -> None:
    add_clade_trait_review_command(comparative_subparsers)
    add_phylogenetic_residual_review_command(comparative_subparsers)
    add_phylogenetic_anova_review_command(comparative_subparsers)
    add_trait_outlier_review_command(comparative_subparsers)
    add_trait_imputation_review_command(comparative_subparsers)


def run_comparative_trait_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    clade_traits_exit_code = run_clade_trait_review_command(args)
    if clade_traits_exit_code is not None:
        return clade_traits_exit_code

    residual_exit_code = run_phylogenetic_residual_review_command(args)
    if residual_exit_code is not None:
        return residual_exit_code

    anova_exit_code = run_phylogenetic_anova_review_command(args)
    if anova_exit_code is not None:
        return anova_exit_code

    outlier_exit_code = run_trait_outlier_review_command(args)
    if outlier_exit_code is not None:
        return outlier_exit_code

    return run_trait_imputation_review_command(args)


__all__ = [
    "add_comparative_trait_commands",
    "run_comparative_trait_command",
]
