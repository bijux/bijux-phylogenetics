from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
)

BRANCH_LENGTH_PRIOR_FAMILIES = (
    "exponential",
    "fixed",
    "gamma",
    "lognormal",
)


@dataclass(frozen=True, slots=True)
class BranchLengthPriorModel:
    """One validated branch-length prior parameterization."""

    family: str
    rate: float | None = None
    shape: float | None = None
    scale: float | None = None
    log_mean: float | None = None
    log_standard_deviation: float | None = None
    fixed_value: float | None = None
    fixed_tolerance: float | None = None

    def parameter_values(self) -> dict[str, float]:
        parameter_values: dict[str, float] = {}
        if self.rate is not None:
            parameter_values["rate"] = self.rate
        if self.shape is not None:
            parameter_values["shape"] = self.shape
        if self.scale is not None:
            parameter_values["scale"] = self.scale
        if self.log_mean is not None:
            parameter_values["log_mean"] = self.log_mean
        if self.log_standard_deviation is not None:
            parameter_values["log_standard_deviation"] = self.log_standard_deviation
        if self.fixed_value is not None:
            parameter_values["fixed_value"] = self.fixed_value
        if self.fixed_tolerance is not None:
            parameter_values["fixed_tolerance"] = self.fixed_tolerance
        return parameter_values


@dataclass(frozen=True, slots=True)
class BranchLengthPriorBranchRow:
    """One branch-level prior density contribution."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    branch_length: float
    log_density: float


@dataclass(slots=True)
class BranchLengthPriorEvaluationReport:
    """One fixed-tree branch-length prior evaluation report."""

    family: str
    tree_newick: str
    branch_count: int
    parameter_values: dict[str, float]
    total_log_prior: float
    branch_rows: list[BranchLengthPriorBranchRow]


def build_exponential_branch_length_prior(rate: float) -> BranchLengthPriorModel:
    """Build one exponential branch-length prior."""
    return BranchLengthPriorModel(
        family="exponential",
        rate=_validate_positive_finite_parameter(
            parameter_name="rate",
            value=rate,
            family="exponential",
        ),
    )


def build_gamma_branch_length_prior(
    *,
    shape: float,
    scale: float,
) -> BranchLengthPriorModel:
    """Build one gamma branch-length prior."""
    return BranchLengthPriorModel(
        family="gamma",
        shape=_validate_positive_finite_parameter(
            parameter_name="shape",
            value=shape,
            family="gamma",
        ),
        scale=_validate_positive_finite_parameter(
            parameter_name="scale",
            value=scale,
            family="gamma",
        ),
    )


def build_lognormal_branch_length_prior(
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> BranchLengthPriorModel:
    """Build one lognormal branch-length prior."""
    return BranchLengthPriorModel(
        family="lognormal",
        log_mean=_validate_finite_parameter(
            parameter_name="log_mean",
            value=log_mean,
            family="lognormal",
        ),
        log_standard_deviation=_validate_positive_finite_parameter(
            parameter_name="log_standard_deviation",
            value=log_standard_deviation,
            family="lognormal",
        ),
    )


def build_fixed_branch_length_prior(
    *,
    fixed_value: float,
    fixed_tolerance: float = 1e-12,
) -> BranchLengthPriorModel:
    """Build one fixed-value branch-length prior."""
    validated_fixed_value = _validate_nonnegative_finite_parameter(
        parameter_name="fixed_value",
        value=fixed_value,
        family="fixed",
    )
    validated_fixed_tolerance = _validate_nonnegative_finite_parameter(
        parameter_name="fixed_tolerance",
        value=fixed_tolerance,
        family="fixed",
    )
    return BranchLengthPriorModel(
        family="fixed",
        fixed_value=validated_fixed_value,
        fixed_tolerance=validated_fixed_tolerance,
    )


def evaluate_branch_length_log_prior(
    branch_length: float,
    prior_model: BranchLengthPriorModel,
) -> float:
    """Evaluate one branch-length log prior density."""
    validated_branch_length = _validate_branch_length(branch_length)
    if prior_model.family == "exponential":
        rate = require_present(
            prior_model.rate,
            owner_name="branch-length prior evaluation",
            field_name="rate",
        )
        return math.log(rate) - (rate * validated_branch_length)
    if prior_model.family == "gamma":
        shape = require_present(
            prior_model.shape,
            owner_name="branch-length prior evaluation",
            field_name="shape",
        )
        scale = require_present(
            prior_model.scale,
            owner_name="branch-length prior evaluation",
            field_name="scale",
        )
        return _gamma_log_density(
            validated_branch_length,
            shape=shape,
            scale=scale,
        )
    if prior_model.family == "lognormal":
        log_mean = require_present(
            prior_model.log_mean,
            owner_name="branch-length prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_model.log_standard_deviation,
            owner_name="branch-length prior evaluation",
            field_name="log_standard_deviation",
        )
        return _lognormal_log_density(
            validated_branch_length,
            log_mean=log_mean,
            log_standard_deviation=log_standard_deviation,
        )
    if prior_model.family == "fixed":
        fixed_value = require_present(
            prior_model.fixed_value,
            owner_name="branch-length prior evaluation",
            field_name="fixed_value",
        )
        fixed_tolerance = require_present(
            prior_model.fixed_tolerance,
            owner_name="branch-length prior evaluation",
            field_name="fixed_tolerance",
        )
        return (
            0.0
            if math.isclose(
                validated_branch_length,
                fixed_value,
                rel_tol=0.0,
                abs_tol=fixed_tolerance,
            )
            else -math.inf
        )
    raise PhylogeneticsError(
        "branch-length prior family is unsupported",
        code="branch_length_prior_family_invalid",
        details={
            "family": prior_model.family,
            "allowed_families": list(BRANCH_LENGTH_PRIOR_FAMILIES),
        },
    )


def evaluate_tree_branch_length_log_prior(
    tree: PhyloTree,
    prior_model: BranchLengthPriorModel,
) -> BranchLengthPriorEvaluationReport:
    """Evaluate one branch-length prior across every explicit branch in one tree."""
    branch_rows: list[BranchLengthPriorBranchRow] = []
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise InvalidBranchLengthError(
                "branch-length prior evaluation requires stable branch identifiers"
            )
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                "branch-length prior evaluation requires explicit branch lengths on every edge"
            )
        log_density = evaluate_branch_length_log_prior(child.branch_length, prior_model)
        branch_rows.append(
            BranchLengthPriorBranchRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                branch_length=float(child.branch_length),
                log_density=log_density,
            )
        )
    return BranchLengthPriorEvaluationReport(
        family=prior_model.family,
        tree_newick=dumps_newick(tree),
        branch_count=len(branch_rows),
        parameter_values=prior_model.parameter_values(),
        total_log_prior=sum(row.log_density for row in branch_rows),
        branch_rows=branch_rows,
    )


def _gamma_log_density(
    branch_length: float,
    *,
    shape: float,
    scale: float,
) -> float:
    if branch_length == 0.0:
        if shape < 1.0:
            return math.inf
        if shape == 1.0:
            return -math.log(scale)
        return -math.inf
    return (
        (shape - 1.0) * math.log(branch_length)
        - (branch_length / scale)
        - math.lgamma(shape)
        - (shape * math.log(scale))
    )


def _lognormal_log_density(
    branch_length: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    if branch_length == 0.0:
        return -math.inf
    centered_log_length = math.log(branch_length) - log_mean
    return (
        -math.log(branch_length)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (
            centered_log_length
            * centered_log_length
            / (2.0 * log_standard_deviation * log_standard_deviation)
        )
    )


def _validate_branch_length(branch_length: float) -> float:
    if math.isnan(branch_length) or math.isinf(branch_length) or branch_length < 0.0:
        raise InvalidBranchLengthError(
            "branch-length prior evaluation requires finite non-negative branch lengths"
        )
    return float(branch_length)


def _validate_positive_finite_parameter(
    *,
    parameter_name: str,
    value: float,
    family: str,
) -> float:
    validated_value = _validate_finite_parameter(
        parameter_name=parameter_name,
        value=value,
        family=family,
    )
    if validated_value <= 0.0:
        raise PhylogeneticsError(
            f"{family} branch-length prior requires '{parameter_name}' to be positive",
            code="branch_length_prior_parameter_invalid",
            details={
                "family": family,
                "parameter_name": parameter_name,
                "value": value,
            },
        )
    return validated_value


def _validate_nonnegative_finite_parameter(
    *,
    parameter_name: str,
    value: float,
    family: str,
) -> float:
    validated_value = _validate_finite_parameter(
        parameter_name=parameter_name,
        value=value,
        family=family,
    )
    if validated_value < 0.0:
        raise PhylogeneticsError(
            f"{family} branch-length prior requires '{parameter_name}' to be non-negative",
            code="branch_length_prior_parameter_invalid",
            details={
                "family": family,
                "parameter_name": parameter_name,
                "value": value,
            },
        )
    return validated_value


def _validate_finite_parameter(
    *,
    parameter_name: str,
    value: float,
    family: str,
) -> float:
    if math.isnan(value) or math.isinf(value):
        raise PhylogeneticsError(
            f"{family} branch-length prior requires '{parameter_name}' to be finite",
            code="branch_length_prior_parameter_invalid",
            details={
                "family": family,
                "parameter_name": parameter_name,
                "value": value,
            },
        )
    return float(value)


__all__ = [
    "BRANCH_LENGTH_PRIOR_FAMILIES",
    "BranchLengthPriorBranchRow",
    "BranchLengthPriorEvaluationReport",
    "BranchLengthPriorModel",
    "build_exponential_branch_length_prior",
    "build_fixed_branch_length_prior",
    "build_gamma_branch_length_prior",
    "build_lognormal_branch_length_prior",
    "evaluate_branch_length_log_prior",
    "evaluate_tree_branch_length_log_prior",
]
