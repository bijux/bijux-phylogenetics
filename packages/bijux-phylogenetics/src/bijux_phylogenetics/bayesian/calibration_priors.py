from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.bayesian.beast import validate_fossil_calibration_table
from bijux_phylogenetics.bayesian.beast.models import ValidatedCalibration
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    PhylogeneticsError,
    UnrootedTreeError,
)

CALIBRATION_PRIOR_FAMILIES = (
    "fixed",
    "uniform",
    "normal",
    "lognormal",
    "offset-exponential",
    "offset-lognormal",
)

_NORMAL_QUANTILE_97_5 = 1.959963984540054
_FIXED_CALIBRATION_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class CalibrationPriorDefinition:
    """One validated calibration prior resolved onto one rooted tree node."""

    calibration_id: str
    requested_distribution: str
    family: str
    target_kind: str
    target_label: str
    descendant_taxa: list[str]
    node_id: str
    node_kind: str
    minimum_age: float | None
    maximum_age: float | None
    translated: bool
    translation_note: str | None
    fixed_age: float | None = None
    fixed_tolerance: float | None = None
    mean_age: float | None = None
    standard_deviation: float | None = None
    log_mean: float | None = None
    log_standard_deviation: float | None = None
    offset_age: float | None = None
    exponential_mean: float | None = None

    def parameter_values(self) -> dict[str, float]:
        """Return the explicit density parameters for this calibration prior."""
        parameter_values: dict[str, float] = {}
        if self.fixed_age is not None:
            parameter_values["fixed_age"] = self.fixed_age
        if self.fixed_tolerance is not None:
            parameter_values["fixed_tolerance"] = self.fixed_tolerance
        if self.mean_age is not None:
            parameter_values["mean_age"] = self.mean_age
        if self.standard_deviation is not None:
            parameter_values["standard_deviation"] = self.standard_deviation
        if self.log_mean is not None:
            parameter_values["log_mean"] = self.log_mean
        if self.log_standard_deviation is not None:
            parameter_values["log_standard_deviation"] = self.log_standard_deviation
        if self.offset_age is not None:
            parameter_values["offset_age"] = self.offset_age
        if self.exponential_mean is not None:
            parameter_values["exponential_mean"] = self.exponential_mean
        return parameter_values


@dataclass(frozen=True, slots=True)
class CalibrationPriorRow:
    """One calibration-prior contribution row."""

    calibration_id: str
    requested_distribution: str
    family: str
    target_kind: str
    target_label: str
    descendant_taxa: list[str]
    node_id: str
    node_kind: str
    node_age: float
    minimum_age: float | None
    maximum_age: float | None
    translated: bool
    translation_note: str | None
    parameter_values: dict[str, float]
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class CalibrationPriorEvaluationReport:
    """One rooted ultrametric calibration-prior evaluation report."""

    tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    calibration_count: int
    translated_calibration_count: int
    total_log_prior: float
    ultrametric_tolerance: float
    calibration_rows: list[CalibrationPriorRow]


def load_calibration_prior_definitions(
    tree_path: Path,
    calibration_path: Path,
) -> list[CalibrationPriorDefinition]:
    """Load and validate calibration priors resolved onto one rooted tree."""
    tree = load_tree(tree_path)
    if not _tree_is_rooted(tree):
        raise UnrootedTreeError(
            "calibration prior evaluation requires a rooted tree",
            code="calibration_prior_requires_rooted_tree",
            details={"rooted": tree.rooted},
        )
    validation_report = validate_fossil_calibration_table(tree_path, calibration_path)
    node_by_descendant_taxa = {
        tuple(node.descendant_taxa): node for node in tree.iter_nodes(order="preorder")
    }
    resolved_calibrations: list[tuple[ValidatedCalibration, TreeNode]] = []
    node_by_calibration_id: dict[str, TreeNode] = {}
    node_kind_by_calibration_id: dict[str, str] = {}
    for calibration in validation_report.calibrations:
        if not calibration.valid:
            continue
        node = node_by_descendant_taxa.get(tuple(sorted(calibration.taxa)))
        if node is None or node.node_id is None:
            raise PhylogeneticsError(
                "valid calibration could not be resolved onto one rooted tree node",
                code="calibration_prior_unresolved_calibration",
                details={"calibration_id": calibration.calibration_id},
            )
        resolved_calibrations.append((calibration, node))
        node_by_calibration_id[calibration.calibration_id] = node
        node_kind_by_calibration_id[calibration.calibration_id] = (
            "root" if node is tree.root else ("tip" if node.is_leaf() else "internal")
        )
    _require_feasible_calibration_age_constraints(resolved_calibrations)
    prior_definitions: list[CalibrationPriorDefinition] = []
    for calibration in validation_report.calibrations:
        if not calibration.valid:
            continue
        prior_definitions.append(
            _build_calibration_prior_definition(
                calibration=calibration,
                node_id=node_by_calibration_id[calibration.calibration_id].node_id
                or "",
                node_kind=node_kind_by_calibration_id[calibration.calibration_id],
            )
        )
    if not prior_definitions:
        raise PhylogeneticsError(
            "calibration prior evaluation requires at least one valid calibration",
            code="calibration_prior_missing_valid_calibrations",
            details={
                "tree_path": str(tree_path),
                "calibration_path": str(calibration_path),
            },
        )
    return sorted(prior_definitions, key=lambda row: row.calibration_id)


def evaluate_calibration_tree_log_prior(
    tree: PhyloTree,
    prior_definitions: list[CalibrationPriorDefinition],
    *,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> CalibrationPriorEvaluationReport:
    """Evaluate one rooted ultrametric tree under one resolved calibration set."""
    _validate_rooted_ultrametric_tree(
        tree,
        prior_name="calibration prior",
        ultrametric_tolerance=ultrametric_tolerance,
    )
    age_by_node_id = _compute_node_age_by_id(
        tree,
        ultrametric_tolerance=ultrametric_tolerance,
    )
    node_by_descendant_taxa = {
        tuple(node.descendant_taxa): node for node in tree.iter_nodes(order="preorder")
    }
    calibration_rows: list[CalibrationPriorRow] = []
    for prior_definition in prior_definitions:
        node = node_by_descendant_taxa.get(tuple(prior_definition.descendant_taxa))
        if node is None or node.node_id is None:
            raise PhylogeneticsError(
                "calibration prior target is absent from the evaluation tree",
                code="calibration_prior_target_missing_from_tree",
                details={
                    "calibration_id": prior_definition.calibration_id,
                    "descendant_taxa": prior_definition.descendant_taxa,
                },
            )
        node_age = age_by_node_id[node.node_id]
        log_prior_contribution = _evaluate_calibration_log_prior(
            node_age=node_age,
            prior_definition=prior_definition,
        )
        calibration_rows.append(
            CalibrationPriorRow(
                calibration_id=prior_definition.calibration_id,
                requested_distribution=prior_definition.requested_distribution,
                family=prior_definition.family,
                target_kind=prior_definition.target_kind,
                target_label=prior_definition.target_label,
                descendant_taxa=prior_definition.descendant_taxa,
                node_id=node.node_id,
                node_kind="root" if node is tree.root else prior_definition.node_kind,
                node_age=float(format(node_age, ".15g")),
                minimum_age=prior_definition.minimum_age,
                maximum_age=prior_definition.maximum_age,
                translated=prior_definition.translated,
                translation_note=prior_definition.translation_note,
                parameter_values=prior_definition.parameter_values(),
                log_prior_contribution=float(format(log_prior_contribution, ".15g")),
            )
        )
    total_log_prior = sum(row.log_prior_contribution for row in calibration_rows)
    return CalibrationPriorEvaluationReport(
        tree_newick=dumps_newick(tree),
        taxa=sorted(tree.tip_names),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        calibration_count=len(calibration_rows),
        translated_calibration_count=sum(
            1 for row in calibration_rows if row.translated
        ),
        total_log_prior=float(format(total_log_prior, ".15g")),
        ultrametric_tolerance=ultrametric_tolerance,
        calibration_rows=calibration_rows,
    )


def _build_calibration_prior_definition(
    *,
    calibration: ValidatedCalibration,
    node_id: str,
    node_kind: str,
) -> CalibrationPriorDefinition:
    requested_distribution = calibration.distribution.strip().casefold() or "uniform"
    minimum_age = calibration.minimum_age
    maximum_age = calibration.maximum_age
    if minimum_age is not None and maximum_age is not None:
        if math.isclose(
            minimum_age,
            maximum_age,
            rel_tol=0.0,
            abs_tol=_FIXED_CALIBRATION_TOLERANCE,
        ):
            return CalibrationPriorDefinition(
                calibration_id=calibration.calibration_id,
                requested_distribution=requested_distribution,
                family="fixed",
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=sorted(calibration.taxa),
                node_id=node_id,
                node_kind=node_kind,
                minimum_age=minimum_age,
                maximum_age=maximum_age,
                translated=(requested_distribution != "fixed"),
                translation_note=(
                    "equal calibration bounds collapse to one fixed calibration prior"
                    if requested_distribution != "fixed"
                    else None
                ),
                fixed_age=float(format(minimum_age, ".15g")),
                fixed_tolerance=_FIXED_CALIBRATION_TOLERANCE,
            )
        if requested_distribution == "uniform":
            return CalibrationPriorDefinition(
                calibration_id=calibration.calibration_id,
                requested_distribution=requested_distribution,
                family="uniform",
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=sorted(calibration.taxa),
                node_id=node_id,
                node_kind=node_kind,
                minimum_age=minimum_age,
                maximum_age=maximum_age,
                translated=False,
                translation_note=None,
            )
        if requested_distribution == "normal":
            mean_age = (minimum_age + maximum_age) / 2.0
            standard_deviation = (maximum_age - minimum_age) / (
                2.0 * _NORMAL_QUANTILE_97_5
            )
            return CalibrationPriorDefinition(
                calibration_id=calibration.calibration_id,
                requested_distribution=requested_distribution,
                family="normal",
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=sorted(calibration.taxa),
                node_id=node_id,
                node_kind=node_kind,
                minimum_age=minimum_age,
                maximum_age=maximum_age,
                translated=False,
                translation_note=None,
                mean_age=float(format(mean_age, ".15g")),
                standard_deviation=float(format(standard_deviation, ".15g")),
            )
        if requested_distribution == "lognormal":
            if minimum_age <= 0.0:
                raise PhylogeneticsError(
                    "lognormal calibration priors require strictly positive bounded ages",
                    code="calibration_prior_invalid_lognormal_bounds",
                    details={
                        "calibration_id": calibration.calibration_id,
                        "minimum_age": minimum_age,
                        "maximum_age": maximum_age,
                    },
                )
            log_mean = (math.log(minimum_age) + math.log(maximum_age)) / 2.0
            log_standard_deviation = (math.log(maximum_age) - math.log(minimum_age)) / (
                2.0 * _NORMAL_QUANTILE_97_5
            )
            return CalibrationPriorDefinition(
                calibration_id=calibration.calibration_id,
                requested_distribution=requested_distribution,
                family="lognormal",
                target_kind=calibration.target_kind,
                target_label=calibration.target_label,
                descendant_taxa=sorted(calibration.taxa),
                node_id=node_id,
                node_kind=node_kind,
                minimum_age=minimum_age,
                maximum_age=maximum_age,
                translated=False,
                translation_note=None,
                log_mean=float(format(log_mean, ".15g")),
                log_standard_deviation=float(format(log_standard_deviation, ".15g")),
            )
        raise PhylogeneticsError(
            "bounded calibration priors require one supported distribution",
            code="calibration_prior_unsupported_bounded_distribution",
            details={
                "calibration_id": calibration.calibration_id,
                "distribution": calibration.distribution,
                "allowed_families": list(CALIBRATION_PRIOR_FAMILIES),
            },
        )
    if minimum_age is None:
        raise PhylogeneticsError(
            "calibration prior evaluation requires one minimum age for offset calibrations",
            code="calibration_prior_missing_minimum_age",
            details={"calibration_id": calibration.calibration_id},
        )
    if requested_distribution in {"lognormal", "offset-lognormal"}:
        return CalibrationPriorDefinition(
            calibration_id=calibration.calibration_id,
            requested_distribution=requested_distribution,
            family="offset-lognormal",
            target_kind=calibration.target_kind,
            target_label=calibration.target_label,
            descendant_taxa=sorted(calibration.taxa),
            node_id=node_id,
            node_kind=node_kind,
            minimum_age=minimum_age,
            maximum_age=maximum_age,
            translated=(requested_distribution != "offset-lognormal"),
            translation_note=(
                "lower-bound-only lognormal calibrations use one offset lognormal prior with default broad shape parameters"
                if requested_distribution != "offset-lognormal"
                else None
            ),
            offset_age=float(format(minimum_age, ".15g")),
            log_mean=1.0,
            log_standard_deviation=1.25,
        )
    if requested_distribution in {
        "uniform",
        "normal",
        "offset",
        "offset-exponential",
    }:
        exponential_mean = max(1.0, minimum_age * 0.25)
        return CalibrationPriorDefinition(
            calibration_id=calibration.calibration_id,
            requested_distribution=requested_distribution,
            family="offset-exponential",
            target_kind=calibration.target_kind,
            target_label=calibration.target_label,
            descendant_taxa=sorted(calibration.taxa),
            node_id=node_id,
            node_kind=node_kind,
            minimum_age=minimum_age,
            maximum_age=maximum_age,
            translated=(requested_distribution not in {"offset", "offset-exponential"}),
            translation_note=(
                "lower-bound-only calibrations use one offset exponential prior when no upper bound is available"
                if requested_distribution not in {"offset", "offset-exponential"}
                else None
            ),
            offset_age=float(format(minimum_age, ".15g")),
            exponential_mean=float(format(exponential_mean, ".15g")),
        )
    raise PhylogeneticsError(
        "lower-bound-only calibration priors require one supported offset family",
        code="calibration_prior_unsupported_offset_distribution",
        details={
            "calibration_id": calibration.calibration_id,
            "distribution": calibration.distribution,
            "allowed_families": list(CALIBRATION_PRIOR_FAMILIES),
        },
    )


def _require_feasible_calibration_age_constraints(
    resolved_calibrations: list[tuple[ValidatedCalibration, TreeNode]],
) -> None:
    merged_bounds_by_node_id: dict[str, tuple[float | None, float | None]] = {}
    node_by_id: dict[str, TreeNode] = {}
    calibration_ids_by_node_id: dict[str, list[str]] = {}
    for calibration, node in resolved_calibrations:
        node_id = require_present(
            node.node_id,
            owner_name="calibration prior age-constraint resolution",
            field_name="node_id",
        )
        node_by_id[node_id] = node
        calibration_ids_by_node_id.setdefault(node_id, []).append(
            calibration.calibration_id
        )
        existing_bounds = merged_bounds_by_node_id.get(node_id)
        minimum_age = calibration.minimum_age
        maximum_age = calibration.maximum_age
        if existing_bounds is None:
            merged_bounds_by_node_id[node_id] = (minimum_age, maximum_age)
            continue
        existing_minimum_age, existing_maximum_age = existing_bounds
        merged_minimum_age = (
            minimum_age
            if existing_minimum_age is None
            else (
                existing_minimum_age
                if minimum_age is None
                else max(existing_minimum_age, minimum_age)
            )
        )
        merged_maximum_age = (
            maximum_age
            if existing_maximum_age is None
            else (
                existing_maximum_age
                if maximum_age is None
                else min(existing_maximum_age, maximum_age)
            )
        )
        if (
            merged_minimum_age is not None
            and merged_maximum_age is not None
            and merged_maximum_age < (merged_minimum_age - _FIXED_CALIBRATION_TOLERANCE)
        ):
            raise PhylogeneticsError(
                "calibration prior constraints are infeasible: one node has contradictory age bounds",
                code="calibration_prior_infeasible_bounds",
                details={
                    "node_id": node.node_id,
                    "calibration_ids": calibration_ids_by_node_id[node.node_id],
                },
            )
        merged_bounds_by_node_id[node_id] = (
            merged_minimum_age,
            merged_maximum_age,
        )

    for ancestor_node_id, ancestor_bounds in merged_bounds_by_node_id.items():
        ancestor_minimum_age, ancestor_maximum_age = ancestor_bounds
        if ancestor_maximum_age is None:
            continue
        ancestor_node = node_by_id[ancestor_node_id]
        for descendant_node_id, descendant_bounds in merged_bounds_by_node_id.items():
            if ancestor_node_id == descendant_node_id:
                continue
            descendant_node = node_by_id[descendant_node_id]
            if not _is_ancestor(ancestor_node, descendant_node):
                continue
            descendant_minimum_age, _descendant_maximum_age = descendant_bounds
            if descendant_minimum_age is None:
                continue
            if ancestor_maximum_age <= (
                descendant_minimum_age + _FIXED_CALIBRATION_TOLERANCE
            ):
                raise PhylogeneticsError(
                    "calibration prior constraints are infeasible: ancestor and descendant calibration windows admit no age ordering",
                    code="calibration_prior_infeasible_chronology",
                    details={
                        "ancestor_node_id": ancestor_node_id,
                        "descendant_node_id": descendant_node_id,
                        "ancestor_calibration_ids": calibration_ids_by_node_id[
                            ancestor_node_id
                        ],
                        "descendant_calibration_ids": calibration_ids_by_node_id[
                            descendant_node_id
                        ],
                    },
                )


def _compute_node_age_by_id(
    tree: PhyloTree,
    *,
    ultrametric_tolerance: float,
) -> dict[str, float]:
    depth_by_node_id: dict[str, float] = {}
    root = tree.root
    root_node_id = require_present(
        root.node_id,
        owner_name="calibration prior node-age computation",
        field_name="root.node_id",
    )
    depth_by_node_id[root_node_id] = 0.0
    for parent, child in tree.iter_edges():
        parent_node_id = require_present(
            parent.node_id,
            owner_name="calibration prior node-age computation",
            field_name="parent.node_id",
        )
        child_node_id = require_present(
            child.node_id,
            owner_name="calibration prior node-age computation",
            field_name="child.node_id",
        )
        child_branch_length = require_present(
            child.branch_length,
            owner_name="calibration prior node-age computation",
            field_name="child.branch_length",
        )
        depth_by_node_id[child_node_id] = depth_by_node_id[parent_node_id] + float(
            child_branch_length
        )
    ultrametric_summary = summarize_ultrametric_tip_depths(
        {
            leaf.name or leaf.node_id or "": depth_by_node_id[leaf.node_id or ""]
            for leaf in tree.iter_leaves()
        },
        tolerance=ultrametric_tolerance,
    )
    if not ultrametric_summary.ultrametric:
        raise NonUltrametricTreeError(
            "calibration prior evaluation requires an ultrametric tree",
            code="calibration_prior_requires_ultrametric_tree",
            details={
                "minimum_tip_depth": ultrametric_summary.minimum_tip_depth,
                "maximum_tip_depth": ultrametric_summary.maximum_tip_depth,
                "max_tip_depth_deviation": ultrametric_summary.max_tip_depth_deviation,
                "offending_taxa": ultrametric_summary.offending_taxa,
                "tolerance": ultrametric_tolerance,
            },
        )
    root_age = ultrametric_summary.root_age
    return {
        node_id: float(format(root_age - depth, ".15g"))
        for node_id, depth in depth_by_node_id.items()
    }


def _evaluate_calibration_log_prior(
    *,
    node_age: float,
    prior_definition: CalibrationPriorDefinition,
) -> float:
    if prior_definition.family == "fixed":
        fixed_age = require_present(
            prior_definition.fixed_age,
            owner_name="calibration prior evaluation",
            field_name="fixed_age",
        )
        fixed_tolerance = require_present(
            prior_definition.fixed_tolerance,
            owner_name="calibration prior evaluation",
            field_name="fixed_tolerance",
        )
        return (
            0.0
            if math.isclose(
                node_age,
                fixed_age,
                rel_tol=0.0,
                abs_tol=fixed_tolerance,
            )
            else -math.inf
        )
    if prior_definition.family == "uniform":
        minimum_age = require_present(
            prior_definition.minimum_age,
            owner_name="calibration prior evaluation",
            field_name="minimum_age",
        )
        maximum_age = require_present(
            prior_definition.maximum_age,
            owner_name="calibration prior evaluation",
            field_name="maximum_age",
        )
        if not (minimum_age <= node_age <= maximum_age):
            return -math.inf
        return -math.log(maximum_age - minimum_age)
    if prior_definition.family == "normal":
        minimum_age = require_present(
            prior_definition.minimum_age,
            owner_name="calibration prior evaluation",
            field_name="minimum_age",
        )
        maximum_age = require_present(
            prior_definition.maximum_age,
            owner_name="calibration prior evaluation",
            field_name="maximum_age",
        )
        mean_age = require_present(
            prior_definition.mean_age,
            owner_name="calibration prior evaluation",
            field_name="mean_age",
        )
        standard_deviation = require_present(
            prior_definition.standard_deviation,
            owner_name="calibration prior evaluation",
            field_name="standard_deviation",
        )
        return _truncated_normal_log_density(
            value=node_age,
            mean=mean_age,
            standard_deviation=standard_deviation,
            lower_bound=minimum_age,
            upper_bound=maximum_age,
        )
    if prior_definition.family == "lognormal":
        minimum_age = require_present(
            prior_definition.minimum_age,
            owner_name="calibration prior evaluation",
            field_name="minimum_age",
        )
        maximum_age = require_present(
            prior_definition.maximum_age,
            owner_name="calibration prior evaluation",
            field_name="maximum_age",
        )
        log_mean = require_present(
            prior_definition.log_mean,
            owner_name="calibration prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_definition.log_standard_deviation,
            owner_name="calibration prior evaluation",
            field_name="log_standard_deviation",
        )
        return _truncated_lognormal_log_density(
            value=node_age,
            log_mean=log_mean,
            log_standard_deviation=log_standard_deviation,
            lower_bound=minimum_age,
            upper_bound=maximum_age,
        )
    if prior_definition.family == "offset-exponential":
        offset_age = require_present(
            prior_definition.offset_age,
            owner_name="calibration prior evaluation",
            field_name="offset_age",
        )
        exponential_mean = require_present(
            prior_definition.exponential_mean,
            owner_name="calibration prior evaluation",
            field_name="exponential_mean",
        )
        return _offset_exponential_log_density(
            value=node_age,
            offset_age=offset_age,
            exponential_mean=exponential_mean,
        )
    if prior_definition.family == "offset-lognormal":
        offset_age = require_present(
            prior_definition.offset_age,
            owner_name="calibration prior evaluation",
            field_name="offset_age",
        )
        log_mean = require_present(
            prior_definition.log_mean,
            owner_name="calibration prior evaluation",
            field_name="log_mean",
        )
        log_standard_deviation = require_present(
            prior_definition.log_standard_deviation,
            owner_name="calibration prior evaluation",
            field_name="log_standard_deviation",
        )
        return _offset_lognormal_log_density(
            value=node_age,
            offset_age=offset_age,
            log_mean=log_mean,
            log_standard_deviation=log_standard_deviation,
        )
    raise PhylogeneticsError(
        "calibration prior family is unsupported",
        code="calibration_prior_family_invalid",
        details={
            "family": prior_definition.family,
            "allowed_families": list(CALIBRATION_PRIOR_FAMILIES),
        },
    )


def _truncated_normal_log_density(
    *,
    value: float,
    mean: float,
    standard_deviation: float,
    lower_bound: float,
    upper_bound: float,
) -> float:
    if value < lower_bound or value > upper_bound:
        return -math.inf
    normalization_mass = _normal_cdf(
        upper_bound, mean, standard_deviation
    ) - _normal_cdf(
        lower_bound,
        mean,
        standard_deviation,
    )
    if normalization_mass <= 0.0:
        return -math.inf
    centered_value = (value - mean) / standard_deviation
    return (
        -math.log(standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - 0.5 * centered_value * centered_value
        - math.log(normalization_mass)
    )


def _truncated_lognormal_log_density(
    *,
    value: float,
    log_mean: float,
    log_standard_deviation: float,
    lower_bound: float,
    upper_bound: float,
) -> float:
    if value <= 0.0 or value < lower_bound or value > upper_bound:
        return -math.inf
    normalization_mass = _normal_cdf(
        math.log(upper_bound),
        log_mean,
        log_standard_deviation,
    ) - _normal_cdf(
        math.log(lower_bound),
        log_mean,
        log_standard_deviation,
    )
    if normalization_mass <= 0.0:
        return -math.inf
    log_value = math.log(value)
    centered_value = (log_value - log_mean) / log_standard_deviation
    return (
        -math.log(value)
        - math.log(log_standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - 0.5 * centered_value * centered_value
        - math.log(normalization_mass)
    )


def _offset_exponential_log_density(
    *,
    value: float,
    offset_age: float,
    exponential_mean: float,
) -> float:
    if value < offset_age:
        return -math.inf
    shifted_value = value - offset_age
    return -math.log(exponential_mean) - (shifted_value / exponential_mean)


def _offset_lognormal_log_density(
    *,
    value: float,
    offset_age: float,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    shifted_value = value - offset_age
    if shifted_value <= 0.0:
        return -math.inf
    centered_value = (math.log(shifted_value) - log_mean) / log_standard_deviation
    return (
        -math.log(shifted_value)
        - math.log(log_standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - 0.5 * centered_value * centered_value
    )


def _normal_cdf(value: float, mean: float, standard_deviation: float) -> float:
    return 0.5 * (
        1.0 + math.erf((value - mean) / (standard_deviation * math.sqrt(2.0)))
    )


def _validate_rooted_ultrametric_tree(
    tree: PhyloTree,
    *,
    prior_name: str,
    ultrametric_tolerance: float,
) -> None:
    if not _tree_is_rooted(tree):
        raise UnrootedTreeError(
            f"{prior_name} evaluation requires a rooted tree",
            code="calibration_prior_requires_rooted_tree",
            details={"rooted": tree.rooted},
        )
    for node in tree.iter_nodes(order="preorder"):
        if node.is_leaf():
            continue
        if len(node.children) != 2:
            raise PhylogeneticsError(
                f"{prior_name} evaluation requires a strictly bifurcating tree",
                code="calibration_prior_requires_bifurcating_tree",
                details={"node_id": node.node_id, "child_count": len(node.children)},
            )
    tip_depth_by_label: dict[str, float] = {}
    for leaf in tree.iter_leaves():
        if leaf.name is None:
            raise PhylogeneticsError(
                f"{prior_name} evaluation requires named tips",
                code="calibration_prior_requires_named_tips",
            )
        tip_depth_by_label[leaf.name] = _tip_depth(leaf)
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_depth_by_label,
        tolerance=ultrametric_tolerance,
    )
    if not ultrametric_summary.ultrametric:
        raise NonUltrametricTreeError(
            f"{prior_name} evaluation requires an ultrametric tree",
            code="calibration_prior_requires_ultrametric_tree",
            details={
                "minimum_tip_depth": ultrametric_summary.minimum_tip_depth,
                "maximum_tip_depth": ultrametric_summary.maximum_tip_depth,
                "max_tip_depth_deviation": ultrametric_summary.max_tip_depth_deviation,
                "offending_taxa": ultrametric_summary.offending_taxa,
                "tolerance": ultrametric_tolerance,
            },
        )


def _tip_depth(leaf: TreeNode) -> float:
    depth = 0.0
    node = leaf
    while node.parent is not None:
        if node.branch_length is None:
            raise InvalidBranchLengthError(
                "calibration prior evaluation requires complete branch lengths on every edge"
            )
        if node.branch_length < 0.0:
            raise InvalidBranchLengthError(
                "calibration prior evaluation requires non-negative branch lengths"
            )
        depth += float(node.branch_length)
        node = node.parent
    return depth


def _is_ancestor(ancestor_node: TreeNode, descendant_node: TreeNode) -> bool:
    node = descendant_node.parent
    while node is not None:
        if node is ancestor_node:
            return True
        node = node.parent
    return False


def _tree_is_rooted(tree: PhyloTree) -> bool:
    if tree.rooted is True:
        return True
    return len(tree.root.children) == 2
