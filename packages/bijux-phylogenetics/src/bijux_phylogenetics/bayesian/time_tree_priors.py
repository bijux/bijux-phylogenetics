from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.node_identity import build_ape_internal_node_map
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    PhylogeneticsError,
    UnrootedTreeError,
)

TIME_TREE_PRIOR_CONDITIONING_MODES = ("fixed-tip-count-and-crown-age",)
YULE_TREE_PRIOR_FAMILIES = ("crown-conditioned-yule",)
BIRTH_DEATH_TREE_PRIOR_FAMILIES = ("crown-conditioned-birth-death",)
COALESCENT_TREE_PRIOR_FAMILIES = (
    "constant-population-coalescent",
    "skyline-coalescent",
)


@dataclass(frozen=True, slots=True)
class YuleTreePriorModel:
    """One validated Yule prior parameterization for rooted crown trees."""

    family: str
    speciation_rate: float


@dataclass(frozen=True, slots=True)
class YuleTreePriorIntervalRow:
    """One deterministic crown-tree interval contribution row."""

    interval_index: int
    older_boundary_age: float
    younger_boundary_age: float
    duration: float
    lineage_count: int
    event_count: int
    interval_log_contribution: float


@dataclass(frozen=True, slots=True)
class YuleTreePriorEvaluationReport:
    """One rooted ultrametric Yule prior evaluation report."""

    family: str
    speciation_rate: float
    tree_newick: str
    tip_count: int
    internal_node_count: int
    post_root_speciation_count: int
    root_age: float
    total_branch_length: float
    ultrametric_tolerance: float
    log_prior: float
    interval_rows: list[YuleTreePriorIntervalRow]


@dataclass(frozen=True, slots=True)
class BirthDeathTreePriorModel:
    """One validated reconstructed birth-death prior parameterization."""

    family: str
    speciation_rate: float
    extinction_rate: float
    sampling_fraction: float


@dataclass(frozen=True, slots=True)
class BirthDeathTreePriorBranchingRow:
    """One post-root branching-time contribution under the birth-death prior."""

    branching_event_index: int
    branching_age: float
    log_branching_probability: float
    log_crown_normalization: float
    log_contribution: float


@dataclass(frozen=True, slots=True)
class BirthDeathTreePriorEvaluationReport:
    """One rooted ultrametric birth-death prior evaluation report."""

    family: str
    conditioning_mode: str
    speciation_rate: float
    extinction_rate: float
    sampling_fraction: float
    transformed_birth_rate: float
    transformed_death_rate: float
    tree_newick: str
    tip_count: int
    internal_node_count: int
    post_root_speciation_count: int
    root_age: float
    ultrametric_tolerance: float
    log_prior: float
    branching_rows: list[BirthDeathTreePriorBranchingRow]


@dataclass(frozen=True, slots=True)
class ConstantPopulationCoalescentPriorModel:
    """One validated constant-population coalescent prior parameterization."""

    family: str
    effective_population_size: float


@dataclass(frozen=True, slots=True)
class ConstantPopulationCoalescentIntervalRow:
    """One backward-time coalescent interval contribution row."""

    interval_index: int
    younger_boundary_age: float
    older_boundary_age: float
    duration: float
    lineage_count: int
    coalescent_event_count: int
    waiting_rate: float
    waiting_log_contribution: float
    event_log_contribution: float
    interval_log_contribution: float


@dataclass(frozen=True, slots=True)
class ConstantPopulationCoalescentPriorEvaluationReport:
    """One rooted ultrametric constant-population coalescent prior report."""

    family: str
    effective_population_size: float
    tree_newick: str
    tip_count: int
    internal_node_count: int
    root_age: float
    total_branch_length: float
    ultrametric_tolerance: float
    log_prior: float
    interval_rows: list[ConstantPopulationCoalescentIntervalRow]


@dataclass(frozen=True, slots=True)
class SkylineCoalescentEpoch:
    """One piecewise-constant population-size epoch for a skyline coalescent prior."""

    younger_boundary_age: float
    older_boundary_age: float | None
    effective_population_size: float


@dataclass(frozen=True, slots=True)
class SkylineCoalescentPriorModel:
    """One validated skyline coalescent prior parameterization."""

    family: str
    epochs: list[SkylineCoalescentEpoch]


@dataclass(frozen=True, slots=True)
class SkylineCoalescentSegmentRow:
    """One skyline-overlap segment from a coalescent interval."""

    coalescent_interval_index: int
    skyline_epoch_index: int
    segment_younger_boundary_age: float
    segment_older_boundary_age: float
    duration: float
    lineage_count: int
    coalescent_event_count: int
    effective_population_size: float
    waiting_rate: float
    waiting_log_contribution: float
    event_log_contribution: float
    segment_log_contribution: float


@dataclass(frozen=True, slots=True)
class SkylineCoalescentPriorEvaluationReport:
    """One rooted ultrametric skyline coalescent prior report."""

    family: str
    epoch_count: int
    tree_newick: str
    tip_count: int
    internal_node_count: int
    root_age: float
    total_branch_length: float
    ultrametric_tolerance: float
    log_prior: float
    segment_rows: list[SkylineCoalescentSegmentRow]


@dataclass(frozen=True, slots=True)
class _CoalescentInterval:
    """One backward-time coalescent interval before skyline segmentation."""

    interval_index: int
    younger_boundary_age: float
    older_boundary_age: float
    duration: float
    lineage_count: int
    coalescent_event_count: int


def build_crown_conditioned_yule_tree_prior(
    speciation_rate: float,
) -> YuleTreePriorModel:
    """Build one crown-conditioned pure-birth prior."""
    if not math.isfinite(speciation_rate) or speciation_rate <= 0.0:
        raise PhylogeneticsError(
            "Yule tree prior requires a strictly positive finite speciation rate",
            code="yule_tree_prior_invalid_speciation_rate",
            details={"speciation_rate": speciation_rate},
        )
    return YuleTreePriorModel(
        family="crown-conditioned-yule",
        speciation_rate=speciation_rate,
    )


def build_crown_conditioned_birth_death_tree_prior(
    *,
    speciation_rate: float,
    extinction_rate: float,
    sampling_fraction: float,
) -> BirthDeathTreePriorModel:
    """Build one reconstructed birth-death prior conditioned on crown age and tip count."""
    if not math.isfinite(speciation_rate) or speciation_rate <= 0.0:
        raise PhylogeneticsError(
            "birth-death tree prior requires a strictly positive finite speciation rate",
            code="birth_death_tree_prior_invalid_speciation_rate",
            details={"speciation_rate": speciation_rate},
        )
    if not math.isfinite(extinction_rate) or extinction_rate < 0.0:
        raise PhylogeneticsError(
            "birth-death tree prior requires a non-negative finite extinction rate",
            code="birth_death_tree_prior_invalid_extinction_rate",
            details={"extinction_rate": extinction_rate},
        )
    if not math.isfinite(sampling_fraction) or not (0.0 < sampling_fraction <= 1.0):
        raise PhylogeneticsError(
            "birth-death tree prior requires sampling_fraction in (0, 1]",
            code="birth_death_tree_prior_invalid_sampling_fraction",
            details={"sampling_fraction": sampling_fraction},
        )
    return BirthDeathTreePriorModel(
        family="crown-conditioned-birth-death",
        speciation_rate=speciation_rate,
        extinction_rate=extinction_rate,
        sampling_fraction=sampling_fraction,
    )


def build_constant_population_coalescent_tree_prior(
    *,
    effective_population_size: float,
) -> ConstantPopulationCoalescentPriorModel:
    """Build one constant-population Kingman coalescent prior."""
    if not math.isfinite(effective_population_size) or effective_population_size <= 0.0:
        raise PhylogeneticsError(
            "constant-population coalescent prior requires a strictly positive finite effective population size",
            code="constant_population_coalescent_prior_invalid_effective_population_size",
            details={"effective_population_size": effective_population_size},
        )
    return ConstantPopulationCoalescentPriorModel(
        family="constant-population-coalescent",
        effective_population_size=effective_population_size,
    )


def build_skyline_coalescent_tree_prior(
    epochs: list[SkylineCoalescentEpoch],
) -> SkylineCoalescentPriorModel:
    """Build one skyline coalescent prior with validated piecewise epochs."""
    validated_epochs = _validate_skyline_coalescent_epochs(epochs)
    return SkylineCoalescentPriorModel(
        family="skyline-coalescent",
        epochs=validated_epochs,
    )


def evaluate_yule_tree_log_prior(
    tree: PhyloTree,
    prior_model: YuleTreePriorModel,
    *,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> YuleTreePriorEvaluationReport:
    """Evaluate one rooted ultrametric crown tree under a pure-birth Yule prior."""
    _validate_time_tree_prior_tree(
        tree,
        prior_name="Yule tree prior",
        code_prefix="yule_tree_prior",
    )
    ultrametric_summary = _validated_ultrametric_summary(
        tree,
        prior_name="Yule tree prior",
        code_prefix="yule_tree_prior",
        ultrametric_tolerance=ultrametric_tolerance,
    )
    root_age = ultrametric_summary.root_age
    branch_rows = _build_yule_interval_rows(
        tree,
        speciation_rate=prior_model.speciation_rate,
        root_age=root_age,
    )
    log_prior = sum(row.interval_log_contribution for row in branch_rows)
    return YuleTreePriorEvaluationReport(
        family=prior_model.family,
        speciation_rate=prior_model.speciation_rate,
        tree_newick=tree.to_newick(),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        post_root_speciation_count=max(tree.internal_node_count - 1, 0),
        root_age=float(format(root_age, ".15g")),
        total_branch_length=float(format(tree.total_branch_length(), ".15g")),
        ultrametric_tolerance=ultrametric_tolerance,
        log_prior=float(format(log_prior, ".15g")),
        interval_rows=branch_rows,
    )


def evaluate_birth_death_tree_log_prior(
    tree: PhyloTree,
    prior_model: BirthDeathTreePriorModel,
    *,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> BirthDeathTreePriorEvaluationReport:
    """Evaluate one rooted ultrametric reconstructed tree under a birth-death prior."""
    _validate_time_tree_prior_tree(
        tree,
        prior_name="birth-death tree prior",
        code_prefix="birth_death_tree_prior",
    )
    ultrametric_summary = _validated_ultrametric_summary(
        tree,
        prior_name="birth-death tree prior",
        code_prefix="birth_death_tree_prior",
        ultrametric_tolerance=ultrametric_tolerance,
    )
    root_age = ultrametric_summary.root_age
    branching_rows = _build_birth_death_branching_rows(
        tree,
        speciation_rate=prior_model.speciation_rate,
        extinction_rate=prior_model.extinction_rate,
        sampling_fraction=prior_model.sampling_fraction,
        crown_age=root_age,
    )
    log_prior = (-math.log(tree.tip_count - 1)) + sum(
        row.log_contribution for row in branching_rows
    )
    return BirthDeathTreePriorEvaluationReport(
        family=prior_model.family,
        conditioning_mode="fixed-tip-count-and-crown-age",
        speciation_rate=prior_model.speciation_rate,
        extinction_rate=prior_model.extinction_rate,
        sampling_fraction=prior_model.sampling_fraction,
        transformed_birth_rate=float(
            format(
                prior_model.speciation_rate * prior_model.sampling_fraction,
                ".15g",
            )
        ),
        transformed_death_rate=float(
            format(
                prior_model.extinction_rate
                - prior_model.speciation_rate * (1.0 - prior_model.sampling_fraction),
                ".15g",
            )
        ),
        tree_newick=tree.to_newick(),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        post_root_speciation_count=max(tree.internal_node_count - 1, 0),
        root_age=float(format(root_age, ".15g")),
        ultrametric_tolerance=ultrametric_tolerance,
        log_prior=float(format(log_prior, ".15g")),
        branching_rows=branching_rows,
    )


def evaluate_constant_population_coalescent_tree_log_prior(
    tree: PhyloTree,
    prior_model: ConstantPopulationCoalescentPriorModel,
    *,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> ConstantPopulationCoalescentPriorEvaluationReport:
    """Evaluate one rooted ultrametric tree under a constant-population coalescent prior."""
    _validate_time_tree_prior_tree(
        tree,
        prior_name="constant-population coalescent prior",
        code_prefix="constant_population_coalescent_prior",
    )
    ultrametric_summary = _validated_ultrametric_summary(
        tree,
        prior_name="constant-population coalescent prior",
        code_prefix="constant_population_coalescent_prior",
        ultrametric_tolerance=ultrametric_tolerance,
    )
    root_age = ultrametric_summary.root_age
    interval_rows = _build_constant_population_coalescent_interval_rows(
        tree,
        effective_population_size=prior_model.effective_population_size,
        root_age=root_age,
    )
    log_prior = sum(row.interval_log_contribution for row in interval_rows)
    return ConstantPopulationCoalescentPriorEvaluationReport(
        family=prior_model.family,
        effective_population_size=prior_model.effective_population_size,
        tree_newick=tree.to_newick(),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        root_age=float(format(root_age, ".15g")),
        total_branch_length=float(format(tree.total_branch_length(), ".15g")),
        ultrametric_tolerance=ultrametric_tolerance,
        log_prior=float(format(log_prior, ".15g")),
        interval_rows=interval_rows,
    )


def evaluate_skyline_coalescent_tree_log_prior(
    tree: PhyloTree,
    prior_model: SkylineCoalescentPriorModel,
    *,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> SkylineCoalescentPriorEvaluationReport:
    """Evaluate one rooted ultrametric tree under a skyline coalescent prior."""
    _validate_time_tree_prior_tree(
        tree,
        prior_name="skyline coalescent prior",
        code_prefix="skyline_coalescent_prior",
    )
    ultrametric_summary = _validated_ultrametric_summary(
        tree,
        prior_name="skyline coalescent prior",
        code_prefix="skyline_coalescent_prior",
        ultrametric_tolerance=ultrametric_tolerance,
    )
    root_age = ultrametric_summary.root_age
    segment_rows = _build_skyline_coalescent_segment_rows(
        tree,
        epochs=prior_model.epochs,
        root_age=root_age,
    )
    log_prior = sum(row.segment_log_contribution for row in segment_rows)
    return SkylineCoalescentPriorEvaluationReport(
        family=prior_model.family,
        epoch_count=len(prior_model.epochs),
        tree_newick=tree.to_newick(),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        root_age=float(format(root_age, ".15g")),
        total_branch_length=float(format(tree.total_branch_length(), ".15g")),
        ultrametric_tolerance=ultrametric_tolerance,
        log_prior=float(format(log_prior, ".15g")),
        segment_rows=segment_rows,
    )


def _validate_time_tree_prior_tree(
    tree: PhyloTree,
    *,
    prior_name: str,
    code_prefix: str,
) -> None:
    if not _tree_is_rooted(tree):
        raise UnrootedTreeError(
            f"{prior_name} requires a rooted tree",
            code=f"{code_prefix}_requires_rooted_tree",
        )
    if tree.tip_count < 2:
        raise PhylogeneticsError(
            f"{prior_name} requires at least two tips",
            code=f"{code_prefix}_requires_two_or_more_tips",
        )
    for _parent, child in tree.iter_edges():
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                f"{prior_name} requires complete branch lengths"
            )
        if child.branch_length < 0.0:
            raise InvalidBranchLengthError(
                f"{prior_name} requires non-negative branch lengths"
            )
    if not _tree_is_strictly_bifurcating(tree):
        raise PhylogeneticsError(
            f"{prior_name} requires a strictly bifurcating tree",
            code=f"{code_prefix}_requires_strictly_bifurcating_tree",
        )


def _tree_is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(
        len(node.children) == 2 for node in tree.iter_nodes() if not node.is_leaf()
    )


def _tree_is_rooted(tree: PhyloTree) -> bool:
    if tree.rooted is True:
        return True
    return len(tree.root.children) == 2


def _validated_ultrametric_summary(
    tree: PhyloTree,
    *,
    prior_name: str,
    code_prefix: str,
    ultrametric_tolerance: float,
):
    tip_depth_by_label = _tip_depth_by_label(tree)
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_depth_by_label,
        tolerance=ultrametric_tolerance,
    )
    if not ultrametric_summary.ultrametric:
        raise NonUltrametricTreeError(
            f"{prior_name} requires an ultrametric tree",
            code=f"{code_prefix}_requires_ultrametric_tree",
            details={
                "minimum_tip_depth": ultrametric_summary.minimum_tip_depth,
                "maximum_tip_depth": ultrametric_summary.maximum_tip_depth,
                "max_tip_depth_deviation": ultrametric_summary.max_tip_depth_deviation,
                "offending_taxa": list(ultrametric_summary.offending_taxa),
                "tolerance": ultrametric_summary.tolerance,
            },
        )
    return ultrametric_summary


def _validate_skyline_coalescent_epochs(
    epochs: list[SkylineCoalescentEpoch],
) -> list[SkylineCoalescentEpoch]:
    if not epochs:
        raise PhylogeneticsError(
            "skyline coalescent prior requires at least one epoch",
            code="skyline_coalescent_prior_requires_one_or_more_epochs",
        )
    validated_epochs: list[SkylineCoalescentEpoch] = []
    expected_younger_boundary = 0.0
    for index, epoch in enumerate(epochs, start=1):
        if (
            not math.isfinite(epoch.younger_boundary_age)
            or epoch.younger_boundary_age < 0.0
        ):
            raise PhylogeneticsError(
                "skyline coalescent prior requires non-negative finite younger epoch boundaries",
                code="skyline_coalescent_prior_invalid_younger_boundary_age",
                details={
                    "epoch_index": index,
                    "younger_boundary_age": epoch.younger_boundary_age,
                },
            )
        if not math.isclose(
            epoch.younger_boundary_age,
            expected_younger_boundary,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise PhylogeneticsError(
                "skyline coalescent prior requires contiguous epochs from age zero",
                code="skyline_coalescent_prior_requires_contiguous_epochs",
                details={
                    "epoch_index": index,
                    "expected_younger_boundary_age": expected_younger_boundary,
                    "younger_boundary_age": epoch.younger_boundary_age,
                },
            )
        if (
            not math.isfinite(epoch.effective_population_size)
            or epoch.effective_population_size <= 0.0
        ):
            raise PhylogeneticsError(
                "skyline coalescent prior requires strictly positive finite epoch effective population sizes",
                code="skyline_coalescent_prior_invalid_effective_population_size",
                details={
                    "epoch_index": index,
                    "effective_population_size": epoch.effective_population_size,
                },
            )
        if epoch.older_boundary_age is None:
            if index != len(epochs):
                raise PhylogeneticsError(
                    "skyline coalescent prior requires only the final epoch to be open ended",
                    code="skyline_coalescent_prior_nonterminal_open_epoch",
                    details={"epoch_index": index},
                )
            older_boundary_age = None
        else:
            if (
                not math.isfinite(epoch.older_boundary_age)
                or epoch.older_boundary_age <= epoch.younger_boundary_age
            ):
                raise PhylogeneticsError(
                    "skyline coalescent prior requires each finite epoch to end after it starts",
                    code="skyline_coalescent_prior_invalid_older_boundary_age",
                    details={
                        "epoch_index": index,
                        "younger_boundary_age": epoch.younger_boundary_age,
                        "older_boundary_age": epoch.older_boundary_age,
                    },
                )
            older_boundary_age = float(format(epoch.older_boundary_age, ".15g"))
            expected_younger_boundary = older_boundary_age
        validated_epochs.append(
            SkylineCoalescentEpoch(
                younger_boundary_age=float(format(epoch.younger_boundary_age, ".15g")),
                older_boundary_age=older_boundary_age,
                effective_population_size=float(
                    format(epoch.effective_population_size, ".15g")
                ),
            )
        )
    if validated_epochs[-1].older_boundary_age is not None:
        raise PhylogeneticsError(
            "skyline coalescent prior requires the final epoch to be open ended",
            code="skyline_coalescent_prior_requires_open_final_epoch",
        )
    return validated_epochs


def _tip_depth_by_label(tree: PhyloTree) -> dict[str, float]:
    depths = _node_depth_lookup(tree)
    return {
        node.name: depths[node.node_id or ""]
        for node in tree.iter_leaves()
        if node.name is not None
    }


def _node_depth_lookup(tree: PhyloTree) -> dict[str, float]:
    root_id = tree.root.node_id or ""
    depths: dict[str, float] = {root_id: 0.0}

    def visit(node: TreeNode) -> None:
        node_id = node.node_id or ""
        base_depth = depths[node_id]
        for child in node.children:
            depths[child.node_id or ""] = base_depth + float(child.branch_length or 0.0)
            visit(child)

    visit(tree.root)
    return depths


def _branching_ages_excluding_root(tree: PhyloTree, *, root_age: float) -> list[float]:
    internal_nodes = build_ape_internal_node_map(tree)
    node_depths = _node_depth_lookup(tree)
    return [
        float(format(root_age - node_depths[node.node_id or ""], ".15g"))
        for node in internal_nodes.values()
        if node is not tree.root
    ]


def _branching_ages_including_root(tree: PhyloTree, *, root_age: float) -> list[float]:
    internal_nodes = build_ape_internal_node_map(tree)
    node_depths = _node_depth_lookup(tree)
    return [
        float(format(root_age - node_depths[node.node_id or ""], ".15g"))
        for node in internal_nodes.values()
    ]


def _build_coalescent_intervals(
    tree: PhyloTree,
    *,
    root_age: float,
) -> list[_CoalescentInterval]:
    branching_ages = _branching_ages_including_root(tree, root_age=root_age)
    events_by_age = Counter(branching_ages)
    older_boundaries = sorted(events_by_age)
    younger_boundary = 0.0
    lineage_count = tree.tip_count
    intervals: list[_CoalescentInterval] = []
    for interval_index, older_boundary in enumerate(older_boundaries, start=1):
        duration = older_boundary - younger_boundary
        coalescent_event_count = events_by_age[older_boundary]
        intervals.append(
            _CoalescentInterval(
                interval_index=interval_index,
                younger_boundary_age=float(format(younger_boundary, ".15g")),
                older_boundary_age=float(format(older_boundary, ".15g")),
                duration=float(format(duration, ".15g")),
                lineage_count=lineage_count,
                coalescent_event_count=coalescent_event_count,
            )
        )
        lineage_count -= coalescent_event_count
        younger_boundary = older_boundary
    return intervals


def _build_yule_interval_rows(
    tree: PhyloTree,
    *,
    speciation_rate: float,
    root_age: float,
) -> list[YuleTreePriorIntervalRow]:
    branching_ages = _branching_ages_excluding_root(tree, root_age=root_age)
    events_by_age = Counter(branching_ages)
    younger_boundaries = sorted(events_by_age, reverse=True) + [0.0]
    older_boundary = float(format(root_age, ".15g"))
    lineage_count = 2
    rows: list[YuleTreePriorIntervalRow] = []
    for interval_index, younger_boundary in enumerate(younger_boundaries, start=1):
        duration = older_boundary - younger_boundary
        event_count = events_by_age.get(younger_boundary, 0)
        interval_log_contribution = (event_count * math.log(speciation_rate)) - (
            lineage_count * speciation_rate * duration
        )
        rows.append(
            YuleTreePriorIntervalRow(
                interval_index=interval_index,
                older_boundary_age=float(format(older_boundary, ".15g")),
                younger_boundary_age=float(format(younger_boundary, ".15g")),
                duration=float(format(duration, ".15g")),
                lineage_count=lineage_count,
                event_count=event_count,
                interval_log_contribution=float(
                    format(interval_log_contribution, ".15g")
                ),
            )
        )
        lineage_count += event_count
        older_boundary = younger_boundary
    return rows


def _build_birth_death_branching_rows(
    tree: PhyloTree,
    *,
    speciation_rate: float,
    extinction_rate: float,
    sampling_fraction: float,
    crown_age: float,
) -> list[BirthDeathTreePriorBranchingRow]:
    crown_normalization = _birth_death_q(
        crown_age,
        speciation_rate=speciation_rate,
        extinction_rate=extinction_rate,
        sampling_fraction=sampling_fraction,
    )
    if crown_normalization <= 0.0:
        raise PhylogeneticsError(
            "birth-death tree prior encountered a non-positive crown normalization term",
            code="birth_death_tree_prior_invalid_crown_normalization",
            details={
                "crown_age": crown_age,
                "speciation_rate": speciation_rate,
                "extinction_rate": extinction_rate,
                "sampling_fraction": sampling_fraction,
                "crown_normalization": crown_normalization,
            },
        )
    branching_ages = sorted(
        _branching_ages_excluding_root(tree, root_age=crown_age),
        reverse=True,
    )
    rows: list[BirthDeathTreePriorBranchingRow] = []
    log_crown_normalization = math.log(crown_normalization)
    for index, branching_age in enumerate(branching_ages, start=1):
        branching_probability = _birth_death_p1(
            branching_age,
            speciation_rate=speciation_rate,
            extinction_rate=extinction_rate,
            sampling_fraction=sampling_fraction,
        )
        if branching_probability <= 0.0:
            raise PhylogeneticsError(
                "birth-death tree prior encountered a non-positive branching probability term",
                code="birth_death_tree_prior_invalid_branching_probability",
                details={
                    "branching_age": branching_age,
                    "speciation_rate": speciation_rate,
                    "extinction_rate": extinction_rate,
                    "sampling_fraction": sampling_fraction,
                    "branching_probability": branching_probability,
                },
            )
        log_branching_probability = math.log(branching_probability)
        rows.append(
            BirthDeathTreePriorBranchingRow(
                branching_event_index=index,
                branching_age=float(format(branching_age, ".15g")),
                log_branching_probability=float(
                    format(log_branching_probability, ".15g")
                ),
                log_crown_normalization=float(format(log_crown_normalization, ".15g")),
                log_contribution=float(
                    format(log_branching_probability - log_crown_normalization, ".15g")
                ),
            )
        )
    return rows


def _build_constant_population_coalescent_interval_rows(
    tree: PhyloTree,
    *,
    effective_population_size: float,
    root_age: float,
) -> list[ConstantPopulationCoalescentIntervalRow]:
    coalescent_intervals = _build_coalescent_intervals(tree, root_age=root_age)
    rows: list[ConstantPopulationCoalescentIntervalRow] = []
    for interval in coalescent_intervals:
        waiting_rate = (
            math.comb(interval.lineage_count, 2) / effective_population_size
            if interval.lineage_count >= 2
            else 0.0
        )
        waiting_log_contribution = -(waiting_rate * interval.duration)
        event_log_contribution = 0.0
        remaining_lineage_count = interval.lineage_count
        for _ in range(interval.coalescent_event_count):
            if remaining_lineage_count < 2:
                raise PhylogeneticsError(
                    "constant-population coalescent prior encountered an invalid lineage count",
                    code="constant_population_coalescent_prior_invalid_lineage_count",
                    details={
                        "lineage_count": remaining_lineage_count,
                        "coalescent_event_count": interval.coalescent_event_count,
                        "older_boundary_age": interval.older_boundary_age,
                    },
                )
            event_log_contribution += math.log(
                math.comb(remaining_lineage_count, 2) / effective_population_size
            )
            remaining_lineage_count -= 1
        interval_log_contribution = waiting_log_contribution + event_log_contribution
        rows.append(
            ConstantPopulationCoalescentIntervalRow(
                interval_index=interval.interval_index,
                younger_boundary_age=interval.younger_boundary_age,
                older_boundary_age=interval.older_boundary_age,
                duration=interval.duration,
                lineage_count=interval.lineage_count,
                coalescent_event_count=interval.coalescent_event_count,
                waiting_rate=float(format(waiting_rate, ".15g")),
                waiting_log_contribution=float(
                    format(waiting_log_contribution, ".15g")
                ),
                event_log_contribution=float(format(event_log_contribution, ".15g")),
                interval_log_contribution=float(
                    format(interval_log_contribution, ".15g")
                ),
            )
        )
    return rows


def _build_skyline_coalescent_segment_rows(
    tree: PhyloTree,
    *,
    epochs: list[SkylineCoalescentEpoch],
    root_age: float,
) -> list[SkylineCoalescentSegmentRow]:
    coalescent_intervals = _build_coalescent_intervals(tree, root_age=root_age)
    rows: list[SkylineCoalescentSegmentRow] = []
    for interval in coalescent_intervals:
        segment_start = interval.younger_boundary_age
        while segment_start < interval.older_boundary_age - 1e-15:
            skyline_epoch_index, epoch = _locate_skyline_epoch(
                epochs,
                age=segment_start,
            )
            epoch_end = (
                interval.older_boundary_age
                if epoch.older_boundary_age is None
                else min(interval.older_boundary_age, epoch.older_boundary_age)
            )
            duration = epoch_end - segment_start
            waiting_rate = (
                math.comb(interval.lineage_count, 2) / epoch.effective_population_size
            )
            waiting_log_contribution = -(waiting_rate * duration)
            ends_with_event = math.isclose(
                epoch_end,
                interval.older_boundary_age,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            event_log_contribution = (
                math.log(waiting_rate) * interval.coalescent_event_count
                if ends_with_event
                else 0.0
            )
            segment_log_contribution = waiting_log_contribution + event_log_contribution
            rows.append(
                SkylineCoalescentSegmentRow(
                    coalescent_interval_index=interval.interval_index,
                    skyline_epoch_index=skyline_epoch_index,
                    segment_younger_boundary_age=float(format(segment_start, ".15g")),
                    segment_older_boundary_age=float(format(epoch_end, ".15g")),
                    duration=float(format(duration, ".15g")),
                    lineage_count=interval.lineage_count,
                    coalescent_event_count=(
                        interval.coalescent_event_count if ends_with_event else 0
                    ),
                    effective_population_size=epoch.effective_population_size,
                    waiting_rate=float(format(waiting_rate, ".15g")),
                    waiting_log_contribution=float(
                        format(waiting_log_contribution, ".15g")
                    ),
                    event_log_contribution=float(
                        format(event_log_contribution, ".15g")
                    ),
                    segment_log_contribution=float(
                        format(segment_log_contribution, ".15g")
                    ),
                )
            )
            segment_start = epoch_end
    return rows


def _locate_skyline_epoch(
    epochs: list[SkylineCoalescentEpoch],
    *,
    age: float,
) -> tuple[int, SkylineCoalescentEpoch]:
    for index, epoch in enumerate(epochs, start=1):
        older_boundary_age = epoch.older_boundary_age
        if older_boundary_age is None:
            if age >= epoch.younger_boundary_age:
                return index, epoch
        elif epoch.younger_boundary_age <= age < older_boundary_age:
            return index, epoch
    raise PhylogeneticsError(
        "skyline coalescent prior does not cover one or more coalescent ages",
        code="skyline_coalescent_prior_age_outside_epoch_coverage",
        details={"age": age},
    )


def _birth_death_p1(
    time_before_present: float,
    *,
    speciation_rate: float,
    extinction_rate: float,
    sampling_fraction: float,
) -> float:
    if math.isclose(speciation_rate, extinction_rate, abs_tol=1e-15):
        return (
            sampling_fraction
            / (1.0 + (sampling_fraction * speciation_rate * time_before_present)) ** 2
        )
    denominator = _birth_death_denominator(
        time_before_present,
        speciation_rate=speciation_rate,
        extinction_rate=extinction_rate,
        sampling_fraction=sampling_fraction,
    )
    rate_gap = speciation_rate - extinction_rate
    return (
        sampling_fraction
        * (rate_gap**2)
        * math.exp(-(rate_gap * time_before_present))
        / (denominator**2)
    )


def _birth_death_q(
    time_before_present: float,
    *,
    speciation_rate: float,
    extinction_rate: float,
    sampling_fraction: float,
) -> float:
    if math.isclose(speciation_rate, extinction_rate, abs_tol=1e-15):
        return (sampling_fraction * time_before_present) / (
            1.0 + (sampling_fraction * speciation_rate * time_before_present)
        )
    denominator = _birth_death_denominator(
        time_before_present,
        speciation_rate=speciation_rate,
        extinction_rate=extinction_rate,
        sampling_fraction=sampling_fraction,
    )
    rate_gap = speciation_rate - extinction_rate
    return (
        sampling_fraction * (1.0 - math.exp(-(rate_gap * time_before_present)))
    ) / denominator


def _birth_death_denominator(
    time_before_present: float,
    *,
    speciation_rate: float,
    extinction_rate: float,
    sampling_fraction: float,
) -> float:
    return speciation_rate * sampling_fraction + (
        (speciation_rate * (1.0 - sampling_fraction)) - extinction_rate
    ) * math.exp(-(speciation_rate - extinction_rate) * time_before_present)
