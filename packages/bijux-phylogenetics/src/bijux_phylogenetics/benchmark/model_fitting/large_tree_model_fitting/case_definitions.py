from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousModeSearchControls,
)

from .contracts import LargeTreeModelFittingThreshold


@dataclass(frozen=True, slots=True)
class ContinuousCaseDefinition:
    case_id: str
    tier: str
    fit_surface: str
    fit_mode: str
    taxon_count: int
    root_state: float
    sigma: float
    lambda_value: float | None
    seed: int
    search_controls: ContinuousModeSearchControls | None
    timeout_seconds: float | None
    threshold: LargeTreeModelFittingThreshold
    geiger_match_tolerance: float


@dataclass(frozen=True, slots=True)
class DiscreteCaseDefinition:
    case_id: str
    tier: str
    fit_surface: str
    taxon_count: int
    states: tuple[str, ...]
    transition_rate: float
    root_state: str
    seed: int
    threshold: LargeTreeModelFittingThreshold
    geiger_match_tolerance: float


SMALL_TIER_CASES: tuple[ContinuousCaseDefinition | DiscreteCaseDefinition, ...] = (
    ContinuousCaseDefinition(
        case_id="fitcontinuous-pagel-lambda-100-taxa",
        tier="small",
        fit_surface="fitcontinuous-pagel-lambda",
        fit_mode="pagel-lambda",
        taxon_count=100,
        root_state=1.5,
        sigma=0.8,
        lambda_value=0.65,
        seed=100,
        search_controls=ContinuousModeSearchControls(
            coarse_grid_point_count=5,
            fine_grid_point_count=5,
        ),
        timeout_seconds=None,
        threshold=LargeTreeModelFittingThreshold(
            max_runtime_seconds=60.0,
            max_peak_memory_bytes=512 * 1024 * 1024,
            max_optimizer_step_count=20,
        ),
        geiger_match_tolerance=0.4,
    ),
    DiscreteCaseDefinition(
        case_id="fitdiscrete-er-binary-100-taxa",
        tier="small",
        fit_surface="fitdiscrete-er",
        taxon_count=100,
        states=("A", "B"),
        transition_rate=0.35,
        root_state="A",
        seed=100,
        threshold=LargeTreeModelFittingThreshold(
            max_runtime_seconds=60.0,
            max_peak_memory_bytes=512 * 1024 * 1024,
            max_optimizer_step_count=500,
        ),
        geiger_match_tolerance=0.5,
    ),
)

HEAVY_TIER_ONLY_CASES: tuple[ContinuousCaseDefinition, ...] = (
    ContinuousCaseDefinition(
        case_id="fitcontinuous-brownian-512-taxa",
        tier="heavy",
        fit_surface="fitcontinuous-brownian",
        fit_mode="brownian",
        taxon_count=512,
        root_state=1.5,
        sigma=0.8,
        lambda_value=None,
        seed=512,
        search_controls=None,
        timeout_seconds=60.0,
        threshold=LargeTreeModelFittingThreshold(
            max_runtime_seconds=60.0,
            max_peak_memory_bytes=2 * 1024 * 1024 * 1024,
            max_optimizer_step_count=0,
        ),
        geiger_match_tolerance=0.2,
    ),
)


def case_definitions_for_tier(
    tier: str,
) -> tuple[ContinuousCaseDefinition | DiscreteCaseDefinition, ...]:
    if tier == "small":
        return SMALL_TIER_CASES
    if tier == "heavy":
        return (*SMALL_TIER_CASES, *HEAVY_TIER_ONLY_CASES)
    raise ValueError(f"unsupported tier '{tier}'; expected one of: heavy, small")
