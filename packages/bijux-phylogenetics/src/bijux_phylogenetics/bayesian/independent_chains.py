from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
import math

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    BayesianLikelihoodUpdate,
    BayesianPriorUpdate,
    BayesianStateProposal,
    MetropolisHastingsRunReport,
    MetropolisHastingsStepRow,
    run_metropolis_hastings_sampler,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    serialize_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainDefinition:
    """One named seed assignment for one independent Metropolis-Hastings chain."""

    chain_name: str
    seed: int


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainReport:
    """One named independent Metropolis-Hastings chain result."""

    chain_name: str
    chain_report: MetropolisHastingsRunReport
    trace_fingerprint: str
    final_state_fingerprint: str

    @property
    def seed(self) -> int:
        return self.chain_report.seed

    @property
    def acceptance_rate(self) -> float:
        return self.chain_report.acceptance_rate

    @property
    def final_state(self) -> BayesianPhylogeneticState:
        return self.chain_report.final_state

    @property
    def step_rows(self) -> list[MetropolisHastingsStepRow]:
        return self.chain_report.step_rows

    @property
    def sampled_states(self) -> list[BayesianPhylogeneticState]:
        return self.chain_report.sampled_states


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainComparisonRow:
    """One pairwise comparison across two independent Metropolis-Hastings chains."""

    left_chain_name: str
    right_chain_name: str
    identical_trace: bool
    identical_final_state: bool
    acceptance_rate_delta: float
    posterior_log_score_delta: float


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsDiagnosticsReport:
    """One merged diagnostic summary over independent Metropolis-Hastings chains."""

    chain_count: int
    mean_acceptance_rate: float
    minimum_acceptance_rate: float
    maximum_acceptance_rate: float
    comparison_rows: list[IndependentMetropolisHastingsChainComparisonRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsRunReport:
    """One completed collection of independent Metropolis-Hastings chains."""

    chain_reports: list[IndependentMetropolisHastingsChainReport]
    diagnostics: IndependentMetropolisHastingsDiagnosticsReport


def build_independent_metropolis_hastings_chain_definition(
    *,
    chain_name: str,
    seed: int,
) -> IndependentMetropolisHastingsChainDefinition:
    """Build one validated independent chain definition."""
    return IndependentMetropolisHastingsChainDefinition(
        chain_name=_validate_nonblank_name(
            value=chain_name,
            field_name="chain_name",
            owner_name="independent metropolis-hastings chain definition",
        ),
        seed=_validate_integer_seed(
            seed=seed,
            owner_name="independent metropolis-hastings chain definition",
        ),
    )


def build_independent_metropolis_hastings_chain_report(
    *,
    chain_name: str,
    chain_report: MetropolisHastingsRunReport,
) -> IndependentMetropolisHastingsChainReport:
    """Build one validated named independent chain report."""
    validated_chain_name = _validate_nonblank_name(
        value=chain_name,
        field_name="chain_name",
        owner_name="independent metropolis-hastings chain report",
    )
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings chain report requires one MetropolisHastingsRunReport",
            code="independent_metropolis_hastings_chain_report_type_invalid",
        )
    return IndependentMetropolisHastingsChainReport(
        chain_name=validated_chain_name,
        chain_report=chain_report,
        trace_fingerprint=_fingerprint_chain_trace(chain_report),
        final_state_fingerprint=_fingerprint_state(chain_report.final_state),
    )


def build_independent_metropolis_hastings_chain_comparison_row(
    *,
    left_chain_report: IndependentMetropolisHastingsChainReport,
    right_chain_report: IndependentMetropolisHastingsChainReport,
) -> IndependentMetropolisHastingsChainComparisonRow:
    """Build one validated pairwise comparison across two independent chains."""
    if not isinstance(left_chain_report, IndependentMetropolisHastingsChainReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings chain comparison requires one left IndependentMetropolisHastingsChainReport",
            code="independent_metropolis_hastings_left_chain_report_type_invalid",
        )
    if not isinstance(right_chain_report, IndependentMetropolisHastingsChainReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings chain comparison requires one right IndependentMetropolisHastingsChainReport",
            code="independent_metropolis_hastings_right_chain_report_type_invalid",
        )
    if left_chain_report.chain_name == right_chain_report.chain_name:
        raise PhylogeneticsError(
            "independent metropolis-hastings chain comparison requires distinct chain names",
            code="independent_metropolis_hastings_chain_name_duplicate",
            details={"chain_name": left_chain_report.chain_name},
        )
    return IndependentMetropolisHastingsChainComparisonRow(
        left_chain_name=left_chain_report.chain_name,
        right_chain_name=right_chain_report.chain_name,
        identical_trace=(
            left_chain_report.trace_fingerprint == right_chain_report.trace_fingerprint
        ),
        identical_final_state=(
            left_chain_report.final_state_fingerprint
            == right_chain_report.final_state_fingerprint
        ),
        acceptance_rate_delta=round(
            right_chain_report.acceptance_rate - left_chain_report.acceptance_rate,
            15,
        ),
        posterior_log_score_delta=round(
            right_chain_report.final_state.posterior_log_score
            - left_chain_report.final_state.posterior_log_score,
            15,
        ),
    )


def build_independent_metropolis_hastings_diagnostics_report(
    *,
    chain_reports: Sequence[IndependentMetropolisHastingsChainReport],
) -> IndependentMetropolisHastingsDiagnosticsReport:
    """Build merged diagnostics across independent Metropolis-Hastings chains."""
    validated_chain_reports = _validate_independent_chain_reports(chain_reports)
    comparison_rows = [
        build_independent_metropolis_hastings_chain_comparison_row(
            left_chain_report=left_chain_report,
            right_chain_report=right_chain_report,
        )
        for left_chain_report, right_chain_report in itertools.combinations(
            validated_chain_reports,
            2,
        )
    ]
    duplicate_trace_rows = [
        comparison_row
        for comparison_row in comparison_rows
        if comparison_row.identical_trace
    ]
    if duplicate_trace_rows:
        duplicated_chain_names = sorted(
            {comparison_row.left_chain_name for comparison_row in duplicate_trace_rows}
            | {
                comparison_row.right_chain_name
                for comparison_row in duplicate_trace_rows
            }
        )
        raise PhylogeneticsError(
            "independent metropolis-hastings diagnostics rejected duplicated chain traces under multiple chain names",
            code="independent_metropolis_hastings_duplicate_trace_detected",
            details={"duplicated_chain_names": duplicated_chain_names},
        )
    acceptance_rates = [
        chain_report.acceptance_rate for chain_report in validated_chain_reports
    ]
    mean_acceptance_rate = math.fsum(acceptance_rates) / len(acceptance_rates)
    warnings: list[str] = []
    if max(acceptance_rates) - min(acceptance_rates) > 0.25:
        warnings.append(
            "independent metropolis-hastings chains show materially different acceptance rates"
        )
    if any(comparison_row.identical_final_state for comparison_row in comparison_rows):
        warnings.append(
            "one or more independent metropolis-hastings chains ended in the same final state"
        )
    return IndependentMetropolisHastingsDiagnosticsReport(
        chain_count=len(validated_chain_reports),
        mean_acceptance_rate=mean_acceptance_rate,
        minimum_acceptance_rate=min(acceptance_rates),
        maximum_acceptance_rate=max(acceptance_rates),
        comparison_rows=comparison_rows,
        warnings=warnings,
    )


def run_independent_metropolis_hastings_chains(
    *,
    initial_state: BayesianPhylogeneticState,
    propose_state: BayesianStateProposal,
    update_prior_components: BayesianPriorUpdate,
    update_log_likelihood: BayesianLikelihoodUpdate,
    chain_definitions: Sequence[IndependentMetropolisHastingsChainDefinition],
    iteration_count: int,
    sample_every: int = 1,
) -> IndependentMetropolisHastingsRunReport:
    """Run multiple independent Metropolis-Hastings chains with merged diagnostics."""
    validated_chain_definitions = _validate_chain_definitions(chain_definitions)
    chain_reports = [
        build_independent_metropolis_hastings_chain_report(
            chain_name=chain_definition.chain_name,
            chain_report=run_metropolis_hastings_sampler(
                initial_state=initial_state,
                propose_state=propose_state,
                update_prior_components=update_prior_components,
                update_log_likelihood=update_log_likelihood,
                iteration_count=iteration_count,
                sample_every=sample_every,
                seed=chain_definition.seed,
            ),
        )
        for chain_definition in validated_chain_definitions
    ]
    return IndependentMetropolisHastingsRunReport(
        chain_reports=chain_reports,
        diagnostics=build_independent_metropolis_hastings_diagnostics_report(
            chain_reports=chain_reports
        ),
    )


def _validate_chain_definitions(
    chain_definitions: Sequence[IndependentMetropolisHastingsChainDefinition],
) -> list[IndependentMetropolisHastingsChainDefinition]:
    validated_chain_definitions = list(chain_definitions)
    if len(validated_chain_definitions) < 2:
        raise PhylogeneticsError(
            "independent metropolis-hastings runner requires at least two chain definitions",
            code="independent_metropolis_hastings_chain_count_too_small",
        )
    if any(
        not isinstance(
            chain_definition,
            IndependentMetropolisHastingsChainDefinition,
        )
        for chain_definition in validated_chain_definitions
    ):
        raise PhylogeneticsError(
            "independent metropolis-hastings runner requires every chain definition to be one IndependentMetropolisHastingsChainDefinition",
            code="independent_metropolis_hastings_chain_definition_type_invalid",
        )
    chain_names = [
        chain_definition.chain_name for chain_definition in validated_chain_definitions
    ]
    if len(chain_names) != len(set(chain_names)):
        raise PhylogeneticsError(
            "independent metropolis-hastings runner requires unique chain names",
            code="independent_metropolis_hastings_chain_name_duplicate",
            details={"chain_names": chain_names},
        )
    seeds = [chain_definition.seed for chain_definition in validated_chain_definitions]
    if len(seeds) != len(set(seeds)):
        raise PhylogeneticsError(
            "independent metropolis-hastings runner requires distinct seeds across chain definitions",
            code="independent_metropolis_hastings_seed_duplicate",
            details={"seeds": seeds},
        )
    return validated_chain_definitions


def _validate_independent_chain_reports(
    chain_reports: Sequence[IndependentMetropolisHastingsChainReport],
) -> list[IndependentMetropolisHastingsChainReport]:
    validated_chain_reports = list(chain_reports)
    if len(validated_chain_reports) < 2:
        raise PhylogeneticsError(
            "independent metropolis-hastings diagnostics require at least two chain reports",
            code="independent_metropolis_hastings_report_count_too_small",
        )
    if any(
        not isinstance(chain_report, IndependentMetropolisHastingsChainReport)
        for chain_report in validated_chain_reports
    ):
        raise PhylogeneticsError(
            "independent metropolis-hastings diagnostics require every chain report to be one IndependentMetropolisHastingsChainReport",
            code="independent_metropolis_hastings_chain_report_type_invalid",
        )
    chain_names = [chain_report.chain_name for chain_report in validated_chain_reports]
    if len(chain_names) != len(set(chain_names)):
        raise PhylogeneticsError(
            "independent metropolis-hastings diagnostics require unique chain names",
            code="independent_metropolis_hastings_chain_name_duplicate",
            details={"chain_names": chain_names},
        )
    seeds = [chain_report.seed for chain_report in validated_chain_reports]
    if len(seeds) != len(set(seeds)):
        raise PhylogeneticsError(
            "independent metropolis-hastings diagnostics require distinct seeds across chain reports",
            code="independent_metropolis_hastings_seed_duplicate",
            details={"seeds": seeds},
        )
    return validated_chain_reports


def _fingerprint_chain_trace(chain_report: MetropolisHastingsRunReport) -> str:
    return _fingerprint_payload(
        {
            "initial_state": serialize_bayesian_phylogenetic_state(
                chain_report.initial_state
            ),
            "final_state": serialize_bayesian_phylogenetic_state(
                chain_report.final_state
            ),
            "sampled_states": [
                serialize_bayesian_phylogenetic_state(sampled_state)
                for sampled_state in chain_report.sampled_states
            ],
            "step_rows": [
                _serialize_step_row(step_row) for step_row in chain_report.step_rows
            ],
        }
    )


def _fingerprint_state(state: BayesianPhylogeneticState) -> str:
    return _fingerprint_payload(serialize_bayesian_phylogenetic_state(state))


def _fingerprint_payload(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()


def _serialize_step_row(step_row: MetropolisHastingsStepRow) -> dict[str, object]:
    return {
        **asdict(step_row),
        "proposal_changed_fields": list(step_row.proposal_changed_fields),
    }


def _validate_nonblank_name(
    *,
    value: str,
    field_name: str,
    owner_name: str,
) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one nonblank name",
            code="independent_metropolis_hastings_name_invalid",
            details={"field_name": field_name},
        )
    return value.strip()


def _validate_integer_seed(*, seed: int, owner_name: str) -> int:
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise PhylogeneticsError(
            f"{owner_name} requires 'seed' to be one integer",
            code="independent_metropolis_hastings_seed_type_invalid",
        )
    return seed
