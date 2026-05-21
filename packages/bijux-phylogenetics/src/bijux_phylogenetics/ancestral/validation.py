from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError
from bijux_phylogenetics.simulation import (
    simulate_brownian_traits,
    simulate_discrete_traits,
    simulate_ou_traits,
)


@dataclass(slots=True)
class ContinuousValidationNodeObservation:
    node: str
    descendant_taxa: list[str]
    true_value: float
    estimated_value: float
    absolute_error: float


@dataclass(slots=True)
class ContinuousValidationReplicate:
    replicate: int
    simulation_model: str
    reconstruction_model: str
    observations: list[ContinuousValidationNodeObservation]
    mean_absolute_error: float
    root_absolute_error: float


@dataclass(slots=True)
class ContinuousAncestralValidationReport:
    tree_path: Path
    simulation_model: str
    reconstruction_model: str
    replicates: int
    node_count: int
    mean_absolute_error: float
    root_mean_absolute_error: float
    root_mean_squared_error: float
    replicate_rows: list[ContinuousValidationReplicate]


@dataclass(slots=True)
class DiscreteValidationNodeObservation:
    node: str
    descendant_taxa: list[str]
    true_state: str
    estimated_state: str
    true_state_probability: float
    correct: bool
    ambiguous: bool


@dataclass(slots=True)
class DiscreteValidationReplicate:
    replicate: int
    simulation_model: str
    reconstruction_model: str
    observations: list[DiscreteValidationNodeObservation]
    accuracy: float
    mean_true_state_probability: float
    calibration_gap: float


@dataclass(slots=True)
class DiscreteAncestralValidationReport:
    tree_path: Path
    simulation_model: str
    reconstruction_model: str
    replicates: int
    node_count: int
    accuracy: float
    mean_true_state_probability: float
    uncertainty_calibration_gap: float
    replicate_rows: list[DiscreteValidationReplicate]


def validate_continuous_ancestral_reconstruction(
    tree_path: Path,
    *,
    simulation_model: str = "brownian",
    reconstruction_model: str | None = None,
    replicates: int = 5,
    seed: int = 1,
    root_state: float = 0.0,
    sigma: float = 0.5,
    alpha: float = 1.0,
    theta: float = 0.0,
) -> ContinuousAncestralValidationReport:
    """Validate continuous ancestral reconstruction against repeated simulated truth."""
    if reconstruction_model is None:
        reconstruction_model = simulation_model
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    replicate_rows: list[ContinuousValidationReplicate] = []
    root_signature = node_signature(load_tree(tree_path).root)
    for replicate in range(replicates):
        if simulation_model == "brownian":
            simulation = simulate_brownian_traits(
                tree_path,
                root_state=root_state,
                sigma=sigma,
                seed=seed + replicate,
            )
        elif simulation_model == "ou":
            simulation = simulate_ou_traits(
                tree_path,
                root_state=root_state,
                sigma=sigma,
                alpha=alpha,
                theta=theta,
                seed=seed + replicate,
            )
        else:
            raise ValueError(
                f"unsupported continuous simulation model: {simulation_model}"
            )
        trait_path = _write_continuous_trait_table(simulation)
        try:
            reconstruction = reconstruct_continuous_ancestral_states(
                tree_path,
                trait_path,
                trait="value",
                model=reconstruction_model,
                alpha=alpha,
            )
        finally:
            trait_path.unlink(missing_ok=True)
        truth_by_node = {
            row.node: row for row in simulation.node_values if not row.is_tip
        }
        observations = [
            ContinuousValidationNodeObservation(
                node=estimate.node,
                descendant_taxa=estimate.descendant_taxa,
                true_value=truth_by_node[estimate.node].value,
                estimated_value=estimate.estimate,
                absolute_error=abs(
                    estimate.estimate - truth_by_node[estimate.node].value
                ),
            )
            for estimate in reconstruction.estimates
            if not estimate.is_tip
        ]
        root_error = next(
            row.absolute_error for row in observations if row.node == root_signature
        )
        replicate_rows.append(
            ContinuousValidationReplicate(
                replicate=replicate + 1,
                simulation_model=simulation_model,
                reconstruction_model=reconstruction_model,
                observations=observations,
                mean_absolute_error=sum(row.absolute_error for row in observations)
                / len(observations),
                root_absolute_error=root_error,
            )
        )
    all_observations = [
        row for replicate in replicate_rows for row in replicate.observations
    ]
    return ContinuousAncestralValidationReport(
        tree_path=tree_path,
        simulation_model=simulation_model,
        reconstruction_model=reconstruction_model,
        replicates=replicates,
        node_count=len(all_observations),
        mean_absolute_error=sum(row.absolute_error for row in all_observations)
        / len(all_observations),
        root_mean_absolute_error=sum(row.root_absolute_error for row in replicate_rows)
        / len(replicate_rows),
        root_mean_squared_error=(
            sum(row.root_absolute_error**2 for row in replicate_rows)
            / len(replicate_rows)
        )
        ** 0.5,
        replicate_rows=replicate_rows,
    )


def validate_discrete_ancestral_reconstruction(
    tree_path: Path,
    *,
    reconstruction_model: str = "fitch",
    replicates: int = 5,
    seed: int = 1,
    states: list[str] | None = None,
    transition_rate: float = 1.0,
    root_state: str | None = None,
) -> DiscreteAncestralValidationReport:
    """Validate discrete ancestral reconstruction against repeated simulated truth."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    simulation_states = states or ["state_a", "state_b", "state_c"]
    replicate_rows: list[DiscreteValidationReplicate] = []
    for replicate in range(replicates):
        simulation = None
        reconstruction = None
        for attempt in range(32):
            simulation = simulate_discrete_traits(
                tree_path,
                states=simulation_states,
                transition_rate=transition_rate,
                root_state=root_state,
                seed=seed + replicate + attempt,
            )
            trait_path = _write_discrete_trait_table(simulation)
            try:
                reconstruction = reconstruct_discrete_ancestral_states(
                    tree_path,
                    trait_path,
                    trait="state",
                    model=reconstruction_model,
                )
                break
            except AncestralReconstructionError as error:
                if "at least two observed states" not in str(error):
                    raise
            finally:
                trait_path.unlink(missing_ok=True)
        if simulation is None or reconstruction is None:
            raise AncestralReconstructionError(
                "discrete ancestral validation could not generate a usable multi-state simulated replicate"
            )
        truth_by_node = {
            row.node: row for row in simulation.node_states if not row.is_tip
        }
        observations = [
            DiscreteValidationNodeObservation(
                node=estimate.node,
                descendant_taxa=estimate.descendant_taxa,
                true_state=truth_by_node[estimate.node].state,
                estimated_state=estimate.most_likely_state,
                true_state_probability=estimate.state_probabilities.get(
                    truth_by_node[estimate.node].state, 0.0
                ),
                correct=estimate.most_likely_state
                == truth_by_node[estimate.node].state,
                ambiguous=estimate.ambiguous,
            )
            for estimate in reconstruction.estimates
            if not estimate.is_tip
        ]
        accuracy = sum(1 for row in observations if row.correct) / len(observations)
        mean_true_probability = sum(
            row.true_state_probability for row in observations
        ) / len(observations)
        replicate_rows.append(
            DiscreteValidationReplicate(
                replicate=replicate + 1,
                simulation_model=simulation.model,
                reconstruction_model=reconstruction_model,
                observations=observations,
                accuracy=accuracy,
                mean_true_state_probability=mean_true_probability,
                calibration_gap=abs(mean_true_probability - accuracy),
            )
        )
    all_observations = [
        row for replicate in replicate_rows for row in replicate.observations
    ]
    accuracy = sum(1 for row in all_observations if row.correct) / len(all_observations)
    mean_true_state_probability = sum(
        row.true_state_probability for row in all_observations
    ) / len(all_observations)
    return DiscreteAncestralValidationReport(
        tree_path=tree_path,
        simulation_model="symmetric-discrete",
        reconstruction_model=reconstruction_model,
        replicates=replicates,
        node_count=len(all_observations),
        accuracy=accuracy,
        mean_true_state_probability=mean_true_state_probability,
        uncertainty_calibration_gap=abs(mean_true_state_probability - accuracy),
        replicate_rows=replicate_rows,
    )


def _write_continuous_trait_table(simulation) -> Path:
    path = Path(
        tempfile.mkstemp(prefix="bijux-ancestral-continuous-", suffix=".tsv")[1]
    )
    write_taxon_rows(
        path,
        columns=["taxon", "value"],
        rows=[
            {"taxon": row.taxon, "value": format(row.value, ".15g")}
            for row in simulation.traits
        ],
    )
    return path


def _write_discrete_trait_table(simulation) -> Path:
    path = Path(tempfile.mkstemp(prefix="bijux-ancestral-discrete-", suffix=".tsv")[1])
    write_taxon_rows(
        path,
        columns=["taxon", "state"],
        rows=[{"taxon": row.taxon, "state": row.state} for row in simulation.traits],
    )
    return path
