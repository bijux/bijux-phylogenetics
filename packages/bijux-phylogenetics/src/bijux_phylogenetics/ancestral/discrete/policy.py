from __future__ import annotations

import math

import numpy

from .models import DiscreteRerootingMethodCompatibility


def resolve_discrete_model_name(model: str) -> str:
    if model == "meristic":
        raise ValueError(
            "geiger::fitDiscrete(model='meristic') is explicitly excluded this round: "
            "local geiger uses a distinct integer-state meristic contract, and bijux ordered-state Mk support is not claimed as meristic parity"
        )
    aliases = {
        "fitch": "fitch",
        "equal-rates": "equal-rates",
        "er": "equal-rates",
        "symmetric": "symmetric",
        "sym": "symmetric",
        "all-rates-different": "all-rates-different",
        "ard": "all-rates-different",
    }
    resolved = aliases.get(model)
    if resolved is None:
        raise ValueError(f"unsupported discrete ancestral model: {model}")
    return resolved


def rerooting_method_compatibility(
    *,
    model: str,
    state_ordering: str,
    root_prior_mode: str,
) -> DiscreteRerootingMethodCompatibility:
    notes: list[str] = []
    reference_model: str | None = None
    reference_root_prior_mode: str | None = None
    comparable = True
    if model == "fitch":
        comparable = False
        notes.append(
            "phytools::rerootingMethod is a likelihood Mk marginal-probability reference and does not apply to Fitch parsimony reconstructions"
        )
    elif model == "all-rates-different":
        comparable = False
        notes.append(
            "phytools::rerootingMethod is invalid for non-symmetric Q matrices such as all-rates-different models in phytools 2.5.2"
        )
    elif model == "symmetric":
        reference_model = "SYM"
    elif model == "equal-rates":
        reference_model = "ER"
    if state_ordering != "unordered":
        comparable = False
        notes.append(
            "phytools::rerootingMethod does not provide a governed ordered-transition parity surface in this repository"
        )
    if root_prior_mode != "equal":
        comparable = False
        notes.append(
            "phytools::rerootingMethod inherits fitMk's default equal root prior; empirical or fixed root-prior runs remain Bijux sensitivity scenarios without direct rerootingMethod parity"
        )
    if comparable:
        reference_root_prior_mode = "equal"
    return DiscreteRerootingMethodCompatibility(
        comparable=comparable,
        reference_model=reference_model,
        reference_root_prior_mode=reference_root_prior_mode,
        notes=notes,
    )


def resolve_state_order(
    observed_states: list[str],
    *,
    state_ordering: str,
    ordered_states: list[str] | None,
) -> list[str]:
    if state_ordering == "unordered":
        return list(observed_states)
    if ordered_states is None:
        return list(observed_states)
    missing_states = sorted(set(observed_states) - set(ordered_states))
    if missing_states:
        raise ValueError(
            "ordered discrete ancestral reconstruction is missing observed states: "
            + ", ".join(missing_states)
        )
    return list(ordered_states)


def parameter_count(
    state_count: int,
    *,
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> int:
    allowed_pair_count = len(allowed_transition_pairs)
    if model == "equal-rates":
        return 1 if allowed_pair_count > 0 else 0
    if state_ordering == "ordered":
        edge_count = max(
            sum(
                1
                for left_index in range(state_count)
                for right_index in range(left_index + 1, state_count)
                if transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
                and transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
            ),
            0,
        )
        if model == "symmetric":
            return edge_count
        return allowed_pair_count
    if model == "symmetric":
        return max(
            sum(
                1
                for left_index in range(state_count)
                for right_index in range(left_index + 1, state_count)
                if transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
                and transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
            ),
            0,
        )
    return allowed_pair_count


def transition_allowed(
    left_index: int,
    right_index: int,
    *,
    state_count: int,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> bool:
    if left_index == right_index:
        return False
    allowed_by_order = state_ordering == "unordered" or (
        abs(left_index - right_index) == 1
        and max(left_index, right_index) < state_count
    )
    if not allowed_by_order:
        return False
    return (left_index, right_index) in allowed_transition_pairs


def resolve_allowed_transition_pairs(
    state_order: list[str],
    *,
    model: str,
    state_ordering: str,
    allowed_transition_pairs: list[tuple[str, str]] | None,
) -> set[tuple[int, int]]:
    state_to_index = {state: index for index, state in enumerate(state_order)}
    if allowed_transition_pairs is None:
        pairs = {
            (left_index, right_index)
            for left_index in range(len(state_order))
            for right_index in range(len(state_order))
            if left_index != right_index
        }
    else:
        pairs: set[tuple[int, int]] = set()
        for source_state, target_state in allowed_transition_pairs:
            if source_state not in state_to_index:
                raise ValueError(
                    "allowed transition source state is not present in the analyzed state vocabulary: "
                    f"{source_state}"
                )
            if target_state not in state_to_index:
                raise ValueError(
                    "allowed transition target state is not present in the analyzed state vocabulary: "
                    f"{target_state}"
                )
            if source_state == target_state:
                raise ValueError(
                    "allowed transition pairs must connect distinct states"
                )
            pairs.add((state_to_index[source_state], state_to_index[target_state]))
    if model == "symmetric":
        asymmetric_pairs = [
            (left_index, right_index)
            for left_index, right_index in sorted(pairs)
            if (right_index, left_index) not in pairs
        ]
        if asymmetric_pairs:
            raise ValueError(
                "symmetric discrete ancestral reconstruction requires bidirectional allowed transitions"
            )
    filtered_pairs = {
        (left_index, right_index)
        for left_index, right_index in pairs
        if left_index != right_index
        and (state_ordering == "unordered" or abs(left_index - right_index) == 1)
    }
    if not filtered_pairs:
        raise ValueError(
            "discrete ancestral reconstruction requires at least one allowed transition after applying constraints"
        )
    return filtered_pairs


def uniform_root_prior(state_count: int) -> numpy.ndarray:
    return numpy.full(state_count, 1.0 / state_count, dtype=float)


def empirical_root_prior(
    state_order: list[str], state_counts: dict[str, int]
) -> numpy.ndarray:
    return normalize_array(
        numpy.array(
            [float(state_counts.get(state, 0)) for state in state_order],
            dtype=float,
        )
    )


def fixed_root_prior(state_order: list[str], fixed_root_state: str) -> numpy.ndarray:
    if fixed_root_state not in state_order:
        raise ValueError(
            "fixed root state is not available in the analyzed state vocabulary: "
            f"{fixed_root_state}"
        )
    prior = numpy.zeros(len(state_order), dtype=float)
    prior[state_order.index(fixed_root_state)] = 1.0
    return prior


def resolve_root_prior(
    state_order: list[str],
    *,
    state_counts: dict[str, int],
    mode: str,
    fixed_root_state: str | None,
    default_root_prior: numpy.ndarray | None = None,
) -> numpy.ndarray:
    if mode == "equal":
        if fixed_root_state is not None:
            raise ValueError("fixed_root_state requires root_prior_mode 'fixed'")
        if default_root_prior is not None:
            return default_root_prior
        return uniform_root_prior(len(state_order))
    if mode == "empirical":
        if fixed_root_state is not None:
            raise ValueError("fixed_root_state requires root_prior_mode 'fixed'")
        return empirical_root_prior(state_order, state_counts)
    if mode == "fixed":
        if fixed_root_state is None:
            raise ValueError("root_prior_mode 'fixed' requires a fixed_root_state")
        return fixed_root_prior(state_order, fixed_root_state)
    raise ValueError(f"unsupported discrete ancestral root prior mode: {mode}")


def reported_discrete_log_likelihood(
    log_likelihood: float,
    *,
    root_prior_mode: str,
    state_count: int,
) -> float:
    if root_prior_mode != "equal" or state_count < 2:
        return log_likelihood
    return log_likelihood + math.log(state_count)


def normalize_array(values: numpy.ndarray) -> numpy.ndarray:
    total = float(values.sum())
    if total <= 0.0:
        return numpy.full(values.shape[0], 1.0 / values.shape[0], dtype=float)
    return values / total
