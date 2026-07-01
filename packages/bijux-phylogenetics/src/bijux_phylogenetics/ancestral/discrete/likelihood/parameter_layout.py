from __future__ import annotations

from dataclasses import dataclass

from ..policy import transition_allowed


@dataclass(frozen=True, slots=True)
class DiscreteMkParameterArc:
    """One directed transition and the shared parameter slot that controls it."""

    source_index: int
    target_index: int
    source_state: str
    target_state: str
    parameter_index: int


@dataclass(frozen=True, slots=True)
class DiscreteMkParameterLayout:
    """One stable parameter layout for an ER, SYM, or ARD discrete Mk surface."""

    state_order: tuple[str, ...]
    model: str
    state_ordering: str
    allowed_transition_pairs: frozenset[tuple[int, int]]
    arcs: tuple[DiscreteMkParameterArc, ...]

    @property
    def parameter_count(self) -> int:
        if not self.arcs:
            return 0
        return 1 + max(arc.parameter_index for arc in self.arcs)

    def parameter_index_for_pair(self, source_index: int, target_index: int) -> int:
        for arc in self.arcs:
            if arc.source_index == source_index and arc.target_index == target_index:
                return arc.parameter_index
        raise KeyError((source_index, target_index))


def build_discrete_mk_parameter_layout(
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> DiscreteMkParameterLayout:
    """Build the stable parameter-slot order for one discrete Mk surface."""
    arcs: list[DiscreteMkParameterArc] = []
    parameter_index = 0
    state_count = len(state_order)
    if model == "equal-rates":
        for left_index in range(state_count):
            for right_index in range(state_count):
                if not transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                arcs.append(
                    DiscreteMkParameterArc(
                        source_index=left_index,
                        target_index=right_index,
                        source_state=state_order[left_index],
                        target_state=state_order[right_index],
                        parameter_index=0,
                    )
                )
        return DiscreteMkParameterLayout(
            state_order=tuple(state_order),
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=frozenset(allowed_transition_pairs),
            arcs=tuple(arcs),
        )
    if model == "symmetric":
        for left_index in range(state_count):
            for right_index in range(left_index + 1, state_count):
                if not transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                if not transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                arcs.append(
                    DiscreteMkParameterArc(
                        source_index=left_index,
                        target_index=right_index,
                        source_state=state_order[left_index],
                        target_state=state_order[right_index],
                        parameter_index=parameter_index,
                    )
                )
                arcs.append(
                    DiscreteMkParameterArc(
                        source_index=right_index,
                        target_index=left_index,
                        source_state=state_order[right_index],
                        target_state=state_order[left_index],
                        parameter_index=parameter_index,
                    )
                )
                parameter_index += 1
        return DiscreteMkParameterLayout(
            state_order=tuple(state_order),
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=frozenset(allowed_transition_pairs),
            arcs=tuple(arcs),
        )
    for left_index in range(state_count):
        for right_index in range(state_count):
            if not transition_allowed(
                left_index,
                right_index,
                state_count=state_count,
                state_ordering=state_ordering,
                allowed_transition_pairs=allowed_transition_pairs,
            ):
                continue
            arcs.append(
                DiscreteMkParameterArc(
                    source_index=left_index,
                    target_index=right_index,
                    source_state=state_order[left_index],
                    target_state=state_order[right_index],
                    parameter_index=parameter_index,
                )
            )
            parameter_index += 1
    return DiscreteMkParameterLayout(
        state_order=tuple(state_order),
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=frozenset(allowed_transition_pairs),
        arcs=tuple(arcs),
    )


def broadcast_discrete_mk_log_parameters(
    source_layout: DiscreteMkParameterLayout,
    target_layout: DiscreteMkParameterLayout,
    source_log_parameters: tuple[float, ...],
) -> tuple[float, ...]:
    """Project one lower-order Mk parameterization into a richer compatible layout."""
    if source_layout.model == target_layout.model:
        if source_layout.parameter_count != target_layout.parameter_count:
            raise ValueError("identical-model discrete Mk layouts must have equal arity")
        return source_log_parameters
    if source_layout.state_order != target_layout.state_order:
        raise ValueError("cannot broadcast discrete Mk parameters across state orders")
    if source_layout.allowed_transition_pairs != target_layout.allowed_transition_pairs:
        raise ValueError(
            "cannot broadcast discrete Mk parameters across transition constraints"
        )
    if source_layout.model == "equal-rates":
        return (source_log_parameters[0],) * target_layout.parameter_count
    if source_layout.model == "symmetric" and target_layout.model == "all-rates-different":
        target_parameters = [0.0] * target_layout.parameter_count
        for arc in target_layout.arcs:
            left_index = min(arc.source_index, arc.target_index)
            right_index = max(arc.source_index, arc.target_index)
            source_parameter_index = source_layout.parameter_index_for_pair(
                left_index,
                right_index,
            )
            target_parameters[arc.parameter_index] = source_log_parameters[
                source_parameter_index
            ]
        return tuple(target_parameters)
    raise ValueError(
        "unsupported discrete Mk parameter broadcast: "
        f"{source_layout.model} -> {target_layout.model}"
    )
