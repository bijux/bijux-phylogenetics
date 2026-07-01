from __future__ import annotations

from bijux_phylogenetics.ancestral.discrete.likelihood.parameter_layout import (
    broadcast_discrete_mk_log_parameters,
    build_discrete_mk_parameter_layout,
)


def _unordered_pairs() -> set[tuple[int, int]]:
    return {
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 2),
        (2, 0),
        (2, 1),
    }


def test_equal_rates_layout_maps_every_allowed_arc_to_one_parameter() -> None:
    layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="equal-rates",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )

    assert layout.parameter_count == 1
    assert len(layout.arcs) == 6
    assert {arc.parameter_index for arc in layout.arcs} == {0}


def test_symmetric_layout_groups_reverse_arcs_into_shared_slots() -> None:
    layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="symmetric",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )

    assert layout.parameter_count == 3
    assert layout.parameter_index_for_pair(0, 1) == layout.parameter_index_for_pair(1, 0)
    assert layout.parameter_index_for_pair(0, 2) == layout.parameter_index_for_pair(2, 0)
    assert layout.parameter_index_for_pair(1, 2) == layout.parameter_index_for_pair(2, 1)


def test_ard_layout_keeps_each_allowed_direction_in_its_own_slot() -> None:
    layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )

    assert layout.parameter_count == 6
    assert layout.parameter_index_for_pair(0, 1) != layout.parameter_index_for_pair(1, 0)
    assert layout.parameter_index_for_pair(0, 2) != layout.parameter_index_for_pair(2, 0)
    assert layout.parameter_index_for_pair(1, 2) != layout.parameter_index_for_pair(2, 1)


def test_broadcast_equal_rates_parameters_to_symmetric_and_ard_layouts() -> None:
    equal_layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="equal-rates",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )
    symmetric_layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="symmetric",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )
    ard_layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )

    symmetric_parameters = broadcast_discrete_mk_log_parameters(
        equal_layout,
        symmetric_layout,
        (0.4,),
    )
    ard_parameters = broadcast_discrete_mk_log_parameters(
        equal_layout,
        ard_layout,
        (0.4,),
    )

    assert symmetric_parameters == (0.4, 0.4, 0.4)
    assert ard_parameters == (0.4, 0.4, 0.4, 0.4, 0.4, 0.4)


def test_broadcast_symmetric_parameters_to_ard_layout_reuses_undirected_pair_values() -> (
    None
):
    symmetric_layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="symmetric",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )
    ard_layout = build_discrete_mk_parameter_layout(
        state_order=["north", "south", "west"],
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs=_unordered_pairs(),
    )

    ard_parameters = broadcast_discrete_mk_log_parameters(
        symmetric_layout,
        ard_layout,
        (0.1, 0.2, 0.3),
    )

    assert ard_parameters == (0.1, 0.2, 0.1, 0.3, 0.2, 0.3)
