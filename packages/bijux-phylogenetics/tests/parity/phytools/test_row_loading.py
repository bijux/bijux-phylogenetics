from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.parity.phytools.runner.comparison import load_rows_table


def test_load_rows_table_preserves_discrete_state_identity_and_booleans(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "fitmk-rate-matrix.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "source_state\ttarget_state\ttransition_allowed\tstep_distance\trate",
                "0\t1\tTrue\t1\t70.77642178117888",
                "north\tsouth\tfalse\t2\t0.15",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = load_rows_table(rows_path)

    assert rows == [
        {
            "source_state": "0",
            "target_state": "1",
            "transition_allowed": True,
            "step_distance": 1,
            "rate": 70.77642178117888,
        },
        {
            "source_state": "north",
            "target_state": "south",
            "transition_allowed": False,
            "step_distance": 2,
            "rate": 0.15,
        },
    ]


def test_load_rows_table_preserves_stochastic_branch_summary_identity(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "stochastic-map-summary-rows.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "row_kind\tlabel\tmean_value\tlower_95_interval\tupper_95_interval\tpresence_fraction",
                "branch_state_occupancy\tA|B->A:0\t1.0\t1.0\t1.0\t1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = load_rows_table(rows_path)

    assert rows == [
        {
            "row_kind": "branch_state_occupancy",
            "label": "A|B->A:0",
            "mean_value": 1.0,
            "lower_95_interval": 1.0,
            "upper_95_interval": 1.0,
            "presence_fraction": 1.0,
        }
    ]


def test_load_rows_table_preserves_rerooting_state_identity(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "rerooting-method-node-probabilities.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "node\tstate\tprobability",
                "A|B|C\t0\t0.75",
                "A|B|C\tnorth\t0.25",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = load_rows_table(rows_path)

    assert rows == [
        {
            "node": "A|B|C",
            "state": "0",
            "probability": 0.75,
        },
        {
            "node": "A|B|C",
            "state": "north",
            "probability": 0.25,
        },
    ]


def test_load_rows_table_preserves_stochastic_map_row_identity(
    tmp_path: Path,
) -> None:
    rows_path = tmp_path / "stochastic-map-summary-rows.tsv"
    rows_path.write_text(
        "\n".join(
            [
                "row_kind\tlabel\tmean_value\tlower_95_interval\tupper_95_interval\tpresence_fraction",
                "transition_count\tnorth->south\t1.5\t0\t4\t0.75",
                "state_time\t0\t2.25\t1\t3\t1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = load_rows_table(rows_path)

    assert rows == [
        {
            "row_kind": "transition_count",
            "label": "north->south",
            "mean_value": 1.5,
            "lower_95_interval": 0,
            "upper_95_interval": 4,
            "presence_fraction": 0.75,
        },
        {
            "row_kind": "state_time",
            "label": "0",
            "mean_value": 2.25,
            "lower_95_interval": 1,
            "upper_95_interval": 3,
            "presence_fraction": 1.0,
        },
    ]
