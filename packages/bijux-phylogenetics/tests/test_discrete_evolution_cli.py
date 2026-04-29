from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.cli import main


FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_discrete_evolution_validate_and_imbalance_cli_report_findings(capsys) -> None:
    validate_exit = main(
        [
            "discrete-evolution",
            "validate-coding",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography_invalid.tsv")),
            "--trait",
            "region",
            "--allowed-states",
            "north,south,island",
            "--json",
        ]
    )
    validate_payload = json.loads(capsys.readouterr().out)
    assert validate_exit == 0
    assert validate_payload["metrics"]["valid"] is False
    assert validate_payload["metrics"]["issue_count"] == 2

    imbalance_exit = main(
        [
            "discrete-evolution",
            "imbalance",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography_single_state.tsv")),
            "--trait",
            "region",
            "--json",
        ]
    )
    imbalance_payload = json.loads(capsys.readouterr().out)
    assert imbalance_exit == 0
    assert imbalance_payload["metrics"]["warning_count"] == 2
    assert "only one observed state remains after pruning to usable tree taxa" in imbalance_payload["warnings"]


def test_discrete_evolution_model_and_compare_cli_write_tables(tmp_path: Path, capsys) -> None:
    node_table = tmp_path / "node-probabilities.tsv"
    transitions_table = tmp_path / "transitions.tsv"
    comparison_table = tmp_path / "model-comparison.tsv"

    model_exit = main(
        [
            "discrete-evolution",
            "model",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "all-rates-different",
            "--node-table-out",
            str(node_table),
            "--transitions-out",
            str(transitions_table),
            "--json",
        ]
    )
    model_payload = json.loads(capsys.readouterr().out)
    assert model_exit == 0
    assert model_payload["metrics"]["model"] == "all-rates-different"
    assert "state_probabilities" in node_table.read_text(encoding="utf-8")
    assert "source_state" in transitions_table.read_text(encoding="utf-8")

    compare_exit = main(
        [
            "discrete-evolution",
            "compare-models",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--table-out",
            str(comparison_table),
            "--json",
        ]
    )
    compare_payload = json.loads(capsys.readouterr().out)
    assert compare_exit == 0
    assert compare_payload["metrics"]["model_count"] == 2
    assert compare_payload["metrics"]["better_model"] in {"equal-rates", "all-rates-different"}
    assert "left_probabilities" in comparison_table.read_text(encoding="utf-8")
