from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

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
    assert (
        "only one observed state remains after pruning to usable tree taxa"
        in imbalance_payload["warnings"]
    )


def test_discrete_evolution_reference_cli_reports_passing_cases(capsys) -> None:
    exit_code = main(["discrete-evolution", "reference", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["case_count"] == 3
    assert payload["metrics"]["all_passed"] is True


def test_discrete_evolution_model_and_compare_cli_write_tables(
    tmp_path: Path, capsys
) -> None:
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
    assert model_payload["metrics"]["strongly_supported_transition_count"] >= 0
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
    assert compare_payload["metrics"]["better_model"] in {
        "equal-rates",
        "all-rates-different",
    }
    assert "left_probabilities" in comparison_table.read_text(encoding="utf-8")


def test_discrete_evolution_model_cli_accepts_symmetric_ordered_states(capsys) -> None:
    exit_code = main(
        [
            "discrete-evolution",
            "model",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "symmetric",
            "--state-ordering",
            "ordered",
            "--ordered-states",
            "north,south,island",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["data"]["model"] == "symmetric"
    assert payload["data"]["state_ordering"] == "ordered"
    assert payload["metrics"]["state_ordering"] == "ordered"


def test_discrete_evolution_model_cli_rejects_meristic_parity_claim(capsys) -> None:
    try:
        main(
            [
                "discrete-evolution",
                "model",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--model",
                "meristic",
                "--state-ordering",
                "ordered",
                "--ordered-states",
                "north,south,island",
                "--json",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    captured = capsys.readouterr()
    assert "explicitly excluded this round" in captured.err
    assert "integer-state meristic contract" in captured.err


def test_discrete_evolution_render_and_report_cli_write_svg_and_html(
    tmp_path: Path, capsys
) -> None:
    svg_path = tmp_path / "geography.svg"
    report_path = tmp_path / "geography-report.html"

    render_exit = main(
        [
            "discrete-evolution",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--layout",
            "circular",
            "--out",
            str(svg_path),
            "--json",
        ]
    )
    render_payload = json.loads(capsys.readouterr().out)
    assert render_exit == 0
    assert render_payload["metrics"]["layout"] == "circular"
    assert 'class="internal-annotation-label"' in svg_path.read_text(encoding="utf-8")

    report_exit = main(
        [
            "discrete-evolution",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--compare-model",
            "all-rates-different",
            "--out",
            str(report_path),
            "--json",
        ]
    )
    report_payload = json.loads(capsys.readouterr().out)
    assert report_exit == 0
    assert report_payload["metrics"]["report_kind"] == "discrete-state-evolution"
    assert report_payload["metrics"]["state_ordering"] == "unordered"
    assert report_path.exists()
    assert report_path.with_suffix(".svg").exists()


def test_discrete_evolution_stochastic_map_and_summary_cli_write_outputs(
    tmp_path: Path, capsys
) -> None:
    collection_path = tmp_path / "stochastic-maps.json"
    summary_path = tmp_path / "stochastic-summary.tsv"
    state_times_path = tmp_path / "stochastic-state-times.tsv"
    branch_occupancy_path = tmp_path / "stochastic-branch-occupancy.tsv"
    count_matrix_path = tmp_path / "stochastic-count-matrix.tsv"
    aggregate_matrix_path = tmp_path / "stochastic-aggregate-matrix.tsv"
    branch_transition_path = tmp_path / "stochastic-branch-transitions.tsv"
    segments_path = tmp_path / "stochastic-segments.tsv"
    events_path = tmp_path / "stochastic-events.tsv"

    stochastic_exit = main(
        [
            "discrete-evolution",
            "stochastic-map",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "symmetric",
            "--replicates",
            "6",
            "--seed",
            "13",
            "--collection-out",
            str(collection_path),
            "--summary-out",
            str(summary_path),
            "--state-times-out",
            str(state_times_path),
            "--branch-occupancy-out",
            str(branch_occupancy_path),
            "--count-matrix-out",
            str(count_matrix_path),
            "--aggregate-matrix-out",
            str(aggregate_matrix_path),
            "--branch-transition-out",
            str(branch_transition_path),
            "--segments-out",
            str(segments_path),
            "--events-out",
            str(events_path),
            "--json",
        ]
    )
    stochastic_payload = json.loads(capsys.readouterr().out)
    assert stochastic_exit == 0
    assert stochastic_payload["metrics"]["requested_replicate_count"] == 6
    assert stochastic_payload["metrics"]["successful_replicate_count"] == 6
    assert stochastic_payload["metrics"]["simulation_failure_count"] == 0
    assert stochastic_payload["metrics"]["mean_total_transition_count"] >= 0.0
    assert stochastic_payload["metrics"]["branch_state_row_count"] > 0
    assert stochastic_payload["metrics"]["count_matrix_row_count"] == 6
    assert stochastic_payload["metrics"]["branch_transition_row_count"] > 0
    assert stochastic_payload["metrics"]["conditioned_on_node_estimates"] is False
    assert stochastic_payload["metrics"]["parameter_count"] == 3
    assert stochastic_payload["metrics"]["optimizer_converged"] is True
    assert stochastic_payload["metrics"]["overparameterized"] is False
    assert stochastic_payload["metrics"]["preferred_model_by_aic"] == "equal-rates"
    assert stochastic_payload["metrics"]["fit_warning_count"] >= 1
    assert collection_path.exists()
    assert state_times_path.exists()
    assert branch_occupancy_path.exists()
    assert count_matrix_path.exists()
    assert aggregate_matrix_path.exists()
    assert branch_transition_path.exists()
    assert segments_path.exists()
    assert events_path.exists()
    assert "transition\tmean_count" in summary_path.read_text(encoding="utf-8")
    assert "state\tmean_time" in state_times_path.read_text(encoding="utf-8")
    assert (
        "branch_index\tparent_node\tchild_node\tstate\tbranch_length\tmean_time"
        in branch_occupancy_path.read_text(encoding="utf-8")
    )
    assert "replicate_index\ttotal_transition_count" in count_matrix_path.read_text(
        encoding="utf-8"
    )
    assert "source_state" in aggregate_matrix_path.read_text(encoding="utf-8")
    assert (
        "branch_index\tparent_node\tchild_node\ttransition\tmean_count"
        in branch_transition_path.read_text(encoding="utf-8")
    )
    assert "replicate_index\tbranch_index" in segments_path.read_text(encoding="utf-8")
    assert (
        "replicate_index\tbranch_index\tparent_node\tchild_node\tevent_index\tsource_state\ttarget_state"
        in events_path.read_text(encoding="utf-8")
    )

    summarize_exit = main(
        [
            "discrete-evolution",
            "summarize-maps",
            str(collection_path),
            "--summary-out",
            str(summary_path),
            "--state-times-out",
            str(state_times_path),
            "--branch-occupancy-out",
            str(branch_occupancy_path),
            "--json",
        ]
    )
    summarize_payload = json.loads(capsys.readouterr().out)
    assert summarize_exit == 0
    assert summarize_payload["metrics"]["replicate_count"] == 6
    assert summarize_payload["metrics"]["mean_total_transition_count"] >= 0.0
    assert summarize_payload["metrics"]["simulation_failure_count"] == 0
    assert summarize_payload["metrics"]["branch_state_row_count"] > 0

    count_exit = main(
        [
            "discrete-evolution",
            "count-maps",
            str(collection_path),
            "--count-matrix-out",
            str(count_matrix_path),
            "--aggregate-matrix-out",
            str(aggregate_matrix_path),
            "--branch-transition-out",
            str(branch_transition_path),
            "--events-out",
            str(events_path),
            "--json",
        ]
    )
    count_payload = json.loads(capsys.readouterr().out)
    assert count_exit == 0
    assert count_payload["metrics"]["replicate_count"] == 6
    assert count_payload["metrics"]["mean_total_transition_count"] >= 0.0
    assert count_payload["metrics"]["count_matrix_row_count"] == 6
    assert count_payload["metrics"]["branch_transition_row_count"] > 0


def test_discrete_evolution_density_maps_cli_writes_branch_density_outputs(
    tmp_path: Path, capsys
) -> None:
    collection_path = tmp_path / "binary-stochastic-maps.json"
    branch_probability_path = tmp_path / "branch-probabilities.tsv"
    density_branch_path = tmp_path / "density-branches.tsv"
    density_slice_path = tmp_path / "density-slices.tsv"
    density_html_path = tmp_path / "density-report.html"

    stochastic_exit = main(
        [
            "discrete-evolution",
            "stochastic-map",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")),
            "--trait",
            "binary_state",
            "--model",
            "equal-rates",
            "--replicates",
            "6",
            "--seed",
            "17",
            "--collection-out",
            str(collection_path),
            "--json",
        ]
    )
    _ = json.loads(capsys.readouterr().out)
    assert stochastic_exit == 0

    density_exit = main(
        [
            "discrete-evolution",
            "density-maps",
            str(collection_path),
            "--resolution",
            "12",
            "--branch-probabilities-out",
            str(branch_probability_path),
            "--density-branches-out",
            str(density_branch_path),
            "--density-slices-out",
            str(density_slice_path),
            "--out",
            str(density_html_path),
            "--json",
        ]
    )
    density_payload = json.loads(capsys.readouterr().out)

    assert density_exit == 0
    assert density_payload["metrics"]["replicate_count"] == 6
    assert density_payload["metrics"]["resolution"] == 12
    assert density_payload["metrics"]["focal_state"] == "1"
    assert density_payload["metrics"]["branch_probability_row_count"] > 0
    assert density_payload["metrics"]["density_branch_row_count"] > 0
    assert density_payload["metrics"]["density_slice_row_count"] > 0
    assert density_payload["metrics"]["rendered_branch_color_count"] > 0
    assert "mean_probability" in branch_probability_path.read_text(encoding="utf-8")
    assert "mean_posterior_probability" in density_branch_path.read_text(
        encoding="utf-8"
    )
    assert "posterior_probability" in density_slice_path.read_text(encoding="utf-8")
    assert density_html_path.exists()
    assert density_html_path.with_suffix(".svg").exists()
