from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.discrete_evolution import (
    StochasticMapBranchHistory,
    StochasticMapCollectionReport,
    StochasticMapModelFitAudit,
    StochasticMapReplicate,
    StochasticMapStateSegment,
    StochasticMapSummaryReport,
    assess_geographic_state_analysis_readiness,
    audit_discrete_state_coding,
    build_biogeographic_interpretation_report,
    compare_discrete_state_models,
    count_discrete_stochastic_map_transitions,
    detect_state_imbalance_problems,
    estimate_ancestral_geographic_states,
    load_stochastic_map_collection,
    render_discrete_state_evolution_report,
    render_stochastic_map_density_artifact,
    render_tree_with_geographic_states,
    run_discrete_state_transition_model,
    simulate_discrete_stochastic_maps,
    simulate_discrete_stochastic_maps_from_fit_report,
    summarize_discrete_stochastic_map_density,
    summarize_discrete_stochastic_maps,
    validate_discrete_state_coding,
    validate_discrete_transition_reference_examples,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_stochastic_map_aggregate_transition_matrix,
    write_stochastic_map_branch_occupancy_table,
    write_stochastic_map_branch_probability_table,
    write_stochastic_map_branch_transition_count_table,
    write_stochastic_map_collection,
    write_stochastic_map_density_branch_table,
    write_stochastic_map_density_slice_table,
    write_stochastic_map_event_table,
    write_stochastic_map_segment_table,
    write_stochastic_map_state_time_table,
    write_stochastic_map_summary_table,
    write_stochastic_map_transition_count_matrix,
    write_transition_summary_table,
)
from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_phytools_comparative_fixture,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

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


def manual_stochastic_map_collection(
    *,
    maps: list[StochasticMapReplicate],
    tree_path: Path | None = None,
) -> StochasticMapCollectionReport:
    report = StochasticMapCollectionReport(
        tree_path=tree_path or Path("manual-tree.nwk"),
        traits_path=Path("manual-traits.tsv"),
        taxon_column="taxon",
        trait="state",
        model="equal-rates",
        state_ordering="unordered",
        ordered_states=[],
        replicates=len(maps),
        seed=0,
        conditioned_on_node_estimates=False,
        fit_audit=StochasticMapModelFitAudit(
            state_order=["0", "1", "2"],
            allowed_transitions=["0->1", "1->0", "1->2", "2->1"],
            parameter_count=1,
            log_likelihood=0.0,
            aic=0.0,
            aicc=0.0,
            overparameterized=False,
            optimizer_converged=True,
            optimizer_iteration_count=0,
            optimizer_function_evaluation_count=0,
            optimizer_hit_lower_parameter_bound=False,
            optimizer_hit_upper_parameter_bound=False,
            baseline_model=None,
            baseline_aic=None,
            baseline_delta_aic=None,
            preferred_model_by_aic=None,
            warnings=[],
        ),
        warnings=[],
        maps=maps,
        failures=[],
        summary=StochasticMapSummaryReport(
            replicate_count=0,
            mean_total_transition_count=0.0,
            lower_95_total_transition_count=0.0,
            upper_95_total_transition_count=0.0,
            rows=[],
            state_time_rows=[],
            branch_occupancy_rows=[],
            simulation_failure_count=0,
            warnings=[],
        ),
    )
    report.summary = summarize_discrete_stochastic_maps(report)
    return report


def test_validate_discrete_state_coding_reports_unsupported_and_delimited_states() -> (
    None
):
    report = validate_discrete_state_coding(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography_invalid.tsv"),
        trait="region",
        allowed_states=["north", "south", "island"],
    )

    assert report.valid is False
    assert report.observed_states == ["island", "north"]
    assert {issue.code for issue in report.issues} == {
        "unsupported-state-delimiter",
        "unsupported-state-label",
    }


def test_audit_discrete_state_coding_tracks_normalized_and_excluded_rows() -> None:
    traits_path = fixture("example_traits_geography_invalid.tsv")
    report = audit_discrete_state_coding(
        fixture("example_tree.nwk"),
        traits_path,
        trait="region",
        allowed_states=["north", "south", "island"],
    )
    assert report.row_count == 4
    assert report.included_row_count == 2
    by_taxon = {row.taxon: row for row in report.rows}
    assert by_taxon["A"].normalized_state == "north"
    assert by_taxon["A"].included is True
    assert by_taxon["B"].issue_code == "unsupported-state-delimiter"
    assert by_taxon["C"].issue_code == "unsupported-state-label"


def test_audit_discrete_state_coding_applies_coding_map_and_tree_overlap(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "geography-audit.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tN\nB\tsouth\nGhost\tisland\n",
        encoding="utf-8",
    )
    report = audit_discrete_state_coding(
        fixture("example_tree.nwk"),
        traits_path,
        trait="region",
        allowed_states=["north", "south", "island"],
        coding_map={"N": "north"},
    )
    by_taxon = {row.taxon: row for row in report.rows}
    assert by_taxon["A"].normalized_state == "north"
    assert by_taxon["A"].included is True
    assert by_taxon["Ghost"].issue_code == "taxon-not-in-tree"


def test_validate_discrete_state_coding_records_ordered_vocabulary() -> None:
    report = validate_discrete_state_coding(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        state_ordering="ordered",
        ordered_states=["north", "south", "island"],
    )
    assert report.valid is True
    assert report.state_ordering == "ordered"
    assert report.ordered_states == ["north", "south", "island"]


def test_validate_discrete_state_coding_rejects_missing_ordered_state_labels() -> None:
    report = validate_discrete_state_coding(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        state_ordering="ordered",
        ordered_states=["north", "south"],
    )
    assert report.valid is False
    assert any(issue.code == "unordered-state-vocabulary" for issue in report.issues)


def test_detect_state_imbalance_problems_flags_single_state_and_dominance() -> None:
    report = detect_state_imbalance_problems(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography_single_state.tsv"),
        trait="region",
    )

    assert report.state_counts == {"north": 4}
    assert {warning.code for warning in report.warnings} == {
        "single-state-dataset",
        "dominant-state-skew",
    }


def test_run_discrete_state_transition_model_returns_transition_matrix_and_events() -> (
    None
):
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )

    assert report.transition_model.model == "equal-rates"
    assert report.transition_model.state_order == ["island", "north", "south"]
    assert report.transition_summary.branch_count == 6
    assert report.transition_summary.transition_count >= 1
    assert report.transition_summary.strongly_supported_transition_count >= 0
    assert (
        len(report.transition_summary.support_rows)
        == report.transition_summary.branch_count
    )
    assert all(
        sum(row.target_rates.values()) > 0.99
        for row in report.transition_model.transition_matrix
    )
    assert report.transition_model.uncertainty.rows


def test_run_discrete_state_transition_model_reports_instability_and_dominant_bias() -> (
    None
):
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    assert isinstance(report.instability.sparse_states, list)
    assert isinstance(report.instability.zero_support_transitions, list)
    assert report.dominant_state_bias.biased is False


def test_run_discrete_state_transition_model_supports_symmetric_rates() -> None:
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="symmetric",
    )
    assert report.transition_model.model == "symmetric"
    assert report.transition_model.parameter_count == 3


def test_run_discrete_state_transition_model_supports_ordered_states() -> None:
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        state_ordering="ordered",
        ordered_states=["north", "south", "island"],
    )
    assert report.state_ordering == "ordered"
    assert report.ordered_states == ["north", "south", "island"]
    assert report.transition_model.ordered_states == ["north", "south", "island"]


def test_run_discrete_state_transition_model_rejects_meristic_parity_claim() -> None:
    with pytest.raises(ValueError) as excinfo:
        run_discrete_state_transition_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            model="meristic",
            state_ordering="ordered",
            ordered_states=["north", "south", "island"],
        )

    assert "explicitly excluded this round" in str(excinfo.value)
    assert "integer-state meristic contract" in str(excinfo.value)


def test_assess_geographic_state_analysis_readiness_allows_balanced_example() -> None:
    report = assess_geographic_state_analysis_readiness(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
    )
    assert report.valid is True
    assert report.blockers == []


def test_assess_geographic_state_analysis_readiness_blocks_sparse_or_single_state_inputs(
    tmp_path: Path,
) -> None:
    sparse_traits = tmp_path / "sparse-geography.tsv"
    sparse_traits.write_text(
        "taxon\tregion\nA\tnorth\nB\tsouth\nC\teast\nD\twest\n",
        encoding="utf-8",
    )
    sparse_report = assess_geographic_state_analysis_readiness(
        fixture("example_tree.nwk"),
        sparse_traits,
        trait="region",
    )
    assert sparse_report.valid is False
    assert any("too sparse" in blocker for blocker in sparse_report.blockers)

    with pytest.raises(AncestralReconstructionError):
        estimate_ancestral_geographic_states(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography_single_state.tsv"),
            trait="region",
        )


def test_estimate_ancestral_geographic_states_and_compare_models_return_node_differences() -> (
    None
):
    report = estimate_ancestral_geographic_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="all-rates-different",
    )
    comparison = compare_discrete_state_models(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
    )

    assert report.estimates[0].state_probabilities
    assert comparison.better_model in {"equal-rates", "all-rates-different"}
    assert len(comparison.rows) == 2
    assert len(comparison.node_differences) == len(report.estimates)
    assert comparison.sensitive_region_count == len(comparison.sensitive_regions)


def test_write_node_probability_table_and_render_report_outputs_files(
    tmp_path: Path,
) -> None:
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    table_path = tmp_path / "node-probabilities.tsv"
    svg_path = tmp_path / "geography.svg"
    html_path = tmp_path / "discrete-report.html"

    write_node_state_probability_table(table_path, report)
    render_result = render_tree_with_geographic_states(
        fixture("example_tree.nwk"),
        report,
        out_path=svg_path,
    )
    report_result = render_discrete_state_evolution_report(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_geography.tsv"),
        trait="region",
        out_path=html_path,
        model="equal-rates",
        compare_model="all-rates-different",
    )

    assert "state_probabilities" in table_path.read_text(encoding="utf-8")
    assert render_result.rendered_internal_annotation_count == 3
    assert render_result.rendered_internal_pie_count == 3
    html = html_path.read_text(encoding="utf-8")
    assert "discrete-state-evolution" in html
    assert "limitations" in html
    assert report_result.report_kind == "discrete-state-evolution"
    assert (
        report_result.machine_manifest["likelihood_method"]
        == "deterministic-node-probability"
    )
    assert "biogeographic-interpretation" in report_result.machine_manifest["sections"]
    assert "limitations" in report_result.machine_manifest["sections"]
    assert report_result.machine_manifest["limitations"]


def test_build_biogeographic_interpretation_report_separates_results_from_guidance() -> (
    None
):
    report = build_biogeographic_interpretation_report(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        compare_model="all-rates-different",
    )
    assert any(row.label == "root_state" for row in report.computed_results)
    assert report.coding_audit_summary["included_row_count"] > 0
    assert report.interpretation_guidance


def test_write_transition_summary_table_exports_branch_state_changes(
    tmp_path: Path,
) -> None:
    report = run_discrete_state_transition_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    table_path = tmp_path / "transitions.tsv"

    write_transition_summary_table(table_path, report)

    contents = table_path.read_text(encoding="utf-8")
    assert (
        "parent_node\tchild_node\tsource_state\ttarget_state\tchanged\tsupport\tstrongly_supported"
        in contents
    )
    assert "false" in contents or "true" in contents


def test_write_discrete_model_comparison_table_exports_node_probabilities(
    tmp_path: Path,
) -> None:
    comparison = compare_discrete_state_models(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
    )
    table_path = tmp_path / "model-comparison.tsv"

    write_discrete_model_comparison_table(table_path, comparison)

    contents = table_path.read_text(encoding="utf-8")
    assert "left_probabilities" in contents
    assert "right_probabilities" in contents


def test_validate_discrete_transition_reference_examples_matches_expected_cases() -> (
    None
):
    report = validate_discrete_transition_reference_examples()
    assert report.case_count == 3
    assert report.all_passed is True
    assert all(
        observation.max_rate_delta <= report.tolerance
        for observation in report.observations
    )


def test_simulate_discrete_stochastic_maps_reports_uncertainty_and_roundtrips(
    tmp_path: Path,
) -> None:
    report = simulate_discrete_stochastic_maps(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="symmetric",
        replicates=8,
        seed=11,
    )
    summary = summarize_discrete_stochastic_maps(report)
    collection_path = tmp_path / "stochastic-maps.json"
    summary_path = tmp_path / "stochastic-summary.tsv"
    state_times_path = tmp_path / "stochastic-state-times.tsv"
    branch_occupancy_path = tmp_path / "stochastic-branch-occupancy.tsv"
    count_matrix_path = tmp_path / "stochastic-count-matrix.tsv"
    aggregate_matrix_path = tmp_path / "stochastic-aggregate-matrix.tsv"
    branch_transition_path = tmp_path / "stochastic-branch-transitions.tsv"
    segments_path = tmp_path / "stochastic-segments.tsv"
    events_path = tmp_path / "stochastic-events.tsv"
    count_report = count_discrete_stochastic_map_transitions(report)

    write_stochastic_map_collection(collection_path, report)
    write_stochastic_map_summary_table(summary_path, summary)
    write_stochastic_map_state_time_table(state_times_path, summary)
    write_stochastic_map_branch_occupancy_table(branch_occupancy_path, summary)
    write_stochastic_map_transition_count_matrix(count_matrix_path, count_report)
    write_stochastic_map_aggregate_transition_matrix(
        aggregate_matrix_path, count_report
    )
    write_stochastic_map_branch_transition_count_table(
        branch_transition_path,
        count_report,
    )
    write_stochastic_map_segment_table(segments_path, report)
    write_stochastic_map_event_table(events_path, report)
    reloaded = load_stochastic_map_collection(collection_path)

    assert report.summary.replicate_count == 8
    assert report.summary.rows
    assert report.summary.state_time_rows
    assert report.summary.mean_total_transition_count >= 0.0
    assert report.conditioned_on_node_estimates is False
    assert report.failures == []
    assert report.maps[0].state_time_totals
    assert report.maps[0].branch_histories[0].segments
    assert (
        abs(
            sum(report.maps[0].state_time_totals.values())
            - sum(history.branch_length for history in report.maps[0].branch_histories)
        )
        < 1e-9
    )
    assert reloaded.summary.replicate_count == 8
    assert reloaded.maps[0].branch_histories
    assert reloaded.maps[0].branch_histories[0].segments
    assert reloaded.maps[0].state_time_totals
    assert reloaded.summary.state_time_rows
    assert reloaded.summary.branch_occupancy_rows
    assert reloaded.summary.simulation_failure_count == 0
    assert count_report.replicate_count == 8
    assert len(count_report.matrix_rows) == 8
    assert count_report.branch_rows
    assert "transition\tmean_count\tlower_95_interval" in summary_path.read_text(
        encoding="utf-8"
    )
    assert "state\tmean_time\tlower_95_interval" in state_times_path.read_text(
        encoding="utf-8"
    )
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
    assert (
        "replicate_index\tbranch_index\tparent_node\tchild_node\tstate\tstart_time_fraction\tend_time_fraction\tduration"
        in segments_path.read_text(encoding="utf-8")
    )
    assert (
        "replicate_index\tbranch_index\tparent_node\tchild_node\tevent_index\tsource_state\ttarget_state"
        in events_path.read_text(encoding="utf-8")
    )


def test_summarize_discrete_stochastic_maps_handles_no_transition_manual_maps() -> None:
    report = manual_stochastic_map_collection(
        maps=[
            StochasticMapReplicate(
                replicate_index=0,
                root_state="0",
                total_transition_count=0,
                transition_counts={},
                state_time_totals={"0": 2.0},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B",
                        child_node="A",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B",
                                child_node="A",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B",
                        child_node="B",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B",
                                child_node="B",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            ),
            StochasticMapReplicate(
                replicate_index=1,
                root_state="0",
                total_transition_count=0,
                transition_counts={},
                state_time_totals={"0": 2.0},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B",
                        child_node="A",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B",
                                child_node="A",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B",
                        child_node="B",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B",
                                child_node="B",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            ),
        ]
    )
    summary = summarize_discrete_stochastic_maps(report)

    assert summary.mean_total_transition_count == 0.0
    assert summary.lower_95_total_transition_count == 0.0
    assert summary.upper_95_total_transition_count == 0.0
    transition_rows = {row.transition: row for row in summary.rows}
    assert transition_rows["0->1"].mean_count == 0.0
    assert transition_rows["0->1"].presence_fraction == 0.0
    count_report = count_discrete_stochastic_map_transitions(report)
    matrix_row = count_report.matrix_rows[0]
    assert matrix_row.total_transition_count == 0
    assert matrix_row.transition_counts["0->1"] == 0
    branch_transition_rows = {
        (row.parent_node, row.child_node, row.transition): row
        for row in count_report.branch_rows
    }
    assert branch_transition_rows[("A|B", "A", "0->1")].mean_count == 0.0
    assert branch_transition_rows[("A|B", "A", "0->1")].presence_fraction == 0.0
    branch_rows = {
        (row.parent_node, row.child_node, row.state): row
        for row in summary.branch_occupancy_rows
    }
    assert branch_rows[("A|B", "A", "0")].mean_time == 1.0
    assert branch_rows[("A|B", "A", "0")].mean_fraction == 1.0
    assert branch_rows[("A|B", "B", "0")].presence_fraction == 1.0


def test_summarize_discrete_stochastic_maps_tracks_multistate_branch_occupancy() -> (
    None
):
    report = manual_stochastic_map_collection(
        maps=[
            StochasticMapReplicate(
                replicate_index=0,
                root_state="0",
                total_transition_count=1,
                transition_counts={"0->1": 1},
                state_time_totals={"0": 1.0, "1": 1.0, "2": 1.0},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B|C",
                        child_node="A|B",
                        branch_length=2.0,
                        start_state="0",
                        end_state="1",
                        event_count=1,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B|C",
                                child_node="A|B",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=0.5,
                                duration=1.0,
                            ),
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B|C",
                                child_node="A|B",
                                state="1",
                                start_time_fraction=0.5,
                                end_time_fraction=1.0,
                                duration=1.0,
                            ),
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B|C",
                        child_node="C",
                        branch_length=1.0,
                        start_state="2",
                        end_state="2",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B|C",
                                child_node="C",
                                state="2",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            ),
            StochasticMapReplicate(
                replicate_index=1,
                root_state="0",
                total_transition_count=1,
                transition_counts={"0->1": 1},
                state_time_totals={"0": 0.5, "1": 1.5, "2": 1.0},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B|C",
                        child_node="A|B",
                        branch_length=2.0,
                        start_state="0",
                        end_state="1",
                        event_count=1,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B|C",
                                child_node="A|B",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=0.25,
                                duration=0.5,
                            ),
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B|C",
                                child_node="A|B",
                                state="1",
                                start_time_fraction=0.25,
                                end_time_fraction=1.0,
                                duration=1.5,
                            ),
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B|C",
                        child_node="C",
                        branch_length=1.0,
                        start_state="2",
                        end_state="2",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B|C",
                                child_node="C",
                                state="2",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            ),
        ]
    )
    summary = summarize_discrete_stochastic_maps(report)

    transition_rows = {row.transition: row for row in summary.rows}
    assert transition_rows["0->1"].mean_count == 1.0
    state_rows = {row.state: row for row in summary.state_time_rows}
    assert state_rows["1"].mean_time == 1.25
    branch_rows = {
        (row.parent_node, row.child_node, row.state): row
        for row in summary.branch_occupancy_rows
    }
    count_report = count_discrete_stochastic_map_transitions(report)
    count_matrix_rows = {row.replicate_index: row for row in count_report.matrix_rows}
    assert count_matrix_rows[0].transition_counts["0->1"] == 1
    branch_transition_rows = {
        (row.parent_node, row.child_node, row.transition): row
        for row in count_report.branch_rows
    }
    assert branch_transition_rows[("A|B|C", "A|B", "0->1")].mean_count == 1.0
    assert branch_transition_rows[("A|B|C", "C", "1->2")].mean_count == 0.0
    assert branch_rows[("A|B|C", "A|B", "1")].mean_time == 1.25
    assert branch_rows[("A|B|C", "A|B", "1")].mean_fraction == 0.625
    assert branch_rows[("A|B|C", "C", "2")].mean_fraction == 1.0


def test_summarize_discrete_stochastic_map_density_tracks_binary_branch_slices(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "binary-density-tree.nwk"
    tree_path.write_text("(A:1,B:1);\n", encoding="utf-8")
    report = manual_stochastic_map_collection(
        tree_path=tree_path,
        maps=[
            StochasticMapReplicate(
                replicate_index=0,
                root_state="0",
                total_transition_count=1,
                transition_counts={"0->1": 1},
                state_time_totals={"0": 1.5, "1": 0.5},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B",
                        child_node="A",
                        branch_length=1.0,
                        start_state="0",
                        end_state="1",
                        event_count=1,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B",
                                child_node="A",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=0.5,
                                duration=0.5,
                            ),
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B",
                                child_node="A",
                                state="1",
                                start_time_fraction=0.5,
                                end_time_fraction=1.0,
                                duration=0.5,
                            ),
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B",
                        child_node="B",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B",
                                child_node="B",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            ),
            StochasticMapReplicate(
                replicate_index=1,
                root_state="0",
                total_transition_count=0,
                transition_counts={},
                state_time_totals={"0": 1.0, "1": 1.0},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B",
                        child_node="A",
                        branch_length=1.0,
                        start_state="1",
                        end_state="1",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B",
                                child_node="A",
                                state="1",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B",
                        child_node="B",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B",
                                child_node="B",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            ),
        ],
    )

    density_report = summarize_discrete_stochastic_map_density(report, resolution=2)

    assert density_report.focal_state == "1"
    assert density_report.baseline_state == "0"
    branch_state_rows = {
        (row.parent_node, row.child_node, row.state): row
        for row in density_report.branch_state_rows
    }
    assert branch_state_rows[("A|B", "A", "1")].mean_probability == 0.75
    assert branch_state_rows[("A|B", "B", "1")].mean_probability == 0.0
    branch_rows = {
        (row.parent_node, row.child_node): row for row in density_report.branch_rows
    }
    assert branch_rows[("A|B", "A")].mean_posterior_probability == 0.75
    assert branch_rows[("A|B", "A")].uncertainty == 0.5
    assert branch_rows[("A|B", "B")].mean_posterior_probability == 0.0
    slice_rows = [
        row
        for row in density_report.density_rows
        if row.parent_node == "A|B" and row.child_node == "A"
    ]
    assert [row.posterior_probability for row in slice_rows] == [0.5, 1.0]
    assert [row.posterior_uncertainty for row in slice_rows] == [1.0, 0.0]

    branch_probability_path = tmp_path / "branch-probabilities.tsv"
    density_branch_path = tmp_path / "density-branches.tsv"
    density_slice_path = tmp_path / "density-slices.tsv"
    density_html_path = tmp_path / "density-report.html"
    write_stochastic_map_branch_probability_table(
        branch_probability_path, density_report
    )
    write_stochastic_map_density_branch_table(density_branch_path, density_report)
    write_stochastic_map_density_slice_table(density_slice_path, density_report)
    render_result = render_stochastic_map_density_artifact(
        density_report,
        tree_path=tree_path,
        out_path=density_html_path,
        layout="phylogram",
    )
    assert "mean_probability" in branch_probability_path.read_text(encoding="utf-8")
    assert "mean_posterior_probability" in density_branch_path.read_text(
        encoding="utf-8"
    )
    assert "posterior_probability" in density_slice_path.read_text(encoding="utf-8")
    assert density_html_path.exists()
    assert density_html_path.with_suffix(".svg").exists()
    assert render_result.rendered_branch_color_count == 2


def test_summarize_discrete_stochastic_map_density_keeps_multistate_branch_probabilities(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "multistate-density-tree.nwk"
    tree_path.write_text("(A:1,B:1,C:1);\n", encoding="utf-8")
    report = manual_stochastic_map_collection(
        tree_path=tree_path,
        maps=[
            StochasticMapReplicate(
                replicate_index=0,
                root_state="0",
                total_transition_count=0,
                transition_counts={},
                state_time_totals={"0": 1.0, "1": 1.0, "2": 1.0},
                branch_histories=[
                    StochasticMapBranchHistory(
                        branch_index=0,
                        parent_node="A|B|C",
                        child_node="A",
                        branch_length=1.0,
                        start_state="0",
                        end_state="0",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=0,
                                parent_node="A|B|C",
                                child_node="A",
                                state="0",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=1,
                        parent_node="A|B|C",
                        child_node="B",
                        branch_length=1.0,
                        start_state="1",
                        end_state="1",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=1,
                                parent_node="A|B|C",
                                child_node="B",
                                state="1",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                    StochasticMapBranchHistory(
                        branch_index=2,
                        parent_node="A|B|C",
                        child_node="C",
                        branch_length=1.0,
                        start_state="2",
                        end_state="2",
                        event_count=0,
                        events=[],
                        segments=[
                            StochasticMapStateSegment(
                                branch_index=2,
                                parent_node="A|B|C",
                                child_node="C",
                                state="2",
                                start_time_fraction=0.0,
                                end_time_fraction=1.0,
                                duration=1.0,
                            )
                        ],
                    ),
                ],
            )
        ],
    )

    density_report = summarize_discrete_stochastic_map_density(report, resolution=4)

    assert density_report.focal_state is None
    assert density_report.density_rows == []
    assert density_report.branch_rows == []
    assert len(density_report.branch_state_rows) == 9
    assert any(
        "density slices require one explicit focal state" in warning
        for warning in density_report.warnings
    )


@pytest.mark.slow
def test_simulate_discrete_stochastic_maps_from_fit_report_matches_path_surface() -> (
    None
):
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_multistate_twenty_four_taxa"
    )
    fit_report = fit_discrete_mk_model(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
    )

    direct_report = simulate_discrete_stochastic_maps_from_fit_report(
        fit_report,
        replicates=8,
        seed=17,
    )
    path_report = simulate_discrete_stochastic_maps(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="symmetric",
        replicates=8,
        seed=17,
    )

    assert direct_report.summary == path_report.summary
    assert (
        direct_report.maps[0].transition_counts == path_report.maps[0].transition_counts
    )
    assert direct_report.fit_audit.parameter_count == 3
    assert direct_report.fit_audit.optimizer_converged is True
    assert direct_report.fit_audit.preferred_model_by_aic == "equal-rates"


@pytest.mark.slow
def test_simulate_discrete_stochastic_maps_reports_ard_fit_instability_honestly() -> (
    None
):
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_multistate_twenty_four_taxa"
    )
    report = simulate_discrete_stochastic_maps(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
        replicates=8,
        seed=17,
    )

    transition_rows = {row.transition: row for row in report.summary.rows}

    assert report.summary.replicate_count == 8
    assert report.summary.simulation_failure_count == 0
    assert report.fit_audit.parameter_count == 12
    assert report.fit_audit.overparameterized is False
    assert report.fit_audit.optimizer_converged is False
    assert report.fit_audit.optimizer_hit_lower_parameter_bound is True
    assert report.fit_audit.preferred_model_by_aic == "equal-rates"
    assert any("weakly identified" in warning for warning in report.fit_audit.warnings)
    assert any("weakly identified" in warning for warning in report.warnings)
    assert "north->south" in transition_rows
    assert "west->south" in transition_rows


def test_simulate_discrete_stochastic_maps_keeps_ard_binary_transition_directions_distinct() -> (
    None
):
    fixture_entry = get_shared_phytools_comparative_fixture(
        "phytools_discrete_ard_binary_twenty_four_taxa"
    )
    report = simulate_discrete_stochastic_maps(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        model="all-rates-different",
        replicates=16,
        seed=17,
    )

    transition_rows = {row.transition: row for row in report.summary.rows}

    assert "0->1" in transition_rows
    assert "1->0" in transition_rows
    assert transition_rows["0->1"].mean_count != transition_rows["1->0"].mean_count
