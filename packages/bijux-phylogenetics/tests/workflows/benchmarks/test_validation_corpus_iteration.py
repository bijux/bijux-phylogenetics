from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.benchmark import (
    benchmark_alignment_site_scaling,
    benchmark_tree_set_consensus,
)
import bijux_phylogenetics.validation as validation_api
from bijux_phylogenetics.validation import (
    build_broken_benchmark_corpus,
    build_clean_benchmark_corpus,
    build_gamma_invariant_family_coverage,
    build_invariant_mixture_family_coverage,
    build_large_alignment_scaling_benchmark_dashboard,
    build_large_tree_scaling_benchmark_dashboard,
    build_large_tree_set_scaling_benchmark_dashboard,
    build_memory_benchmark_dashboard,
    build_messy_benchmark_corpus,
    build_method_accuracy_dashboard,
    build_method_limitation_registry,
    build_regression_dataset_corpus,
    build_runtime_benchmark_dashboard,
    build_scientific_validation_report,
    build_shared_continuous_fixture_catalog,
    build_shared_discrete_fixture_catalog,
    build_substitution_parameter_prior_coverage,
    build_workflow_practical_limit_dashboard,
    validate_simulation_reproducibility,
    write_validation_corpus_json,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_build_clean_benchmark_corpus_keeps_core_dataset_ready() -> None:
    report = build_clean_benchmark_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 242
    assert report.passed is True
    assert report.passed_case_count == 1
    case = report.cases[0]
    assert case.readiness_decision == "ready_with_warnings"
    assert "comparative" in case.allowed_analyses
    assert "publication" in case.allowed_analyses
    assert case.blockers == []


def test_build_broken_benchmark_corpus_preserves_expected_failure_signatures() -> None:
    report = build_broken_benchmark_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 243
    assert report.passed is True
    observed = {case.name: case for case in report.cases}
    assert observed["duplicate_tip_tree"].observed_code == "duplicate_taxon_error"
    assert (
        observed["invalid_alignment_lengths"].observed_code == "invalid_alignment_error"
    )
    assert observed["dataset_missing_metadata_taxon"].readiness_decision == "blocked"
    assert (
        "metadata table is missing one or more tree taxa"
        in observed["dataset_missing_metadata_taxon"].blockers
    )


def test_build_messy_benchmark_corpus_captures_multi_surface_warning_cases() -> None:
    report = build_messy_benchmark_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 244
    assert report.passed is True
    observed = {case.name: case for case in report.cases}
    first = observed["reordered_alignment_extra_taxa_invalid_dates_and_calibrations"]
    assert (
        "calibration table contains invalid fossil calibration targets or ages"
        in first.blockers
    )
    assert (
        "one or more dataset surfaces silently reorder shared taxa relative to the tree"
        in first.warnings
    )
    second = observed["low_information_alignment_with_trait_mismatch"]
    assert (
        "alignment is not currently safe for core inference workflows"
        in second.blockers
    )
    assert "alignment contains near-duplicate sequences" in second.warnings


def test_build_regression_dataset_corpus_matches_checked_in_summaries() -> None:
    report = build_regression_dataset_corpus(fixtures_root=FIXTURES)

    assert report.goal_id == 245
    assert report.passed is True
    assert report.passed_case_count == 2
    observed = {case.name: case for case in report.cases}
    assert observed["core_inference_ready_dataset"].observed["risk_level"] == "low"
    assert observed["warning_rich_dataset"].observed["warning_count"] == 12


def test_benchmark_alignment_site_scaling_reports_site_axis_observations() -> None:
    report = benchmark_alignment_site_scaling(
        replicates=1, site_counts=[24, 48], sequence_count=4
    )

    assert report.sequence_count == 4
    assert [row.item_count for row in report.observations] == [24, 48]


def test_benchmark_tree_set_consensus_reports_tree_count_scaling() -> None:
    report = benchmark_tree_set_consensus(replicates=1, tree_counts=[4, 8], tip_count=4)

    assert report.tip_count == 4
    assert [row.item_count for row in report.observations] == [4, 8]


def test_build_large_tree_scaling_benchmark_dashboard_tracks_goal_221() -> None:
    report = build_large_tree_scaling_benchmark_dashboard(
        replicates=1,
        tip_counts=[8, 16],
    )

    assert report.goal_id == 221
    assert {workflow.workflow for workflow in report.workflows} == {
        "tree-validation",
        "tree-comparison",
        "tree-rendering",
        "tree-reporting",
    }
    assert all(
        [row.item_count for row in workflow.observations] == [8, 16]
        for workflow in report.workflows
    )


@pytest.mark.slow
def test_build_large_alignment_scaling_benchmark_dashboard_tracks_goal_222() -> None:
    report = build_large_alignment_scaling_benchmark_dashboard(
        replicates=1,
        size_classes=[
            ("sequences-4-sites-16", 4, 16),
            ("sequences-6-sites-24", 6, 24),
        ],
    )

    assert report.goal_id == 222
    assert {workflow.workflow for workflow in report.workflows} == {
        "alignment-diagnostics",
        "alignment-trimming",
        "distance-analysis",
        "alignment-readiness",
    }
    assert all(
        [row.sequence_count for row in workflow.observations] == [4, 6]
        for workflow in report.workflows
    )


def test_build_large_tree_set_scaling_benchmark_dashboard_tracks_goal_223() -> None:
    report = build_large_tree_set_scaling_benchmark_dashboard(
        replicates=1,
        size_classes=[
            ("trees-8-taxa-6", 8, 6),
            ("trees-12-taxa-8", 12, 8),
        ],
    )

    assert report.goal_id == 223
    assert {workflow.workflow for workflow in report.workflows} == {
        "tree-set-consensus",
        "pairwise-rf-diversity",
        "topology-clustering",
        "uncertainty-summaries",
    }
    assert all(
        [row.tree_count for row in workflow.observations] == [8, 12]
        for workflow in report.workflows
    )


@pytest.mark.slow
def test_build_workflow_practical_limit_dashboard_tracks_goal_224() -> None:
    report = build_workflow_practical_limit_dashboard(
        replicates=1,
        tree_tip_counts=[8, 16],
        alignment_size_classes=[
            ("sequences-4-sites-16", 4, 16),
            ("sequences-6-sites-24", 6, 24),
        ],
        tree_set_size_classes=[
            ("trees-8-taxa-6", 8, 6),
            ("trees-12-taxa-8", 12, 8),
        ],
        stress_tiers=["small"],
    )

    assert report.goal_id == 224
    assert any(entry.workflow == "tree-validation" for entry in report.entries)
    assert any(
        entry.workflow == "posterior-tree-set-consensus" for entry in report.entries
    )
    assert (
        max(
            entry.tested_taxon_limit
            for entry in report.entries
            if entry.tested_taxon_limit is not None
        )
        == 256
    )


@pytest.mark.slow
def test_build_method_accuracy_dashboard_summarizes_fixture_and_corpus_pass_rates() -> (
    None
):
    report = build_method_accuracy_dashboard(fixtures_root=FIXTURES)

    assert report.goal_id == 246
    surfaces = {row.surface: row for row in report.rows}
    assert surfaces["level1-reference-validation"].coverage_count > 0
    assert surfaces["clean-benchmark-corpus"].accuracy == 1.0
    assert surfaces["regression-dataset-corpus"].failed_count == 0


@pytest.mark.slow
def test_build_runtime_and_memory_dashboards_cover_sites_and_posterior_samples() -> (
    None
):
    runtime = build_runtime_benchmark_dashboard(replicates=1)
    memory = build_memory_benchmark_dashboard(replicates=1)

    assert runtime.goal_id == 247
    assert memory.goal_id == 248
    assert {row.scaling_axis for row in runtime.rows} >= {
        "sites",
        "posterior_samples",
        "taxa",
    }
    assert {row.scaling_axis for row in memory.rows} >= {
        "sites",
        "posterior_samples",
        "taxa",
    }


def test_build_method_limitation_registry_marks_experimental_and_validated_surfaces() -> (
    None
):
    report = build_method_limitation_registry()

    assert report.goal_id == 250
    statuses = {entry.method: entry.status for entry in report.entries}
    assert statuses["tree-validation"] == "validated"
    assert statuses["bayesian-time-tree"] == "experimental"


@pytest.mark.slow
def test_build_scientific_validation_report_separates_claim_statuses() -> None:
    report = build_scientific_validation_report(fixtures_root=FIXTURES)

    assert report.goal_id == 249
    statuses = {claim.status for claim in report.claims}
    assert statuses >= {"validated", "experimental", "unvalidated", "unsafe"}
    unvalidated_claim = next(
        claim for claim in report.claims if claim.status == "unvalidated"
    )
    assert "goal-252 remains outside this iteration" not in unvalidated_claim.evidence
    assert any(
        evidence.startswith("goal-252 shared continuous fixture catalog")
        for evidence in unvalidated_claim.evidence
    )


def test_validate_simulation_reproducibility_confirms_same_seed_repeatability() -> None:
    report = validate_simulation_reproducibility(fixtures_root=FIXTURES)

    assert report.goal_id == 251
    assert report.passed is True
    assert len(report.cases) == 6
    assert all(case.digest for case in report.cases)


def test_build_invariant_mixture_family_coverage_proves_supported_families() -> None:
    report = build_invariant_mixture_family_coverage()

    assert report.goal_id == 29
    assert report.passed is True
    assert report.family_count == 7
    assert report.model_count == 17
    assert report.passed_family_count == 7
    assert report.passed_model_count == 17
    observed = {entry.family_name: entry for entry in report.families}
    assert set(observed) == {
        "nucleotide",
        "protein-poisson",
        "protein-empirical",
        "codon-ctmc",
        "discrete-mk-unordered",
        "discrete-mk-ordered",
        "discrete-mk-irreversible",
    }

    nucleotide = observed["nucleotide"]
    assert nucleotide.supported_model_names == ["JC69", "K80", "F81", "HKY85", "GTR"]
    assert nucleotide.covered_model_names == nucleotide.supported_model_names
    assert nucleotide.passed_model_count == 5

    gtr = next(case for case in nucleotide.model_cases if case.model_name == "GTR")
    assert gtr.surface_model_name == "GTR"
    assert gtr.invariant_proportion == 0.35
    assert gtr.likelihood_changed is True
    assert gtr.invariant_component_active is True
    assert gtr.passed is True
    assert gtr.invariant_component_likelihood > 0.0
    assert gtr.variable_component_likelihood > 0.0


def test_build_gamma_invariant_family_coverage_proves_supported_families() -> None:
    report = build_gamma_invariant_family_coverage()

    assert report.goal_id == 30
    assert report.passed is True
    assert report.family_count == 7
    assert report.model_count == 17
    assert report.passed_family_count == 7
    assert report.passed_model_count == 17
    observed = {entry.family_name: entry for entry in report.families}
    assert set(observed) == {
        "nucleotide",
        "protein-poisson",
        "protein-empirical",
        "codon-ctmc",
        "discrete-mk-unordered",
        "discrete-mk-ordered",
        "discrete-mk-irreversible",
    }

    empirical = observed["protein-empirical"]
    assert empirical.supported_model_names == ["JTT", "WAG", "LG"]
    assert empirical.covered_model_names == empirical.supported_model_names
    assert empirical.passed_model_count == 3

    wag = next(case for case in empirical.model_cases if case.model_name == "WAG")
    assert wag.alpha == 0.8
    assert wag.category_count == 4
    assert wag.invariant_proportion == 0.35
    assert wag.category_likelihood_count == 4
    assert wag.invariant_component_active is True
    assert wag.variable_component_active is True
    assert wag.boundary_warnings == ["invariant proportion hit lower search boundary"]
    assert wag.passed is True


def test_build_shared_continuous_fixture_catalog_tracks_live_consumed_fixtures() -> (
    None
):
    report = build_shared_continuous_fixture_catalog()

    assert report.goal_id == 252
    assert report.fixture_count == 6
    assert report.live_case_count >= report.fixture_count
    observed = {fixture.fixture_id: fixture for fixture in report.fixtures}
    missing_values = observed["geiger_continuous_missing_values_twenty_four_taxa"]

    assert "trend" in missing_values.supported_model_names
    assert "trend" not in missing_values.covered_model_names
    assert missing_values.covered_model_names == [
        "BM",
        "OU",
        "delta",
        "kappa",
        "lambda",
        "white",
    ]
    assert all(fixture.live_case_count > 0 for fixture in report.fixtures)


def test_build_shared_discrete_fixture_catalog_tracks_live_consumed_fixtures() -> None:
    report = build_shared_discrete_fixture_catalog()

    assert report.goal_id == 253
    assert report.fixture_count == 17
    assert report.live_case_count >= report.fixture_count
    observed = {fixture.fixture_id: fixture for fixture in report.fixtures}
    assert "geiger_discrete_constant_negative_twenty_four_taxa" not in observed
    sym_fixture = observed["geiger_discrete_sym_three_state_twenty_four_taxa"]

    assert sym_fixture.covered_model_names == ["SYM"]
    assert sym_fixture.trait_kind == "three-state-discrete"
    assert sym_fixture.transition_matrix_metadata == {
        "states": ["north", "south", "central"],
        "matrix_shape": "3x3",
        "parameterization": "symmetric",
        "pairwise_rate": 0.45,
    }
    assert all(fixture.live_case_count > 0 for fixture in report.fixtures)


def test_build_substitution_parameter_prior_coverage_reports_supported_targets() -> (
    None
):
    report = build_substitution_parameter_prior_coverage()

    assert report.goal_id == 76
    assert report.passed is True
    assert report.target_count == 7
    assert report.passed_target_count == 7
    observed = {entry.target_name: entry for entry in report.entries}
    assert set(observed) == {
        "kappa",
        "exchangeabilities",
        "base-frequencies",
        "gamma-alpha",
        "invariant-proportion",
        "codon-frequencies",
        "mk-rates",
    }

    codon_frequencies = observed["codon-frequencies"]
    assert codon_frequencies.surface_name == "substitution-parameter-bundle"
    assert codon_frequencies.prior_family == "dirichlet"
    assert codon_frequencies.contribution_count == 1
    assert len(codon_frequencies.covered_names) == 61
    assert codon_frequencies.changed is True
    assert codon_frequencies.passed is True

    mk_rates = observed["mk-rates"]
    assert mk_rates.surface_name == "discrete-trait-rate-prior"
    assert mk_rates.prior_family == "exponential"
    assert mk_rates.covered_names == ["0->1", "1->0"]
    assert mk_rates.changed is True
    assert mk_rates.passed is True


def test_package_root_exports_validation_corpus_surfaces() -> None:
    assert validation_api.build_clean_benchmark_corpus is build_clean_benchmark_corpus
    assert (
        validation_api.build_large_alignment_scaling_benchmark_dashboard
        is build_large_alignment_scaling_benchmark_dashboard
    )
    assert (
        validation_api.build_large_tree_set_scaling_benchmark_dashboard
        is build_large_tree_set_scaling_benchmark_dashboard
    )
    assert (
        validation_api.build_workflow_practical_limit_dashboard
        is build_workflow_practical_limit_dashboard
    )
    assert (
        validation_api.build_large_tree_scaling_benchmark_dashboard
        is build_large_tree_scaling_benchmark_dashboard
    )
    assert (
        validation_api.build_method_limitation_registry
        is build_method_limitation_registry
    )
    assert (
        validation_api.build_shared_continuous_fixture_catalog
        is build_shared_continuous_fixture_catalog
    )
    assert (
        validation_api.build_gamma_invariant_family_coverage
        is build_gamma_invariant_family_coverage
    )
    assert (
        validation_api.build_invariant_mixture_family_coverage
        is build_invariant_mixture_family_coverage
    )
    assert (
        validation_api.build_shared_discrete_fixture_catalog
        is build_shared_discrete_fixture_catalog
    )
    assert (
        validation_api.build_substitution_parameter_prior_coverage
        is build_substitution_parameter_prior_coverage
    )
    assert (
        validation_api.validate_simulation_reproducibility
        is validate_simulation_reproducibility
    )


def test_write_validation_corpus_json_serializes_report_payloads(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "validation.json",
        build_clean_benchmark_corpus(fixtures_root=FIXTURES),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 242' in text
    assert '"corpus": "clean-benchmark-corpus"' in text


def test_write_validation_corpus_json_serializes_large_tree_scaling_dashboard(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "large-tree-scaling.json",
        build_large_tree_scaling_benchmark_dashboard(
            replicates=1,
            tip_counts=[8],
        ),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 221' in text
    assert '"tree-rendering"' in text


def test_write_validation_corpus_json_serializes_large_alignment_scaling_dashboard(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "large-alignment-scaling.json",
        build_large_alignment_scaling_benchmark_dashboard(
            replicates=1,
            size_classes=[("sequences-4-sites-16", 4, 16)],
        ),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 222' in text
    assert '"alignment-trimming"' in text


def test_write_validation_corpus_json_serializes_large_tree_set_scaling_dashboard(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "large-tree-set-scaling.json",
        build_large_tree_set_scaling_benchmark_dashboard(
            replicates=1,
            size_classes=[("trees-8-taxa-6", 8, 6)],
        ),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 223' in text
    assert '"pairwise-rf-diversity"' in text


def test_write_validation_corpus_json_serializes_continuous_fixture_catalog(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "continuous-fixture-catalog.json",
        build_shared_continuous_fixture_catalog(),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 252' in text
    assert '"geiger_continuous_missing_values_twenty_four_taxa"' in text


def test_write_validation_corpus_json_serializes_invariant_family_coverage(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "invariant-mixture-family-coverage.json",
        build_invariant_mixture_family_coverage(),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 29' in text
    assert '"protein-empirical"' in text
    assert '"LG"' in text


def test_write_validation_corpus_json_serializes_gamma_invariant_family_coverage(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "gamma-invariant-family-coverage.json",
        build_gamma_invariant_family_coverage(),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 30' in text
    assert '"protein-empirical"' in text
    assert '"boundary_warnings"' in text


def test_write_validation_corpus_json_serializes_discrete_fixture_catalog(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "discrete-fixture-catalog.json",
        build_shared_discrete_fixture_catalog(),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 253' in text
    assert '"geiger_discrete_sym_three_state_twenty_four_taxa"' in text


def test_write_validation_corpus_json_serializes_substitution_prior_coverage(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "substitution-prior-coverage.json",
        build_substitution_parameter_prior_coverage(),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 76' in text
    assert '"codon-frequencies"' in text
    assert '"mk-rates"' in text


@pytest.mark.slow
def test_write_validation_corpus_json_serializes_workflow_practical_limit_dashboard(
    tmp_path: Path,
) -> None:
    path = write_validation_corpus_json(
        tmp_path / "workflow-practical-limits.json",
        build_workflow_practical_limit_dashboard(
            replicates=1,
            tree_tip_counts=[8],
            alignment_size_classes=[("sequences-4-sites-16", 4, 16)],
            tree_set_size_classes=[("trees-8-taxa-6", 8, 6)],
            stress_tiers=["small"],
        ),
    )

    text = path.read_text(encoding="utf-8")
    assert '"goal_id": 224' in text
    assert '"posterior-tree-set-consensus"' in text
