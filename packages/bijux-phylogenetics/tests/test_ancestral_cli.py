from __future__ import annotations

import json
from pathlib import Path
import shutil

import pytest

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


def test_ancestral_continuous_cli_can_export_table(tmp_path: Path, capsys) -> None:
    table_path = tmp_path / "ancestral.tsv"
    summary_path = tmp_path / "ancestral-summary.tsv"
    uncertainty_path = tmp_path / "ancestral-uncertainty.tsv"
    exclusions_path = tmp_path / "ancestral-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "continuous",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--table-out",
            str(table_path),
            "--summary-out",
            str(summary_path),
            "--uncertainty-out",
            str(uncertainty_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["estimator"] == "ace-pic"
    assert payload["metrics"]["internal_node_count"] == 3
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["unstable_node_count"] >= 0
    assert payload["metrics"]["tree_is_ultrametric"] is True
    assert payload["metrics"]["covariance_near_singular"] is False
    assert payload["metrics"]["covariance_condition_number"] > 0.0
    assert payload["metrics"]["log_likelihood"] is not None
    assert payload["metrics"]["residual_sigma_squared"] > 0.0
    assert payload["data"]["brownian_fit_diagnostics"]["covariance_model"] == (
        "brownian-shared-path"
    )
    assert "estimate\tstandard_error" in table_path.read_text(encoding="utf-8")
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\testimator\talpha"
    )
    assert uncertainty_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\testimate\tstandard_error"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_continuous_cli_supports_fast_anc_estimator(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-summary.tsv"
    uncertainty_path = tmp_path / "ancestral-uncertainty.tsv"

    exit_code = main(
        [
            "ancestral",
            "continuous",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--model",
            "brownian",
            "--estimator",
            "fast-anc",
            "--summary-out",
            str(summary_path),
            "--uncertainty-out",
            str(uncertainty_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["estimator"] == "fast-anc"
    assert "\tfast-anc\t" in summary_path.read_text(encoding="utf-8")
    assert "A|B\t" in uncertainty_path.read_text(encoding="utf-8")


def test_ancestral_continuous_cli_supports_anc_ml_estimator(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-summary.tsv"
    uncertainty_path = tmp_path / "ancestral-uncertainty.tsv"

    exit_code = main(
        [
            "ancestral",
            "continuous",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--model",
            "brownian",
            "--estimator",
            "anc-ml",
            "--summary-out",
            str(summary_path),
            "--uncertainty-out",
            str(uncertainty_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["estimator"] == "anc-ml"
    assert payload["metrics"]["optimizer_name"] == "closed-form-profile-solution"
    assert payload["metrics"]["optimizer_converged"] is True
    assert payload["metrics"]["optimizer_iteration_count"] == 0
    assert payload["metrics"]["optimizer_function_evaluation_count"] == 1
    assert "\tanc-ml\t" in summary_path.read_text(encoding="utf-8")
    assert "A|B|C|D\t" in uncertainty_path.read_text(encoding="utf-8")


def test_ancestral_discrete_cli_reports_sparse_state_warning(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert (
        "one or more discrete states are represented by fewer than two taxa and should be interpreted cautiously"
        in payload["warnings"]
    )
    assert (
        "one or more discrete ancestral nodes remain unstable across candidate states"
        in payload["warnings"]
    )


def test_ancestral_discrete_cli_can_export_probability_review(
    tmp_path: Path, capsys
) -> None:
    table_path = tmp_path / "ancestral-discrete.tsv"
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    probabilities_path = tmp_path / "ancestral-discrete-probabilities.tsv"
    transitions_path = tmp_path / "ancestral-discrete-transitions.tsv"
    exclusions_path = tmp_path / "ancestral-discrete-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "equal-rates",
            "--table-out",
            str(table_path),
            "--summary-out",
            str(summary_path),
            "--probabilities-out",
            str(probabilities_path),
            "--transitions-out",
            str(transitions_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["internal_node_count"] == 3
    assert payload["metrics"]["excluded_taxon_count"] == 0
    assert payload["metrics"]["unstable_node_count"] >= 0
    assert payload["metrics"]["log_likelihood"] is not None
    assert payload["metrics"]["parameter_count"] == 1
    assert payload["metrics"]["aic"] is not None
    assert payload["metrics"]["root_prior_mode"] == "equal"
    assert payload["metrics"]["phytools_rerooting_method_comparable"] is True
    assert payload["metrics"]["transition_rate_count"] == 6
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert probabilities_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_state\tstate_set"
    )
    assert transitions_path.read_text(encoding="utf-8").startswith(
        "source_state\ttarget_state\ttransition_allowed\tstep_distance\trate"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"
    assert "most_likely_state\tstate_set" in table_path.read_text(encoding="utf-8")


def test_ancestral_discrete_cli_supports_fixed_root_prior_policy(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    transitions_path = tmp_path / "ancestral-discrete-transitions.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "equal-rates",
            "--root-prior-mode",
            "fixed",
            "--fixed-root-state",
            "north",
            "--summary-out",
            str(summary_path),
            "--transitions-out",
            str(transitions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["root_prior_mode"] == "fixed"
    assert payload["metrics"]["fixed_root_state"] == "north"
    assert payload["metrics"]["phytools_rerooting_method_comparable"] is False
    assert payload["metrics"]["transition_rate_count"] == 6
    assert (
        "phytools::rerootingMethod inherits fitMk's default equal root prior; empirical or fixed root-prior runs remain Bijux sensitivity scenarios without direct rerootingMethod parity"
        in payload["warnings"]
    )
    assert "\tfixed\tnorth\t" in summary_path.read_text(encoding="utf-8")
    assert transitions_path.read_text(encoding="utf-8").startswith(
        "source_state\ttarget_state\ttransition_allowed\tstep_distance\trate"
    )


def test_ancestral_discrete_cli_reports_sym_fit_diagnostics_and_er_comparison(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    fit_path = tmp_path / "ancestral-discrete-fit.tsv"
    comparison_path = tmp_path / "ancestral-discrete-comparison.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_geography_biased.tsv")),
            "--trait",
            "region",
            "--model",
            "symmetric",
            "--compare-model",
            "equal-rates",
            "--summary-out",
            str(summary_path),
            "--fit-out",
            str(fit_path),
            "--comparison-out",
            str(comparison_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "symmetric"
    assert payload["metrics"]["optimizer_converged"] is True
    assert payload["metrics"]["optimizer_iteration_count"] > 0
    assert payload["metrics"]["optimizer_function_evaluation_count"] > 0
    assert payload["metrics"]["overparameterized"] is False
    assert payload["metrics"]["baseline_model"] == "equal-rates"
    assert payload["metrics"]["baseline_delta_aic"] > 0.0
    assert payload["metrics"]["preferred_model_by_aic"] == "equal-rates"
    assert payload["metrics"]["comparison_node_count"] == 5
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert fit_path.read_text(encoding="utf-8").startswith(
        "model\ttaxon_count\tstate_count\tparameter_count\tlog_likelihood"
    )
    assert comparison_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tleft_model\tright_model"
    )


@pytest.mark.slow
def test_ancestral_discrete_cli_reports_ard_fit_diagnostics_and_weak_fit_warning(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    fit_path = tmp_path / "ancestral-discrete-fit.tsv"
    comparison_path = tmp_path / "ancestral-discrete-comparison.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree_ladderized.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "all-rates-different",
            "--compare-model",
            "equal-rates",
            "--summary-out",
            str(summary_path),
            "--fit-out",
            str(fit_path),
            "--comparison-out",
            str(comparison_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "all-rates-different"
    assert payload["metrics"]["phytools_rerooting_method_comparable"] is False
    assert payload["metrics"]["optimizer_iteration_count"] > 0
    assert payload["metrics"]["overparameterized"] is True
    assert payload["metrics"]["baseline_model"] == "equal-rates"
    assert payload["metrics"]["preferred_model_by_aic"] == "equal-rates"
    assert (
        "one or more discrete rate parameters hit an optimizer bound and should be interpreted as weakly identified"
        in payload["warnings"]
    )
    assert (
        "phytools::rerootingMethod is invalid for non-symmetric Q matrices such as all-rates-different models in phytools 2.5.2"
        in payload["warnings"]
    )
    assert fit_path.read_text(encoding="utf-8").startswith(
        "model\ttaxon_count\tstate_count\tparameter_count\tlog_likelihood"
    )
    assert comparison_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tleft_model\tright_model"
    )


def test_ancestral_discrete_cli_can_export_parsimony_comparison(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-discrete-summary.tsv"
    comparison_path = tmp_path / "ancestral-discrete-comparison.tsv"
    exit_code = main(
        [
            "ancestral",
            "discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--model",
            "fitch",
            "--compare-model",
            "equal-rates",
            "--summary-out",
            str(summary_path),
            "--comparison-out",
            str(comparison_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "fitch"
    assert payload["metrics"]["minimal_change_count"] == 1
    assert payload["metrics"]["ambiguous_internal_node_count"] == 1
    assert payload["metrics"]["comparison_node_count"] == 3
    assert payload["metrics"]["comparison_differing_node_count"] >= 0
    assert "minimal_change_count" in summary_path.read_text(encoding="utf-8")
    assert comparison_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tleft_model\tright_model"
    )


@pytest.mark.slow
def test_ancestral_discrete_reference_cli_reports_passing_cases(capsys) -> None:
    exit_code = main(["ancestral", "discrete-reference", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["case_count"] == 10
    assert payload["metrics"]["external_case_count"] == 6
    assert payload["metrics"]["all_passed"] is True


def test_ancestral_confidence_cli_can_export_single_tree_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-confidence-summary.tsv"
    confidence_path = tmp_path / "ancestral-confidence.tsv"
    exit_code = main(
        [
            "ancestral",
            "confidence",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--summary-out",
            str(summary_path),
            "--confidence-out",
            str(confidence_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["kind"] == "discrete"
    assert payload["metrics"]["source_kind"] == "tree"
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["confidence_row_count"] == 3
    assert payload["metrics"]["high_entropy_count"] == 2
    assert payload["metrics"]["top_uncertain_id"] == "C|D"
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tsource_kind\treconstruction_kind\ttarget_kind"
    )
    assert confidence_path.read_text(encoding="utf-8").startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_state\tstate_set"
    )


def test_ancestral_root_sensitivity_cli_can_export_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-root-sensitivity-summary.tsv"
    assumptions_path = tmp_path / "ancestral-root-sensitivity-assumptions.tsv"
    nodes_path = tmp_path / "ancestral-root-sensitivity-nodes.tsv"
    exit_code = main(
        [
            "ancestral",
            "root-sensitivity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "equal-rates",
            "--fixed-root-state",
            "island",
            "--summary-out",
            str(summary_path),
            "--assumptions-out",
            str(assumptions_path),
            "--nodes-out",
            str(nodes_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["assumption_count"] == 3
    assert payload["metrics"]["compared_node_count"] == 3
    assert payload["metrics"]["state_changed_node_count"] == 2
    assert payload["metrics"]["support_changed_node_count"] == 1
    assert payload["metrics"]["top_sensitive_node"] == "A|B|C|D"
    assert payload["metrics"]["fixed_root_state"] == "island"
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering\tanalyzed_taxon_count"
    )
    assert assumptions_path.read_text(encoding="utf-8").startswith(
        "assumption_id\troot_prior_mode\tfixed_root_state\troot_prior_distribution"
    )
    assert nodes_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tassumption_states\tassumption_confidences"
    )


def test_ancestral_ordered_discrete_cli_can_export_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ordered-discrete-summary.tsv"
    fits_path = tmp_path / "ordered-discrete-fits.tsv"
    nodes_path = tmp_path / "ordered-discrete-nodes.tsv"
    transitions_path = tmp_path / "ordered-discrete-transitions.tsv"
    exit_code = main(
        [
            "ancestral",
            "ordered-discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "equal-rates",
            "--ordered-states",
            "north,south,island",
            "--summary-out",
            str(summary_path),
            "--fits-out",
            str(fits_path),
            "--nodes-out",
            str(nodes_path),
            "--transitions-out",
            str(transitions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["ordered_state_count"] == 3
    assert payload["metrics"]["fit_count"] == 2
    assert payload["metrics"]["restricted_transition_count"] == 2
    assert payload["metrics"]["preferred_ordering"] == "ordered"
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tanalyzed_taxon_count\tstate_count"
    )
    assert fits_path.read_text(encoding="utf-8").startswith(
        "ordering_mode\tmodel\tstate_ordering\tordered_states"
    )
    assert nodes_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tordered_state\tunordered_state"
    )
    assert transitions_path.read_text(encoding="utf-8").startswith(
        "source_state\ttarget_state\tstep_distance\tordered_transition_allowed"
    )


@pytest.mark.slow
def test_ancestral_irreversible_discrete_cli_can_export_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "irreversible-discrete-summary.tsv"
    fits_path = tmp_path / "irreversible-discrete-fits.tsv"
    nodes_path = tmp_path / "irreversible-discrete-nodes.tsv"
    transitions_path = tmp_path / "irreversible-discrete-transitions.tsv"
    exit_code = main(
        [
            "ancestral",
            "irreversible-discrete",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "all-rates-different",
            "--allowed-transitions",
            "north->south,south->island",
            "--summary-out",
            str(summary_path),
            "--fits-out",
            str(fits_path),
            "--nodes-out",
            str(nodes_path),
            "--transitions-out",
            str(transitions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["model"] == "all-rates-different"
    assert payload["metrics"]["allowed_transition_count"] == 2
    assert payload["metrics"]["fit_count"] == 2
    assert payload["metrics"]["forbidden_transition_count"] == 4
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tanalyzed_taxon_count\tconstrained_log_likelihood"
    )
    assert fits_path.read_text(encoding="utf-8").startswith(
        "constraint_mode\tmodel\tanalyzed_taxon_count\tlog_likelihood"
    )
    assert nodes_path.read_text(encoding="utf-8").startswith(
        "node\tdescendant_taxa\tconstrained_state\tunconstrained_state"
    )
    assert transitions_path.read_text(encoding="utf-8").startswith(
        "source_state\ttarget_state\tconstrained_transition_allowed"
    )


def test_ancestral_confidence_cli_can_export_tree_set_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-tree-set-confidence-summary.tsv"
    confidence_path = tmp_path / "ancestral-tree-set-confidence.tsv"
    exit_code = main(
        [
            "ancestral",
            "confidence",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "body_mass",
            "--kind",
            "continuous",
            "--tree-set",
            "--summary-out",
            str(summary_path),
            "--confidence-out",
            str(confidence_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["kind"] == "continuous"
    assert payload["metrics"]["source_kind"] == "tree_set"
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["confidence_row_count"] == 14
    assert payload["metrics"]["top_uncertain_id"] == "A|B|C|D|E|F"
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tsource_kind\treconstruction_kind\ttarget_kind"
    )
    assert confidence_path.read_text(encoding="utf-8").startswith(
        "clade_id\tclade_taxa\ttree_presence_count\ttree_presence_fraction"
    )


def test_ancestral_tree_set_cli_can_export_continuous_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-tree-set-summary.tsv"
    trees_path = tmp_path / "ancestral-tree-set-trees.tsv"
    nodes_path = tmp_path / "ancestral-tree-set-nodes.tsv"
    clades_path = tmp_path / "ancestral-tree-set-clades.tsv"
    exclusions_path = tmp_path / "ancestral-tree-set-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "tree-set",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "body_mass",
            "--kind",
            "continuous",
            "--summary-out",
            str(summary_path),
            "--trees-out",
            str(trees_path),
            "--nodes-out",
            str(nodes_path),
            "--clades-out",
            str(clades_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["kind"] == "continuous"
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["clade_summary_count"] == 14
    assert payload["metrics"]["unstable_clade_count"] == 14
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\talpha"
    )
    assert trees_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert nodes_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tclade_id"
    )
    assert clades_path.read_text(encoding="utf-8").startswith(
        "clade_id\tclade_taxa\ttree_presence_count\ttree_presence_fraction"
    )
    assert exclusions_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_tree_set_cli_reports_discrete_stability_warnings(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "tree-set",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "habitat",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["kind"] == "discrete"
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["clade_summary_count"] == 14
    assert payload["metrics"]["unstable_clade_count"] == 14
    assert (
        "one or more comparable ancestral clades are absent from some retained trees"
        in payload["warnings"]
    )
    assert (
        "one or more discrete ancestral clades change state or support profile across retained trees"
        in payload["warnings"]
    )


def test_ancestral_transitions_cli_can_export_branch_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-transition-summary.tsv"
    branch_path = tmp_path / "ancestral-transition-branches.tsv"
    count_path = tmp_path / "ancestral-transition-counts.tsv"
    exclusion_path = tmp_path / "ancestral-transition-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "transitions",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_ancestral_sparse.tsv")),
            "--trait",
            "habitat",
            "--summary-out",
            str(summary_path),
            "--branches-out",
            str(branch_path),
            "--counts-out",
            str(count_path),
            "--exclusions-out",
            str(exclusion_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_set"] is False
    assert payload["metrics"]["model"] == "fitch"
    assert payload["metrics"]["total_branch_count"] == 6
    assert payload["metrics"]["changed_branch_count"] == 2
    assert payload["metrics"]["uncertain_change_count"] == 2
    assert payload["metrics"]["transition_pair_count"] == 2
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert branch_path.read_text(encoding="utf-8").startswith(
        "parent_node\tchild_node\tchild_descendant_taxa\tbranch_length"
    )
    assert count_path.read_text(encoding="utf-8").startswith(
        "transition\tsource_state\ttarget_state\tcertain_change_count"
    )
    assert exclusion_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_transitions_cli_can_export_tree_set_review(
    tmp_path: Path, capsys
) -> None:
    summary_path = tmp_path / "ancestral-transition-tree-set-summary.tsv"
    tree_path = tmp_path / "ancestral-transition-tree-set-trees.tsv"
    branch_path = tmp_path / "ancestral-transition-tree-set-branches.tsv"
    count_path = tmp_path / "ancestral-transition-tree-set-counts.tsv"
    exclusion_path = tmp_path / "ancestral-transition-tree-set-excluded.tsv"
    exit_code = main(
        [
            "ancestral",
            "transitions",
            str(fixture("example_posterior_tree_set_six_taxa.nwk")),
            str(fixture("example_traits_clade_summary.tsv")),
            "--trait",
            "habitat",
            "--model",
            "equal-rates",
            "--tree-set",
            "--summary-out",
            str(summary_path),
            "--trees-out",
            str(tree_path),
            "--branches-out",
            str(branch_path),
            "--counts-out",
            str(count_path),
            "--exclusions-out",
            str(exclusion_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_set"] is True
    assert payload["metrics"]["model"] == "equal-rates"
    assert payload["metrics"]["kept_tree_count"] == 5
    assert payload["metrics"]["transition_pair_count"] >= 2
    assert (
        payload["metrics"]["topology_sensitive_transition_pair_count"]
        + payload["metrics"]["uncertainty_sensitive_transition_pair_count"]
        >= 1
    )
    assert summary_path.read_text(encoding="utf-8").startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert tree_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert branch_path.read_text(encoding="utf-8").startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tparent_node"
    )
    assert count_path.read_text(encoding="utf-8").startswith(
        "transition\tsource_state\ttarget_state\ttree_presence_count"
    )
    assert exclusion_path.read_text(encoding="utf-8") == "taxon\treason\n"


def test_ancestral_render_cli_writes_svg_with_internal_annotations(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "ancestral.svg"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["format"] == "svg"
    assert payload["metrics"]["rendered_internal_annotation_count"] == 3
    assert 'class="internal-annotation-label"' in output.read_text(encoding="utf-8")


def test_ancestral_render_cli_can_export_discrete_pie_html(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "ancestral-discrete.html"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--discrete-node-style",
            "pies",
            "--branch-coloring",
            "state",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["format"] == "html"
    assert payload["metrics"]["rendered_internal_pie_count"] == 3
    assert payload["metrics"]["rendered_branch_color_count"] >= 4
    assert "<figure><svg" in output.read_text(encoding="utf-8")
    assert output.with_suffix(".svg").exists()


def test_ancestral_render_cli_can_export_png(tmp_path: Path, capsys) -> None:
    if shutil.which("rsvg-convert") is None and shutil.which("sips") is None:
        return
    output = tmp_path / "ancestral-discrete.png"
    exit_code = main(
        [
            "ancestral",
            "render",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--discrete-node-style",
            "pies",
            "--branch-coloring",
            "state",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["format"] == "png"
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert output.with_suffix(".svg").exists()


def test_ancestral_report_cli_writes_html_and_svg(tmp_path: Path, capsys) -> None:
    output = tmp_path / "ancestral-report.html"
    exit_code = main(
        [
            "ancestral",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "ancestral-state"
    assert output.exists()
    assert output.with_suffix(".svg").exists()


def test_ancestral_report_cli_can_export_full_review_package(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "ancestral-report-package"
    exit_code = main(
        [
            "ancestral",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--kind",
            "discrete",
            "--model",
            "equal-rates",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "ancestral-report-package"
    assert payload["metrics"]["reconstruction_kind"] == "discrete"
    assert payload["metrics"]["artifact_count"] == 13
    assert payload["metrics"]["methods_summary_warning_count"] >= 0
    assert output_dir.exists()
    assert (output_dir / "ancestral-report.html").exists()
    assert (output_dir / "ancestral-methods-summary.md").exists()
    assert (output_dir / "reviewer-audit-checklist.tsv").exists()
    assert (output_dir / "ancestral-figure.svg").exists()
    assert (output_dir / "summary.tsv").exists()
    assert (output_dir / "node-table.tsv").exists()
    assert (output_dir / "uncertainty-table.tsv").exists()
    assert (output_dir / "transition-counts.tsv").exists()
    assert (output_dir / "transition-branches.tsv").exists()
    assert (output_dir / "exclusions.tsv").exists()


def test_ancestral_sensitivity_cli_reports_available_comparisons(capsys) -> None:
    exit_code = main(
        [
            "ancestral",
            "sensitivity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--compare-model",
            "ou",
            "--compare-tree",
            str(fixture("example_tree_topology_diff.nwk")),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["has_model_sensitivity"] is True
    assert payload["metrics"]["has_tree_sensitivity"] is True


def test_ancestral_package_cli_writes_publication_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "ancestral-package"
    exit_code = main(
        [
            "ancestral",
            "package",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--kind",
            "continuous",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["artifact_count"] == 12
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["internal_state_visible"] is True
    assert payload["metrics"]["uncertainty_visible"] is True
    assert (output_dir / "figure-manifest.json").exists()
    assert (output_dir / "figure-reproducibility.manifest.json").exists()
    assert (output_dir / "ancestral-figure.png").exists()
    assert (output_dir / "ancestral-figure.html").exists()
    assert (output_dir / "ancestral-figure-review.html").exists()
    assert (output_dir / "node-uncertainty-review.tsv").exists()


def test_ancestral_discrete_cli_rejects_ordered_fitch_model(capsys) -> None:
    try:
        main(
            [
                "ancestral",
                "report",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--kind",
                "discrete",
                "--state-ordering",
                "ordered",
                "--out",
                "ignored.html",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    captured = capsys.readouterr()
    assert (
        "ordered ancestral discrete reconstruction requires a likelihood model"
        in captured.err
    )


def test_ancestral_discrete_cli_rejects_meristic_parity_claim(capsys) -> None:
    try:
        main(
            [
                "ancestral",
                "report",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--kind",
                "discrete",
                "--model",
                "meristic",
                "--state-ordering",
                "ordered",
                "--ordered-states",
                "north,south,island",
                "--out",
                "ignored.html",
            ]
        )
    except SystemExit as error:
        assert error.code == 2
    captured = capsys.readouterr()
    assert "explicitly excluded this round" in captured.err
    assert "ordered-state Mk support is not claimed as meristic parity" in (
        captured.err
    )
