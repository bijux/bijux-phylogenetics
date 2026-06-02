from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from bijux_phylogenetics.ancestral.discrete import (
    discrete_ancestral_exclusions,
    reconstruct_discrete_ancestral_states,
    summarize_discrete_ancestral_report,
    write_discrete_ancestral_exclusion_table,
    write_discrete_ancestral_fit_table,
    write_discrete_ancestral_probability_table,
    write_discrete_ancestral_summary_table,
    write_discrete_ancestral_transition_table,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_trait_table_fixture,
    get_shared_tree_fixture,
)

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


def test_summarize_discrete_ancestral_report_tracks_root_and_exclusions(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "traits.tsv"
    table_path.write_text(
        "taxon\thabitat\nA\tforest\nB\tforest\nC\ttundra\nD\t\n",
        encoding="utf-8",
    )
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        table_path,
        trait="habitat",
        model="equal-rates",
    )
    summary = summarize_discrete_ancestral_report(report)
    exclusions = discrete_ancestral_exclusions(report)
    assert summary.model == "equal-rates"
    assert summary.analyzed_taxon_count == 3
    assert summary.excluded_taxon_count == 1
    assert summary.internal_node_count == 2
    assert summary.root_node == "A|B|C"
    assert summary.root_most_likely_state in {"forest", "tundra"}
    assert 0.0 <= summary.root_confidence <= 1.0
    assert exclusions == [
        type(exclusions[0])(taxon="D", reason="missing_discrete_trait_state")
    ]


def test_fitch_parsimony_reports_minimal_changes_and_ambiguous_nodes(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "traits.tsv"
    table_path.write_text(
        "taxon\thabitat\nA\tforest\nB\tdesert\nC\tforest\nD\tdesert\n",
        encoding="utf-8",
    )
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        table_path,
        trait="habitat",
        model="fitch",
    )
    summary = summarize_discrete_ancestral_report(report)
    assert report.minimal_change_count == 2
    assert report.parsimonious_root_state_count == 2
    assert summary.ambiguous_internal_node_count == 3
    assert summary.minimal_change_count == 2
    assert summary.parsimonious_root_state_count == 2


def test_write_discrete_ancestral_review_tables(tmp_path: Path) -> None:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("multistate_discrete_match")
    report = reconstruct_discrete_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="region",
        model="equal-rates",
    )
    summary_path = tmp_path / "discrete-ancestral-summary.tsv"
    probability_path = tmp_path / "discrete-ancestral-probabilities.tsv"
    exclusion_path = tmp_path / "discrete-ancestral-excluded.tsv"
    write_discrete_ancestral_summary_table(summary_path, report)
    write_discrete_ancestral_probability_table(probability_path, report)
    write_discrete_ancestral_exclusion_table(exclusion_path, report)
    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    probability_rows = probability_path.read_text(encoding="utf-8").splitlines()
    exclusion_rows = exclusion_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel\tstate_ordering")
    assert len(summary_rows) == 2
    assert probability_rows[0].startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_state\tstate_set"
    )
    assert len(probability_rows) == 4
    assert exclusion_rows == ["taxon\treason"]


def test_write_discrete_ancestral_er_fit_tables(tmp_path: Path) -> None:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("multistate_discrete_match")
    report = reconstruct_discrete_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="region",
        model="equal-rates",
        root_prior_mode="empirical",
    )
    summary_path = tmp_path / "discrete-ancestral-summary.tsv"
    transition_path = tmp_path / "discrete-ancestral-transitions.tsv"
    write_discrete_ancestral_summary_table(summary_path, report)
    write_discrete_ancestral_transition_table(transition_path, report)
    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    transition_rows = transition_path.read_text(encoding="utf-8").splitlines()
    assert "root_prior_mode" in summary_rows[0]
    assert "log_likelihood" in summary_rows[0]
    assert "parameter_count" in summary_rows[0]
    assert "aic" in summary_rows[0]
    assert "\tempirical\t" in summary_rows[1]
    assert transition_rows[0] == (
        "source_state\ttarget_state\ttransition_allowed\tstep_distance\trate"
    )
    assert len(transition_rows) == 7


def test_write_discrete_ancestral_sym_fit_table_tracks_optimizer_and_baseline(
    tmp_path: Path,
) -> None:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_six_taxon")
    trait_fixture = get_shared_trait_table_fixture("ace_discrete_sym_six_taxon")
    report = reconstruct_discrete_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="region",
        model="symmetric",
    )
    fit_path = tmp_path / "discrete-ancestral-fit.tsv"
    write_discrete_ancestral_fit_table(fit_path, report)
    fit_rows = fit_path.read_text(encoding="utf-8").splitlines()
    assert fit_rows[0].startswith(
        "model\ttaxon_count\tstate_count\tparameter_count\tlog_likelihood"
    )
    assert "\tsymmetric\t6\t3\t3\t" in f"\t{fit_rows[1]}\t"
    assert "\tnelder-mead\t" in f"\t{fit_rows[1]}\t"
    assert report.optimizer_diagnostics is not None
    assert report.optimizer_diagnostics.optimizer_name == "nelder-mead"
    assert report.baseline_comparison is not None
    assert report.baseline_comparison.baseline_model == "equal-rates"
    assert report.baseline_comparison.delta_aic > 0.0
    assert report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    transition_pairs = {
        (row.source_state, row.target_state): row.rate
        for row in report.transition_rate_rows
    }
    assert all(
        math.isfinite(rate) and rate >= 0.0 for rate in transition_pairs.values()
    )
    assert transition_pairs[("island", "north")] == pytest.approx(
        transition_pairs[("north", "island")],
        rel=1e-9,
        abs=1e-12,
    )
    assert transition_pairs[("island", "south")] == pytest.approx(
        transition_pairs[("south", "island")],
        rel=1e-9,
        abs=1e-12,
    )
    assert transition_pairs[("north", "south")] == pytest.approx(
        transition_pairs[("south", "north")],
        rel=1e-9,
        abs=1e-12,
    )


def test_reconstruct_discrete_ancestral_states_flags_overparameterized_sym_case() -> (
    None
):
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography_four_state.tsv"),
        trait="region",
        model="symmetric",
    )
    summary = summarize_discrete_ancestral_report(report)
    assert report.overparameterized is True
    assert summary.overparameterized is True
    assert report.optimizer_diagnostics is not None
    assert (
        "the discrete likelihood fit is likely overparameterized relative to the analyzed taxon count"
        in report.warnings
    )


@pytest.mark.slow
def test_write_discrete_ancestral_ard_fit_table_tracks_directional_rates_and_weak_fit_warnings(
    tmp_path: Path,
) -> None:
    tree_fixture = get_shared_tree_fixture("pectinate_rooted_non_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("ace_discrete_ard_pectinate")
    report = reconstruct_discrete_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="region",
        model="all-rates-different",
    )
    fit_path = tmp_path / "discrete-ancestral-fit.tsv"
    write_discrete_ancestral_fit_table(fit_path, report)
    fit_rows = fit_path.read_text(encoding="utf-8").splitlines()
    assert "\tall-rates-different\t4\t3\t6\t" in f"\t{fit_rows[1]}\t"
    assert "\ttrue\t" in f"\t{fit_rows[1].lower()}\t"
    assert report.overparameterized is True
    assert report.optimizer_diagnostics is not None
    assert report.optimizer_diagnostics.hit_lower_parameter_bound is True
    assert report.baseline_comparison is not None
    assert report.baseline_comparison.preferred_model_by_aic == "equal-rates"
    transition_pairs = {
        (row.source_state, row.target_state): row.rate
        for row in report.transition_rate_rows
    }
    assert (
        transition_pairs[("island", "south")] != transition_pairs[("south", "island")]
    )
    assert (
        "one or more discrete rate parameters hit an optimizer bound and should be interpreted as weakly identified"
        in report.warnings
    )


def test_fitch_summary_table_includes_parsimony_counts(tmp_path: Path) -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_ancestral_sparse.tsv"),
        trait="habitat",
        model="fitch",
    )
    summary_path = tmp_path / "discrete-ancestral-summary.tsv"
    write_discrete_ancestral_summary_table(summary_path, report)
    rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering\troot_prior_mode\tfixed_root_state\tanalyzed_taxon_count"
    )
    assert "ambiguous_internal_node_count" in rows[0]
    assert "minimal_change_count" in rows[0]
    assert "\t1\tA|B|C|D\tforest\t1.0\t" in rows[1]


@pytest.mark.slow
def test_reconstruct_discrete_ancestral_states_matches_checked_ace_reference() -> None:
    reference = json.loads(
        fixture("discrete_ancestral_ace_reference.json").read_text(encoding="utf-8")
    )
    repository_root = Path(__file__).resolve().parents[3]
    for case in reference["cases"]:
        report = reconstruct_discrete_ancestral_states(
            repository_root / case["tree"],
            repository_root / case["table"],
            trait=case["trait"],
            taxon_column=case["taxon_column"],
            model=case["model"],
        )
        observed_rows = {
            estimate.node: estimate.state_probabilities
            for estimate in report.estimates
            if not estimate.is_tip
        }
        for row in case["rows"]:
            for state, expected_probability in row["probabilities"].items():
                assert math.isclose(
                    observed_rows[row["node"]][state],
                    expected_probability,
                    rel_tol=0.0,
                    abs_tol=case["abs_tolerance"],
                )


@pytest.mark.slow
def test_reconstruct_discrete_ancestral_states_tracks_live_ape_ace_when_available() -> (
    None
):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available")
    repository_root = Path(__file__).resolve().parents[3]
    r_library = repository_root / "artifacts" / "r-lib"
    environment = dict(os.environ)
    if r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    package_check = subprocess.run(
        [
            rscript,
            "-e",
            (
                "cat(requireNamespace('ape', quietly=TRUE), '\\n');"
                "cat(requireNamespace('jsonlite', quietly=TRUE), '\\n')"
            ),
        ],
        capture_output=True,
        check=False,
        cwd=repository_root,
        env=environment,
        text=True,
    )
    if package_check.returncode != 0:
        pytest.skip("R package availability could not be checked")
    package_flags = [
        line.strip() for line in package_check.stdout.splitlines() if line.strip()
    ]
    if package_flags != ["TRUE", "TRUE"]:
        pytest.skip("ape and jsonlite are required for live ace validation")
    reference = json.loads(
        fixture("discrete_ancestral_ace_reference.json").read_text(encoding="utf-8")
    )
    r_models = {
        "equal-rates": "ER",
        "symmetric": "SYM",
        "all-rates-different": "ARD",
    }
    for case in reference["cases"]:
        expr = (
            "library(ape);"
            "library(jsonlite);"
            f"tree <- read.tree('{case['tree']}');"
            f"dat <- read.delim('{case['table']}', stringsAsFactors=FALSE);"
            f"values <- setNames(dat${case['trait']}, dat${case['taxon_column']});"
            f"ace_result <- ace(values[tree$tip.label], tree, type='discrete', model='{r_models[case['model']]}', CI=TRUE);"
            "leaf_descendants <- function(phy, node) {"
            "  if (node <= length(phy$tip.label)) return(phy$tip.label[[node]]);"
            "  children <- phy$edge[phy$edge[, 1] == node, 2];"
            "  sort(unique(unlist(lapply(children, function(child) leaf_descendants(phy, child)))))"
            "};"
            "node_signature <- function(phy, node) paste(leaf_descendants(phy, node), collapse='|');"
            "internal_nodes <- seq(length(tree$tip.label) + 1, length(tree$tip.label) + tree$Nnode);"
            "state_names <- colnames(ace_result$lik.anc);"
            "rows <- lapply(seq_along(internal_nodes), function(index) list(node=node_signature(tree, internal_nodes[[index]]), probabilities=as.list(setNames(as.numeric(ace_result$lik.anc[index, ]), state_names))));"
            "cat(toJSON(rows, auto_unbox=TRUE))"
        )
        live_fit = subprocess.run(
            [rscript, "-e", expr],
            capture_output=True,
            check=False,
            cwd=repository_root,
            env=environment,
            text=True,
        )
        if live_fit.returncode != 0:
            pytest.skip("live ape::ace execution failed in this environment")
        expected_rows = {
            row["node"]: row["probabilities"] for row in json.loads(live_fit.stdout)
        }
        report = reconstruct_discrete_ancestral_states(
            repository_root / case["tree"],
            repository_root / case["table"],
            trait=case["trait"],
            taxon_column=case["taxon_column"],
            model=case["model"],
        )
        observed_rows = {
            estimate.node: estimate.state_probabilities
            for estimate in report.estimates
            if not estimate.is_tip
        }
        for node, expected_probabilities in expected_rows.items():
            for state, expected_probability in expected_probabilities.items():
                assert math.isclose(
                    observed_rows[node][state],
                    expected_probability,
                    rel_tol=0.0,
                    abs_tol=next(
                        entry["abs_tolerance"]
                        for entry in reference["cases"]
                        if entry["case_id"] == case["case_id"]
                    ),
                )


def test_reconstruct_discrete_ancestral_states_tracks_live_ape_ace_on_shared_binary_fixture_when_available() -> (
    None
):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is not available")
    repository_root = Path(__file__).resolve().parents[3]
    r_library = repository_root / "artifacts" / "r-lib"
    environment = dict(os.environ)
    if r_library.is_dir():
        environment["R_LIBS_USER"] = str(r_library)
    package_check = subprocess.run(
        [
            rscript,
            "-e",
            (
                "cat(requireNamespace('ape', quietly=TRUE), '\\n');"
                "cat(requireNamespace('jsonlite', quietly=TRUE), '\\n')"
            ),
        ],
        capture_output=True,
        check=False,
        cwd=repository_root,
        env=environment,
        text=True,
    )
    if package_check.returncode != 0:
        pytest.skip("R package availability could not be checked")
    package_flags = [
        line.strip() for line in package_check.stdout.splitlines() if line.strip()
    ]
    if package_flags != ["TRUE", "TRUE"]:
        pytest.skip("ape and jsonlite are required for live ace validation")
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("binary_discrete_match")
    live_fit = subprocess.run(
        [
            rscript,
            "-e",
            (
                "library(ape);"
                "library(jsonlite);"
                f"tree <- read.tree('{tree_fixture.path}');"
                f"dat <- read.delim('{trait_fixture.path}', stringsAsFactors=FALSE);"
                "values <- setNames(dat$presence, dat$taxon);"
                "ace_result <- ace(values[tree$tip.label], tree, type='discrete', model='ER', CI=TRUE);"
                "leaf_descendants <- function(phy, node) {"
                "  if (node <= length(phy$tip.label)) return(phy$tip.label[[node]]);"
                "  children <- phy$edge[phy$edge[, 1] == node, 2];"
                "  sort(unique(unlist(lapply(children, function(child) leaf_descendants(phy, child)))))"
                "};"
                "node_signature <- function(phy, node) paste(leaf_descendants(phy, node), collapse='|');"
                "internal_nodes <- seq(length(tree$tip.label) + 1, length(tree$tip.label) + tree$Nnode);"
                "state_names <- colnames(ace_result$lik.anc);"
                "rows <- lapply(seq_along(internal_nodes), function(index) list(node=node_signature(tree, internal_nodes[[index]]), probabilities=as.list(setNames(as.numeric(ace_result$lik.anc[index, ]), state_names))));"
                "cat(toJSON(rows, auto_unbox=TRUE))"
            ),
        ],
        capture_output=True,
        check=False,
        cwd=repository_root,
        env=environment,
        text=True,
    )
    if live_fit.returncode != 0:
        pytest.skip("live ape::ace execution failed in this environment")
    expected_rows = {
        row["node"]: row["probabilities"] for row in json.loads(live_fit.stdout)
    }
    report = reconstruct_discrete_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="presence",
        model="equal-rates",
    )
    observed_rows = {
        estimate.node: estimate.state_probabilities
        for estimate in report.estimates
        if not estimate.is_tip
    }
    assert observed_rows.keys() == expected_rows.keys()
    for node, expected_probabilities in expected_rows.items():
        for state, expected_probability in expected_probabilities.items():
            assert math.isclose(
                observed_rows[node][state],
                expected_probability,
                rel_tol=0.0,
                abs_tol=5e-5,
            )
