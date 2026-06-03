from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
    summarize_continuous_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
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


def test_summarize_continuous_ancestral_report_tracks_root_and_exclusions() -> None:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("continuous_tree_mismatch")
    report = reconstruct_continuous_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="value",
    )
    summary = summarize_continuous_ancestral_report(report)
    assert summary.trait == "value"
    assert summary.model == "brownian"
    assert summary.analyzed_taxon_count == 3
    assert summary.excluded_taxon_count == 1
    assert summary.missing_tip_taxon_count == 1
    assert summary.non_numeric_tip_taxon_count == 0
    assert summary.internal_node_count == 2
    assert summary.root_node == "A|B|C"
    assert summary.root_lower_95_interval < summary.root_estimate
    assert summary.root_upper_95_interval > summary.root_estimate
    assert summary.tree_is_ultrametric is True
    assert summary.covariance_near_singular is False
    assert summary.covariance_condition_number is not None
    assert summary.log_likelihood is not None
    assert summary.residual_sigma_squared is not None


def test_write_continuous_ancestral_review_tables(tmp_path: Path) -> None:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("continuous_tree_mismatch")
    report = reconstruct_continuous_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="value",
    )
    summary_path = tmp_path / "continuous-ancestral-summary.tsv"
    uncertainty_path = tmp_path / "continuous-ancestral-uncertainty.tsv"
    exclusion_path = tmp_path / "continuous-ancestral-excluded.tsv"
    write_continuous_ancestral_summary_table(summary_path, report)
    write_continuous_ancestral_uncertainty_table(uncertainty_path, report)
    write_continuous_ancestral_exclusion_table(exclusion_path, report)
    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    uncertainty_rows = uncertainty_path.read_text(encoding="utf-8").splitlines()
    exclusion_rows = exclusion_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel\testimator\talpha")
    assert len(summary_rows) == 2
    assert "log_likelihood" in summary_rows[0]
    assert "covariance_condition_number" in summary_rows[0]
    assert "optimizer_name" in summary_rows[0]
    assert "optimizer_converged" in summary_rows[0]
    assert uncertainty_rows[0].startswith(
        "node\tnode_name\tdescendant_taxa\testimate\tstandard_error"
    )
    assert len(uncertainty_rows) == 3
    assert exclusion_rows == [
        "taxon\treason",
        "D\tmissing_trait_value",
    ]


def test_summarize_continuous_ancestral_report_tracks_fast_anc_estimator() -> None:
    tree_fixture = get_shared_tree_fixture("balanced_rooted_ultrametric")
    trait_fixture = get_shared_trait_table_fixture("continuous_tree_mismatch")
    report = reconstruct_continuous_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="value",
        estimator="fast-anc",
    )
    summary = summarize_continuous_ancestral_report(report)

    assert report.estimator == "fast-anc"
    assert summary.estimator == "fast-anc"
    assert summary.tree_is_ultrametric is True
    assert summary.root_standard_error > 0.0


def test_summarize_continuous_ancestral_report_tracks_anc_ml_optimizer() -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        estimator="anc-ml",
    )
    summary = summarize_continuous_ancestral_report(report)

    assert report.estimator == "anc-ml"
    assert summary.estimator == "anc-ml"
    assert summary.optimizer_name == "closed-form-profile-solution"
    assert summary.optimizer_converged is True
    assert summary.optimizer_iteration_count == 0
    assert summary.optimizer_function_evaluation_count == 1
    assert summary.optimizer_convergence_status == "profile-solved"


@pytest.mark.slow
def test_reconstruct_continuous_ancestral_states_matches_primate_reference_fixture() -> (
    None
):
    repository_root = Path(__file__).resolve().parents[3]
    parity_fixture = (
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-pgls-and-signal"
        / "evidence-009"
        / "results"
        / "ancestral-mode-parity.json"
    )
    reference = json.loads(parity_fixture.read_text(encoding="utf-8"))
    report = reconstruct_continuous_ancestral_states(
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_primate.csv",
        trait="longevity",
        taxon_column="species",
        model="brownian",
    )
    expected_rows = {
        row["node"]: row["estimate"]
        for row in reference["r_ancestral_reconstruction"]["brownian"]["rows"]
    }
    observed_rows = {
        estimate.node: estimate.estimate
        for estimate in report.estimates
        if not estimate.is_tip
    }
    assert (
        len(observed_rows)
        == reference["r_ancestral_reconstruction"]["brownian"]["node_count"]
    )
    assert observed_rows.keys() == expected_rows.keys()
    for node, expected_estimate in expected_rows.items():
        assert math.isclose(
            observed_rows[node],
            expected_estimate,
            rel_tol=0.0,
            abs_tol=5e-5,
        )


@pytest.mark.slow
def test_reconstruct_continuous_ancestral_states_tracks_live_ape_ace_when_available() -> (
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
    live_fit = subprocess.run(
        [
            rscript,
            "-e",
            (
                "library(ape);"
                "library(jsonlite);"
                "tree <- read.tree('evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk');"
                "dat <- read.csv('evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv', stringsAsFactors=FALSE);"
                "values <- dat$longevity[match(tree$tip.label, dat$species)];"
                "ace_result <- ace(values, tree, type='continuous', method='pic');"
                "leaf_descendants <- function(phy, node) {"
                "  if (node <= length(phy$tip.label)) return(phy$tip.label[[node]]);"
                "  children <- phy$edge[phy$edge[, 1] == node, 2];"
                "  sort(unique(unlist(lapply(children, function(child) leaf_descendants(phy, child)))))"
                "};"
                "node_signature <- function(phy, node) paste(leaf_descendants(phy, node), collapse='|');"
                "internal_nodes <- seq(length(tree$tip.label) + 1, length(tree$tip.label) + tree$Nnode);"
                "rows <- lapply(seq_along(internal_nodes), function(index) list(node=node_signature(tree, internal_nodes[[index]]), estimate=unname(signif(as.numeric(ace_result$ace[[index]]), 15))));"
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
        row["node"]: row["estimate"] for row in json.loads(live_fit.stdout)
    }
    report = reconstruct_continuous_ancestral_states(
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_primate.csv",
        trait="longevity",
        taxon_column="species",
        model="brownian",
    )
    observed_rows = {
        estimate.node: estimate.estimate
        for estimate in report.estimates
        if not estimate.is_tip
    }
    assert report.brownian_fit_diagnostics is not None
    assert math.isfinite(report.brownian_fit_diagnostics.log_likelihood)
    assert observed_rows.keys() == expected_rows.keys()
    for node, expected_estimate in expected_rows.items():
        assert math.isclose(
            observed_rows[node],
            expected_estimate,
            rel_tol=0.0,
            abs_tol=5e-5,
        )


def test_reconstruct_continuous_ancestral_states_tracks_live_ape_ace_on_shared_fixture_when_available() -> (
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
    trait_fixture = get_shared_trait_table_fixture("categorical_predictor_match")
    live_fit = subprocess.run(
        [
            rscript,
            "-e",
            (
                "library(ape);"
                "library(jsonlite);"
                f"tree <- read.tree('{tree_fixture.path}');"
                f"dat <- read.delim('{trait_fixture.path}', stringsAsFactors=FALSE);"
                "values <- dat$response[match(tree$tip.label, dat$taxon)];"
                "ace_result <- ace(values, tree, type='continuous', method='pic');"
                "leaf_descendants <- function(phy, node) {"
                "  if (node <= length(phy$tip.label)) return(phy$tip.label[[node]]);"
                "  children <- phy$edge[phy$edge[, 1] == node, 2];"
                "  sort(unique(unlist(lapply(children, function(child) leaf_descendants(phy, child)))))"
                "};"
                "node_signature <- function(phy, node) paste(leaf_descendants(phy, node), collapse='|');"
                "internal_nodes <- seq(length(tree$tip.label) + 1, length(tree$tip.label) + tree$Nnode);"
                "rows <- lapply(seq_along(internal_nodes), function(index) list(node=node_signature(tree, internal_nodes[[index]]), estimate=unname(signif(as.numeric(ace_result$ace[[index]]), 15))));"
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
        row["node"]: row["estimate"] for row in json.loads(live_fit.stdout)
    }
    report = reconstruct_continuous_ancestral_states(
        tree_fixture.path,
        trait_fixture.path,
        trait="response",
    )
    observed_rows = {
        estimate.node: estimate.estimate
        for estimate in report.estimates
        if not estimate.is_tip
    }
    assert report.brownian_fit_diagnostics is not None
    assert report.brownian_fit_diagnostics.tree_is_ultrametric is True
    assert observed_rows.keys() == expected_rows.keys()
    for node, expected_estimate in expected_rows.items():
        assert math.isclose(
            observed_rows[node],
            expected_estimate,
            rel_tol=0.0,
            abs_tol=5e-5,
        )
