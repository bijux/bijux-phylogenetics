from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.sensitivity import build_ancestral_sensitivity_report

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


def test_build_continuous_ancestral_sensitivity_report_summarizes_model_tree_and_pruning() -> (
    None
):
    report = build_ancestral_sensitivity_report(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_comparative.tsv"),
        trait="response",
        reconstruction_kind="continuous",
        model="brownian",
        compare_model="ou",
        compare_tree_path=fixture("example_tree_topology_diff.nwk"),
        drop_taxa=["D"],
    )
    assert report.model_sensitivity is not None
    assert report.tree_sensitivity is not None
    assert report.pruning_sensitivity is not None


def test_build_discrete_ancestral_sensitivity_report_summarizes_coding_changes() -> (
    None
):
    report = build_ancestral_sensitivity_report(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_comparative.tsv"),
        trait="habitat",
        reconstruction_kind="discrete",
        model="fitch",
        compare_model="equal-rates",
        coding_map={"forest": "temperate", "tundra": "cold"},
    )
    assert report.model_sensitivity is not None
    assert report.trait_coding_sensitivity is not None
