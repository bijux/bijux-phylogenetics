from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.ancestral.report_package import build_ancestral_report_package


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


def test_build_continuous_ancestral_report_package_writes_review_bundle(
    tmp_path: Path,
) -> None:
    result = build_ancestral_report_package(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_comparative.tsv"),
        trait="response",
        reconstruction_kind="continuous",
        out_dir=tmp_path / "continuous-report",
        model="brownian",
    )

    assert result.report_path.exists()
    assert result.figure_path.exists()
    assert result.figure_png_path.exists()
    assert result.figure_html_path.exists()
    assert result.summary_table_path.exists()
    assert result.node_table_path.exists()
    assert result.uncertainty_table_path.exists()
    assert result.transition_count_table_path.exists()
    assert result.transition_branch_table_path.exists()
    assert result.exclusion_table_path.exists()
    assert result.manifest_path.exists()
    assert result.figure_png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")

    html = result.report_path.read_text(encoding="utf-8")
    assert "Bijux Ancestral Reconstruction Report" in html
    assert "Tree Visualization" in html
    assert "Transition Review" in html

    summary_rows = result.summary_table_path.read_text(encoding="utf-8").splitlines()
    node_rows = result.node_table_path.read_text(encoding="utf-8").splitlines()
    uncertainty_rows = result.uncertainty_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    transition_count_rows = result.transition_count_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    transition_branch_rows = result.transition_branch_table_path.read_text(
        encoding="utf-8"
    ).splitlines()

    assert summary_rows[0].startswith("trait\ttaxon_column\tmodel\talpha")
    assert node_rows[0].startswith("node\tnode_name\tis_tip\tdescendant_taxa")
    assert uncertainty_rows[0].startswith(
        "node\tnode_name\tdescendant_taxa\testimate\tstandard_error"
    )
    assert transition_count_rows[0].startswith(
        "direction\tbranch_count\tbranch_fraction"
    )
    assert transition_branch_rows[0].startswith(
        "parent_node\tchild_node\tchild_descendant_taxa"
    )

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["report_kind"] == "ancestral_report_package"
    assert manifest["reconstruction_kind"] == "continuous"
    assert manifest["metrics"]["analyzed_taxon_count"] == 4


def test_build_discrete_ancestral_report_package_writes_probabilities_and_transitions(
    tmp_path: Path,
) -> None:
    result = build_ancestral_report_package(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_geography.tsv"),
        trait="region",
        reconstruction_kind="discrete",
        out_dir=tmp_path / "discrete-report",
        model="equal-rates",
    )

    assert result.report_path.exists()
    assert result.figure_path.exists()
    assert result.transition_report is not None

    uncertainty_rows = result.uncertainty_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    transition_count_rows = result.transition_count_table_path.read_text(
        encoding="utf-8"
    ).splitlines()
    transition_branch_rows = result.transition_branch_table_path.read_text(
        encoding="utf-8"
    ).splitlines()

    assert uncertainty_rows[0].startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_state\tstate_set"
    )
    assert "state_probabilities" in uncertainty_rows[0]
    assert "confidence" in uncertainty_rows[0]
    assert transition_count_rows[0].startswith(
        "transition\tsource_state\ttarget_state"
    )
    assert transition_branch_rows[0].startswith(
        "parent_node\tchild_node\tchild_descendant_taxa"
    )
    assert 'class="internal-pie-slice"' in result.figure_path.read_text(
        encoding="utf-8"
    )
