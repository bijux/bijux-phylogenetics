from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.package import build_ancestral_figure_package

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


def test_build_ancestral_figure_package_writes_publication_artifacts(
    tmp_path: Path,
) -> None:
    result = build_ancestral_figure_package(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_comparative.tsv"),
        trait="response",
        reconstruction_kind="continuous",
        out_dir=tmp_path / "package",
        model="brownian",
    )
    assert result.figure_path.exists()
    assert result.figure_png_path.exists()
    assert result.figure_html_path.exists()
    assert result.node_table_path.exists()
    assert result.uncertainty_table_path.exists()
    assert result.legend_path.exists()
    assert result.model_description_path.exists()
    assert result.caption_path.exists()
    assert result.manifest_path.exists()
    assert "uncertainty" in result.legend_path.read_text(encoding="utf-8").lower()
    assert "model" in result.model_description_path.read_text(encoding="utf-8").lower()
    assert result.figure_png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert "<figure><svg" in result.figure_html_path.read_text(encoding="utf-8")


def test_build_ancestral_figure_package_uses_discrete_probability_ledger(
    tmp_path: Path,
) -> None:
    result = build_ancestral_figure_package(
        tree_path=fixture("example_tree.nwk"),
        traits_path=fixture("example_traits_geography.tsv"),
        trait="region",
        reconstruction_kind="discrete",
        out_dir=tmp_path / "discrete-package",
        model="equal-rates",
    )
    uncertainty_rows = result.uncertainty_table_path.read_text(encoding="utf-8")
    assert "state_probabilities" in uncertainty_rows
    assert "most_likely_state" in uncertainty_rows
    assert 'class="internal-pie-slice"' in result.figure_path.read_text(encoding="utf-8")
