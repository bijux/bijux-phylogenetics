from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.ancestral.presentation import build_ancestral_figure_package

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
    assert result.review_path.exists()
    assert result.node_table_path.exists()
    assert result.uncertainty_table_path.exists()
    assert result.node_review_path.exists()
    assert result.legend_path.exists()
    assert result.model_description_path.exists()
    assert result.caption_path.exists()
    assert result.manifest_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert "uncertainty" in result.legend_path.read_text(encoding="utf-8").lower()
    assert "model" in result.model_description_path.read_text(encoding="utf-8").lower()
    assert result.figure_png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert "<figure><svg" in result.figure_html_path.read_text(encoding="utf-8")
    assert result.audit.publication_ready is True
    assert result.audit.internal_state_visible is True
    assert result.audit.uncertainty_visible is True
    assert result.review_path.read_text(encoding="utf-8").count("Reviewer Summary") == 1
    assert "uncertainty_label" in result.node_review_path.read_text(encoding="utf-8")
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert reproducibility["report_kind"] == "ancestral_figure_package"
    assert reproducibility["model"]["kind"] == "continuous"


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
    assert 'class="internal-pie-slice"' in result.figure_path.read_text(
        encoding="utf-8"
    )
    assert result.audit.publication_ready is True
    assert result.audit.rendered_internal_pie_count >= 1
    assert result.audit.rendered_internal_annotation_count >= 1
    assert result.node_review_path.exists()
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    assert reproducibility["model"]["kind"] == "discrete"
