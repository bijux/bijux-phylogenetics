from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative import build_comparative_model_figure_package

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


def test_build_comparative_model_figure_package_writes_publication_bundle(
    tmp_path: Path,
) -> None:
    result = build_comparative_model_figure_package(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="signal_strong",
        out_dir=tmp_path / "comparative-model-figure-package",
    )

    assert result.criteria_figure_path.exists()
    assert result.likelihood_figure_path.exists()
    assert result.parameter_figure_path.exists()
    assert result.fit_figure_path.exists()
    assert result.criteria_table_path.exists()
    assert result.likelihood_table_path.exists()
    assert result.parameter_table_path.exists()
    assert result.fit_table_path.exists()
    assert result.legend_path.exists()
    assert result.caption_path.exists()
    assert result.review_path.exists()
    assert result.manifest_path.exists()
    assert result.reproducibility_manifest_path.exists()
    assert result.audit.publication_ready is True
    assert result.audit.criteria_surface_visible is True
    assert result.audit.likelihood_surface_visible is True
    assert result.audit.parameter_surface_visible is True
    assert result.audit.fit_surface_visible is True
    assert result.audit.support_distinct is True
    assert result.audit.selected_model == "brownian"
    assert result.audit.finite_aicc_model_count == 2
    assert result.audit.plotted_model_count == 2
    assert result.audit.rendered_parameter_count == 5
    assert result.audit.rendered_fit_row_count == 2
    assert result.audit.aicc_delta is not None
    assert result.audit.aicc_delta > 2.0

    html = result.review_path.read_text(encoding="utf-8")
    assert "Bijux Comparative Model Comparison Review" in html
    assert "Information Criteria" in html
    assert "Likelihood" in html
    assert "Parameters" in html
    assert "Fit Summary" in html

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    assert manifest["report_kind"] == "comparative_model_figure_package"
    assert manifest["metrics"]["publication_ready"] is True
    assert manifest["metrics"]["selected_model"] == "brownian"
    assert manifest["reproducibility_manifest_path"] == str(
        result.reproducibility_manifest_path
    )
    assert reproducibility["report_kind"] == "comparative_model_figure_package"
    assert reproducibility["model"]["name"] == "brownian"


@pytest.mark.slow
def test_build_comparative_model_figure_package_blocks_ambiguous_support(
    tmp_path: Path,
) -> None:
    result = build_comparative_model_figure_package(
        fixture("example_tree_phytools_ultrametric_one_hundred_twenty_eight_taxa.nwk"),
        fixture("example_traits_phytools_signal_one_hundred_twenty_eight_taxa.tsv"),
        trait="signal_strong",
        out_dir=tmp_path / "comparative-model-figure-package-ambiguous",
    )

    assert result.audit.publication_ready is False
    assert result.audit.criteria_surface_visible is True
    assert result.audit.likelihood_surface_visible is True
    assert result.audit.parameter_surface_visible is True
    assert result.audit.fit_surface_visible is True
    assert result.audit.support_distinct is False
    assert result.audit.selected_model == "ou"
    assert result.audit.finite_aicc_model_count == 2
    assert result.audit.aicc_delta is not None
    assert result.audit.aicc_delta < 2.0
    assert any(
        "publication threshold of 2.0" in limitation
        for limitation in result.audit.limitations
    )
    assert result.review_path.exists()
    assert result.manifest_path.exists()
    reproducibility = json.loads(
        result.reproducibility_manifest_path.read_text(encoding="utf-8")
    )
    assert reproducibility["model"]["name"] == "ou"
