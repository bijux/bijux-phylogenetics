from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral import (
    build_ancestral_methods_summary_text,
    reconstruct_continuous_ancestral_states,
    reconstruct_discrete_ancestral_states,
    write_ancestral_methods_summary_text,
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


def test_build_continuous_ancestral_methods_summary_text_reports_model_and_uncertainty() -> (
    None
):
    reconstruction = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
        model="brownian",
    )

    text = build_ancestral_methods_summary_text(
        reconstruction_kind="continuous",
        reconstruction=reconstruction,
    )

    assert "Ancestral Reconstruction Methods Summary" in text
    assert "- continuous model: `brownian`" in text
    assert "- estimator: `ace-pic`" in text
    assert "- reconstructed internal node count: `3`" in text
    assert "- node uncertainty is reported as standard error plus 95% interval" in text
    assert "- reconstruction warnings:" in text


def test_write_discrete_ancestral_methods_summary_text_writes_markdown(
    tmp_path: Path,
) -> None:
    reconstruction = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    output_path = tmp_path / "ancestral-methods-summary.md"

    result = write_ancestral_methods_summary_text(
        output_path,
        reconstruction_kind="discrete",
        reconstruction=reconstruction,
    )

    assert result.output_path == output_path
    assert result.reconstruction_kind == "discrete"
    assert result.model == "equal-rates"
    assert result.analyzed_taxon_count == 4
    assert "Ancestral Reconstruction Methods Summary" in result.text
    assert "- discrete model: `equal-rates`" in result.text
    assert "- root prior mode: `equal`" in result.text
    assert (
        "- node uncertainty is reported as marginal state probabilities" in result.text
    )
    assert output_path.read_text(encoding="utf-8") == result.text
