from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.reports import write_tree_inference_methods_summary_text

FIXTURES = Path(__file__).parent / "fixtures" / "expected" / "fasta_to_tree"


def workflow_fixture(dataset: str) -> Path:
    return FIXTURES / dataset / f"{dataset}.manifest.json"


def test_write_tree_inference_methods_summary_text_reports_fasta_to_tree_steps(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "tree-inference-methods-summary.md"

    result = write_tree_inference_methods_summary_text(
        output_path,
        workflow_manifest_path=workflow_fixture("pleistocene-bear-cytb-fragments"),
    )

    assert result.output_path == output_path
    assert result.selected_model == "HKY+F"
    assert result.bootstrap_replicates == 1000
    assert result.trimmed_alignment_length == 1140
    assert result.supported_node_count == 2
    assert "Tree Inference Methods Summary" in result.text
    assert "- alignment mode: `auto`" in result.text
    assert "- trimming mode: `gap-threshold`" in result.text
    assert "- selected substitution model: `HKY+F`" in result.text
    assert "- governing information criterion: `BIC`" in result.text
    assert "- bootstrap replicates: `1000`" in result.text
    assert (
        "- no outgroup rooting, midpoint rerooting, or posterior summarization step is recorded in this fasta-to-tree workflow manifest"
        in result.text
    )
    assert output_path.read_text(encoding="utf-8") == result.text
