from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    CandidateTreeQuartetScoreReport,
    compute_candidate_tree_quartet_score,
    write_candidate_tree_quartet_score_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_tree_gateway_exports_quartet_score_surface() -> None:
    assert trees_api.CandidateTreeQuartetScoreReport is CandidateTreeQuartetScoreReport
    assert (
        trees_api.compute_candidate_tree_quartet_score
        is compute_candidate_tree_quartet_score
    )
    assert (
        trees_api.write_candidate_tree_quartet_score_table
        is write_candidate_tree_quartet_score_table
    )


def test_compute_candidate_tree_quartet_score_prefers_known_better_candidate() -> None:
    gene_tree_set_path = fixture("quartet_concordance_gene_trees_4_taxa.nwk")
    higher = compute_candidate_tree_quartet_score(
        fixture("quartet_score_candidate_high_4_taxa.nwk"),
        gene_tree_set_path,
    )
    lower = compute_candidate_tree_quartet_score(
        fixture("quartet_score_candidate_low_4_taxa.nwk"),
        gene_tree_set_path,
    )

    assert higher.quartet_score == 2
    assert higher.normalized_quartet_score == 0.5
    assert lower.quartet_score == 1
    assert lower.normalized_quartet_score == 0.25
    assert higher.quartet_score > lower.quartet_score
    assert higher.rows[0].branch_id == "A|B::C|D"
    assert lower.rows[0].branch_id == "A|C::B|D"


def test_compute_candidate_tree_quartet_score_requires_exact_taxon_match(
    tmp_path: Path,
) -> None:
    candidate_tree = tmp_path / "quartet-score-candidate.nwk"
    gene_tree_set = tmp_path / "quartet-score-gene-tree-set.nwk"
    candidate_tree.write_text("((A,B),(C,D));\n", encoding="utf-8")
    gene_tree_set.write_text("((A,B),(C,E));\n", encoding="utf-8")

    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        compute_candidate_tree_quartet_score(
            candidate_tree,
            gene_tree_set,
        )
