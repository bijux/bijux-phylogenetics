from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_nucleotide_sh_like_branch_support,
    evaluate_nucleotide_sh_like_branch_support_from_alignment,
    validate_nucleotide_sh_like_branch_support_replicate_count,
    write_nucleotide_sh_like_branch_support_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"
REFERENCE_TREE_NEWICK = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_sh_like_branch_support_ranks_strong_and_weak_local_signal() -> None:
    reference_tree = loads_newick(REFERENCE_TREE_NEWICK)
    strong_report = evaluate_nucleotide_sh_like_branch_support(
        reference_tree,
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        resampling_replicate_count=16,
        resampling_seed=7,
    )
    weak_report = evaluate_nucleotide_sh_like_branch_support(
        reference_tree,
        fixture("alignments", "jc69_likelihood_equal_best_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        resampling_replicate_count=16,
        resampling_seed=7,
    )

    assert strong_report.algorithm == "nucleotide-sh-like-branch-support"
    assert strong_report.model_name == "JC69"
    assert strong_report.branch_count == 2
    assert strong_report.site_count == 12
    assert strong_report.pattern_count == 2
    assert len(strong_report.local_topology_rows) == 6
    assert len(strong_report.resampling_rows) == 32
    assert all(
        row.caution_label == strong_report.caution_label
        for row in strong_report.branch_support_rows
    )
    assert all(row.support_fraction == 1.0 for row in strong_report.branch_support_rows)
    assert all(
        row.reference_is_observed_best for row in strong_report.branch_support_rows
    )
    assert all(
        strong_row.support_fraction > weak_row.support_fraction
        for strong_row, weak_row in zip(
            strong_report.branch_support_rows,
            weak_report.branch_support_rows,
            strict=True,
        )
    )
    assert all(
        row.caution_label == weak_report.caution_label
        for row in weak_report.branch_support_rows
    )
    assert all(
        row.support_fraction == 0.5625 for row in weak_report.branch_support_rows
    )


def test_sh_like_branch_support_writes_governed_outputs(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "reference-tree.nwk"
    tree_path.write_text(REFERENCE_TREE_NEWICK + "\n", encoding="utf-8")
    report = evaluate_nucleotide_sh_like_branch_support_from_alignment(
        tree_path,
        fixture("alignments", "jc69_likelihood_equal_best_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        resampling_replicate_count=8,
        resampling_seed=7,
    )

    outputs = write_nucleotide_sh_like_branch_support_artifacts(
        tmp_path / "sh-like-branch-support-run",
        report,
    )

    assert set(outputs) == {
        "reference_tree_path",
        "branch_support_path",
        "local_topology_path",
        "resampling_path",
        "run_json_path",
    }
    assert (
        outputs["branch_support_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tnode_label\tdescendant_taxa\talternative_arrangement_count\treference_log_likelihood\tbest_alternative_tree_id\tbest_alternative_tree_label\tbest_alternative_topology_fingerprint\tbest_alternative_log_likelihood\tobserved_delta_log_likelihood\treference_is_observed_best\tsupport_replicate_count\tsupport_fraction\tsupport_percent\tcaution_label\n"
        )
    )
    assert (
        outputs["local_topology_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tnode_label\tdescendant_taxa\tcandidate_tree_id\tcandidate_tree_label\tlocal_arrangement_kind\ttopology_fingerprint\tobserved_log_likelihood\tobserved_delta_log_likelihood\tobserved_best_local_arrangement\ttree_newick\n"
        )
    )
    assert (
        outputs["resampling_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "branch_id\tdescendant_taxa\treplicate_index\treference_resampled_log_likelihood\tbest_local_tree_id\tbest_local_tree_label\tbest_local_resampled_log_likelihood\tbest_alternative_tree_id\tbest_alternative_tree_label\tbest_alternative_resampled_log_likelihood\treference_delta_log_likelihood\treference_matches_or_beats_alternatives\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-sh-like-branch-support"
    assert payload["model_name"] == "JC69"
    assert payload["branch_count"] == 2
    assert payload["resampling_replicate_count"] == 8
    assert len(payload["branch_support_rows"]) == 2
    assert len(payload["local_topology_rows"]) == 6
    assert len(payload["resampling_rows"]) == 16


def test_validate_sh_like_branch_support_rejects_zero_replicates() -> None:
    try:
        validate_nucleotide_sh_like_branch_support_replicate_count(0)
    except ValueError as error:
        assert str(error) == "resampling_replicate_count must be at least one"
    else:
        raise AssertionError(
            "SH-like branch support must require at least one replicate"
        )
