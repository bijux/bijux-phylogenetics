from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative import compute_phylogenetic_signal_test
from bijux_phylogenetics.core import copy_example_inputs, example_resource_root


def test_example_resource_root_exposes_packaged_example_files() -> None:
    resource_root = example_resource_root()

    assert resource_root.joinpath("alignments", "example_alignment.fasta").exists()
    assert resource_root.joinpath("trees", "example_tree.nwk").exists()
    assert resource_root.joinpath("metadata", "example_traits.tsv").exists()


def test_copy_example_inputs_materializes_comparative_ready_examples(
    tmp_path: Path,
) -> None:
    copied = copy_example_inputs(tmp_path / "examples")
    report = compute_phylogenetic_signal_test(
        copied["tree"],
        copied["traits"],
        trait="value",
        permutations=19,
        seed=7,
    )

    assert set(copied) == {"alignment", "alt_tree", "metadata", "traits", "tree"}
    assert copied["alignment"].exists()
    assert copied["tree"].exists()
    assert copied["traits"].exists()
    assert report.taxon_count == 4
    assert report.input_audit.pruned_missing_value_taxa == []
