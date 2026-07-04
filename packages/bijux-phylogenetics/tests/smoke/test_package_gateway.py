from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics as package
from bijux_phylogenetics.core import copy_example_inputs, example_resource_root

pytestmark = pytest.mark.smoke


def test_package_gateway_loads_curated_public_modules() -> None:
    assert package.__version__
    assert package.__all__ == [
        "__version__",
        "ancestral",
        "api",
        "bayesian",
        "biogeography",
        "comparative",
        "datasets",
        "distance",
        "evidence",
        "parsimony",
        "parity",
        "phylo",
        "trees",
    ]

    for module_name in package.__all__[1:]:
        module = getattr(package, module_name)
        assert module.__name__ == f"bijux_phylogenetics.{module_name}"

    assert set(package.__all__).issubset(set(dir(package)))


def test_copy_example_inputs_creates_writable_package_fixture_set(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "example-inputs"
    copied = copy_example_inputs(destination)

    assert set(copied) == {"alignment", "alt_tree", "metadata", "traits", "tree"}

    source_root = example_resource_root()
    for path in copied.values():
        assert path.exists()
        assert path.is_file()
        assert path.is_relative_to(destination)
        assert not path.is_relative_to(source_root)

    assert copied["tree"].read_text(encoding="utf-8").strip().endswith(";")
