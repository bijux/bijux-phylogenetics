from __future__ import annotations

import json

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.core import copy_example_inputs

pytestmark = pytest.mark.smoke


def test_inspect_command_reports_example_tree_metrics(tmp_path, capsys) -> None:
    inputs = copy_example_inputs(tmp_path / "example-inputs")

    exit_code = main(["inspect", str(inputs["tree"]), "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tip_count"] == 4
    assert payload["metrics"]["is_binary"] is True


def test_alignment_validate_input_reports_example_alignment_metrics(
    tmp_path, capsys
) -> None:
    inputs = copy_example_inputs(tmp_path / "example-inputs")

    exit_code = main(
        ["alignment", "validate-input", str(inputs["alignment"]), "--json"]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["sequence_count"] == 4
