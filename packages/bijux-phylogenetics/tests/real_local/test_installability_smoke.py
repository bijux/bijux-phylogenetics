from __future__ import annotations

import json
from pathlib import Path

import pytest

from ..support.installability_smoke import (
    PACKAGE_ROOT,
    assert_distribution_contains_packaged_resources,
    build_installable_distributions,
    copy_installed_example_inputs,
    create_clean_virtualenv,
    install_distribution,
    run_installed_cli,
)

pytestmark = [pytest.mark.real_local, pytest.mark.slow]


def test_installable_distributions_run_core_cli_smoke_commands(
    tmp_path: Path,
) -> None:
    wheel_path, sdist_path = build_installable_distributions(tmp_path / "dist")

    for distribution_path in (wheel_path, sdist_path):
        assert_distribution_contains_packaged_resources(distribution_path)

        install_root = tmp_path / distribution_path.name.replace(".", "-")
        venv_root = install_root / "venv"
        work_root = install_root / "work"
        work_root.mkdir(parents=True, exist_ok=True)
        venv_python = create_clean_virtualenv(venv_root)
        install_distribution(venv_python, distribution_path)

        help_result = run_installed_cli(venv_root, ["--help"], cwd=work_root)
        assert help_result.returncode == 0
        assert "bijux-phylogenetics" in help_result.stdout

        copied_inputs = copy_installed_example_inputs(
            venv_python,
            work_root / "example-inputs",
        )
        assert set(copied_inputs) == {
            "alignment",
            "alt_tree",
            "metadata",
            "traits",
            "tree",
        }
        for path in copied_inputs.values():
            assert path.exists()
            assert not path.is_relative_to(PACKAGE_ROOT)

        validation_result = run_installed_cli(
            venv_root,
            [
                "alignment",
                "validate-input",
                str(copied_inputs["alignment"]),
                "--json",
            ],
            cwd=work_root,
        )
        validation_payload = json.loads(validation_result.stdout)
        assert validation_payload["status"] == "ok"
        assert validation_payload["metrics"]["sequence_count"] == 4

        tree_report_path = work_root / "tree-report.html"
        tree_result = run_installed_cli(
            venv_root,
            [
                "report",
                "tree",
                str(copied_inputs["tree"]),
                "--out",
                str(tree_report_path),
                "--json",
            ],
            cwd=work_root,
        )
        tree_payload = json.loads(tree_result.stdout)
        assert tree_payload["status"] == "ok"
        assert tree_report_path.exists()

        signal_result = run_installed_cli(
            venv_root,
            [
                "comparative",
                "signal",
                str(copied_inputs["tree"]),
                str(copied_inputs["traits"]),
                "--trait",
                "value",
                "--permutations",
                "19",
                "--seed",
                "7",
                "--json",
            ],
            cwd=work_root,
        )
        signal_payload = json.loads(signal_result.stdout)
        assert signal_payload["status"] == "ok"
        assert signal_payload["metrics"]["taxon_count"] == 4
        assert signal_payload["metrics"]["signal_seed"] == 7
