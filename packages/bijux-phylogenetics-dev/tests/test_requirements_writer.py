from __future__ import annotations

from pathlib import Path

from pytest import MonkeyPatch

from bijux_phylogenetics_dev.quality.requirements_writer import main


def test_requirements_writer_merges_package_and_companion_dev_groups(
    tmp_path: Path, monkeypatch: MonkeyPatch
) -> None:
    package_dir = tmp_path / "packages"
    app_dir = package_dir / "bijux-phylogenetics"
    dev_dir = package_dir / "bijux-phylogenetics-dev"
    app_dir.mkdir(parents=True)
    dev_dir.mkdir(parents=True)

    (app_dir / "pyproject.toml").write_text(
        """
[project]
dependencies = ["biopython>=1.87,<2.0"]

[project.optional-dependencies]
dev = ["pytest>=9.0.3,<10.0", "ruff>=0.13.0,<1.0"]
""".strip()
        + "\n",
        encoding="utf-8",
    )
    (dev_dir / "pyproject.toml").write_text(
        """
[project]
dependencies = ["PyYAML>=6.0,<7.0"]

[project.optional-dependencies]
dev = ["ruff>=0.13.0,<1.0", "mypy>=1.18.2,<3.0"]
""".strip()
        + "\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "requirements.txt"

    monkeypatch.chdir(tmp_path)
    exit_code = main(
        [
            "--pyproject",
            str(app_dir / "pyproject.toml"),
            "--group",
            "dev",
            "--optional-group",
            "dev",
            "--output",
            str(output_path),
        ]
    )

    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "biopython>=1.87,<2.0",
        "pytest>=9.0.3,<10.0",
        "ruff>=0.13.0,<1.0",
        "PyYAML>=6.0,<7.0",
        "mypy>=1.18.2,<3.0",
    ]
