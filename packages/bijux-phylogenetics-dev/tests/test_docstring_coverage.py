from __future__ import annotations

from pathlib import Path

from pytest import CaptureFixture, MonkeyPatch

from bijux_phylogenetics_dev.quality.docstring_coverage import (
    analyze_file,
    load_config,
    main,
)


def test_docstring_coverage_honors_package_policy(tmp_path: Path) -> None:
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(
        """
[tool.interrogate]
fail-under = 75
ignore-init-method = true
ignore-module = true
ignore-private = true
ignore-semiprivate = true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    source = tmp_path / "module.py"
    source.write_text(
        '''
"""module docstring"""

def documented() -> None:
    """public function"""


def undocumented() -> None:
    pass


def _helper() -> None:
    pass


class Example:
    """class docstring"""

    def __init__(self) -> None:
        pass

    def method(self) -> None:
        """documented method"""
'''.strip()
        + "\n",
        encoding="utf-8",
    )

    config = load_config(pyproject, None)
    coverage = analyze_file(source, config)

    assert coverage.documented == 3
    assert coverage.total == 4
    assert coverage.percentage == 75.0


def test_docstring_coverage_main_reports_failure(
    tmp_path: Path, capsys: CaptureFixture[str], monkeypatch: MonkeyPatch
) -> None:
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(
        """
[tool.interrogate]
fail-under = 90
ignore-module = true
ignore-private = true
ignore-semiprivate = true
""".strip()
        + "\n",
        encoding="utf-8",
    )
    source_dir = tmp_path / "src"
    source_dir.mkdir()
    (source_dir / "sample.py").write_text(
        '''
def documented() -> None:
    """present"""


def undocumented() -> None:
    pass
'''.strip()
        + "\n",
        encoding="utf-8",
    )

    monkeypatch.chdir(tmp_path)
    exit_code = main(["--pyproject", "pyproject.toml", "--fail-under", "90", "src"])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert "FILE: src/sample.py | 50.0%" in captured.out
    assert "RESULT: FAILED (50.0%)" in captured.out
