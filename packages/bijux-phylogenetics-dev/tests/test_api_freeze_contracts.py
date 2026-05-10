from __future__ import annotations

from pathlib import Path

from _pytest.capture import CaptureFixture

from bijux_phylogenetics_dev.api.freeze_contracts import run


def test_api_freeze_skips_cleanly_when_repository_has_no_checked_in_schemas(
    tmp_path: Path,
    capsys: CaptureFixture[str],
) -> None:
    repo_root = tmp_path / "demo-repo"
    repo_root.mkdir()

    assert run(repo_root) == 0
    captured = capsys.readouterr()
    assert captured.out.strip() == "No checked-in OpenAPI schemas found; skipping."
    assert captured.err == ""
