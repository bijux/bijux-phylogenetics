from __future__ import annotations

import pytest

from bijux_phylogenetics.command_line.bootstrap import main


def test_comparative_covariance_audit_cli_help_lists_supported_analyses(
    capsys,
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["comparative", "covariance-audit", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--analysis {pgls,brownian-trait,ou-trait}" in captured.out
    assert "--summary-out SUMMARY_OUT" in captured.out
