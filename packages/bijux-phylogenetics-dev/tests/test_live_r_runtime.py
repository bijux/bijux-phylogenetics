from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics_dev.quality import live_r_runtime


def test_ensure_live_r_test_packages_reports_existing_runtime(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        live_r_runtime,
        "_resolve_rscript_executable",
        lambda executable: executable or "/usr/local/bin/Rscript",
    )
    monkeypatch.setattr(
        live_r_runtime,
        "query_live_r_test_packages",
        lambda repo_root, **_: {
            package_name: True
            for package_name in live_r_runtime.DEFAULT_LIVE_R_TEST_PACKAGES
        },
    )
    install_calls: list[tuple[str, ...]] = []
    monkeypatch.setattr(
        live_r_runtime,
        "install_live_r_test_packages",
        lambda repo_root, *, package_names, **__: install_calls.append(package_names),
    )

    report = live_r_runtime.ensure_live_r_test_packages(tmp_path)

    assert report["status"] == "ready"
    assert report["installed_now"] == []
    assert report["already_available"] == list(
        live_r_runtime.DEFAULT_LIVE_R_TEST_PACKAGES
    )
    assert install_calls == []


def test_ensure_live_r_test_packages_installs_missing_runtime_packages(
    monkeypatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        live_r_runtime,
        "_resolve_rscript_executable",
        lambda executable: executable or "/usr/local/bin/Rscript",
    )
    package_state = {
        package_name: package_name != "phylolm"
        for package_name in live_r_runtime.DEFAULT_LIVE_R_TEST_PACKAGES
    }

    def fake_query(repo_root: Path, **_: object) -> dict[str, bool]:
        return dict(package_state)

    install_calls: list[tuple[str, ...]] = []

    def fake_install(
        repo_root: Path,
        *,
        package_names: tuple[str, ...],
        **__: object,
    ) -> None:
        install_calls.append(package_names)
        for package_name in package_names:
            package_state[package_name] = True

    monkeypatch.setattr(
        live_r_runtime,
        "query_live_r_test_packages",
        fake_query,
    )
    monkeypatch.setattr(
        live_r_runtime,
        "install_live_r_test_packages",
        fake_install,
    )

    report = live_r_runtime.ensure_live_r_test_packages(tmp_path)

    assert report["status"] == "ready"
    assert report["installed_now"] == ["phylolm"]
    assert "phylolm" not in report["already_available"]
    assert install_calls == [("phylolm",)]
