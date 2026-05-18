from __future__ import annotations

import json
import math
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_comparative_brownian_cli_reports_parameters(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "brownian",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert math.isclose(payload["metrics"]["root_state"], 2.8055555543209874)
    assert math.isclose(payload["metrics"]["rate"], 4.774305191647407)


def test_comparative_compare_models_cli_reports_selected_model(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "compare-models",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["better_model"] == "brownian"


def test_comparative_model_comparison_package_cli_writes_review_bundle(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "comparative-model-figure-package"

    exit_code = main(
        [
            "comparative",
            "model-comparison-package",
            str(fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk")),
            str(fixture("example_traits_phytools_signal_twenty_four_taxa.tsv")),
            "--trait",
            "signal_strong",
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["selected_model"] == "brownian"
    assert payload["metrics"]["support_distinct"] is True
    assert payload["metrics"]["plotted_model_count"] == 2
    assert payload["metrics"]["rendered_parameter_count"] == 5
    assert payload["metrics"]["rendered_fit_row_count"] == 2
    assert len(payload["outputs"]) == 13
    assert (output_dir / "model-comparison-criteria.svg").exists()
    assert (output_dir / "model-comparison-likelihood.svg").exists()
    assert (output_dir / "model-comparison-parameters.svg").exists()
    assert (output_dir / "model-comparison-fit-summary.svg").exists()
    assert (output_dir / "model-comparison-review.html").exists()
    assert (output_dir / "model-comparison-package.manifest.json").exists()
    assert (output_dir / "figure-reproducibility.manifest.json").exists()


def test_comparative_validate_reference_cli_reports_pass(capsys) -> None:
    exit_code = main(["comparative", "validate-reference", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["all_passed"] is True
    assert payload["metrics"]["interval_audit_passed"] is True
    assert payload["metrics"]["identifiability_audit_passed"] is True


def test_comparative_sensitivity_cli_reports_influential_taxa(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "sensitivity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--trait",
            "response",
            "--model",
            "brownian",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["influential_taxa"] == 3


def test_comparative_maturity_cli_reports_residual_surfaces(capsys) -> None:
    exit_code = main(
        [
            "comparative",
            "maturity",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_comparative.tsv")),
            "--response",
            "response",
            "--predictors",
            "predictor_one",
            "--lambda-value",
            "1.0",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["residual_surface_count"] == 3
    assert payload["metrics"]["reference_validation_passed"] is True
