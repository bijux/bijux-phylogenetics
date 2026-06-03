from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.benchmark.macroevolution import (
    RealDatasetMacroevolutionAlignmentReviewRow,
    RealDatasetMacroevolutionBenchmarkBundle,
    RealDatasetMacroevolutionBenchmarkDemoResult,
    RealDatasetMacroevolutionBenchmarkReport,
    RealDatasetMacroevolutionModelRow,
    RealDatasetMacroevolutionParityRow,
    RealDatasetMacroevolutionSummaryRow,
)
import bijux_phylogenetics.command_line as command_line_api
from bijux_phylogenetics.command_line import main
import bijux_phylogenetics.command_line.demo as demo_command_module
from bijux_phylogenetics.datasets.central_european_seashore_flora import (
    CentralEuropeanSeashoreFloraDataset,
    CentralEuropeanSeashoreFloraDatasetExportResult,
)


def test_cli_benchmark_real_dataset_macroevolution_reports_metrics(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    dataset = CentralEuropeanSeashoreFloraDataset(
        dataset_id="central_european_seashore_flora",
        label="Central European seashore flora dataset",
        dataset_root=Path("/tmp/central-european-seashore-flora"),
        tree_path=Path("/tmp/central-european-seashore-flora/tree.nwk"),
        traits_path=Path("/tmp/central-european-seashore-flora/traits.csv"),
        reference_output_root=Path("/tmp/central-european-seashore-flora/expected"),
        taxon_column="species",
        taxon_count=42,
        continuous_traits=("seed_mass",),
        categorical_traits=("lifeform",),
        workflow_continuous_trait="seed_mass",
        workflow_pgls_predictor="plant_height",
        workflow_discrete_trait="lifeform",
        workflow_clade_trait="lifeform",
        source_summary="governed test fixture",
    )
    report = RealDatasetMacroevolutionBenchmarkReport(
        dataset=dataset,
        provenance_citation="fixture",
        provenance_doi="10.5061/dryad.0st06f0",
        summary_rows=[
            RealDatasetMacroevolutionSummaryRow(
                surface_id="seed-mass-native-model-table",
                trait="seed_mass",
                trait_kind="continuous",
                review_scope="native-model-table",
                bijux_selected_model="ornstein-uhlenbeck",
                geiger_selected_model="ornstein-uhlenbeck",
                selection_matches_geiger=True,
                bijux_selected_model_akaike_weight=0.8,
                geiger_selected_model_akaike_weight=0.8,
                stable_conclusion_supported=False,
                aligned_taxa_count=42,
                dropped_tree_taxon_count=0,
                dropped_trait_taxon_count=0,
                dropped_missing_value_taxon_count=0,
                biological_interpretation="fixture",
                notes=["fixture"],
            )
        ],
        model_rows=[
            RealDatasetMacroevolutionModelRow(
                surface_id="seed-mass-native-model-table",
                trait="seed_mass",
                trait_kind="continuous",
                model="ornstein-uhlenbeck",
                bijux_rank=1,
                geiger_rank=1,
                bijux_selected=True,
                geiger_selected=True,
                bijux_parameter_count=3,
                geiger_parameter_count=3,
                bijux_log_likelihood=-10.0,
                geiger_log_likelihood=-10.0,
                bijux_aic=26.0,
                geiger_aic=26.0,
                bijux_aicc=27.0,
                geiger_aicc=27.0,
                bijux_akaike_weight=0.8,
                geiger_akaike_weight=0.8,
                bijux_parameter_name="alpha",
                geiger_parameter_name="alpha",
                bijux_parameter_value=0.02,
                geiger_parameter_value=0.02,
                bijux_rate=1.0,
                geiger_rate=1.0,
                bijux_root_state=5.0,
                geiger_root_state=5.0,
                notes=["fixture"],
            )
        ],
        alignment_review_rows=[
            RealDatasetMacroevolutionAlignmentReviewRow(
                surface_id="seed-mass-alignment-review",
                trait="seed_mass",
                model="ornstein-uhlenbeck",
                original_tree_taxa=42,
                original_trait_taxa=42,
                aligned_taxa_count=40,
                dropped_tree_taxa=["Triglochin_maritimum"],
                dropped_trait_taxa=["unmatched_review_taxon"],
                dropped_missing_value_taxa=["Juncus_maritimus"],
                geiger_overlap_taxa=41,
                geiger_usable_taxa=40,
                notes=["fixture"],
            )
        ],
        parity_rows=[
            RealDatasetMacroevolutionParityRow(
                surface_id="seed-mass-native-model-table",
                trait="seed_mass",
                model="ornstein-uhlenbeck",
                comparison_scope="native-model-table",
                bijux_log_likelihood=-10.0,
                geiger_log_likelihood=-10.0,
                absolute_log_likelihood_delta=0.0,
                bijux_aicc=27.0,
                geiger_aicc=27.0,
                absolute_aicc_delta=0.0,
                bijux_parameter_name="alpha",
                geiger_parameter_name="alpha",
                bijux_parameter_value=0.02,
                geiger_parameter_value=0.02,
                absolute_parameter_delta=0.0,
                within_log_likelihood_tolerance=True,
                within_aicc_tolerance=True,
                within_parameter_tolerance=True,
                notes=["fixture"],
            )
        ],
        limitations=["fixture"],
    )
    monkeypatch.setattr(
        command_line_api,
        "benchmark_real_dataset_macroevolution",
        lambda: report,
    )

    exit_code = main(["benchmark", "real-dataset-macroevolution", "--json"])
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["command"] == "benchmark"
    assert payload["metrics"]["summary_row_count"] == 1
    assert payload["metrics"]["model_row_count"] == 1
    assert payload["metrics"]["alignment_review_row_count"] == 1
    assert payload["metrics"]["parity_row_count"] == 1
    assert payload["data"]["dataset"]["dataset_id"] == "central_european_seashore_flora"


def test_cli_demo_real_dataset_macroevolution_reports_metrics(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    expected_root = tmp_path / "expected"
    expected_root.mkdir()
    (expected_root / "expected.tsv").write_text("ok\n", encoding="utf-8")
    dataset = CentralEuropeanSeashoreFloraDataset(
        dataset_id="central_european_seashore_flora",
        label="Central European seashore flora dataset",
        dataset_root=tmp_path / "dataset-source",
        tree_path=tmp_path / "dataset-source" / "tree.nwk",
        traits_path=tmp_path / "dataset-source" / "traits.csv",
        reference_output_root=expected_root,
        taxon_column="species",
        taxon_count=42,
        continuous_traits=("seed_mass",),
        categorical_traits=("lifeform",),
        workflow_continuous_trait="seed_mass",
        workflow_pgls_predictor="plant_height",
        workflow_discrete_trait="lifeform",
        workflow_clade_trait="lifeform",
        source_summary="fixture",
    )
    export = CentralEuropeanSeashoreFloraDatasetExportResult(
        output_root=tmp_path / "dataset",
        readme_path=tmp_path / "dataset" / "README.md",
        tree_path=tmp_path / "dataset" / "tree.nwk",
        traits_path=tmp_path / "dataset" / "traits.csv",
        expected_output_root=expected_root,
    )
    bundle = RealDatasetMacroevolutionBenchmarkBundle(
        output_root=tmp_path / "benchmark",
        review_traits_path=tmp_path / "benchmark" / "alignment-review-traits.csv",
        summary_path=tmp_path / "benchmark" / "benchmark-summary.tsv",
        model_table_path=tmp_path / "benchmark" / "model-table.tsv",
        alignment_review_path=tmp_path / "benchmark" / "alignment-review.tsv",
        parity_table_path=tmp_path / "benchmark" / "geiger-parity.tsv",
        geiger_reference_path=tmp_path / "benchmark" / "geiger-reference.tsv",
    )
    for path in [
        export.readme_path,
        export.tree_path,
        export.traits_path,
        bundle.review_traits_path,
        bundle.summary_path,
        bundle.model_table_path,
        bundle.alignment_review_path,
        bundle.parity_table_path,
        bundle.geiger_reference_path,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n", encoding="utf-8")
    overview_path = tmp_path / "README.md"
    overview_path.write_text("fixture\n", encoding="utf-8")
    result = RealDatasetMacroevolutionBenchmarkDemoResult(
        output_root=tmp_path,
        dataset=dataset,
        dataset_export=export,
        benchmark_bundle=bundle,
        overview_path=overview_path,
    )
    monkeypatch.setattr(
        demo_command_module,
        "run_real_dataset_macroevolution_benchmark_demo",
        lambda destination: result,
    )

    exit_code = main(
        ["demo", "real-dataset-macroevolution", "--out", str(tmp_path), "--json"]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["command"] == "demo"
    assert payload["metrics"]["dataset_taxon_count"] == 42
    assert payload["metrics"]["summary_row_count"] == 4
    assert payload["metrics"]["native_model_row_count"] == 8
    assert payload["metrics"]["alignment_review_row_count"] == 2
    assert payload["metrics"]["parity_row_count"] == 10
