from __future__ import annotations

import json
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


def test_cli_alignment_distance_quality_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-quality",
            str(fixture("example_alignment_distance_saturated.fasta")),
            "--model",
            "jukes-cantor",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["decision"] == "risky"
    assert payload["metrics"]["saturated_pair_count"] > 0


def test_cli_alignment_distance_saturation_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-saturation",
            str(fixture("example_alignment_distance_saturated.fasta")),
            "--model",
            "jc69",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["warning_pair_count"] == 2
    assert payload["metrics"]["blocking_warning_count"] == 2
    assert payload["metrics"]["blocks_tree_inference"] is True
    assert [
        (
            row["left_identifier"],
            row["right_identifier"],
            row["warning_kind"],
        )
        for row in payload["data"]["warning_rows"]
    ] == [
        ("A", "B", "undefined-corrected-distance"),
        ("B", "C", "infinite-corrected-distance"),
    ]


def test_cli_alignment_distance_ultrametricity_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-ultrametricity",
            str(fixture("example_alignment_distance.fasta")),
            "--model",
            "p-distance",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tested_triple_count"] == 4
    assert payload["metrics"]["violating_triple_count"] == 4
    assert payload["metrics"]["max_violation"] == 0.125
    assert payload["metrics"]["tolerance"] == 1e-6
    assert payload["metrics"]["ultrametric"] is False


def test_cli_alignment_distance_additivity_writes_governed_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "distance-additivity"
    exit_code = main(
        [
            "alignment",
            "distance-additivity",
            str(fixture("example_alignment_distance.fasta")),
            "--model",
            "p-distance",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tested_quartet_count"] == 1
    assert payload["metrics"]["violating_quartet_count"] == 1
    assert payload["metrics"]["max_violation"] == 0.25
    assert payload["metrics"]["additive"] is False
    table_lines = (
        (out_dir / "four_point_violations.tsv").read_text(encoding="utf-8").splitlines()
    )
    assert table_lines[1] == "A|B|C|D\t0.25\t1\t1.25\tA,B|C,D\t0.25"
    assert (out_dir / "run.json").exists()


def test_cli_alignment_distance_assumptions_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-assumptions",
            str(fixture("example_alignment_distance.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["ultrametric_compatible"] is False
    assert payload["metrics"]["upgma_violation_count"] > 0


def test_cli_distance_ultrametricity_json_output(capsys) -> None:
    exit_code = main(
        [
            "distance",
            "ultrametricity",
            str(fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tested_triple_count"] == 4
    assert payload["metrics"]["violating_triple_count"] == 4
    assert payload["metrics"]["max_violation"] == 3.0
    assert payload["metrics"]["tolerance"] == 1e-6
    assert payload["metrics"]["ultrametric"] is False


def test_cli_distance_additivity_writes_governed_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "distance-additivity"
    exit_code = main(
        [
            "distance",
            "additivity",
            str(fixture("example_distance_matrix_equal_row_sum_nonultrametric.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["tested_quartet_count"] == 1
    assert payload["metrics"]["violating_quartet_count"] == 1
    assert payload["metrics"]["max_violation"] == 6.0
    assert payload["metrics"]["additive"] is False
    table_lines = (
        (out_dir / "four_point_violations.tsv").read_text(encoding="utf-8").splitlines()
    )
    assert table_lines[1] == "A|B|C|D\t2\t8\t14\tA,B|C,D\t6"
    assert (out_dir / "run.json").exists()


def test_cli_distance_build_tree_supports_wpgma(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_wpgma_uneven_cluster.tsv")),
            "--method",
            "wpgma",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((((A:0.5,D:0.5)Inner1:2.25,E:2.75)Inner2:0.75,C:3.5)Inner3:0.5625,B:4.0625)Inner4;\n"
    )
    assert payload["metrics"]["method"] == "wpgma"


def test_cli_distance_build_tree_reports_missing_distance_imputation(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_missing_pair_four_taxon.tsv")),
            "--method",
            "neighbor-joining",
            "--missing-distance-policy",
            "triangle-bound",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").endswith(";\n")
    assert payload["metrics"]["missing_distance_policy"] == "triangle-bound"
    assert payload["metrics"]["imputed_pair_count"] == 1
    assert payload["data"]["missing_distance_policy_report"]["missing_pairs"] == ["A/C"]
    assert (
        payload["data"]["missing_distance_policy_report"]["imputed_rows"][0][
            "imputed_distance"
        ]
        == 5.0
    )


def test_cli_distance_build_tree_supports_single_linkage(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_single_linkage_chain.tsv")),
            "--method",
            "single-linkage",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((((A:0.5,B:0.5)Inner1:0.5,C:1)Inner2:0.5,D:1.5)Inner3:0.5,E:2)Inner4;\n"
    )
    assert payload["metrics"]["method"] == "single-linkage"


def test_cli_distance_build_tree_supports_complete_linkage(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(
                fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv")
            ),
            "--method",
            "complete-linkage",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "(((A:1,D:1)Inner2:2.5,E:3.5)Inner3:2,(B:1,C:1)Inner1:4.5)Inner4;\n"
    )
    assert payload["metrics"]["method"] == "complete-linkage"


def test_cli_alignment_build_tree_reports_missing_distance_imputation(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "alignment-distance-tree.nwk"
    exit_code = main(
        [
            "alignment",
            "build-tree",
            str(fixture("example_alignment_distance_missing_pair.fasta")),
            "--method",
            "neighbor-joining",
            "--model",
            "p-distance",
            "--missing-distance-policy",
            "nearest-valid",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8").endswith(";\n")
    assert payload["metrics"]["missing_distance_policy"] == "nearest-valid"
    assert payload["metrics"]["imputed_pair_count"] == 1
    assert payload["data"]["missing_distance_policy_report"]["missing_pairs"] == ["A/C"]
    assert (
        payload["data"]["missing_distance_policy_report"]["imputed_rows"][0][
            "imputed_distance"
        ]
        == 0.5
    )


def test_cli_distance_minimum_evolution_writes_fitted_tree(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "minimum-evolution-tree.nwk"
    exit_code = main(
        [
            "distance",
            "minimum-evolution",
            str(fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv")),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:1,B:2):3,C:4,(D:5,E:6):7);\n"
    )
    assert payload["metrics"]["criterion"] == "minimum-evolution"
    assert payload["metrics"]["minimum_evolution_score"] == 28.0


def test_cli_distance_fitch_margoliash_writes_fitted_tree(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "fitch-margoliash-tree.nwk"
    exit_code = main(
        [
            "distance",
            "fitch-margoliash",
            str(fixture("example_distance_matrix_fitch_margoliash_five_taxon.tsv")),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:1.186387219954,B:0.813612780046):3.461532690732,C:1.538467309268,(D:4.498139020967,E:8.501860979033):6.564266585263);\n"
    )
    assert payload["metrics"]["criterion"] == "fitch-margoliash"
    assert payload["metrics"]["weighting_power"] == 2.0
    assert payload["metrics"]["residual_sum_squares"] == 4.562310999055
    assert payload["metrics"]["weighted_residual_sum_squares"] == 0.020472720566


def test_cli_distance_ordinary_least_squares_writes_fitted_tree(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "ordinary-least-squares-tree.nwk"
    exit_code = main(
        [
            "distance",
            "ordinary-least-squares",
            str(fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv")),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:1,B:2):3,C:4,(D:5,E:6):7);\n"
    )
    assert payload["metrics"]["criterion"] == "ordinary-least-squares"
    assert payload["metrics"]["residual_sum_squares"] == 0.0
    assert payload["metrics"]["matrix_rank"] == 7
    assert payload["metrics"]["condition_number"] == 4.08308918215
    assert payload["metrics"]["negative_branch_count"] == 0


def test_cli_distance_nonnegative_least_squares_writes_fitted_tree(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "nonnegative-least-squares-tree.nwk"
    exit_code = main(
        [
            "distance",
            "nonnegative-least-squares",
            str(
                fixture(
                    "example_distance_matrix_ordinary_least_squares_negative_branch_five_taxon.tsv"
                )
            ),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:0.777777777778,B:3.777777777778):0,C:6.444444444444,(D:2.666666666667,E:7.333333333333):7.333333333333);\n"
    )
    assert payload["metrics"]["criterion"] == "nonnegative-least-squares"
    assert payload["metrics"]["residual_sum_squares"] == 57.777777777778
    assert payload["metrics"]["condition_number"] == 4.08308918215
    assert payload["metrics"]["active_constraint_count"] == 1


def test_cli_distance_patristic_residuals_writes_ranked_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "patristic-residuals"
    exit_code = main(
        [
            "distance",
            "patristic-residuals",
            str(fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv")),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    residual_lines = (
        (out_dir / "distance_residuals.tsv").read_text(encoding="utf-8").splitlines()
    )
    run_payload = json.loads((out_dir / "run.json").read_text(encoding="utf-8"))
    assert exit_code == 0
    assert residual_lines[:4] == [
        "left_identifier\tright_identifier\tobserved_distance\tfitted_distance\tresidual\tabsolute_residual\trank",
        "A\tD\t16\t36\t-20\t20\t1",
        "A\tC\t8\t27\t-19\t19\t2",
        "A\tE\t17\t36\t-19\t19\t3",
    ]
    assert payload["metrics"]["criterion"] == "patristic-residuals"
    assert payload["metrics"]["pair_count"] == 10
    assert payload["metrics"]["residual_sum_squares"] == 2626.0
    assert payload["metrics"]["max_absolute_residual"] == 20.0
    assert run_payload["pair_count"] == 10
    assert run_payload["residual_sum_squares"] == 2626.0
    assert run_payload["rows"][0]["left_identifier"] == "A"
    assert run_payload["rows"][0]["right_identifier"] == "D"
    assert run_payload["rows"][0]["rank"] == 1


def test_cli_distance_taxon_influence_writes_ranked_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "taxon-influence"
    exit_code = main(
        [
            "distance",
            "taxon-influence",
            str(
                fixture(
                    "example_distance_matrix_taxon_influence_missing_noisy_five_taxon.tsv"
                )
            ),
            str(fixture("example_tree_minimum_evolution_five_taxon.nwk")),
            "--method",
            "neighbor-joining",
            "--missing-distance-policy",
            "mean-impute",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    influence_lines = (
        (out_dir / "taxon_influence.tsv").read_text(encoding="utf-8").splitlines()
    )
    run_payload = json.loads((out_dir / "run.json").read_text(encoding="utf-8"))
    assert exit_code == 0
    assert influence_lines[:3] == [
        "influence_rank\ttaxon\tretained_taxa\traw_missing_pair_count\tbaseline_residual_sum_squares\tleave_one_out_residual_sum_squares\tresidual_sum_squares_improvement\tbaseline_rooted_robinson_foulds_distance\tleave_one_out_rooted_robinson_foulds_distance\trooted_robinson_foulds_improvement\tbaseline_rooted_normalized_robinson_foulds\tleave_one_out_rooted_normalized_robinson_foulds\trooted_normalized_robinson_foulds_improvement\ttopology_improved\tresidual_improved",
        "1\tB\tA|C|D|E\t0\t20.6342592593\t1.21\t19.4242592593\t4\t2\t2\t1\t1\t0\ttrue\ttrue",
        "2\tA\tB|C|D|E\t1\t20.6342592593\t2.25\t18.3842592593\t4\t2\t2\t1\t1\t0\ttrue\ttrue",
    ]
    assert payload["metrics"]["criterion"] == "distance-taxon-influence"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["baseline_residual_sum_squares"] == 20.634259259263
    assert payload["metrics"]["baseline_rooted_robinson_foulds_distance"] == 4
    assert payload["metrics"]["top_taxon"] == "B"
    assert run_payload["baseline_residual_sum_squares"] == 20.634259259263
    assert run_payload["rows"][0]["taxon"] == "B"
    assert run_payload["rows"][0]["raw_missing_pair_count"] == 0
    assert run_payload["rows"][1]["taxon"] == "A"
    assert run_payload["rows"][1]["raw_missing_pair_count"] == 1


def test_cli_distance_taxon_jackknife_writes_rebuilt_tree_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "taxon-jackknife"
    exit_code = main(
        [
            "distance",
            "taxon-jackknife",
            str(
                fixture(
                    "example_distance_matrix_taxon_influence_missing_noisy_five_taxon.tsv"
                )
            ),
            "--method",
            "neighbor-joining",
            "--missing-distance-policy",
            "mean-impute",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    jackknife_lines = (
        (out_dir / "taxon_jackknife.tsv").read_text(encoding="utf-8").splitlines()
    )
    run_payload = json.loads((out_dir / "run.json").read_text(encoding="utf-8"))
    assert exit_code == 0
    assert (out_dir / "baseline_tree.nwk").read_text(encoding="utf-8") == (
        "(((A:0.722222222222223,C:7.27777777777778)Inner1:5.58333333333333,E:5.58333333333333)Inner2:1.91666666666667,B:-3.75,D:4.75)Inner3;\n"
    )
    assert jackknife_lines[:3] == [
        "removed_taxon\tretained_taxa\trooted_robinson_foulds_distance\trooted_normalized_robinson_foulds\tpruned_baseline_residual_sum_squares\trebuilt_residual_sum_squares\tresidual_sum_squares_change\taffected_clades\ttopology_changed\tpruned_baseline_tree_newick\trebuilt_tree_newick",
        "A\tB|C|D|E\t2\t1\t6.87962962963\t2.25\t4.62962962963\tB|C;C|E\ttrue\t(B:-3.75,(C:12.8611111111111,E:5.58333333333333)Inner2:1.91666666666667,D:4.75)Inner3;\t((B:-3.25,C:13.25)Inner1:0.75,D:4.25,E:6.75)Inner2;",
        "B\tA|C|D|E\t1\t0.333333333333\t21.902808642\t1.21\t20.692808642\tA|C|E\ttrue\t(((A:0.722222222222223,C:7.27777777777778)Inner1:5.58333333333333,E:5.58333333333333)Inner2:1.91666666666667,D:4.75)Inner3;\t((A:1.95,C:6.05)Inner1:7.95,D:5.55,E:5.45)Inner2;",
    ]
    assert payload["metrics"]["criterion"] == "distance-taxon-jackknife"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["baseline_residual_sum_squares"] == 20.634259259263
    assert payload["metrics"]["topology_changed_taxon_count"] == 3
    assert run_payload["baseline_residual_sum_squares"] == 20.634259259263
    assert run_payload["rows"][0]["removed_taxon"] == "A"
    assert run_payload["rows"][0]["affected_clades"] == ["B|C", "C|E"]
    assert run_payload["rows"][2]["removed_taxon"] == "C"
    assert run_payload["rows"][2]["affected_clades"] == []


def test_cli_distance_method_comparison_writes_tree_rf_and_warning_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "method-comparison"
    exit_code = main(
        [
            "distance",
            "method-comparison",
            str(fixture("example_distance_matrix_bionj_noisy.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    score_lines = (
        (out_dir / "method_scores.tsv").read_text(encoding="utf-8").splitlines()
    )
    rf_lines = (out_dir / "rf_matrix.tsv").read_text(encoding="utf-8").splitlines()
    warning_lines = (
        (out_dir / "assumption_warnings.tsv").read_text(encoding="utf-8").splitlines()
    )
    run_payload = json.loads((out_dir / "run.json").read_text(encoding="utf-8"))
    assert exit_code == 0
    assert (out_dir / "neighbor-joining.nwk").read_text(encoding="utf-8") == (
        "((A:1.125,(B:1,C:2)Inner1:5.375)Inner2:5.375,D:4.625,E:-2.625)Inner3;\n"
    )
    assert (out_dir / "bionj.nwk").read_text(encoding="utf-8") == (
        "((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;\n"
    )
    assert score_lines[:3] == [
        "method\tpatristic_residual_sum_squares\tbalanced_minimum_evolution_score\tordinary_least_squares_residual_sum_squares\tordinary_least_squares_negative_branch_count\ttree_newick",
        "neighbor-joining\t168.09375\t28.25\t165.75\t1\t((A:1.125,(B:1,C:2)Inner1:5.375)Inner2:5.375,D:4.625,E:-2.625)Inner3;",
        "bionj\t194.315733333\t28.25\t165.75\t1\t((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;",
    ]
    assert rf_lines == [
        "method\tneighbor-joining\tbionj\tupgma\twpgma",
        "neighbor-joining\t0\t0\t3\t3",
        "bionj\t0\t0\t3\t3",
        "upgma\t3\t3\t0\t0",
        "wpgma\t3\t3\t0\t0",
    ]
    assert warning_lines[:4] == [
        "warning_rank\tscope\tmethod\twarning",
        "1\tmatrix\t\tdistance matrix violates triangle inequality for one or more taxon triples",
        "2\tmatrix\t\tpairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated",
        "3\tmethod\tbionj\tbionj remains a distance-summary method rather than a full likelihood inference",
    ]
    assert payload["metrics"]["criterion"] == "distance-method-comparison"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["method_count"] == 4
    assert payload["metrics"]["warning_count"] == 8
    assert run_payload["compared_methods"] == [
        "neighbor-joining",
        "bionj",
        "upgma",
        "wpgma",
    ]
    assert run_payload["rf_rows"][0]["left_method"] == "neighbor-joining"
    assert run_payload["rf_rows"][0]["right_method"] == "bionj"
    assert run_payload["rows"][2]["method"] == "upgma"
    assert run_payload["warning_rows"][-1]["method"] == "wpgma"


def test_cli_distance_bme_nni_search_writes_governed_artifacts(
    tmp_path: Path, capsys
) -> None:
    out_dir = tmp_path / "bme-nni-search"
    exit_code = main(
        [
            "distance",
            "bme-nni-search",
            str(
                fixture(
                    "example_distance_matrix_balanced_minimum_evolution_nni_five_taxon.tsv"
                )
            ),
            "--start-method",
            "bionj",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert (out_dir / "start_tree.nwk").read_text(encoding="utf-8") == (
        "(((A,E)Inner1,D)Inner2,B,C)Inner3;\n"
    )
    assert (out_dir / "final_tree.nwk").read_text(encoding="utf-8") == (
        "(((A,D)Inner1,E)Inner2,B,C)Inner3;\n"
    )
    assert (out_dir / "search_trace.tsv").read_text(encoding="utf-8") == (
        "event_index\tevent_kind\titeration\tscore_before\tscore_after\tscore_delta\ttree_before_newick\ttree_after_newick\tpivot_branch_id\tsibling_clade_id\texchanged_clade_id\tstopping_reason\n"
        "1\tstart\t0\t\t34\t\t\t(((A,E)Inner1,D)Inner2,B,C)Inner3;\t\t\t\t\n"
        "2\taccepted-move\t1\t34\t33.75\t-0.25\t(((A,E)Inner1,D)Inner2,B,C)Inner3;\t(((A,D)Inner1,E)Inner2,B,C)Inner3;\tA|E\tD\tE\t\n"
        "3\tfinal\t1\t\t33.75\t\t\t(((A,D)Inner1,E)Inner2,B,C)Inner3;\t\t\t\tno-improving-neighbor\n"
    )
    run_payload = json.loads((out_dir / "run.json").read_text(encoding="utf-8"))
    assert run_payload["algorithm"] == "balanced-minimum-evolution-nni-search"
    assert run_payload["start_method"] == "bionj"
    assert run_payload["accepted_move_count"] == 1
    assert run_payload["final_score"] == 33.75
    assert payload["metrics"]["algorithm"] == "balanced-minimum-evolution-nni-search"
    assert payload["metrics"]["start_method"] == "bionj"
    assert payload["metrics"]["accepted_move_count"] == 1
    assert payload["metrics"]["start_score"] == 34.0
    assert payload["metrics"]["final_score"] == 33.75
    assert payload["metrics"]["stopping_reason"] == "no-improving-neighbor"


def test_cli_alignment_bootstrap_tree_writes_outputs(tmp_path: Path, capsys) -> None:
    support_path = tmp_path / "support.tsv"
    tree_set_path = tmp_path / "bootstrap.trees"
    draws_path = tmp_path / "bootstrap-draws.tsv"
    exit_code = main(
        [
            "alignment",
            "bootstrap-tree",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "neighbor-joining",
            "--replicates",
            "5",
            "--seed",
            "5",
            "--support-out",
            str(support_path),
            "--tree-set-out",
            str(tree_set_path),
            "--draws-out",
            str(draws_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert support_path.exists()
    assert tree_set_path.exists()
    assert draws_path.exists()
    assert support_path.read_text(encoding="utf-8").splitlines() == [
        "clade\ttree_count\tfrequency",
        "A|B\t5\t1",
    ]
    assert tree_set_path.read_text(encoding="utf-8").splitlines() == [
        "((A:0.0625,B:0.0625)Inner1:0.5625,C:0.0625,D:0.0625)Inner2;",
        "((A:0.0625,B:0.0625)Inner1:0.5625,C:0.0625,D:0.0625)Inner2;",
        "((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;",
        "((A:0,B:0)Inner1:1,C:0,D:0)Inner2;",
        "((A:0,B:0)Inner1:0.5,C:0,D:0)Inner2;",
    ]
    assert draws_path.read_text(encoding="utf-8").splitlines() == [
        "replicate_index\tsampled_site_indices\ttree_newick",
        "1\t4,5,0,7,3,0,2,1\t((A:0.0625,B:0.0625)Inner1:0.5625,C:0.0625,D:0.0625)Inner2;",
        "2\t5,7,3,6,1,3,0,3\t((A:0.0625,B:0.0625)Inner1:0.5625,C:0.0625,D:0.0625)Inner2;",
        "3\t6,4,2,6,2,1,2,7\t((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;",
        "4\t2,2,0,0,3,3,2,2\t((A:0,B:0)Inner1:1,C:0,D:0)Inner2;",
        "5\t4,5,3,3,2,3,6,4\t((A:0,B:0)Inner1:0.5,C:0,D:0)Inner2;",
    ]
    assert payload["metrics"]["replicate_count"] == 5


def test_cli_alignment_distance_support_summary_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-support-summary",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "neighbor-joining",
            "--replicates",
            "5",
            "--seed",
            "3",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["replicates"] == 5
    assert payload["metrics"]["clade_count"] > 0


def test_cli_alignment_distance_models_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-models",
            str(fixture("example_alignment_distance.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    model_rows = payload["data"]["rows"]
    assert payload["metrics"]["model_count"] == len(model_rows)
    assert {row["model"] for row in model_rows} >= {
        "p-distance",
        "jukes-cantor",
        "kimura-2-parameter",
        "felsenstein-81",
        "tamura-nei-93",
    }


def test_cli_alignment_distance_gap_sensitivity_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-gap-sensitivity",
            str(fixture("example_alignment_distance_gaps.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["changed_pair_count"] > 0


def test_cli_alignment_distance_maturity_json_output(capsys) -> None:
    exit_code = main(
        [
            "alignment",
            "distance-maturity",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "neighbor-joining",
            "--replicates",
            "5",
            "--seed",
            "3",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["check_count"] > 0
    assert payload["metrics"]["decision"] in {
        "production_candidate",
        "validated_with_limits",
    }


def test_cli_alignment_build_tree_supports_bionj_json(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "alignment",
            "build-tree",
            str(fixture("example_alignment_distance.fasta")),
            "--method",
            "bionj",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;\n"
    )
    assert payload["metrics"]["method"] == "bionj"


def test_cli_distance_build_tree_supports_bionj(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "distance-tree.nwk"
    exit_code = main(
        [
            "distance",
            "build-tree",
            str(fixture("example_distance_matrix_bionj_noisy.tsv")),
            "--method",
            "bionj",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert output_path.read_text(encoding="utf-8") == (
        "((A:2,(B:1,C:2)Inner1:5.66666666666667)Inner2:4.5,D:6.02,E:-4.02)Inner3;\n"
    )
    assert payload["metrics"]["method"] == "bionj"


def test_cli_distance_reference_json_output(capsys) -> None:
    exit_code = main(["distance", "reference", "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["all_passed"] is True
    observations = payload["data"]["observations"]
    assert payload["metrics"]["case_count"] == len(observations)
    assert {row["case"] for row in observations} >= {
        "dna-p-distance",
        "dna-jukes-cantor",
        "dna-kimura-2-parameter",
        "dna-felsenstein-81",
        "dna-tamura-nei-93",
        "protein-p-distance",
        "ambiguity-partial-match",
        "ambiguity-strict-mismatch",
        "ambiguity-report-only",
    }


def test_cli_distance_assumptions_json_output(capsys) -> None:
    exit_code = main(
        [
            "distance",
            "assumptions",
            str(fixture("example_distance_matrix_nonultrametric.tsv")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["ultrametric_compatible"] is False
    assert payload["metrics"]["upgma_violation_count"] > 0


def test_cli_distance_quality_json_output(capsys) -> None:
    exit_code = main(
        ["distance", "quality", str(fixture("example_distance_matrix.tsv")), "--json"]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["saturation_audit_scale"] == "unit-interval-like"
    assert payload["metrics"]["low_information_pair_count"] == 3
