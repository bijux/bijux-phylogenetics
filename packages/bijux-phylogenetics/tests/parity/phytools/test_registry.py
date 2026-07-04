from __future__ import annotations

from bijux_phylogenetics.parity import list_phytools_parity_cases


def test_list_phytools_parity_cases_returns_governed_registry() -> None:
    cases = list_phytools_parity_cases()

    assert [case.case_id for case in cases] == [
        "phylosig-lambda-non-ultrametric-strong-signal-twenty-four-taxa",
        "phylosig-lambda-weak-signal-twenty-four-taxa",
        "phylosig-k-strong-signal-twenty-four-taxa",
        "phylosig-k-weak-signal-twenty-four-taxa",
        "fitmk-er-binary-twenty-four-taxa",
        "fitmk-er-multistate-twenty-four-taxa",
        "fitmk-er-binary-missing-twenty-four-taxa",
        "fitmk-er-multistate-missing-twenty-four-taxa",
        "fitmk-sym-multistate-twenty-four-taxa",
        "fitmk-sym-multistate-missing-twenty-four-taxa",
        "fitmk-ard-binary-twenty-four-taxa",
        "fitmk-ard-multistate-twenty-four-taxa",
        "fitmk-ard-binary-missing-twenty-four-taxa",
        "fitmk-ard-multistate-missing-twenty-four-taxa",
        "simmap-er-binary-twenty-four-taxa",
        "simmap-er-multistate-twenty-four-taxa",
        "simmap-er-binary-missing-twenty-four-taxa",
        "simmap-sym-multistate-twenty-four-taxa",
        "simmap-sym-multistate-missing-twenty-four-taxa",
        "simmap-ard-binary-twenty-four-taxa",
        "simmap-ard-multistate-twenty-four-taxa",
        "simmap-ard-binary-missing-twenty-four-taxa",
        "simmap-ard-multistate-missing-twenty-four-taxa",
        "count-simmap-er-binary-twenty-four-taxa",
        "count-simmap-er-multistate-twenty-four-taxa",
        "count-simmap-sym-multistate-twenty-four-taxa",
        "count-simmap-er-binary-missing-twenty-four-taxa",
        "describe-simmap-er-binary-twenty-four-taxa",
        "describe-simmap-er-multistate-twenty-four-taxa",
        "describe-simmap-sym-multistate-twenty-four-taxa",
        "describe-simmap-er-binary-missing-twenty-four-taxa",
        "density-map-er-binary-twenty-four-taxa",
        "density-map-er-binary-missing-twenty-four-taxa",
        "sim-history-binary-no-change-example-tree",
        "sim-history-binary-high-rate-example-tree",
        "sim-history-multistate-no-change-six-taxa",
        "sim-history-multistate-high-rate-six-taxa",
        "sim-history-binary-root-prior-example-tree",
        "rerooting-er-binary-twenty-four-taxa",
        "rerooting-er-multistate-twenty-four-taxa",
        "rerooting-er-binary-missing-twenty-four-taxa",
        "rerooting-sym-multistate-twenty-four-taxa",
        "rerooting-sym-multistate-missing-twenty-four-taxa",
        "fast-anc-strong-signal-twenty-four-taxa",
        "fast-anc-weak-signal-twenty-four-taxa",
        "fast-anc-non-ultrametric-strong-signal-twenty-four-taxa",
        "fast-anc-missing-values-twenty-four-taxa",
        "anc-ml-strong-signal-twenty-four-taxa",
        "anc-ml-weak-signal-twenty-four-taxa",
        "anc-ml-non-ultrametric-strong-signal-twenty-four-taxa",
        "anc-ml-missing-values-twenty-four-taxa",
        "fastbm-example-tree-low-variance",
        "fastbm-example-tree-root-shift-high-variance",
        "fastbm-six-taxa-root-shift",
        "simcorrs-example-tree-low-correlation",
        "simcorrs-example-tree-negative-correlation-root-shift",
        "simcorrs-six-taxa-three-trait",
        "pgls-sey-brownian-continuous-four-taxa",
        "pgls-sey-brownian-categorical-eight-taxa",
        "pgls-sey-brownian-interaction-eight-taxa",
        "phyl-resid-bm-allometry-six-taxa",
        "phyl-resid-lambda-allometry-six-taxa",
        "phyl-resid-lambda-missing-six-taxa",
        "phyl-anova-group-effect-six-taxa",
        "phyl-anova-group-effect-missing-six-taxa",
    ]
    assert cases[0].function_name == "phytools::phylosig(method='lambda')"
    assert cases[1].function_name == "phytools::phylosig(method='lambda')"
    assert cases[2].function_name == "phytools::phylosig(method='K')"
    assert cases[3].function_name == "phytools::phylosig(method='K')"
    assert cases[4].function_name == "phytools::fitMk(model='ER')"
    assert cases[7].function_name == "phytools::fitMk(model='ER')"
    assert cases[8].function_name == "phytools::fitMk(model='SYM')"
    assert cases[9].function_name == "phytools::fitMk(model='SYM')"
    assert cases[10].function_name == "phytools::fitMk(model='ARD')"
    assert cases[13].function_name == "phytools::fitMk(model='ARD')"
    assert cases[14].function_name == "phytools::make.simmap(model='ER')"
    assert cases[16].function_name == "phytools::make.simmap(model='ER')"
    assert cases[17].function_name == "phytools::make.simmap(model='SYM')"
    assert cases[18].function_name == "phytools::make.simmap(model='SYM')"
    assert cases[19].function_name == "phytools::make.simmap(model='ARD')"
    assert cases[22].function_name == "phytools::make.simmap(model='ARD')"
    assert cases[23].function_name == "phytools::countSimmap"
    assert cases[26].function_name == "phytools::countSimmap"
    assert cases[27].function_name == "phytools::describe.simmap"
    assert cases[30].function_name == "phytools::describe.simmap"
    assert cases[31].function_name == "phytools::densityMap"
    assert cases[32].function_name == "phytools::densityMap"
    assert cases[33].function_name == "phytools::sim.history"
    assert cases[36].function_name == "phytools::sim.history"
    assert cases[37].function_name == "phytools::sim.history"
    assert cases[38].function_name == "phytools::rerootingMethod"
    assert cases[42].function_name == "phytools::rerootingMethod"
    assert cases[43].function_name == "phytools::fastAnc"
    assert cases[46].function_name == "phytools::fastAnc"
    assert cases[47].function_name == "phytools::anc.ML"
    assert cases[50].function_name == "phytools::anc.ML"
    assert cases[51].function_name == "phytools::fastBM"
    assert cases[53].function_name == "phytools::fastBM"
    assert cases[54].function_name == "phytools::sim.corrs"
    assert cases[56].function_name == "phytools::sim.corrs"
    assert (
        cases[0].fixture_id
        == "phytools_continuous_strong_signal_non_ultrametric_twenty_four_taxa"
    )
    assert cases[2].permutation_count == 199
    assert cases[3].permutation_seed == 17
    assert cases[7].row_field_tolerances == {"rate": 1e-5}
    assert cases[8].row_field_tolerances == {"rate": 1e-4}
    assert cases[10].row_field_tolerances == {"rate": 1e-3}
    assert cases[14].stochastic_map_replicate_count == 128
    assert cases[16].stochastic_map_seed == 17
    assert cases[17].row_field_tolerances == {
        "mean_value": 2.5,
        "lower_95_interval": 5.0,
        "upper_95_interval": 5.0,
        "presence_fraction": 0.25,
    }
    assert cases[20].compare_rows is False
    assert cases[22].compare_rows is False
    assert cases[23].row_field_tolerances == {
        "mean_value": 1.5,
        "lower_95_interval": 3.0,
        "upper_95_interval": 3.0,
        "presence_fraction": 0.2,
    }
    assert cases[25].row_field_tolerances == {
        "mean_value": 2.5,
        "lower_95_interval": 5.0,
        "upper_95_interval": 5.0,
        "presence_fraction": 0.25,
    }
    assert cases[27].row_field_tolerances == {
        "mean_value": 1.5,
        "lower_95_interval": 3.0,
        "upper_95_interval": 3.0,
        "presence_fraction": 0.2,
    }
    assert cases[29].row_field_tolerances == {
        "mean_value": 2.5,
        "lower_95_interval": 5.0,
        "upper_95_interval": 5.0,
        "presence_fraction": 0.25,
    }
    assert cases[31].row_field_tolerances == {
        "mean_posterior_probability": 2e-3,
        "minimum_posterior_probability": 2e-3,
        "maximum_posterior_probability": 2e-3,
        "uncertainty": 2e-3,
    }
    assert cases[33].row_field_tolerances == {
        "mean_value": 0.25,
        "lower_95_interval": 0.35,
        "upper_95_interval": 0.35,
        "presence_fraction": 0.2,
    }
    assert cases[36].row_field_tolerances == {
        "mean_value": 1.5,
        "lower_95_interval": 3.0,
        "upper_95_interval": 3.0,
        "presence_fraction": 0.2,
    }
    assert cases[37].simulation_root_state_probabilities == {"0": 0.2, "1": 0.8}
    assert cases[37].field_tolerances == {
        "mean_total_transition_count": 1.0,
        "lower_95_total_transition_count": 2.0,
        "upper_95_total_transition_count": 2.0,
    }
    assert cases[38].row_field_tolerances == {"probability": 1e-5}
    assert cases[41].row_field_tolerances == {"probability": 5e-5}
    assert cases[42].row_field_tolerances == {"probability": 5e-5}
    assert cases[50].row_field_tolerances == {
        "estimate": 1e-8,
        "standard_error": 5e-8,
        "lower_95_interval": 5e-8,
        "upper_95_interval": 5e-8,
    }
    assert cases[46].row_field_tolerances == {
        "estimate": 1e-8,
        "standard_error": 1e-8,
    }
    assert cases[51].continuous_sigma_squared == 0.25
    assert cases[52].continuous_root_state == 2.5
    assert cases[53].continuous_replicate_count == 512
    assert cases[54].continuous_trait_names == ("trait_alpha", "trait_beta")
    assert cases[55].continuous_root_states == (2.0, -1.0)
    assert cases[56].continuous_replicate_count == 512
    assert cases[57].function_name == "phytools::pgls.SEy"
    assert cases[59].function_name == "phytools::pgls.SEy"
    assert cases[57].comparative_formula == "response ~ predictor_one"
    assert cases[58].comparative_formula == "response ~ habitat + diet"
    assert cases[59].comparative_formula == "response ~ habitat * diet"
    assert cases[59].comparative_lambda_value == 1.0
    assert cases[60].function_name == "phytools::phyl.resid(method='BM')"
    assert cases[62].function_name == "phytools::phyl.resid(method='lambda')"
    assert cases[60].comparative_predictors == ("body_mass",)
    assert cases[60].comparative_lambda_value == 1.0
    assert cases[61].field_tolerances == {
        "lambda_value": 5e-4,
        "log_likelihood": 5e-4,
    }
    assert cases[63].function_name == "phytools::phylANOVA"
    assert cases[64].comparative_predictors == ("habitat",)
    assert cases[63].permutation_count == 199
    assert cases[64].permutation_seed == 17


def test_list_phytools_parity_cases_resolves_existing_fixture_inputs() -> None:
    cases = list_phytools_parity_cases()

    missing_inputs = [
        (case.case_id, str(path))
        for case in cases
        for path in case.input_fixtures
        if not path.exists()
    ]

    assert missing_inputs == []
