from __future__ import annotations

from bijux_phylogenetics.parity import list_geiger_parity_cases


def test_list_geiger_parity_cases_returns_governed_registry() -> None:
    cases = list_geiger_parity_cases()

    assert [case.case_id for case in cases] == [
        "fitcontinuous-bm-example-tree",
        "fitcontinuous-bm-brownian-sigma-recovery",
        "fitcontinuous-bm-missing-values-review",
        "fitcontinuous-white-strong-signal-review",
        "fitcontinuous-white-weak-signal-review",
        "fitcontinuous-white-missing-values-review",
        "fitcontinuous-lambda-strong-signal-review",
        "fitcontinuous-lambda-weak-signal-review",
        "fitcontinuous-lambda-missing-values-review",
        "fitcontinuous-lambda-bounded-control-review",
        "fitcontinuous-kappa-strong-signal-review",
        "fitcontinuous-kappa-weak-signal-review",
        "fitcontinuous-kappa-missing-values-review",
        "fitcontinuous-delta-strong-signal-review",
        "fitcontinuous-delta-weak-signal-review",
        "fitcontinuous-delta-missing-values-review",
        "fitcontinuous-ou-ou-parameter-recovery",
        "fitcontinuous-ou-missing-values-review",
        "fitcontinuous-ou-lower-boundary-review",
        "fitcontinuous-ou-bounded-control-review",
        "fitcontinuous-eb-early-burst-rate-recovery",
        "fitcontinuous-eb-lower-boundary-review",
        "fitcontinuous-model-comparison-brownian-review",
        "fitcontinuous-model-comparison-ou-review",
        "fitcontinuous-model-comparison-early-burst-review",
        "fitcontinuous-model-comparison-white-review",
        "fitdiscrete-er-binary-twenty-four-taxa",
        "fitdiscrete-er-multistate-missing-twenty-four-taxa",
        "fitdiscrete-er-multistate-tip-intersection-review",
        "fitdiscrete-lambda-strong-signal-review",
        "fitdiscrete-lambda-weak-signal-review",
        "fitdiscrete-lambda-missing-values-review",
        "fitdiscrete-kappa-strong-signal-review",
        "fitdiscrete-kappa-weak-signal-review",
        "fitdiscrete-kappa-missing-values-review",
        "fitdiscrete-delta-late-change-boundary-review",
        "fitdiscrete-delta-earliest-change-review",
        "fitdiscrete-delta-missing-values-review",
        "fitdiscrete-early-burst-early-change-review",
        "fitdiscrete-early-burst-weak-signal-review",
        "fitdiscrete-early-burst-late-change-review",
        "fitdiscrete-early-burst-missing-values-review",
        "fitdiscrete-sym-three-state-twenty-four-taxa",
        "fitdiscrete-sym-four-state-twenty-four-taxa",
        "fitdiscrete-sym-multistate-missing-twenty-four-taxa",
        "fitdiscrete-ard-binary-twenty-four-taxa",
        "fitdiscrete-ard-four-state-weak-identification-review",
        "fitdiscrete-ard-multistate-missing-twenty-four-taxa",
    ]
    assert cases[0].function_name == "geiger::fitContinuous(model='BM')"
    assert cases[3].function_name == "geiger::fitContinuous(model='white')"
    assert cases[6].function_name == "geiger::fitContinuous(model='lambda')"
    assert cases[10].function_name == "geiger::fitContinuous(model='kappa')"
    assert cases[13].function_name == "geiger::fitContinuous(model='delta')"
    assert cases[16].function_name == "geiger::fitContinuous(model='OU')"
    assert cases[20].function_name == "geiger::fitContinuous(model='EB')"
    assert cases[22].function_name == "geiger::fitContinuous(model comparison)"
    assert cases[26].function_name == "geiger::fitDiscrete(model='ER')"
    assert (
        cases[29].function_name == "geiger::fitDiscrete(model='ER', transform='lambda')"
    )
    assert (
        cases[32].function_name == "geiger::fitDiscrete(model='ER', transform='kappa')"
    )
    assert (
        cases[34].function_name == "geiger::fitDiscrete(model='SYM', transform='kappa')"
    )
    assert (
        cases[35].function_name == "geiger::fitDiscrete(model='ER', transform='delta')"
    )
    assert (
        cases[36].function_name == "geiger::fitDiscrete(model='SYM', transform='delta')"
    )
    assert cases[38].function_name == "geiger::fitDiscrete(model='ER', transform='EB')"
    assert cases[42].function_name == "geiger::fitDiscrete(model='SYM')"
    assert cases[45].function_name == "geiger::fitDiscrete(model='ARD')"
    assert cases[1].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[2].fixture_id == "geiger_continuous_missing_values_twenty_four_taxa"
    assert cases[3].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[4].fixture_id == "geiger_continuous_white_noise_twenty_four_taxa"
    assert cases[5].fixture_id == "geiger_continuous_missing_values_twenty_four_taxa"
    assert cases[6].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[7].fixture_id == "geiger_continuous_white_noise_twenty_four_taxa"
    assert cases[8].fixture_id == "geiger_continuous_missing_values_twenty_four_taxa"
    assert cases[9].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[10].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[11].fixture_id == "geiger_continuous_white_noise_twenty_four_taxa"
    assert cases[12].fixture_id == "geiger_continuous_missing_values_twenty_four_taxa"
    assert cases[16].fixture_id == "geiger_continuous_ou_known_truth_twenty_four_taxa"
    assert cases[17].fixture_id == "geiger_continuous_missing_values_twenty_four_taxa"
    assert (
        cases[18].fixture_id
        == "geiger_continuous_nonultrametric_control_twenty_four_taxa"
    )
    assert cases[19].fixture_id == "geiger_continuous_ou_known_truth_twenty_four_taxa"
    assert (
        cases[20].fixture_id
        == "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
    )
    assert cases[21].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[22].fixture_id == "geiger_continuous_brownian_signal_twenty_four_taxa"
    assert cases[23].fixture_id == "geiger_continuous_ou_known_truth_twenty_four_taxa"
    assert (
        cases[24].fixture_id
        == "geiger_continuous_early_burst_known_truth_twenty_four_taxa"
    )
    assert cases[25].fixture_id == "geiger_continuous_white_noise_twenty_four_taxa"
    assert cases[26].fixture_id == "geiger_discrete_er_binary_twenty_four_taxa"
    assert (
        cases[27].fixture_id == "geiger_discrete_missing_three_state_twenty_four_taxa"
    )
    assert (
        cases[28].fixture_id == "geiger_discrete_mismatch_four_state_twenty_four_taxa"
    )
    assert cases[29].fixture_id == "geiger_discrete_er_binary_twenty_four_taxa"
    assert (
        cases[30].fixture_id == "geiger_discrete_transform_weak_signal_twenty_four_taxa"
    )
    assert (
        cases[31].fixture_id == "geiger_discrete_lambda_missing_binary_twenty_four_taxa"
    )
    assert (
        cases[32].fixture_id
        == "geiger_discrete_kappa_branch_sensitive_twenty_four_taxa"
    )
    assert cases[33].fixture_id == "geiger_discrete_kappa_weak_signal_twenty_four_taxa"
    assert (
        cases[34].fixture_id
        == "geiger_discrete_kappa_missing_three_state_twenty_four_taxa"
    )
    assert (
        cases[35].fixture_id
        == "geiger_discrete_delta_late_change_binary_twenty_four_taxa"
    )
    assert cases[36].fixture_id == "geiger_discrete_sym_three_state_twenty_four_taxa"
    assert (
        cases[37].fixture_id == "geiger_discrete_missing_three_state_twenty_four_taxa"
    )
    assert (
        cases[38].fixture_id
        == "geiger_discrete_early_burst_early_change_twenty_four_taxa"
    )
    assert (
        cases[39].fixture_id
        == "geiger_discrete_early_burst_weak_signal_twenty_four_taxa"
    )
    assert (
        cases[40].fixture_id
        == "geiger_discrete_early_burst_late_change_twenty_four_taxa"
    )
    assert (
        cases[41].fixture_id
        == "geiger_discrete_early_burst_missing_binary_twenty_four_taxa"
    )
    assert cases[42].fixture_id == "geiger_discrete_sym_three_state_twenty_four_taxa"
    assert cases[43].fixture_id == "geiger_discrete_sym_four_state_twenty_four_taxa"
    assert (
        cases[44].fixture_id == "geiger_discrete_missing_three_state_twenty_four_taxa"
    )
    assert cases[45].fixture_id == "geiger_discrete_ard_binary_twenty_four_taxa"
    assert cases[46].fixture_id == "geiger_discrete_ard_four_state_twenty_four_taxa"
    assert (
        cases[47].fixture_id == "geiger_discrete_missing_three_state_twenty_four_taxa"
    )
    assert cases[2].comparison_fields[:7] == (
        "taxon_count",
        "trait_name",
        "model_name",
        "excluded_taxon_count",
        "excluded_taxa",
        "missing_value_taxa",
        "non_numeric_taxa",
    )
    assert cases[1].optimizer_settings is not None
    assert cases[3].optimizer_settings["bijux_optimizer_name"] == (
        "closed-form-profile-solution"
    )
    assert cases[6].optimizer_settings["bijux_optimizer_name"] == (
        "governed-two-stage-grid-search"
    )
    assert cases[9].optimizer_settings["bijux_optimizer_name"] == (
        "governed-two-stage-grid-search"
    )
    assert cases[10].optimizer_settings["bijux_optimizer_name"] == (
        "governed-two-stage-grid-search"
    )
    assert cases[13].optimizer_settings["bijux_optimizer_name"] == (
        "governed-two-stage-grid-search"
    )
    assert cases[16].optimizer_settings["bijux_optimizer_name"] == (
        "governed-two-stage-grid-search"
    )
    assert "aicc" in cases[3].comparison_fields
    assert "aicc" in cases[4].comparison_fields
    assert "excluded_taxa" in cases[5].comparison_fields
    assert "aicc" in cases[6].comparison_fields
    assert "hit_lower_parameter_boundary" in cases[7].comparison_fields
    assert "excluded_taxa" in cases[8].comparison_fields
    assert "optimizer_settings" in cases[9].comparison_fields
    assert "aicc" in cases[10].comparison_fields
    assert "hit_upper_parameter_boundary" in cases[11].comparison_fields
    assert "excluded_taxa" in cases[12].comparison_fields
    assert "hit_lower_parameter_boundary" in cases[18].comparison_fields
    assert "optimizer_settings" in cases[19].comparison_fields
    assert "aicc" in cases[20].comparison_fields
    assert "hit_lower_parameter_boundary" in cases[21].comparison_fields
    assert cases[22].candidate_model_names == (
        "BM",
        "white",
        "lambda",
        "kappa",
        "delta",
        "OU",
        "EB",
    )
    assert "selected_model" in cases[22].comparison_fields
    assert "runner_up_aicc_delta" in cases[22].comparison_fields
    assert cases[26].comparison_fields[:6] == (
        "taxon_count",
        "trait_name",
        "model_name",
        "observed_state_count",
        "state_order",
        "excluded_taxon_count",
    )
    assert cases[26].optimizer_settings is not None
    assert cases[26].optimizer_settings["bijux_optimizer_name"] == (
        "golden-section-search"
    )
    assert cases[29].optimizer_settings is not None
    assert cases[29].optimizer_settings["bijux_optimizer_name"] == (
        "bounded-coarse-and-golden-search"
    )
    assert cases[32].optimizer_settings is not None
    assert cases[32].optimizer_settings["bijux_optimizer_name"] == (
        "bounded-coarse-and-golden-search"
    )
    assert cases[35].optimizer_settings is not None
    assert cases[35].optimizer_settings["bijux_optimizer_name"] == (
        "bounded-coarse-and-golden-search"
    )
    assert cases[38].optimizer_settings is not None
    assert cases[38].optimizer_settings["bijux_optimizer_name"] == (
        "bounded-coarse-and-golden-search"
    )
    assert cases[42].optimizer_settings is not None
    assert cases[42].optimizer_settings["bijux_optimizer_name"] == "nelder-mead"
    assert cases[45].optimizer_settings is not None
    assert cases[45].optimizer_settings["bijux_optimizer_name"] == "nelder-mead"
    assert all(path.is_file() for case in cases for path in case.input_fixtures)
