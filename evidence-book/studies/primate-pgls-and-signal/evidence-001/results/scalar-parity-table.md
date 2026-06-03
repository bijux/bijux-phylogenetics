# Scalar Parity Table

| Row | Family | Metric | Kind | Verdict | R | Bijux |
| --- | --- | --- | --- | --- | --- | --- |
| reload-object-count | workflow-contracts | object_name_count | exact_answer | matched | 2 | 2 |
| reload-primate-row-count | workflow-contracts | primate_row_count | exact_answer | matched | 75 | 75 |
| reload-tree-tip-count | workflow-contracts | tree_tip_count | exact_answer | matched | 75 | 75 |
| ou-alpha-1-branch-count | transformed-tree-workflows | ou_alpha_1_branch_count | exact_answer | matched | 148 | 148 |
| ou-alpha-1-total-branch-length | transformed-tree-workflows | ou_alpha_1_total_branch_length | exact_answer | mismatch_unexplained | 11.8428761446323 | 11.8429 |
| ou-alpha-10-total-branch-length | transformed-tree-workflows | ou_alpha_10_total_branch_length | exact_answer | mismatch_unexplained | 3.20350911964606 | 3.2035 |
| early-burst-2-total-branch-length | transformed-tree-workflows | early_burst_2_total_branch_length | exact_answer | mismatch_unexplained | 87.5076762141263 | 4.1181 |
| late-burst-minus-2-total-branch-length | transformed-tree-workflows | late_burst_minus_2_total_branch_length | exact_answer | mismatch_unexplained | 4.11805690444875 | 87.5077 |
| brownian-root-state | continuous-model-fitting | brownian_root_state | tolerance | matched_with_tolerance | 262.35002249971 | 262.35 |
| brownian-rate | continuous-model-fitting | brownian_rate | tolerance | matched_with_tolerance | 88740.5017326116 | 88740.4378 |
| brownian-log-likelihood | continuous-model-fitting | brownian_log_likelihood | tolerance | matched_with_tolerance | -476.42343794489 | -476.4234 |
| ou-alpha | continuous-model-fitting | ou_alpha | tolerance | matched_with_tolerance | 5.3343754665625 | 5.3356 |
| ou-log-likelihood | continuous-model-fitting | ou_log_likelihood | tolerance | matched_with_tolerance | -461.253095886538 | -461.2531 |
| early-burst-rate-change | continuous-model-fitting | early_burst_rate_change | tolerance | mismatch_unexplained | 10.6718757865625 | 0.0 |
| early-burst-log-likelihood | continuous-model-fitting | early_burst_log_likelihood | tolerance | mismatch_unexplained | -461.253095658179 | -476.4234 |
| brownian-ou-lrt-statistic | likelihood-ratio-tests | brownian_vs_ornstein_uhlenbeck_statistic | tolerance | matched_with_tolerance | 30.340684116704 | 30.3406 |
| brownian-eb-lrt-statistic | likelihood-ratio-tests | brownian_vs_early_burst_statistic | tolerance | mismatch_unexplained | 30.340684573422 | 0.0 |
| ou-eb-lrt-statistic | likelihood-ratio-tests | ornstein_uhlenbeck_vs_early_burst_statistic | tolerance | matched_with_tolerance | 4.56717998531531e-07 | 0.0 |
| ancestral-brownian-node-count | ancestral-reconstruction | brownian_node_count | exact_answer | matched | 74 | 74 |
| ancestral-brownian-first-five | ancestral-reconstruction | brownian_first_five_estimates | exact_answer | mismatch_unexplained | [262.350022499711, 295.392237823336, 369.773311672215, 411.2760353975, 334.804088053515] | [262.35, 295.3922, 369.7733, 411.276, 334.8041] |
| ancestral-brownian-recent-five | ancestral-reconstruction | brownian_recent_five_estimates | exact_answer | mismatch_unexplained | [192.371916302131, 189.874069323138, 171.574858834618, 174, 214.8] | [192.3719, 189.8741, 171.5749, 174.0, 214.8] |
| ancestral-eb-node-count | ancestral-reconstruction | early_burst_node_count | exact_answer | matched | 74 | 74 |
| ancestral-eb-first-five | ancestral-reconstruction | early_burst_first_five_estimates | exact_answer | mismatch_unexplained | [282.721003635748, 336.828890951664, 355.392585511741, 386.473721532351, 341.736336661997] | [262.35, 268.7938, 324.9443, 387.4761, 353.8869] |
| ancestral-eb-recent-five | ancestral-reconstruction | early_burst_recent_five_estimates | exact_answer | mismatch_unexplained | [190.301546647852, 187.945676006467, 171.753821601142, 174, 214.8] | [205.799, 202.3003, 182.0545, 177.8185, 213.5094] |
| baseline-intercept | baseline-regression | intercept | tolerance | matched_with_tolerance | 263.916162624481 | 263.916162626577 |
| baseline-slope | baseline-regression | social_group_size | tolerance | matched_with_tolerance | 3.59115553051772 | 3.59115553047129 |
| baseline-log-likelihood | baseline-regression | log_likelihood | tolerance | matched_with_tolerance | -456.310146562924 | -456.310146563131 |
| baseline-r-squared | baseline-regression | r_squared | tolerance | matched | 0.257169890828586 | 0.257169890828586 |
| fixed-reference-lambda-value | phylogenetic-regression | fixed_reference_lambda_value | tolerance | matched | 0.768656451898616 | 0.768656451898616 |
| fixed-reference-pgls-intercept | phylogenetic-regression | fixed_reference_intercept | tolerance | matched_with_tolerance | 249.9584349618 | 249.958434911554 |
| fixed-reference-pgls-slope | phylogenetic-regression | fixed_reference_social_group_size | tolerance | matched_with_tolerance | 1.66214823716133 | 1.66214825503763 |
| fixed-reference-pgls-intercept-standard-error | phylogenetic-regression | fixed_reference_intercept_standard_error | tolerance | matched_with_tolerance | 49.7282720332178 | 49.7282709322746 |
| fixed-reference-pgls-slope-standard-error | phylogenetic-regression | fixed_reference_social_group_size_standard_error | tolerance | matched_with_tolerance | 0.679751504116516 | 0.679751503543842 |
| fixed-reference-pgls-intercept-p-value | phylogenetic-regression | fixed_reference_intercept_p_value | tolerance | matched_with_tolerance | 3.44325081102684e-06 | 3.44324934209439e-06 |
| fixed-reference-pgls-slope-p-value | phylogenetic-regression | fixed_reference_social_group_size_p_value | tolerance | matched_with_tolerance | 0.0168886063098727 | 0.0168886050849197 |
| fixed-reference-pgls-log-likelihood | phylogenetic-regression | fixed_reference_log_likelihood | tolerance | matched_with_tolerance | -444.069148162235 | -444.069148162264 |
| fixed-reference-pgls-aic | phylogenetic-regression | fixed_reference_aic | tolerance | matched_with_tolerance | 894.13829632447 | 894.138296324528 |
| estimated-lambda-value | phylogenetic-regression | lambda_value | tolerance | matched_with_tolerance | 0.768656451898616 | 0.76865646856019 |
| estimated-pgls-intercept | phylogenetic-regression | intercept | tolerance | matched_with_tolerance | 249.9584349618 | 249.95843497104 |
| estimated-pgls-slope | phylogenetic-regression | social_group_size | tolerance | matched_with_tolerance | 1.66214823716133 | 1.66214823566609 |
| estimated-pgls-intercept-standard-error | phylogenetic-regression | intercept_standard_error | tolerance | matched_with_tolerance | 49.7282720332178 | 49.7282721257287 |
| estimated-pgls-slope-standard-error | phylogenetic-regression | social_group_size_standard_error | tolerance | matched_with_tolerance | 0.679751504116516 | 0.67975150416116 |
| estimated-pgls-log-likelihood | phylogenetic-regression | log_likelihood | tolerance | matched_with_tolerance | -444.069148162235 | -444.069148162264 |
| estimated-pgls-aic | phylogenetic-regression | aic | tolerance | matched_with_tolerance | 896.13829632447 | 896.138296324527 |
| estimated-pgls-intercept-p-value | phylogenetic-regression | intercept_p_value | tolerance | matched | 3.44325081102684e-06 | 3.44325093304398e-06 |
| estimated-pgls-slope-p-value | phylogenetic-regression | social_group_size_p_value | tolerance | matched_with_tolerance | 0.0168886063098727 | 0.0168886064118261 |
| estimated-pgls-slope-significance | phylogenetic-regression | social_group_size_significant_under_0.05 | scientific_equivalence | matched_with_tolerance | True | True |
| signal-estimated-lambda | phylogenetic-signal | estimated_lambda | tolerance | matched_with_tolerance | 0.80273442404245 | 0.805 |
| signal-likelihood-ratio | phylogenetic-signal | likelihood_ratio | tolerance | matched_with_tolerance | 41.0304186930398 | 41.0296165956062 |
| signal-reject-lambda-zero | phylogenetic-signal | p_value_below_0.05 | scientific_equivalence | matched_with_tolerance | True | True |
| baseline-diagnostic-qq-correlation | diagnostics | qq_correlation | tolerance | matched_with_tolerance | 0.974950137688634 | 0.974935270488255 |
| baseline-diagnostic-fitted-correlation | diagnostics | abs_residual_fitted_correlation | tolerance | matched_with_tolerance | 0.0234397354786275 | 0.0234397354757235 |
| baseline-diagnostic-outlier-pressure | diagnostics | outlier_count_abs_z_ge_2_close | scientific_equivalence | matched_with_tolerance | 5 | 5 |
| estimated-diagnostic-qq-correlation | diagnostics | qq_correlation | tolerance | matched_with_tolerance | 0.980223544138575 | 0.977685033649863 |
| estimated-diagnostic-fitted-correlation | diagnostics | abs_residual_fitted_correlation | tolerance | mismatch_unexplained | 0.205562548325984 | 0.109972191549302 |
| estimated-diagnostic-outlier-pressure | diagnostics | outlier_count_abs_z_ge_2_close | scientific_equivalence | matched_with_tolerance | 3 | 4 |

## Verdict Counts

- `matched`: `9`
- `matched_with_tolerance`: `35`
- `mismatch_unexplained`: `12`
