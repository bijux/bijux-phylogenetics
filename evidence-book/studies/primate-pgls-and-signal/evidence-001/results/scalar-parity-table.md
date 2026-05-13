# Scalar Parity Table

| Row | Family | Metric | Kind | Verdict | R | Bijux |
| --- | --- | --- | --- | --- | --- | --- |
| reload-object-count | workflow-contracts | object_name_count | exact_answer | matched | 2 | 2 |
| reload-primate-row-count | workflow-contracts | primate_row_count | exact_answer | matched | 75 | 75 |
| reload-tree-tip-count | workflow-contracts | tree_tip_count | exact_answer | matched | 75 | 75 |
| ou-alpha-1-branch-count | transformed-tree-workflows | ou_alpha_1_branch_count | exact_answer | matched | 148 | 148 |
| ou-alpha-1-total-branch-length | transformed-tree-workflows | ou_alpha_1_total_branch_length | exact_answer | matched | 11.8429 | 11.8429 |
| ou-alpha-10-total-branch-length | transformed-tree-workflows | ou_alpha_10_total_branch_length | exact_answer | matched | 3.2035 | 3.2035 |
| early-burst-2-total-branch-length | transformed-tree-workflows | early_burst_2_total_branch_length | exact_answer | matched | 87.5077 | 87.5077 |
| late-burst-minus-2-total-branch-length | transformed-tree-workflows | late_burst_minus_2_total_branch_length | exact_answer | matched | 4.1181 | 4.1181 |
| brownian-root-state | continuous-model-fitting | brownian_root_state | tolerance | matched | 262.35 | 262.35 |
| brownian-rate | continuous-model-fitting | brownian_rate | tolerance | matched_with_tolerance | 88740.5017 | 88740.4378 |
| brownian-log-likelihood | continuous-model-fitting | brownian_log_likelihood | tolerance | matched | -476.4234 | -476.4234 |
| ou-alpha | continuous-model-fitting | ou_alpha | tolerance | matched | 5.3344 | 5.3344 |
| ou-log-likelihood | continuous-model-fitting | ou_log_likelihood | tolerance | matched | -461.2531 | -461.2531 |
| early-burst-rate-change | continuous-model-fitting | early_burst_rate_change | tolerance | matched | 10.6719 | 10.6719 |
| early-burst-log-likelihood | continuous-model-fitting | early_burst_log_likelihood | tolerance | matched | -461.2531 | -461.2531 |
| brownian-ou-lrt-statistic | likelihood-ratio-tests | brownian_vs_ornstein_uhlenbeck_statistic | tolerance | matched | 30.3407 | 30.3407 |
| brownian-eb-lrt-statistic | likelihood-ratio-tests | brownian_vs_early_burst_statistic | tolerance | matched_with_tolerance | 30.3407 | 30.3406 |
| ou-eb-lrt-statistic | likelihood-ratio-tests | ornstein_uhlenbeck_vs_early_burst_statistic | tolerance | matched_with_tolerance | 4.5672e-07 | 0.0 |
| ancestral-brownian-node-count | ancestral-reconstruction | brownian_node_count | exact_answer | matched | 74 | 74 |
| ancestral-brownian-first-five | ancestral-reconstruction | brownian_first_five_estimates | exact_answer | matched | [262.35, 295.3922, 369.7733, 411.276, 334.8041] | [262.35, 295.3922, 369.7733, 411.276, 334.8041] |
| ancestral-brownian-recent-five | ancestral-reconstruction | brownian_recent_five_estimates | exact_answer | matched | [192.3719, 189.8741, 171.5749, 174, 214.8] | [192.3719, 189.8741, 171.5749, 174.0, 214.8] |
| ancestral-eb-node-count | ancestral-reconstruction | early_burst_node_count | exact_answer | matched | 74 | 74 |
| ancestral-eb-first-five | ancestral-reconstruction | early_burst_first_five_estimates | exact_answer | mismatch_unexplained | [282.721, 336.8289, 355.3926, 386.4737, 341.7363] | [260.1962, 285.3532, 371.7088, 415.4058, 334.0808] |
| ancestral-eb-recent-five | ancestral-reconstruction | early_burst_recent_five_estimates | exact_answer | mismatch_unexplained | [190.3015, 187.9457, 171.7538, 174, 214.8] | [192.7689, 190.2801, 171.5387, 174.0, 214.8] |
| baseline-intercept | baseline-regression | intercept | tolerance | mismatch_explained | 263.9162 | 263.916162626577 |
| baseline-slope | baseline-regression | social_group_size | tolerance | mismatch_explained | 3.5912 | 3.59115553047129 |
| baseline-log-likelihood | baseline-regression | log_likelihood | tolerance | mismatch_explained | -456.3101 | -456.310146563131 |
| baseline-r-squared | baseline-regression | r_squared | tolerance | mismatch_explained | 0.2572 | 0.257169890828586 |
| estimated-lambda-value | phylogenetic-regression | lambda_value | tolerance | matched_with_tolerance | 0.7687 | 0.805 |
| estimated-pgls-intercept | phylogenetic-regression | intercept | tolerance | matched_with_tolerance | 249.9584 | 250.105945368603 |
| estimated-pgls-slope | phylogenetic-regression | social_group_size | tolerance | matched_with_tolerance | 1.6621 | 1.61904779942812 |
| estimated-pgls-log-likelihood | phylogenetic-regression | log_likelihood | tolerance | matched_with_tolerance | -444.0691 | -444.155346705612 |
| estimated-pgls-slope-significance | phylogenetic-regression | social_group_size_significant_under_0.05 | scientific_equivalence | matched_with_tolerance | True | True |
| signal-estimated-lambda | phylogenetic-signal | estimated_lambda | tolerance | matched_with_tolerance | 0.8027 | 0.805 |
| signal-likelihood-ratio | phylogenetic-signal | likelihood_ratio | tolerance | matched_with_tolerance | 41.0304 | 41.0296165956062 |
| signal-reject-lambda-zero | phylogenetic-signal | p_value_below_0.05 | scientific_equivalence | matched_with_tolerance | True | True |
| baseline-diagnostic-qq-correlation | diagnostics | qq_correlation | tolerance | matched_with_tolerance | 0.975 | 0.974935270488255 |
| baseline-diagnostic-fitted-correlation | diagnostics | abs_residual_fitted_correlation | tolerance | matched_with_tolerance | 0.0234 | 0.0234397354757235 |
| baseline-diagnostic-outlier-pressure | diagnostics | outlier_count_abs_z_ge_2_close | scientific_equivalence | matched_with_tolerance | 5 | 5 |
| estimated-diagnostic-qq-correlation | diagnostics | qq_correlation | tolerance | matched_with_tolerance | 0.9802 | 0.977832436449786 |
| estimated-diagnostic-fitted-correlation | diagnostics | abs_residual_fitted_correlation | tolerance | mismatch_unexplained | 0.2056 | 0.112199761717797 |
| estimated-diagnostic-outlier-pressure | diagnostics | outlier_count_abs_z_ge_2_close | scientific_equivalence | matched_with_tolerance | 3 | 4 |

## Verdict Counts

- `matched`: `19`
- `matched_with_tolerance`: `16`
- `mismatch_explained`: `4`
- `mismatch_unexplained`: `3`
