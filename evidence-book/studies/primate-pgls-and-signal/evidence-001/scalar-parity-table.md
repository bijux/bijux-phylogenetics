# Scalar Parity Table

| Row | Family | Metric | Kind | Verdict | R | Bijux |
| --- | --- | --- | --- | --- | --- | --- |
| reload-object-count | workflow-contracts | object_name_count | exact_answer | matched | 2 | 2 |
| reload-primate-row-count | workflow-contracts | primate_row_count | exact_answer | matched | 75 | 75 |
| reload-tree-tip-count | workflow-contracts | tree_tip_count | exact_answer | matched | 75 | 75 |
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

- `matched`: `3`
- `matched_with_tolerance`: `13`
- `mismatch_explained`: `4`
- `mismatch_unexplained`: `1`
