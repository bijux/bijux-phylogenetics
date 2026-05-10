args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: primate_pgls_and_signal_reference_r.R <repo_root> <out_path>")
}

repo_root <- normalizePath(args[[1]], mustWork = TRUE)
out_path <- normalizePath(args[[2]], mustWork = FALSE)
source_csv <- file.path(
  repo_root,
  "evidence-book",
  "studies",
  "primate-longevity-signal",
  "evidence-001",
  "reference_primate.csv"
)
source_tree <- file.path(
  repo_root,
  "evidence-book",
  "studies",
  "primate-longevity-signal",
  "evidence-001",
  "reference_trimmed_primatetree.nwk"
)

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(jsonlite))

primate <- read.csv(source_csv, stringsAsFactors = FALSE)
primate$species <- as.character(primate$species)
primatetree <- read.tree(source_tree)

diagnostic_summary <- function(fitted_values, residuals) {
  residual_mean <- mean(residuals)
  residual_variance <- if (length(residuals) > 1) var(residuals) else 0
  residual_sd <- sqrt(max(residual_variance, 0))
  standardized <- if (residual_sd == 0) {
    rep(0, length(residuals))
  } else {
    (residuals - residual_mean) / residual_sd
  }
  abs_residuals <- abs(residuals)
  abs_residual_fitted_correlation <- suppressWarnings(cor(fitted_values, abs_residuals))
  if (is.na(abs_residual_fitted_correlation)) {
    abs_residual_fitted_correlation <- 0
  }
  qq_correlation <- suppressWarnings(
    cor(
      qnorm(ppoints(length(residuals))),
      sort(residuals)
    )
  )
  if (is.na(qq_correlation)) {
    qq_correlation <- 0
  }
  list(
    residual_mean = unname(signif(residual_mean, 15)),
    residual_variance = unname(signif(residual_variance, 15)),
    residual_sd = unname(signif(residual_sd, 15)),
    max_abs_z_residual = unname(signif(max(abs(standardized)), 15)),
    abs_residual_fitted_correlation = unname(signif(abs_residual_fitted_correlation, 15)),
    qq_correlation = unname(signif(qq_correlation, 15)),
    outlier_count_abs_z_ge_2 = as.integer(sum(abs(standardized) >= 2))
  )
}

r_squared <- function(observed, fitted_values) {
  total <- sum((observed - mean(observed)) ^ 2)
  residual <- sum((observed - fitted_values) ^ 2)
  if (total == 0) {
    return(1)
  }
  unname(signif(1 - (residual / total), 15))
}

gls1 <- gls(longevity ~ social_group_size, data = primate, method = "ML")
pgls1 <- gls(
  longevity ~ social_group_size,
  data = primate,
  correlation = corPagel(0, primatetree, fixed = TRUE, form = ~species),
  method = "ML"
)
pgls2 <- gls(
  longevity ~ social_group_size,
  data = primate,
  correlation = corPagel(0, primatetree, fixed = FALSE, form = ~species),
  method = "ML"
)
pgls3 <- gls(
  longevity ~ 1,
  data = primate,
  correlation = corPagel(0, primatetree, fixed = FALSE, form = ~species),
  method = "ML"
)
pgls4 <- gls(
  longevity ~ 1,
  data = primate,
  correlation = corPagel(0, primatetree, fixed = TRUE, form = ~species),
  method = "ML"
)

baseline_t_table <- summary(gls1)$tTable
estimated_t_table <- summary(pgls2)$tTable
signal_lambda <- unname(coef(pgls3$modelStruct$corStruct, unconstrained = FALSE)[["lambda"]])
signal_log_likelihood <- as.numeric(logLik(pgls3))
signal_null_log_likelihood <- as.numeric(logLik(pgls4))
signal_likelihood_ratio <- -2 * (signal_null_log_likelihood - signal_log_likelihood)
signal_p_value <- pchisq(signal_likelihood_ratio, df = 1, lower.tail = FALSE)

payload <- list(
  schema_version = 1,
  study_id = "primate-pgls-and-signal",
  source_contract = list(
    object_names = c("primate", "primatetree"),
    object_name_count = 2L,
    row_count = nrow(primate),
    tip_count = length(primatetree$tip.label),
    species_tip_match = setequal(primate$species, primatetree$tip.label)
  ),
  baseline_gls = list(
    coefficients = list(
      intercept = unname(signif(coef(gls1)[["(Intercept)"]], 15)),
      social_group_size = unname(signif(coef(gls1)[["social_group_size"]], 15))
    ),
    p_values = list(
      intercept = unname(signif(baseline_t_table["(Intercept)", "p-value"], 15)),
      social_group_size = unname(signif(baseline_t_table["social_group_size", "p-value"], 15))
    ),
    log_likelihood = unname(signif(as.numeric(logLik(gls1)), 15)),
    r_squared = r_squared(primate$longevity, fitted(gls1)),
    diagnostics = diagnostic_summary(fitted(gls1), residuals(gls1, type = "n"))
  ),
  fixed_lambda_gls_matches_baseline = isTRUE(all.equal(
    as.numeric(coef(gls1)),
    as.numeric(coef(pgls1)),
    tolerance = 1e-12
  )) && isTRUE(all.equal(
    as.numeric(logLik(gls1)),
    as.numeric(logLik(pgls1)),
    tolerance = 1e-12
  )),
  estimated_lambda_pgls = list(
    lambda_value = unname(signif(unname(coef(pgls2$modelStruct$corStruct, unconstrained = FALSE)[["lambda"]]), 15)),
    coefficients = list(
      intercept = unname(signif(coef(pgls2)[["(Intercept)"]], 15)),
      social_group_size = unname(signif(coef(pgls2)[["social_group_size"]], 15))
    ),
    p_values = list(
      intercept = unname(signif(estimated_t_table["(Intercept)", "p-value"], 15)),
      social_group_size = unname(signif(estimated_t_table["social_group_size", "p-value"], 15))
    ),
    log_likelihood = unname(signif(as.numeric(logLik(pgls2)), 15)),
    r_squared = r_squared(primate$longevity, fitted(pgls2)),
    diagnostics = diagnostic_summary(fitted(pgls2), residuals(pgls2, type = "n"))
  ),
  signal_test = list(
    estimated_lambda = unname(signif(signal_lambda, 15)),
    estimated_log_likelihood = unname(signif(signal_log_likelihood, 15)),
    null_log_likelihood = unname(signif(signal_null_log_likelihood, 15)),
    likelihood_ratio = unname(signif(signal_likelihood_ratio, 15)),
    p_value = unname(signif(signal_p_value, 15))
  ),
  coverage_boundaries = list(
    uncovered_fragments = c(
      "transformed-tree-workflows",
      "continuous-model-comparison",
      "ancestral-mode-comparison",
      "mode-linked-intercept-models"
    ),
    notes = c(
      "The lecture EB surfaces remain outside the current canonical runtime parity surface.",
      "Transformed-tree and ancestral comparison fragments stay explicitly indexed until governed runtime coverage exists."
    )
  )
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(toJSON(payload, auto_unbox = TRUE, pretty = TRUE), con = out_path)
