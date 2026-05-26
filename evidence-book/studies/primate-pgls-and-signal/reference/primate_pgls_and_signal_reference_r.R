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
  "datasets",
  "reference_primate.csv"
)
source_tree <- file.path(
  repo_root,
  "evidence-book",
  "studies",
  "primate-longevity-signal",
  "datasets",
  "reference_trimmed_primatetree.nwk"
)

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(jsonlite))

rounded <- function(value) {
  unname(signif(as.numeric(value), 15))
}

leaf_descendants <- function(phy, node) {
  if (node <= length(phy$tip.label)) {
    return(phy$tip.label[[node]])
  }
  children <- phy$edge[phy$edge[, 1] == node, 2]
  sort(unique(unlist(lapply(children, function(child) leaf_descendants(phy, child)))))
}

node_signature <- function(phy, node) {
  paste(leaf_descendants(phy, node), collapse = "|")
}

node_depths <- function(phy) {
  total_nodes <- length(phy$tip.label) + phy$Nnode
  depths <- rep(0, total_nodes)
  repeat {
    updated <- FALSE
    for (edge_index in seq_len(nrow(phy$edge))) {
      parent <- phy$edge[edge_index, 1]
      child <- phy$edge[edge_index, 2]
      edge_length <- phy$edge.length[[edge_index]]
      candidate <- depths[[parent]] + edge_length
      if (!isTRUE(all.equal(depths[[child]], candidate))) {
        depths[[child]] <- candidate
        updated <- TRUE
      }
    }
    if (!updated) {
      break
    }
  }
  names(depths) <- seq_len(total_nodes)
  depths
}

tree_total_depth <- function(phy) {
  max(node_depths(phy)[seq_along(phy$tip.label)])
}

rescale_ou_tree <- function(phy, alpha, sigsq = 1) {
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  transformed <- phy
  depths <- node_depths(phy)
  total_depth <- tree_total_depth(phy)
  transformed$edge.length <- apply(
    phy$edge,
    1,
    function(edge) {
      parent_depth <- depths[[as.character(edge[[1]])]]
      child_depth <- depths[[as.character(edge[[2]])]]
      term <- function(depth) {
        (1 / (2 * alpha)) *
          exp(-2 * alpha * (total_depth - depth)) *
          (1 - exp(-2 * alpha * depth))
      }
      rounded(max(0, (term(child_depth) - term(parent_depth)) * sigsq))
    }
  )
  transformed
}

rescale_eb_tree <- function(phy, a, sigsq = 1) {
  transformed <- phy
  depths <- node_depths(phy)
  transformed$edge.length <- apply(
    phy$edge,
    1,
    function(edge) {
      parent_depth <- depths[[as.character(edge[[1]])]]
      child_depth <- depths[[as.character(edge[[2]])]]
      if (isTRUE(all.equal(a, 0))) {
        return(rounded((child_depth - parent_depth) * sigsq))
      }
      rounded(max(0, ((exp(a * child_depth) - exp(a * parent_depth)) / a) * sigsq))
    }
  )
  transformed
}

tree_branch_rows <- function(phy) {
  depths <- node_depths(phy)
  rows <- lapply(
    seq_len(nrow(phy$edge)),
    function(edge_index) {
      parent <- phy$edge[edge_index, 1]
      child <- phy$edge[edge_index, 2]
      descendant_taxa <- leaf_descendants(phy, child)
      list(
        node = node_signature(phy, child),
        descendant_taxa = descendant_taxa,
        branch_length = rounded(phy$edge.length[[edge_index]]),
        parent_depth = rounded(depths[[as.character(parent)]]),
        child_depth = rounded(depths[[as.character(child)]])
      )
    }
  )
  ordered <- rows[order(vapply(rows, function(row) row$node, character(1)))]
  list(
    branch_count = length(ordered),
    total_branch_length = rounded(sum(phy$edge.length)),
    branch_rows = ordered
  )
}

brownian_intercept_fit <- function(phy, named_values, parameter_count) {
  tip_values <- named_values[phy$tip.label]
  covariance <- vcv.phylo(phy)[phy$tip.label, phy$tip.label]
  inverse_covariance <- solve(covariance)
  ones <- rep(1, length(tip_values))
  theta <- as.numeric(
    (t(ones) %*% inverse_covariance %*% tip_values) /
      (t(ones) %*% inverse_covariance %*% ones)
  )
  fitted_values <- rep(theta, length(tip_values))
  residuals <- as.numeric(tip_values - theta)
  sigma_squared <- as.numeric(
    (t(residuals) %*% inverse_covariance %*% residuals) / length(tip_values)
  )
  log_likelihood <- as.numeric(
    -0.5 * (
      length(tip_values) * log(2 * pi * sigma_squared) +
        determinant(covariance, logarithm = TRUE)$modulus[[1]] +
        length(tip_values)
    )
  )
  aic <- as.numeric((2 * parameter_count) - (2 * log_likelihood))
  list(
    tip_values = as.numeric(tip_values),
    root_state = rounded(theta),
    rate = rounded(sigma_squared),
    log_likelihood = rounded(log_likelihood),
    aic = rounded(aic),
    parameter_count = as.integer(parameter_count)
  )
}

mode_grid <- function(lower, upper) {
  coarse <- seq(lower, upper, length.out = 81)
  coarse
}

best_mode_fit <- function(phy, named_values, mode, bounds) {
  coarse <- mode_grid(bounds[[1]], bounds[[2]])
  best_parameter <- coarse[[1]]
  best_tree <- if (mode == "ornstein-uhlenbeck") {
    rescale_ou_tree(phy, alpha = best_parameter)
  } else {
    rescale_eb_tree(phy, a = best_parameter)
  }
  best_fit <- brownian_intercept_fit(best_tree, named_values, parameter_count = 3)
  best_index <- 1
  for (index in seq_along(coarse)[-1]) {
    parameter <- coarse[[index]]
    transformed_tree <- if (mode == "ornstein-uhlenbeck") {
      rescale_ou_tree(phy, alpha = parameter)
    } else {
      rescale_eb_tree(phy, a = parameter)
    }
    fit <- brownian_intercept_fit(transformed_tree, named_values, parameter_count = 3)
    if (fit$log_likelihood > best_fit$log_likelihood) {
      best_parameter <- parameter
      best_tree <- transformed_tree
      best_fit <- fit
      best_index <- index
    }
  }
  left <- coarse[[max(1, best_index - 1)]]
  right <- coarse[[min(length(coarse), best_index + 1)]]
  for (parameter in seq(left, right, length.out = 81)) {
    if (isTRUE(all.equal(parameter, best_parameter))) {
      next
    }
    transformed_tree <- if (mode == "ornstein-uhlenbeck") {
      rescale_ou_tree(phy, alpha = parameter)
    } else {
      rescale_eb_tree(phy, a = parameter)
    }
    fit <- brownian_intercept_fit(transformed_tree, named_values, parameter_count = 3)
    if (fit$log_likelihood > best_fit$log_likelihood) {
      best_parameter <- parameter
      best_tree <- transformed_tree
      best_fit <- fit
    }
  }
  list(
    parameter_value = rounded(best_parameter),
    tree = best_tree,
    fit = best_fit
  )
}

likelihood_ratio <- function(left_fit, right_fit) {
  statistic <- max(0, -2 * (left_fit$log_likelihood - right_fit$log_likelihood))
  list(
    statistic = rounded(statistic),
    p_value = rounded(pchisq(statistic, df = 1, lower.tail = FALSE))
  )
}

ancestral_rows <- function(phy, named_values) {
  ace_result <- ace(named_values[phy$tip.label], phy, type = "continuous", method = "pic")
  internal_nodes <- seq(length(phy$tip.label) + 1, length(phy$tip.label) + phy$Nnode)
  rows <- lapply(
    seq_along(internal_nodes),
    function(index) {
      node <- internal_nodes[[index]]
      list(
        node_index = as.integer(index),
        node = node_signature(phy, node),
        estimate = rounded(ace_result$ace[[index]])
      )
    }
  )
  list(
    node_count = length(rows),
    first_five_estimates = unname(vapply(rows[1:5], function(row) row$estimate, numeric(1))),
    recent_five_estimates = unname(vapply(rows[70:74], function(row) row$estimate, numeric(1))),
    rows = rows
  )
}

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
    residual_mean = rounded(residual_mean),
    residual_variance = rounded(residual_variance),
    residual_sd = rounded(residual_sd),
    max_abs_z_residual = rounded(max(abs(standardized))),
    abs_residual_fitted_correlation = rounded(abs_residual_fitted_correlation),
    qq_correlation = rounded(qq_correlation),
    outlier_count_abs_z_ge_2 = as.integer(sum(abs(standardized) >= 2))
  )
}

r_squared <- function(observed, fitted_values) {
  total <- sum((observed - mean(observed)) ^ 2)
  residual <- sum((observed - fitted_values) ^ 2)
  if (total == 0) {
    return(1)
  }
  rounded(1 - (residual / total))
}

gls_payload <- function(fit, observed, lambda_value = NULL) {
  coefficient_table <- summary(fit)$tTable
  payload <- list(
    coefficients = list(
      intercept = rounded(coef(fit)[["(Intercept)"]]),
      social_group_size = rounded(coef(fit)[["social_group_size"]])
    ),
    standard_errors = list(
      intercept = rounded(coefficient_table["(Intercept)", "Std.Error"]),
      social_group_size = rounded(coefficient_table["social_group_size", "Std.Error"])
    ),
    p_values = list(
      intercept = rounded(coefficient_table["(Intercept)", "p-value"]),
      social_group_size = rounded(coefficient_table["social_group_size", "p-value"])
    ),
    log_likelihood = rounded(as.numeric(logLik(fit))),
    aic = rounded(AIC(fit)),
    r_squared = r_squared(observed, fitted(fit)),
    diagnostics = diagnostic_summary(fitted(fit), residuals(fit, type = "n"))
  )
  if (!is.null(lambda_value)) {
    payload$lambda_value <- rounded(lambda_value)
  }
  payload
}

primate <- read.csv(source_csv, stringsAsFactors = FALSE)
primate$species <- as.character(primate$species)
primatetree <- read.tree(source_tree)
primateLL <- as.numeric(primate$longevity[match(primatetree$tip.label, primate$species)])
names(primateLL) <- primatetree$tip.label

ou1_tree <- rescale_ou_tree(primatetree, alpha = 1)
ou10_tree <- rescale_ou_tree(primatetree, alpha = 10)
eb2_tree <- rescale_eb_tree(primatetree, a = 2)
lb_minus2_tree <- rescale_eb_tree(primatetree, a = -2)
bm_fit <- brownian_intercept_fit(primatetree, primateLL, parameter_count = 2)
ou_fit <- best_mode_fit(primatetree, primateLL, mode = "ornstein-uhlenbeck", bounds = c(0.000001, 10))
eb_fit <- best_mode_fit(primatetree, primateLL, mode = "early-burst", bounds = c(0.000001, 50))
bm_ancestral <- ancestral_rows(primatetree, primateLL)
eb_ancestral <- ancestral_rows(eb_fit$tree, primateLL)

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
estimated_lambda_value <- unname(
  coef(pgls2$modelStruct$corStruct, unconstrained = FALSE)[["lambda"]]
)
pgls2_fixed <- gls(
  longevity ~ social_group_size,
  data = primate,
  correlation = corPagel(
    estimated_lambda_value,
    primatetree,
    fixed = TRUE,
    form = ~species
  ),
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
  tree_rescaling = list(
    ou_alpha_1 = tree_branch_rows(ou1_tree),
    ou_alpha_10 = tree_branch_rows(ou10_tree),
    early_burst_2 = tree_branch_rows(eb2_tree),
    late_burst_minus_2 = tree_branch_rows(lb_minus2_tree)
  ),
  continuous_mode_fits = list(
    brownian = bm_fit,
    ornstein_uhlenbeck = list(
      alpha = ou_fit$parameter_value,
      root_state = ou_fit$fit$root_state,
      rate = ou_fit$fit$rate,
      log_likelihood = ou_fit$fit$log_likelihood,
      aic = ou_fit$fit$aic
    ),
    early_burst = list(
      rate_change = eb_fit$parameter_value,
      root_state = eb_fit$fit$root_state,
      rate = eb_fit$fit$rate,
      log_likelihood = eb_fit$fit$log_likelihood,
      aic = eb_fit$fit$aic
    )
  ),
  likelihood_ratio_tests = list(
    brownian_vs_ornstein_uhlenbeck = likelihood_ratio(bm_fit, ou_fit$fit),
    brownian_vs_early_burst = likelihood_ratio(bm_fit, eb_fit$fit),
    ornstein_uhlenbeck_vs_early_burst = likelihood_ratio(ou_fit$fit, eb_fit$fit)
  ),
  ancestral_reconstruction = list(
    brownian = bm_ancestral,
    early_burst = list(
      rate_change = eb_fit$parameter_value,
      node_count = eb_ancestral$node_count,
      first_five_estimates = eb_ancestral$first_five_estimates,
      recent_five_estimates = eb_ancestral$recent_five_estimates,
      rows = eb_ancestral$rows
    )
  ),
  baseline_gls = gls_payload(gls1, primate$longevity),
  fixed_lambda_gls_matches_baseline = isTRUE(all.equal(
    as.numeric(coef(gls1)),
    as.numeric(coef(pgls1)),
    tolerance = 1e-12
  )) && isTRUE(all.equal(
    as.numeric(logLik(gls1)),
    as.numeric(logLik(pgls1)),
    tolerance = 1e-12
  )),
  fixed_reference_lambda_pgls = gls_payload(
    pgls2_fixed,
    primate$longevity,
    lambda_value = estimated_lambda_value
  ),
  estimated_lambda_pgls = gls_payload(
    pgls2,
    primate$longevity,
    lambda_value = estimated_lambda_value
  ),
  signal_test = list(
    estimated_lambda = rounded(signal_lambda),
    estimated_log_likelihood = rounded(signal_log_likelihood),
    null_log_likelihood = rounded(signal_null_log_likelihood),
    likelihood_ratio = rounded(signal_likelihood_ratio),
    p_value = rounded(signal_p_value)
  ),
  coverage_boundaries = list(
    uncovered_fragments = c("mode-linked-intercept-models"),
    notes = c(
      "The lecture corBlomberg likelihood sweep remains outside the current canonical runtime parity surface.",
      "The governed evidence closes transformed-tree, fitContinuous, likelihood-ratio, and ancestral-state parity without overstating the remaining intercept-mode boundary."
    )
  )
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(
  toJSON(
    payload,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA
  ),
  con = out_path
)
