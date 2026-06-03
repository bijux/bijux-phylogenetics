args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("expected case-json path and output-root path")
}

case_path <- args[[1]]
output_root <- args[[2]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required for phytools parity execution")
}

write_payload <- function(path, payload) {
  writeLines(
    jsonlite::toJSON(
      payload,
      auto_unbox = TRUE,
      digits = 16,
      null = "null",
      pretty = TRUE
    ),
    con = path
  )
}

normalize_table_cell <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("")
  }
  if (length(value) > 1) {
    return(jsonlite::toJSON(unname(value), auto_unbox = TRUE))
  }
  value
}

rows_to_data_frame <- function(rows) {
  if (is.data.frame(rows)) {
    return(rows)
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  if (is.list(rows) && all(vapply(rows, is.list, logical(1)))) {
    normalized_rows <- lapply(rows, function(row) {
      normalized_row <- lapply(row, normalize_table_cell)
      as.data.frame(normalized_row, stringsAsFactors = FALSE)
    })
    return(do.call(rbind.data.frame, c(normalized_rows, list(stringsAsFactors = FALSE))))
  }
  as.data.frame(rows, stringsAsFactors = FALSE)
}

write_table <- function(path, rows) {
  utils::write.table(
    rows_to_data_frame(rows),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

format_table_value <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("")
  }
  if (is.list(value)) {
    return(jsonlite::toJSON(unname(value), auto_unbox = TRUE))
  }
  if (length(value) > 1) {
    return(jsonlite::toJSON(unname(value), auto_unbox = TRUE))
  }
  value
}

execution_path <- file.path(output_root, "reference-execution.json")
summary_path <- file.path(output_root, "reference-summary.json")
summary_table_path <- file.path(output_root, "reference-summary.tsv")
fitmk_rows_path <- file.path(output_root, "fitmk-rate-matrix.tsv")
stochastic_map_rows_path <- file.path(output_root, "stochastic-map-summary-rows.tsv")
rerooting_rows_path <- file.path(output_root, "rerooting-method-node-probabilities.tsv")
fast_anc_rows_path <- file.path(output_root, "fast-anc-node-estimates.tsv")
anc_ml_rows_path <- file.path(output_root, "anc-ml-node-estimates.tsv")
fastbm_rows_path <- file.path(output_root, "fastbm-summary-rows.tsv")
simcorrs_rows_path <- file.path(output_root, "simcorrs-summary-rows.tsv")
pgls_rows_path <- file.path(output_root, "pgls-summary-rows.tsv")
phyl_resid_rows_path <- file.path(output_root, "phyl-resid-summary-rows.tsv")
phyl_anova_rows_path <- file.path(output_root, "phyl-anova-summary-rows.tsv")
case_payload <- jsonlite::fromJSON(case_path)
r_version <- as.character(getRversion())

if (!requireNamespace("phytools", quietly = TRUE)) {
  write_payload(
    execution_path,
    list(
      status = "unavailable",
      mismatch_reason = "phytools_package_unavailable",
      r_version = r_version,
      phytools_version = NULL
    )
  )
  quit(save = "no", status = 0)
}

library(ape)
library(phytools)

input_fixtures <- unname(unlist(case_payload$input_fixtures))
tree <- ape::read.tree(input_fixtures[[1]])
traits_path <- if (length(input_fixtures) >= 2) input_fixtures[[2]] else NULL
trait_name <- case_payload$trait_name
taxon_column <- case_payload$taxon_column
trait_table <- NULL
taxon_names <- NULL
is_discrete_operation <- case_payload$operation %in% c(
  "discrete-fit-mk",
  "discrete-stochastic-map",
  "discrete-stochastic-map-count",
  "discrete-stochastic-map-description",
  "discrete-stochastic-map-density",
  "discrete-ancestral-rerooting"
)
if (!is.null(traits_path)) {
  trait_table <- utils::read.table(
    traits_path,
    header = TRUE,
    sep = if (grepl("\\.csv$", traits_path, ignore.case = TRUE)) "," else "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  taxon_names <- trait_table[[taxon_column]]
}
if (is_discrete_operation && !is.null(traits_path)) {
  raw_trait_values <- trait_table[[trait_name]]
  trait_values <- stats::setNames(
    ifelse(is.na(raw_trait_values), "", trimws(as.character(raw_trait_values))),
    taxon_names
  )
  kept_taxa <- names(trait_values)[
    nzchar(trait_values) &
      names(trait_values) %in% tree$tip.label
  ]
  trait_values <- trait_values[kept_taxa]
} else if (!is.null(traits_path)) {
  trait_values <- stats::setNames(as.numeric(trait_table[[trait_name]]), taxon_names)
  kept_taxa <- names(trait_values)[
    !is.na(trait_values) &
      names(trait_values) %in% tree$tip.label
  ]
  trait_values <- trait_values[kept_taxa]
} else {
  trait_values <- NULL
  kept_taxa <- tree$tip.label
}
excluded_taxa <- sort(setdiff(tree$tip.label, kept_taxa))
if (length(excluded_taxa) > 0) {
  tree <- ape::drop.tip(tree, excluded_taxa)
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

build_lambda_summary <- function(tree, trait_values, trait_name) {
  fit <- phytools::phylosig(
    tree,
    trait_values,
    method = "lambda",
    test = TRUE
  )
  list(
    taxon_count = length(trait_values),
    trait_name = trait_name,
    lambda_value = as.numeric(fit$lambda),
    log_likelihood = as.numeric(fit$logL),
    null_log_likelihood = as.numeric(fit$logL0),
    p_value = as.numeric(fit$P)
  )
}

build_k_summary <- function(tree, trait_values, trait_name) {
  requested_permutation_count <- case_payload$permutation_count
  permutation_seed <- case_payload$permutation_seed
  if (!is.null(permutation_seed)) {
    set.seed(as.integer(permutation_seed))
  }
  fit <- phytools::phylosig(
    tree,
    trait_values,
    method = "K",
    test = TRUE,
    nsim = as.integer(requested_permutation_count)
  )
  simulated_values <- as.numeric(fit$sim.K)
  observed_k <- as.numeric(fit$K)
  matching_index <- which(abs(simulated_values - observed_k) < 1e-12)[1]
  if (!is.na(matching_index)) {
    simulated_values <- simulated_values[-matching_index]
  }
  list(
    taxon_count = length(trait_values),
    trait_name = trait_name,
    k = observed_k,
    p_value = as.numeric(fit$P),
    permutation_count = as.integer(requested_permutation_count),
    permutation_seed = as.integer(permutation_seed),
    null_distribution_count = length(simulated_values),
    simulated_k_minimum = as.numeric(min(simulated_values)),
    simulated_k_mean = as.numeric(mean(simulated_values)),
    simulated_k_maximum = as.numeric(max(simulated_values))
  )
}

fit_aicc <- function(aic, sample_size, parameter_count) {
  denominator <- sample_size - parameter_count - 1
  if (denominator <= 0) {
    return(Inf)
  }
  aic + ((2 * parameter_count * (parameter_count + 1)) / denominator)
}

build_fitmk_result <- function(tree, trait_values, trait_name, excluded_taxa, discrete_model) {
  phytools_model <- switch(
    discrete_model,
    "equal-rates" = "ER",
    "symmetric" = "SYM",
    "all-rates-different" = "ARD",
    stop(paste("unsupported fitMk parity model:", discrete_model))
  )
  fit <- phytools::fitMk(
    tree,
    trait_values,
    model = phytools_model
  )
  q_matrix <- matrix(
    0,
    nrow = nrow(fit$index.matrix),
    ncol = ncol(fit$index.matrix),
    dimnames = list(fit$states, fit$states)
  )
  for (row_index in seq_len(nrow(fit$index.matrix))) {
    for (column_index in seq_len(ncol(fit$index.matrix))) {
      rate_index <- fit$index.matrix[row_index, column_index]
      if (!is.na(rate_index)) {
        q_matrix[row_index, column_index] <- fit$rates[[rate_index]]
      }
    }
  }
  diag(q_matrix) <- -rowSums(q_matrix)
  rate_rows <- list()
  for (source_state in rownames(q_matrix)) {
    for (target_state in colnames(q_matrix)) {
      if (identical(source_state, target_state)) {
        next
      }
      rate_rows[[length(rate_rows) + 1]] <- list(
        source_state = source_state,
        target_state = target_state,
        transition_allowed = TRUE,
        step_distance = 1,
        rate = unname(as.numeric(q_matrix[source_state, target_state]))
      )
    }
  }
  rate_rows <- rate_rows[order(
    vapply(rate_rows, function(row) row$source_state, character(1)),
    vapply(rate_rows, function(row) row$target_state, character(1))
  )]
  parameter_count <- length(unname(fit$rates))
  log_likelihood <- unname(as.numeric(stats::logLik(fit)))
  aic <- unname(as.numeric(stats::AIC(fit)))
  list(
    summary = list(
      taxon_count = length(trait_values),
      trait_name = trait_name,
      excluded_taxon_count = length(excluded_taxa),
      excluded_taxa = unname(as.list(excluded_taxa)),
      model = discrete_model,
      state_count = length(unique(unname(trait_values))),
      parameter_count = parameter_count,
      log_likelihood = log_likelihood,
      aic = aic,
      aicc = fit_aicc(aic, length(trait_values), parameter_count),
      overparameterized = parameter_count >= length(trait_values),
      baseline_model = if (identical(discrete_model, "equal-rates")) NULL else "equal-rates",
      preferred_model_by_aic = if (identical(discrete_model, "equal-rates")) {
        NULL
      } else {
        baseline_fit <- phytools::fitMk(
          tree,
          trait_values,
          model = "ER"
        )
        if (aic <= unname(as.numeric(stats::AIC(baseline_fit)))) {
          discrete_model
        } else {
          "equal-rates"
        }
      }
    ),
    rows = rate_rows
  )
}

normalize_pgls_term_component <- function(component, data) {
  if (identical(component, "(Intercept)")) {
    return("intercept")
  }
  if (component %in% names(data)) {
    return(component)
  }
  for (column_name in names(data)) {
    if (!is.factor(data[[column_name]])) {
      next
    }
    levels_to_encode <- levels(data[[column_name]])
    if (length(levels_to_encode) <= 1) {
      next
    }
    for (level_name in levels_to_encode[-1]) {
      if (identical(component, paste0(column_name, level_name))) {
        return(paste0(column_name, "[", level_name, "]"))
      }
    }
  }
  component
}

normalize_pgls_column_name <- function(column_name, data) {
  if (identical(column_name, "(Intercept)")) {
    return("intercept")
  }
  components <- strsplit(column_name, ":", fixed = TRUE)[[1]]
  paste(
    vapply(
      components,
      normalize_pgls_term_component,
      character(1),
      data = data
    ),
    collapse = ":"
  )
}

build_pgls_result <- function(tree, trait_table, taxon_column, formula_string, lambda_value) {
  formula_object <- stats::as.formula(formula_string)
  formula_variables <- all.vars(formula_object)
  pgls_data <- trait_table
  rownames(pgls_data) <- pgls_data[[taxon_column]]
  pgls_data[[taxon_column]] <- NULL
  if (length(formula_variables) > 0) {
    complete_rows <- stats::complete.cases(pgls_data[, formula_variables, drop = FALSE])
    pgls_data <- pgls_data[complete_rows, , drop = FALSE]
  }
  pgls_data <- pgls_data[rownames(pgls_data) %in% tree$tip.label, , drop = FALSE]
  kept_taxa <- rownames(pgls_data)
  excluded_taxa <- sort(setdiff(tree$tip.label, kept_taxa))
  if (length(excluded_taxa) > 0) {
    tree <- ape::drop.tip(tree, excluded_taxa)
  }
  for (column_name in names(pgls_data)) {
    if (!is.numeric(pgls_data[[column_name]])) {
      pgls_data[[column_name]] <- factor(pgls_data[[column_name]])
    }
  }
  fit <- phytools::pgls.SEy(
    formula_object,
    data = pgls_data,
    tree = tree,
    method = "ML"
  )
  coefficient_table <- as.data.frame(summary(fit)$tTable, stringsAsFactors = FALSE)
  coefficient_table$coefficient_name <- vapply(
    rownames(coefficient_table),
    normalize_pgls_column_name,
    character(1),
    data = pgls_data
  )
  model_matrix <- stats::model.matrix(formula_object, data = pgls_data)
  normalized_column_names <- vapply(
    colnames(model_matrix),
    normalize_pgls_column_name,
    character(1),
    data = pgls_data
  )
  colnames(model_matrix) <- normalized_column_names
  term_labels <- attr(stats::terms(formula_object), "term.labels")
  predictor_components <- unique(unlist(strsplit(term_labels, ":", fixed = TRUE)))
  categorical_predictor_count <- sum(
    predictor_components %in% names(pgls_data) &
      vapply(
        predictor_components,
        function(component_name) is.factor(pgls_data[[component_name]]),
        logical(1)
      )
  )
  rows <- list()
  for (row_index in seq_len(nrow(coefficient_table))) {
    coefficient_name <- coefficient_table$coefficient_name[[row_index]]
    rows[[length(rows) + 1]] <- list(
      row_kind = "coefficient_estimate",
      label = coefficient_name,
      value = unname(as.numeric(coefficient_table$Value[[row_index]]))
    )
    rows[[length(rows) + 1]] <- list(
      row_kind = "coefficient_standard_error",
      label = coefficient_name,
      value = unname(as.numeric(coefficient_table$Std.Error[[row_index]]))
    )
    rows[[length(rows) + 1]] <- list(
      row_kind = "coefficient_p_value",
      label = coefficient_name,
      value = unname(as.numeric(coefficient_table[["p-value"]][[row_index]]))
    )
  }
  for (taxon_name in rownames(model_matrix)) {
    for (column_name in colnames(model_matrix)) {
      rows[[length(rows) + 1]] <- list(
        row_kind = "model_matrix",
        label = paste0(taxon_name, ":", column_name),
        value = unname(as.numeric(model_matrix[taxon_name, column_name]))
      )
    }
  }
  rows <- rows[order(
    vapply(rows, function(row) row$row_kind, character(1)),
    vapply(rows, function(row) row$label, character(1))
  )]
  list(
    summary = list(
      taxon_count = nrow(pgls_data),
      trait_name = formula_variables[[1]],
      formula = formula_string,
      analysis_taxon_count = nrow(pgls_data),
      coefficient_count = nrow(coefficient_table),
      model_matrix_row_count = nrow(model_matrix),
      model_matrix_column_count = ncol(model_matrix),
      categorical_predictor_count = categorical_predictor_count,
      interaction_term_count = sum(grepl(":", term_labels, fixed = TRUE)),
      lambda_value = as.numeric(lambda_value),
      lambda_estimation_mode = "fixed",
      log_likelihood = unname(as.numeric(stats::logLik(fit))),
      aic = unname(as.numeric(stats::AIC(fit)))
    ),
    rows = rows
  )
}

build_phyl_resid_result <- function(tree, trait_table, taxon_column, trait_name, predictor_name, lambda_value) {
  taxa_in_tree <- trait_table[[taxon_column]] %in% tree$tip.label
  complete_rows <- stats::complete.cases(
    trait_table[, c(taxon_column, predictor_name, trait_name), drop = FALSE]
  )
  residual_data <- trait_table[taxa_in_tree & complete_rows, , drop = FALSE]
  kept_taxa <- residual_data[[taxon_column]]
  missing_tree_taxa <- setdiff(tree$tip.label, trait_table[[taxon_column]])
  missing_value_taxa <- trait_table[[taxon_column]][
    trait_table[[taxon_column]] %in% tree$tip.label &
      !complete_rows
  ]
  absent_from_tree_taxa <- trait_table[[taxon_column]][!taxa_in_tree]
  excluded_taxa <- sort(unique(c(
    missing_tree_taxa,
    missing_value_taxa,
    absent_from_tree_taxa
  )))
  pruned_tree <- tree
  if (length(excluded_taxa) > 0) {
    tree_taxa_to_drop <- intersect(pruned_tree$tip.label, excluded_taxa)
    if (length(tree_taxa_to_drop) > 0) {
      pruned_tree <- ape::drop.tip(pruned_tree, tree_taxa_to_drop)
    }
  }
  x <- stats::setNames(
    as.numeric(residual_data[[predictor_name]]),
    residual_data[[taxon_column]]
  )
  y <- stats::setNames(
    as.numeric(residual_data[[trait_name]]),
    residual_data[[taxon_column]]
  )
  method_name <- if (
    !is.null(lambda_value) && abs(as.numeric(lambda_value) - 1.0) < 1e-12
  ) {
    "BM"
  } else {
    "lambda"
  }
  fit <- phytools::phyl.resid(pruned_tree, x, y, method = method_name)
  fitted_values <- as.vector(cbind(1, x[pruned_tree$tip.label]) %*% fit$beta)
  names(fitted_values) <- pruned_tree$tip.label
  residual_values <- as.vector(fit$resid)
  names(residual_values) <- pruned_tree$tip.label
  rows <- list(
    list(
      row_kind = "coefficient_estimate",
      label = "intercept",
      value = unname(as.numeric(fit$beta[1, 1]))
    ),
    list(
      row_kind = "coefficient_estimate",
      label = predictor_name,
      value = unname(as.numeric(fit$beta[2, 1]))
    )
  )
  for (taxon_name in residual_data[[taxon_column]]) {
    rows[[length(rows) + 1]] <- list(
      row_kind = "taxon_value",
      label = taxon_name,
      observed_value = unname(as.numeric(y[[taxon_name]])),
      fitted_value = unname(as.numeric(fitted_values[[taxon_name]])),
      residual = unname(as.numeric(residual_values[[taxon_name]]))
    )
  }
  rows <- rows[order(
    vapply(rows, function(row) row$row_kind, character(1)),
    vapply(rows, function(row) row$label, character(1))
  )]
  summary <- list(
    taxon_count = nrow(residual_data),
    trait_name = trait_name,
    predictor_name = predictor_name,
    method = if (identical(method_name, "BM")) "brownian" else "lambda",
    excluded_taxon_count = length(excluded_taxa),
    excluded_taxa = unname(as.list(excluded_taxa))
  )
  if (identical(method_name, "lambda")) {
    summary$lambda_value <- unname(as.numeric(fit$lambda[[1]]))
    summary$log_likelihood <- unname(as.numeric(fit$logL[[1]]))
  }
  list(summary = summary, rows = rows)
}

build_phyl_anova_result <- function(tree, trait_table, taxon_column, trait_name, group_column, nsim, seed) {
  taxa_in_tree <- trait_table[[taxon_column]] %in% tree$tip.label
  complete_rows <- stats::complete.cases(
    trait_table[, c(taxon_column, group_column, trait_name), drop = FALSE]
  )
  anova_data <- trait_table[taxa_in_tree & complete_rows, , drop = FALSE]
  analyzed_taxa <- anova_data[[taxon_column]]
  missing_tree_taxa <- setdiff(tree$tip.label, trait_table[[taxon_column]])
  missing_value_taxa <- trait_table[[taxon_column]][
    trait_table[[taxon_column]] %in% tree$tip.label &
      !complete_rows
  ]
  absent_from_tree_taxa <- trait_table[[taxon_column]][!taxa_in_tree]
  excluded_taxa <- sort(unique(c(
    missing_tree_taxa,
    missing_value_taxa,
    absent_from_tree_taxa
  )))
  pruned_tree <- tree
  tree_taxa_to_drop <- setdiff(pruned_tree$tip.label, analyzed_taxa)
  if (length(tree_taxa_to_drop) > 0) {
    pruned_tree <- ape::drop.tip(pruned_tree, tree_taxa_to_drop)
  }
  x <- stats::setNames(as.character(anova_data[[group_column]]), analyzed_taxa)
  y <- stats::setNames(as.numeric(anova_data[[trait_name]]), analyzed_taxa)
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  fit <- phytools::phylANOVA(
    pruned_tree,
    x,
    y,
    nsim = as.integer(nsim),
    posthoc = TRUE,
    p.adj = "holm"
  )
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  x_factor <- as.factor(x[pruned_tree$tip.label])
  y_values <- y[pruned_tree$tip.label]
  sigma_squared <- mean(pic(y_values, multi2di(pruned_tree, random = FALSE))^2)
  group_levels <- levels(x_factor)
  group_count <- length(group_levels)
  observed_t <- tTests(x_factor, y_values)
  simulated_values <- fastBM(pruned_tree, sig2 = sigma_squared, nsim = (nsim - 1))
  null_t <- array(
    NA,
    dim = c(group_count, group_count, nsim),
    dimnames = list(group_levels, group_levels, NULL)
  )
  null_t[, , 1] <- observed_t
  for (simulation_index in 2:nsim) {
    null_t[, , simulation_index] <- tTests(
      x_factor,
      simulated_values[, simulation_index - 1]
    )
  }
  uncorrected_p <- matrix(
    NA,
    group_count,
    group_count,
    dimnames = list(group_levels, group_levels)
  )
  for (left_index in seq_len(group_count)) {
    for (right_index in left_index:group_count) {
      p_value <- sum(
        abs(null_t[left_index, right_index, ]) >=
          abs(observed_t[left_index, right_index])
      ) / nsim
      uncorrected_p[left_index, right_index] <- p_value
      uncorrected_p[right_index, left_index] <- p_value
    }
  }
  group_rows <- list()
  group_sizes <- stats::setNames(numeric(), character())
  for (group_name in group_levels) {
    group_taxa <- analyzed_taxa[x[analyzed_taxa] == group_name]
    group_values <- as.numeric(y[group_taxa])
    group_sizes[[group_name]] <- length(group_taxa)
    group_rows[[length(group_rows) + 1]] <- list(
      row_kind = "group_summary",
      label = group_name,
      taxon_count = length(group_taxa),
      taxa = paste(group_taxa, collapse = ","),
      mean = as.numeric(mean(group_values)),
      variance = as.numeric(stats::var(group_values)),
      minimum = as.numeric(min(group_values)),
      maximum = as.numeric(max(group_values))
    )
  }
  pairwise_rows <- list()
  for (left_index in seq_len(group_count - 1)) {
    for (right_index in (left_index + 1):group_count) {
      left_group <- group_levels[[left_index]]
      right_group <- group_levels[[right_index]]
      pairwise_rows[[length(pairwise_rows) + 1]] <- list(
        row_kind = "pairwise_comparison",
        label = paste0(left_group, "|", right_group),
        left_taxon_count = group_sizes[[left_group]],
        right_taxon_count = group_sizes[[right_group]],
        observed_t_statistic = as.numeric(observed_t[left_group, right_group]),
        uncorrected_p_value = as.numeric(uncorrected_p[left_group, right_group]),
        adjusted_p_value = as.numeric(fit$Pt[left_group, right_group])
      )
    }
  }
  rows <- c(group_rows, pairwise_rows)
  rows <- rows[order(
    vapply(rows, function(row) row$row_kind, character(1)),
    vapply(rows, function(row) row$label, character(1))
  )]
  list(
    summary = list(
      taxon_count = length(analyzed_taxa),
      trait_name = trait_name,
      group_column = group_column,
      excluded_taxon_count = length(excluded_taxa),
      excluded_taxa = unname(as.list(excluded_taxa)),
      group_count = group_count,
      simulation_count = as.integer(nsim),
      seed = as.integer(seed),
      pairwise_adjustment_method = fit$method,
      brownian_sigma_squared = as.numeric(sigma_squared),
      sum_of_squares_between = as.numeric(fit$"Sum Sq"[[1]]),
      sum_of_squares_within = as.numeric(fit$"Sum Sq"[[2]]),
      mean_square_between = as.numeric(fit$"Mean Sq"[[1]]),
      mean_square_within = as.numeric(fit$"Mean Sq"[[2]]),
      f_statistic = as.numeric(fit$F),
      p_value = as.numeric(fit$Pf),
      low_sample_group_count = sum(unname(group_sizes) < 3)
    ),
    rows = rows
  )
}

build_simmap_summary_rows <- function(fit, description, include_branch_occupancy) {
  count_table <- as.data.frame(description$count, stringsAsFactors = FALSE)
  time_table <- as.data.frame(description$times, stringsAsFactors = FALSE)
  transition_columns <- setdiff(colnames(count_table), "N")
  state_columns <- setdiff(colnames(time_table), "total")
  summary_rows <- list()
  for (transition_label in transition_columns) {
    values <- as.numeric(count_table[[transition_label]])
    summary_rows[[length(summary_rows) + 1]] <- list(
      row_kind = "transition_count",
      label = gsub(",", "->", transition_label, fixed = TRUE),
      mean_value = as.numeric(mean(values)),
      lower_95_interval = as.numeric(stats::quantile(values, probs = 0.025, names = FALSE)),
      upper_95_interval = as.numeric(stats::quantile(values, probs = 0.975, names = FALSE)),
      presence_fraction = as.numeric(mean(values > 0))
    )
  }
  for (state_label in state_columns) {
    values <- as.numeric(time_table[[state_label]])
    summary_rows[[length(summary_rows) + 1]] <- list(
      row_kind = "state_time",
      label = state_label,
      mean_value = as.numeric(mean(values)),
      lower_95_interval = as.numeric(stats::quantile(values, probs = 0.025, names = FALSE)),
      upper_95_interval = as.numeric(stats::quantile(values, probs = 0.975, names = FALSE)),
      presence_fraction = 1.0
    )
  }
  if (include_branch_occupancy) {
    reference_tree <- fit[[1]]
    state_order <- colnames(reference_tree$mapped.edge)
    for (edge_index in seq_len(nrow(reference_tree$edge))) {
      parent_node <- node_signature(reference_tree, reference_tree$edge[edge_index, 1])
      child_node <- node_signature(reference_tree, reference_tree$edge[edge_index, 2])
      for (state_label in state_order) {
        values <- vapply(
          fit,
          function(map_tree) {
            if (!(state_label %in% colnames(map_tree$mapped.edge))) {
              return(0.0)
            }
            unname(as.numeric(map_tree$mapped.edge[edge_index, state_label]))
          },
          numeric(1)
        )
        summary_rows[[length(summary_rows) + 1]] <- list(
          row_kind = "branch_state_occupancy",
          label = paste0(parent_node, "->", child_node, ":", state_label),
          mean_value = as.numeric(mean(values)),
          lower_95_interval = as.numeric(stats::quantile(values, probs = 0.025, names = FALSE)),
          upper_95_interval = as.numeric(stats::quantile(values, probs = 0.975, names = FALSE)),
          presence_fraction = as.numeric(mean(values > 0))
        )
      }
    }
  }
  summary_rows[order(
    vapply(summary_rows, function(row) row$row_kind, character(1)),
    vapply(summary_rows, function(row) row$label, character(1))
  )]
}

build_countsimmap_summary_rows <- function(count_matrix) {
  count_table <- as.data.frame(count_matrix, stringsAsFactors = FALSE)
  transition_columns <- setdiff(colnames(count_table), "N")
  summary_rows <- list()
  for (transition_label in transition_columns) {
    values <- as.numeric(count_table[[transition_label]])
    summary_rows[[length(summary_rows) + 1]] <- list(
      row_kind = "transition_count",
      label = gsub(",", "->", transition_label, fixed = TRUE),
      mean_value = as.numeric(mean(values)),
      lower_95_interval = as.numeric(stats::quantile(values, probs = 0.025, names = FALSE)),
      upper_95_interval = as.numeric(stats::quantile(values, probs = 0.975, names = FALSE)),
      presence_fraction = as.numeric(mean(values > 0))
    )
  }
  summary_rows[order(
    vapply(summary_rows, function(row) row$row_kind, character(1)),
    vapply(summary_rows, function(row) row$label, character(1))
  )]
}

binary_entropy <- function(probability) {
  bounded <- min(max(as.numeric(probability), 0.0), 1.0)
  if (bounded <= 0.0 || bounded >= 1.0) {
    return(0.0)
  }
  as.numeric(-bounded * log(bounded, 2) - (1.0 - bounded) * log(1.0 - bounded, 2))
}

build_densitymap_branch_rows <- function(density_map) {
  reference_tree <- density_map$tree
  branch_rows <- list()
  for (edge_index in seq_len(nrow(reference_tree$edge))) {
    branch_segments <- as.numeric(reference_tree$maps[[edge_index]])
    branch_probabilities <- as.numeric(names(reference_tree$maps[[edge_index]])) / 1000.0
    branch_length <- sum(branch_segments)
    if (branch_length <= 0.0) {
      next
    }
    parent_node <- node_signature(reference_tree, reference_tree$edge[edge_index, 1])
    child_node <- node_signature(reference_tree, reference_tree$edge[edge_index, 2])
    branch_rows[[length(branch_rows) + 1]] <- list(
      label = paste0(parent_node, "->", child_node),
      mean_posterior_probability = as.numeric(sum(branch_segments * branch_probabilities) / branch_length),
      minimum_posterior_probability = as.numeric(min(branch_probabilities)),
      maximum_posterior_probability = as.numeric(max(branch_probabilities)),
      uncertainty = as.numeric(sum(branch_segments * vapply(branch_probabilities, binary_entropy, numeric(1))) / branch_length),
      slice_count = length(branch_probabilities)
    )
  }
  branch_rows[order(vapply(branch_rows, function(row) row$label, character(1)))]
}

build_sim_history_summary_rows <- function(simmap_result) {
  description <- phytools::describe.simmap(simmap_result, plot = FALSE)
  summary_rows <- build_simmap_summary_rows(
    simmap_result,
    description,
    include_branch_occupancy = FALSE
  )
  reference_tree <- if (inherits(simmap_result, "multiSimmap")) simmap_result[[1]] else simmap_result
  simmap_list <- if (inherits(simmap_result, "simmap")) list(simmap_result) else simmap_result
  state_order <- sort(unique(unname(unlist(lapply(simmap_list, function(map_tree) map_tree$states)))))
  tip_rows <- list()
  for (taxon in reference_tree$tip.label) {
    for (state_label in state_order) {
      values <- vapply(
        simmap_list,
        function(map_tree) as.numeric(map_tree$states[[taxon]] == state_label),
        numeric(1)
      )
      tip_rows[[length(tip_rows) + 1]] <- list(
        row_kind = "tip_state_frequency",
        label = paste0(taxon, ":", state_label),
        mean_value = as.numeric(mean(values)),
        lower_95_interval = as.numeric(stats::quantile(values, probs = 0.025, names = FALSE)),
        upper_95_interval = as.numeric(stats::quantile(values, probs = 0.975, names = FALSE)),
        presence_fraction = as.numeric(mean(values))
      )
    }
  }
  rows <- c(summary_rows, tip_rows)
  rows[order(
    vapply(rows, function(row) row$row_kind, character(1)),
    vapply(rows, function(row) row$label, character(1))
  )]
}

build_sim_history_result <- function(tree, trait_name) {
  states <- unname(unlist(case_payload$simulation_states))
  rate_rows <- rows_to_data_frame(case_payload$simulation_rate_rows)
  q_matrix <- matrix(
    0,
    nrow = length(states),
    ncol = length(states),
    dimnames = list(states, states)
  )
  for (row_index in seq_len(nrow(rate_rows))) {
    source_state <- as.character(rate_rows[row_index, "source_state"])
    target_state <- as.character(rate_rows[row_index, "target_state"])
    rate <- as.numeric(rate_rows[row_index, "rate"])
    q_matrix[source_state, target_state] <- rate
  }
  diag(q_matrix) <- -rowSums(q_matrix)
  anc <- case_payload$simulation_root_state
  if (is.null(anc) && !is.null(case_payload$simulation_root_state_probabilities)) {
    anc <- unlist(case_payload$simulation_root_state_probabilities)[states]
  }
  requested_replicate_count <- as.integer(case_payload$simulation_replicate_count)
  simulation_seed <- as.integer(case_payload$simulation_seed)
  set.seed(simulation_seed)
  simmap_result <- phytools::sim.history(
    tree,
    q_matrix,
    anc = anc,
    nsim = requested_replicate_count,
    direction = "row_to_column",
    message = FALSE
  )
  description <- phytools::describe.simmap(simmap_result, plot = FALSE)
  count_table <- as.data.frame(description$count, stringsAsFactors = FALSE)
  total_transition_counts <- as.numeric(count_table$N)
  list(
    summary = list(
      taxon_count = length(tree$tip.label),
      trait_name = trait_name,
      branch_count = nrow(tree$edge),
      state_count = length(states),
      requested_replicate_count = requested_replicate_count,
      successful_replicate_count = requested_replicate_count,
      fixed_root_state = case_payload$simulation_root_state,
      root_prior_probabilities = case_payload$simulation_root_state_probabilities,
      seed = simulation_seed,
      mean_total_transition_count = as.numeric(mean(total_transition_counts)),
      lower_95_total_transition_count = as.numeric(stats::quantile(total_transition_counts, probs = 0.025, names = FALSE)),
      upper_95_total_transition_count = as.numeric(stats::quantile(total_transition_counts, probs = 0.975, names = FALSE))
    ),
    rows = build_sim_history_summary_rows(simmap_result)
  )
}

normalize_fastbm_matrix <- function(simulated, tip_labels) {
  if (is.null(dim(simulated))) {
    matrix_result <- matrix(
      as.numeric(simulated),
      nrow = length(tip_labels),
      ncol = 1,
      dimnames = list(names(simulated), "simulation_1")
    )
  } else {
    matrix_result <- as.matrix(simulated)
  }
  if (!is.null(rownames(matrix_result)) && identical(rownames(matrix_result), tip_labels)) {
    return(matrix_result)
  }
  if (!is.null(colnames(matrix_result)) && identical(colnames(matrix_result), tip_labels)) {
    return(t(matrix_result))
  }
  rownames(matrix_result) <- tip_labels
  matrix_result
}

build_fastbm_summary_rows <- function(simulated_matrix) {
  rows <- list()
  tip_frame <- t(simulated_matrix)
  for (taxon in rownames(simulated_matrix)) {
    values <- as.numeric(simulated_matrix[taxon, ])
    rows[[length(rows) + 1]] <- list(
      row_kind = "tip_distribution",
      label = taxon,
      mean_value = as.numeric(mean(values)),
      standard_deviation = as.numeric(stats::sd(values)),
      minimum = as.numeric(min(values)),
      median = as.numeric(stats::median(values)),
      maximum = as.numeric(max(values)),
      covariance = "",
      correlation = ""
    )
  }
  tip_labels <- rownames(simulated_matrix)
  for (left_index in seq_along(tip_labels)) {
    for (right_index in left_index:length(tip_labels)) {
      left_taxon <- tip_labels[[left_index]]
      right_taxon <- tip_labels[[right_index]]
      left_values <- tip_frame[, left_taxon]
      right_values <- tip_frame[, right_taxon]
      correlation <- if (stats::sd(left_values) == 0 || stats::sd(right_values) == 0) {
        0.0
      } else {
        as.numeric(stats::cor(left_values, right_values))
      }
      rows[[length(rows) + 1]] <- list(
        row_kind = "tip_covariance",
        label = paste0(left_taxon, "|", right_taxon),
        mean_value = "",
        standard_deviation = "",
        minimum = "",
        median = "",
        maximum = "",
        covariance = as.numeric(stats::cov(left_values, right_values)),
        correlation = correlation
      )
    }
  }
  rows[order(
    vapply(rows, function(row) row$row_kind, character(1)),
    vapply(rows, function(row) row$label, character(1))
  )]
}

build_fastbm_result <- function(tree) {
  requested_replicate_count <- as.integer(case_payload$continuous_replicate_count)
  simulation_seed <- as.integer(case_payload$continuous_seed)
  root_state <- as.numeric(case_payload$continuous_root_state)
  sigma_squared <- as.numeric(case_payload$continuous_sigma_squared)
  set.seed(simulation_seed)
  simulated <- phytools::fastBM(
    tree,
    a = root_state,
    sig2 = sigma_squared,
    nsim = requested_replicate_count,
    internal = FALSE
  )
  simulated_matrix <- normalize_fastbm_matrix(simulated, tree$tip.label)
  list(
    summary = list(
      taxon_count = length(tree$tip.label),
      branch_count = nrow(tree$edge),
      requested_replicate_count = requested_replicate_count,
      successful_replicate_count = requested_replicate_count,
      seed = simulation_seed,
      root_state = root_state,
      sigma_squared = sigma_squared
    ),
    rows = build_fastbm_summary_rows(simulated_matrix)
  )
}

normalize_simcorrs_matrix <- function(simulated_matrix, tip_labels, trait_names) {
  matrix_result <- as.matrix(simulated_matrix)
  rownames(matrix_result) <- tip_labels
  if (is.null(colnames(matrix_result)) || length(colnames(matrix_result)) != length(trait_names)) {
    colnames(matrix_result) <- trait_names
  }
  matrix_result[tip_labels, trait_names, drop = FALSE]
}

build_simcorrs_summary_rows <- function(simulated_matrices, trait_names, root_states, covariance_matrix) {
  rows <- list()
  for (index in seq_along(trait_names)) {
    rows[[length(rows) + 1]] <- list(
      row_kind = "root_state",
      label = trait_names[[index]],
      mean_value = as.numeric(root_states[[index]]),
      standard_deviation = "",
      minimum = "",
      median = "",
      maximum = "",
      covariance = "",
      correlation = ""
    )
  }
  for (left_index in seq_along(trait_names)) {
    for (right_index in left_index:length(trait_names)) {
      covariance <- as.numeric(covariance_matrix[left_index, right_index])
      left_variance <- as.numeric(covariance_matrix[left_index, left_index])
      right_variance <- as.numeric(covariance_matrix[right_index, right_index])
      correlation <- if (left_variance <= 0 || right_variance <= 0) {
        0.0
      } else {
        covariance / sqrt(left_variance * right_variance)
      }
      rows[[length(rows) + 1]] <- list(
        row_kind = "evolutionary_covariance",
        label = paste0(trait_names[[left_index]], "|", trait_names[[right_index]]),
        mean_value = "",
        standard_deviation = "",
        minimum = "",
        median = "",
        maximum = "",
        covariance = covariance,
        correlation = as.numeric(correlation)
      )
    }
  }
  dimension_labels <- list()
  for (taxon in rownames(simulated_matrices[[1]])) {
    for (trait_name in trait_names) {
      dimension_labels[[length(dimension_labels) + 1]] <- c(taxon, trait_name)
    }
  }
  values_by_dimension <- lapply(dimension_labels, function(label_parts) {
    vapply(
      simulated_matrices,
      function(simulated_matrix) as.numeric(simulated_matrix[label_parts[[1]], label_parts[[2]]]),
      numeric(1)
    )
  })
  for (index in seq_along(dimension_labels)) {
    label_parts <- dimension_labels[[index]]
    values <- values_by_dimension[[index]]
    rows[[length(rows) + 1]] <- list(
      row_kind = "tip_distribution",
      label = paste0(label_parts[[1]], "|", label_parts[[2]]),
      mean_value = as.numeric(mean(values)),
      standard_deviation = as.numeric(stats::sd(values)),
      minimum = as.numeric(min(values)),
      median = as.numeric(stats::median(values)),
      maximum = as.numeric(max(values)),
      covariance = "",
      correlation = ""
    )
  }
  for (left_index in seq_along(dimension_labels)) {
    left_label <- dimension_labels[[left_index]]
    left_values <- values_by_dimension[[left_index]]
    for (right_index in left_index:length(dimension_labels)) {
      right_label <- dimension_labels[[right_index]]
      right_values <- values_by_dimension[[right_index]]
      correlation <- if (stats::sd(left_values) == 0 || stats::sd(right_values) == 0) {
        0.0
      } else {
        as.numeric(stats::cor(left_values, right_values))
      }
      rows[[length(rows) + 1]] <- list(
        row_kind = "tip_covariance",
        label = paste0(
          left_label[[1]], "|", left_label[[2]], "||",
          right_label[[1]], "|", right_label[[2]]
        ),
        mean_value = "",
        standard_deviation = "",
        minimum = "",
        median = "",
        maximum = "",
        covariance = as.numeric(stats::cov(left_values, right_values)),
        correlation = correlation
      )
    }
  }
  rows[order(
    vapply(rows, function(row) row$row_kind, character(1)),
    vapply(rows, function(row) row$label, character(1))
  )]
}

build_simcorrs_result <- function(tree) {
  requested_replicate_count <- as.integer(case_payload$continuous_replicate_count)
  simulation_seed <- as.integer(case_payload$continuous_seed)
  trait_names <- as.character(unname(case_payload$continuous_trait_names))
  root_states <- as.numeric(unname(case_payload$continuous_root_states))
  covariance_values <- as.numeric(unlist(case_payload$continuous_covariance_matrix))
  covariance_matrix <- matrix(
    covariance_values,
    nrow = length(trait_names),
    byrow = TRUE
  )
  rownames(covariance_matrix) <- trait_names
  colnames(covariance_matrix) <- trait_names
  set.seed(simulation_seed)
  simulated_matrices <- lapply(seq_len(requested_replicate_count), function(index) {
    normalize_simcorrs_matrix(
      phytools::sim.corrs(
        tree,
        vcv = covariance_matrix,
        anc = root_states,
        internal = FALSE
      ),
      tree$tip.label,
      trait_names
    )
  })
  list(
    summary = list(
      taxon_count = length(tree$tip.label),
      branch_count = nrow(tree$edge),
      trait_count = length(trait_names),
      requested_replicate_count = requested_replicate_count,
      successful_replicate_count = requested_replicate_count,
      seed = simulation_seed
    ),
    rows = build_simcorrs_summary_rows(
      simulated_matrices,
      trait_names,
      root_states,
      covariance_matrix
    )
  )
}

build_make_simmap_result <- function(tree, trait_values, trait_name, excluded_taxa, discrete_model) {
  phytools_model <- switch(
    discrete_model,
    "equal-rates" = "ER",
    "symmetric" = "SYM",
    "all-rates-different" = "ARD",
    stop(paste("unsupported make.simmap parity model:", discrete_model))
  )
  requested_replicate_count <- as.integer(case_payload$stochastic_map_replicate_count)
  stochastic_map_seed <- as.integer(case_payload$stochastic_map_seed)
  fitmk_fit <- phytools::fitMk(
    tree,
    trait_values,
    model = phytools_model
  )
  parameter_count <- length(unname(fitmk_fit$rates))
  log_likelihood <- unname(as.numeric(stats::logLik(fitmk_fit)))
  aic <- unname(as.numeric(stats::AIC(fitmk_fit)))
  baseline_fit <- NULL
  if (!identical(discrete_model, "equal-rates")) {
    baseline_fit <- phytools::fitMk(
      tree,
      trait_values,
      model = "ER"
    )
  }
  set.seed(stochastic_map_seed)
  fit <- phytools::make.simmap(
    tree,
    trait_values,
    model = phytools_model,
    nsim = requested_replicate_count,
    pi = "equal",
    message = FALSE
  )
  description <- phytools::describe.simmap(fit, plot = FALSE)
  count_table <- as.data.frame(description$count, stringsAsFactors = FALSE)
  total_transition_counts <- as.numeric(count_table$N)
  list(
    summary = list(
      taxon_count = length(trait_values),
      trait_name = trait_name,
      excluded_taxon_count = length(excluded_taxa),
      excluded_taxa = unname(as.list(excluded_taxa)),
      model = discrete_model,
      state_count = length(unique(unname(trait_values))),
      parameter_count = parameter_count,
      log_likelihood = log_likelihood,
      aic = aic,
      aicc = fit_aicc(aic, length(trait_values), parameter_count),
      overparameterized = parameter_count >= length(trait_values),
      baseline_model = if (is.null(baseline_fit)) NULL else "equal-rates",
      preferred_model_by_aic = if (is.null(baseline_fit)) {
        NULL
      } else if (aic <= unname(as.numeric(stats::AIC(baseline_fit)))) {
        discrete_model
      } else {
        "equal-rates"
      },
      requested_replicate_count = requested_replicate_count,
      successful_replicate_count = requested_replicate_count,
      simulation_failure_count = 0L,
      conditioned_on_node_estimates = FALSE,
      seed = stochastic_map_seed,
      mean_total_transition_count = as.numeric(mean(total_transition_counts)),
      lower_95_total_transition_count = as.numeric(stats::quantile(total_transition_counts, probs = 0.025, names = FALSE)),
      upper_95_total_transition_count = as.numeric(stats::quantile(total_transition_counts, probs = 0.975, names = FALSE))
    ),
    rows = build_simmap_summary_rows(
      fit,
      description,
      include_branch_occupancy = FALSE
    )
  )
}

build_describe_simmap_result <- function(tree, trait_values, trait_name, excluded_taxa, discrete_model) {
  simmap_result <- build_make_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    discrete_model
  )
  phytools_model <- switch(
    discrete_model,
    "equal-rates" = "ER",
    "symmetric" = "SYM",
    "all-rates-different" = "ARD",
    stop(paste("unsupported describe.simmap parity model:", discrete_model))
  )
  requested_replicate_count <- as.integer(case_payload$stochastic_map_replicate_count)
  stochastic_map_seed <- as.integer(case_payload$stochastic_map_seed)
  set.seed(stochastic_map_seed)
  fit <- phytools::make.simmap(
    tree,
    trait_values,
    model = phytools_model,
    nsim = requested_replicate_count,
    pi = "equal",
    message = FALSE
  )
  description <- phytools::describe.simmap(fit, plot = FALSE)
  simmap_result$summary$branch_count <- nrow(fit[[1]]$edge)
  simmap_result$rows <- build_simmap_summary_rows(
    fit,
    description,
    include_branch_occupancy = TRUE
  )
  simmap_result
}

build_count_simmap_result <- function(tree, trait_values, trait_name, excluded_taxa, discrete_model) {
  simmap_result <- build_make_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    discrete_model
  )
  phytools_model <- switch(
    discrete_model,
    "equal-rates" = "ER",
    "symmetric" = "SYM",
    "all-rates-different" = "ARD",
    stop(paste("unsupported countSimmap parity model:", discrete_model))
  )
  requested_replicate_count <- as.integer(case_payload$stochastic_map_replicate_count)
  stochastic_map_seed <- as.integer(case_payload$stochastic_map_seed)
  set.seed(stochastic_map_seed)
  fit <- phytools::make.simmap(
    tree,
    trait_values,
    model = phytools_model,
    nsim = requested_replicate_count,
    pi = "equal",
    message = FALSE
  )
  count_matrix <- phytools::countSimmap(fit, message = FALSE)
  total_transition_counts <- as.numeric(count_matrix[, "N"])
  simmap_result$summary$mean_total_transition_count <- as.numeric(mean(total_transition_counts))
  simmap_result$summary$lower_95_total_transition_count <- as.numeric(stats::quantile(total_transition_counts, probs = 0.025, names = FALSE))
  simmap_result$summary$upper_95_total_transition_count <- as.numeric(stats::quantile(total_transition_counts, probs = 0.975, names = FALSE))
  simmap_result$rows <- build_countsimmap_summary_rows(count_matrix)
  simmap_result
}

build_density_map_result <- function(tree, trait_values, trait_name, excluded_taxa, discrete_model) {
  simmap_result <- build_make_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    discrete_model
  )
  if (!identical(discrete_model, "equal-rates")) {
    stop("densityMap parity is only governed for equal-rates binary stochastic maps")
  }
  requested_replicate_count <- as.integer(case_payload$stochastic_map_replicate_count)
  stochastic_map_seed <- as.integer(case_payload$stochastic_map_seed)
  density_resolution <- as.integer(case_payload$density_resolution)
  set.seed(stochastic_map_seed)
  fit <- phytools::make.simmap(
    tree,
    trait_values,
    model = "ER",
    nsim = requested_replicate_count,
    pi = "equal",
    message = FALSE
  )
  capture.output({
    density_map <- phytools::densityMap(
      fit,
      plot = FALSE,
      res = density_resolution
    )
  })
  simmap_result$summary$branch_count <- nrow(density_map$tree$edge)
  simmap_result$summary$focal_state <- density_map$states[[2]]
  simmap_result$summary$baseline_state <- density_map$states[[1]]
  simmap_result$summary$resolution <- density_resolution
  simmap_result$summary$total_tree_depth <- as.numeric(max(nodeHeights(density_map$tree)))
  simmap_result$rows <- build_densitymap_branch_rows(density_map)
  simmap_result
}

build_rerooting_method_result <- function(tree, trait_values, trait_name, excluded_taxa, discrete_model) {
  phytools_model <- switch(
    discrete_model,
    "equal-rates" = "ER",
    "symmetric" = "SYM",
    stop(paste("unsupported rerootingMethod parity model:", discrete_model))
  )
  fit <- phytools::rerootingMethod(
    tree,
    trait_values,
    model = phytools_model
  )
  node_rows <- list()
  marginal <- fit$marginal.anc
  for (row_index in seq_len(nrow(marginal))) {
    node <- as.integer(rownames(marginal)[[row_index]])
    for (state in colnames(marginal)) {
      node_rows[[length(node_rows) + 1]] <- list(
        node = node_signature(tree, node),
        state = state,
        probability = unname(as.numeric(marginal[row_index, state]))
      )
    }
  }
  node_rows <- node_rows[order(
    vapply(node_rows, function(row) row$node, character(1)),
    vapply(node_rows, function(row) row$state, character(1))
  )]
  list(
    summary = list(
      taxon_count = length(trait_values),
      trait_name = trait_name,
      excluded_taxon_count = length(excluded_taxa),
      excluded_taxa = unname(as.list(excluded_taxa)),
      model = discrete_model,
      state_count = length(unique(unname(trait_values))),
      internal_node_count = nrow(marginal),
      root_prior_mode = "equal",
      phytools_rerooting_method_comparable = TRUE
    ),
    rows = node_rows
  )
}

build_fast_anc_result <- function(tree, trait_values, trait_name, excluded_taxa) {
  fit <- phytools::fastAnc(
    tree,
    trait_values,
    vars = TRUE,
    CI = TRUE
  )
  internal_nodes <- seq(
    length(tree$tip.label) + 1,
    length(tree$tip.label) + tree$Nnode
  )
  node_rows <- lapply(seq_along(internal_nodes), function(index) {
    node <- internal_nodes[[index]]
    list(
      node = node_signature(tree, node),
      estimate = unname(as.numeric(fit$ace[[index]])),
      standard_error = sqrt(unname(as.numeric(fit$var[[index]]))),
      lower_95_interval = unname(as.numeric(fit$CI95[[index, 1]])),
      upper_95_interval = unname(as.numeric(fit$CI95[[index, 2]]))
    )
  })
  node_rows <- node_rows[order(vapply(node_rows, function(row) row$node, character(1)))]
  list(
    summary = list(
      taxon_count = length(trait_values),
      trait_name = trait_name,
      internal_node_count = length(node_rows),
      excluded_taxon_count = length(excluded_taxa),
      excluded_taxa = unname(as.list(excluded_taxa)),
      tree_is_ultrametric = isTRUE(ape::is.ultrametric(tree))
    ),
    rows = node_rows
  )
}

build_anc_ml_result <- function(tree, trait_values, trait_name, excluded_taxa) {
  fit <- phytools::anc.ML(
    tree,
    trait_values,
    model = "BM",
    CI = TRUE,
    vars = TRUE
  )
  internal_nodes <- seq(
    length(tree$tip.label) + 1,
    length(tree$tip.label) + tree$Nnode
  )
  node_rows <- lapply(seq_along(internal_nodes), function(index) {
    node <- internal_nodes[[index]]
    list(
      node = node_signature(tree, node),
      estimate = unname(as.numeric(fit$ace[[index]])),
      standard_error = sqrt(unname(as.numeric(fit$var[[index]]))),
      lower_95_interval = unname(as.numeric(fit$CI95[[index, 1]])),
      upper_95_interval = unname(as.numeric(fit$CI95[[index, 2]]))
    )
  })
  node_rows <- node_rows[order(vapply(node_rows, function(row) row$node, character(1)))]
  list(
    summary = list(
      taxon_count = length(trait_values),
      trait_name = trait_name,
      internal_node_count = length(node_rows),
      excluded_taxon_count = length(excluded_taxa),
      excluded_taxa = unname(as.list(excluded_taxa)),
      tree_is_ultrametric = isTRUE(ape::is.ultrametric(tree)),
      sigma_squared = as.numeric(fit$sig2),
      log_likelihood = as.numeric(fit$logLik)
    ),
    rows = node_rows
  )
}

result_payload <- switch(
  case_payload$operation,
  "phylogenetic-signal-lambda" = list(
    summary = build_lambda_summary(
      tree,
      trait_values,
      trait_name
    ),
    rows = NULL
  ),
  "phylogenetic-signal-k" = list(
    summary = build_k_summary(
      tree,
      trait_values,
      trait_name
    ),
    rows = NULL
  ),
  "discrete-fit-mk" = build_fitmk_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
  ),
  "discrete-stochastic-map" = build_make_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
  ),
  "discrete-stochastic-map-count" = build_count_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
  ),
  "discrete-stochastic-map-description" = build_describe_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
  ),
  "discrete-stochastic-map-density" = build_density_map_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
  ),
  "simulate-discrete-history" = build_sim_history_result(
    tree,
    trait_name
  ),
  "simulate-continuous-brownian" = build_fastbm_result(
    tree
  ),
  "simulate-continuous-correlated-brownian" = build_simcorrs_result(
    tree
  ),
  "comparative-pgls-brownian" = build_pgls_result(
    tree,
    trait_table,
    taxon_column,
    case_payload$comparative_formula,
    case_payload$comparative_lambda_value
  ),
  "phylogenetic-residuals" = build_phyl_resid_result(
    tree,
    trait_table,
    taxon_column,
    trait_name,
    case_payload$comparative_predictors[[1]],
    case_payload$comparative_lambda_value
  ),
  "phylogenetic-anova" = build_phyl_anova_result(
    tree,
    trait_table,
    taxon_column,
    trait_name,
    case_payload$comparative_predictors[[1]],
    case_payload$permutation_count,
    case_payload$permutation_seed
  ),
  "discrete-ancestral-rerooting" = build_rerooting_method_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
  ),
  "continuous-ancestral-fast-anc" = build_fast_anc_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa
  ),
  "continuous-ancestral-anc-ml" = build_anc_ml_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa
  ),
  stop(paste("unsupported phytools parity operation:", case_payload$operation))
)
summary_payload <- result_payload$summary

write_payload(
  execution_path,
  list(
    status = "ok",
    r_version = r_version,
    phytools_version = as.character(utils::packageVersion("phytools"))
  )
)
write_payload(summary_path, summary_payload)
write_table(
  summary_table_path,
  lapply(names(summary_payload), function(key) {
    list(metric = key, value = format_table_value(summary_payload[[key]]))
  })
)
if (!is.null(result_payload$rows)) {
  if (identical(case_payload$operation, "discrete-fit-mk")) {
    write_table(fitmk_rows_path, result_payload$rows)
  } else if (case_payload$operation %in% c(
    "discrete-stochastic-map",
    "discrete-stochastic-map-count",
    "discrete-stochastic-map-description",
    "discrete-stochastic-map-density",
    "simulate-discrete-history"
  )) {
    write_table(stochastic_map_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "simulate-continuous-brownian")) {
    write_table(fastbm_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "simulate-continuous-correlated-brownian")) {
    write_table(simcorrs_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "comparative-pgls-brownian")) {
    write_table(pgls_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "phylogenetic-residuals")) {
    write_table(phyl_resid_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "phylogenetic-anova")) {
    write_table(phyl_anova_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "discrete-ancestral-rerooting")) {
    write_table(rerooting_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "continuous-ancestral-fast-anc")) {
    write_table(fast_anc_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "continuous-ancestral-anc-ml")) {
    write_table(anc_ml_rows_path, result_payload$rows)
  }
}
