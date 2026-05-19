args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("expected case-json path and output-root path")
}

case_path <- args[[1]]
output_root <- args[[2]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required for geiger parity execution")
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

write_table <- function(path, rows) {
  if (length(rows) == 0) {
    utils::write.table(
      data.frame(parameter = character(), value = numeric()),
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    return(invisible(NULL))
  }
  normalized_rows <- lapply(rows, function(row) {
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  utils::write.table(
    do.call(rbind.data.frame, c(normalized_rows, list(stringsAsFactors = FALSE))),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

execution_path <- file.path(output_root, "reference-execution.json")
summary_path <- file.path(output_root, "reference-summary.json")
rows_path <- file.path(output_root, "reference-parameters.tsv")
case_payload <- jsonlite::fromJSON(case_path)
r_version <- as.character(getRversion())

if (!requireNamespace("geiger", quietly = TRUE)) {
  write_payload(
    execution_path,
    list(
      status = "unavailable",
      mismatch_reason = "geiger_package_unavailable",
      r_version = r_version,
      geiger_version = NULL
    )
  )
  quit(save = "no", status = 0)
}

library(ape)
library(geiger)

fit_aicc <- function(aic, sample_size, parameter_count) {
  denominator <- sample_size - parameter_count - 1
  if (denominator <= 0) {
    return(NULL)
  }
  aic + ((2 * parameter_count * (parameter_count + 1)) / denominator)
}

json_array <- function(values) {
  if (length(values) == 0) {
    return(list())
  }
  as.list(unname(values))
}

public_parameter_bounds <- function(case_payload) {
  if (identical(case_payload$model_name, "lambda") && !is.null(case_payload$lambda_bounds)) {
    return(as.numeric(unlist(case_payload$lambda_bounds)))
  }
  if (identical(case_payload$model_name, "kappa") && !is.null(case_payload$kappa_bounds)) {
    return(as.numeric(unlist(case_payload$kappa_bounds)))
  }
  if (identical(case_payload$model_name, "delta") && !is.null(case_payload$delta_bounds)) {
    return(as.numeric(unlist(case_payload$delta_bounds)))
  }
  if (identical(case_payload$model_name, "OU") && !is.null(case_payload$ou_bounds)) {
    return(as.numeric(unlist(case_payload$ou_bounds)))
  }
  if (identical(case_payload$model_name, "EB") && !is.null(case_payload$early_burst_bounds)) {
    return(as.numeric(unlist(case_payload$early_burst_bounds)))
  }
  NULL
}

hit_parameter_boundary <- function(value, bounds) {
  if (is.null(value) || is.null(bounds) || length(bounds) != 2) {
    return(list(hit_lower = FALSE, hit_upper = FALSE))
  }
  tolerance <- max(1e-9, abs(bounds[[2]] - bounds[[1]]) * 1e-6)
  list(
    hit_lower = isTRUE(all.equal(as.numeric(value), bounds[[1]], tolerance = tolerance)),
    hit_upper = isTRUE(all.equal(as.numeric(value), bounds[[2]], tolerance = tolerance))
  )
}

fitcontinuous_bounds <- function(case_payload) {
  if (identical(case_payload$model_name, "lambda") && !is.null(case_payload$lambda_bounds)) {
    bounds <- as.numeric(unlist(case_payload$lambda_bounds))
    return(list(lambda = bounds))
  }
  if (identical(case_payload$model_name, "kappa") && !is.null(case_payload$kappa_bounds)) {
    bounds <- as.numeric(unlist(case_payload$kappa_bounds))
    return(list(kappa = bounds))
  }
  if (identical(case_payload$model_name, "delta") && !is.null(case_payload$delta_bounds)) {
    bounds <- as.numeric(unlist(case_payload$delta_bounds))
    return(list(delta = bounds))
  }
  if (identical(case_payload$model_name, "OU") && !is.null(case_payload$ou_bounds)) {
    bounds <- as.numeric(unlist(case_payload$ou_bounds))
    return(list(alpha = bounds))
  }
  if (identical(case_payload$model_name, "EB") && !is.null(case_payload$early_burst_bounds)) {
    public_bounds <- sort(as.numeric(unlist(case_payload$early_burst_bounds)))
    lower_public <- public_bounds[[1]]
    upper_public <- public_bounds[[2]]
    upper_a <- if (lower_public <= 0) -1e-06 else -lower_public
    return(list(a = c(-upper_public, upper_a)))
  }
  list()
}

fitcontinuous_control <- function(case_payload) {
  if (is.null(case_payload$reference_control)) {
    return(list())
  }
  case_payload$reference_control
}

matrix_attempt_rows <- function(result_matrix) {
  if (!(is.matrix(result_matrix) || is.data.frame(result_matrix))) {
    return(NULL)
  }
  row_names <- rownames(result_matrix)
  rows <- lapply(seq_len(nrow(result_matrix)), function(index) {
    attempt <- as.list(result_matrix[index, , drop = FALSE])
    for (column_name in names(attempt)) {
      value <- attempt[[column_name]]
      if (is.factor(value)) {
        value <- as.character(value)
      }
      if (is.numeric(value) || is.integer(value)) {
        attempt[[column_name]] <- as.numeric(value)
      } else if (is.logical(value)) {
        attempt[[column_name]] <- isTRUE(value)
      } else {
        attempt[[column_name]] <- value
      }
    }
    attempt$attempt_index <- index
    if (!is.null(row_names) && length(row_names) >= index) {
      attempt$method <- as.character(row_names[[index]])
    }
    if (!is.null(attempt$lnL)) {
      attempt$log_likelihood <- as.numeric(attempt$lnL)
    }
    attempt
  })
  json_array(rows)
}

normalize_optimizer_result <- function(fit) {
  result <- list()
  if (!is.null(fit$opt$method)) {
    result$best_method <- as.character(fit$opt$method)
  }
  if (is.matrix(fit$res) || is.data.frame(fit$res)) {
    result$attempt_rows <- matrix_attempt_rows(fit$res)
    result$attempt_count <- nrow(fit$res)
    if (!is.null(rownames(fit$res))) {
      result$attempted_methods <- sort(unique(as.character(rownames(fit$res))))
    }
    if ("convergence" %in% colnames(fit$res)) {
      convergence_codes <- as.integer(fit$res[, "convergence"])
      result$converged_attempt_count <- sum(convergence_codes == 0, na.rm = TRUE)
      best_method <- result$best_method
      if (!is.null(best_method) && best_method %in% rownames(fit$res)) {
        result$convergence_code <- as.integer(fit$res[best_method[[1]], "convergence"])
      } else if (length(convergence_codes) > 0) {
        result$convergence_code <- convergence_codes[[1]]
      }
    }
    if ("lnL" %in% colnames(fit$res)) {
      result$best_log_likelihood <- max(as.numeric(fit$res[, "lnL"]), na.rm = TRUE)
    }
  } else if (is.list(fit$res)) {
    if (!is.null(fit$res$convergence)) {
      result$convergence_code <- as.integer(fit$res$convergence)
    }
    if (!is.null(fit$res$message)) {
      result$message <- as.character(fit$res$message)
    }
    if (!is.null(fit$res$counts)) {
      counts <- unclass(fit$res$counts)
      result$counts <- as.list(counts)
    }
    if (!is.null(fit$res$value)) {
      result$objective_value <- as.numeric(fit$res$value)
    }
  }
  if (length(result) == 0) {
    return(NULL)
  }
  result
}

normalize_fitdiscrete_result <- function(fit) {
  result <- list()
  if (!is.null(fit$opt$method)) {
    result$best_method <- as.character(fit$opt$method)
  }
  if (!is.null(fit$opt$k)) {
    result$parameter_count <- as.integer(fit$opt$k)
  }
  if (is.matrix(fit$res) || is.data.frame(fit$res)) {
    result$attempt_rows <- matrix_attempt_rows(fit$res)
    result$attempt_count <- nrow(fit$res)
    convergence_codes <- as.integer(fit$res[, ncol(fit$res)])
    result$converged_attempt_count <- sum(convergence_codes == 0, na.rm = TRUE)
    if (length(convergence_codes) > 0) {
      result$convergence_code <- convergence_codes[[1]]
    }
    result$best_log_likelihood <- max(as.numeric(fit$res[, 2]), na.rm = TRUE)
  }
  if (length(result) == 0) {
    return(NULL)
  }
  result
}

decode_transition_parameter <- function(parameter_name, state_count) {
  code <- substring(parameter_name, 2)
  split_candidates <- seq.int(1, nchar(code) - 1)
  for (split_index in split_candidates) {
    left_index <- as.integer(substring(code, 1, split_index))
    right_index <- as.integer(substring(code, split_index + 1, nchar(code)))
    if (!is.na(left_index) &&
        !is.na(right_index) &&
        left_index >= 1 &&
        left_index <= state_count &&
        right_index >= 1 &&
        right_index <= state_count) {
      return(c(left_index, right_index))
    }
  }
  stop(sprintf("could not decode fitDiscrete transition parameter '%s'", parameter_name))
}

fitdiscrete_rate_rows <- function(opt, state_levels) {
  parameter_names <- grep("^q[0-9]+[0-9]+$", names(opt), value = TRUE)
  rows <- lapply(parameter_names, function(parameter_name) {
    decoded <- decode_transition_parameter(parameter_name, length(state_levels))
    list(
      source_state = state_levels[[decoded[[1]]]],
      target_state = state_levels[[decoded[[2]]]],
      transition_allowed = TRUE,
      step_distance = 1,
      rate = as.numeric(opt[[parameter_name]])
    )
  })
  rows[order(
    vapply(rows, function(row) row$source_state, character(1)),
    vapply(rows, function(row) row$target_state, character(1))
  )]
}

parameter_surface <- function(model_name, opt) {
  if (identical(model_name, "OU")) {
    return(list(parameter_name = "alpha", parameter_value = as.numeric(opt$alpha)))
  }
  if (identical(model_name, "EB")) {
    return(list(parameter_name = "rate_change", parameter_value = -as.numeric(opt$a)))
  }
  if (identical(model_name, "lambda")) {
    return(list(parameter_name = "lambda", parameter_value = as.numeric(opt$lambda)))
  }
  if (identical(model_name, "kappa")) {
    return(list(parameter_name = "kappa", parameter_value = as.numeric(opt$kappa)))
  }
  if (identical(model_name, "delta")) {
    return(list(parameter_name = "delta", parameter_value = as.numeric(opt$delta)))
  }
  list(parameter_name = NULL, parameter_value = NULL)
}

public_mode_name <- function(model_name) {
  if (identical(model_name, "BM")) {
    return("brownian")
  }
  if (identical(model_name, "white")) {
    return("white-noise")
  }
  if (identical(model_name, "lambda")) {
    return("pagel-lambda")
  }
  if (identical(model_name, "kappa")) {
    return("pagel-kappa")
  }
  if (identical(model_name, "delta")) {
    return("pagel-delta")
  }
  if (identical(model_name, "OU")) {
    return("ornstein-uhlenbeck")
  }
  if (identical(model_name, "EB")) {
    return("early-burst")
  }
  stop(sprintf("unsupported fitContinuous model name: %s", model_name))
}

fitcontinuous_parameter_count <- function(model_name) {
  if (identical(model_name, "BM") || identical(model_name, "white")) {
    return(2)
  }
  3
}

build_fitcontinuous_payload <- function(tree, trait_values, excluded_taxa, case_payload) {
  fit <- do.call(
    geiger::fitContinuous,
    list(
      phy = tree,
      dat = trait_values,
      model = case_payload$model_name,
      bounds = fitcontinuous_bounds(case_payload),
      control = fitcontinuous_control(case_payload)
    )
  )
  opt <- fit$opt
  parameter_surface_result <- parameter_surface(case_payload$model_name, opt)
  bounds <- public_parameter_bounds(case_payload)
  boundary_hits <- hit_parameter_boundary(
    parameter_surface_result$parameter_value,
    bounds
  )
  parameter_count <- if (is.null(parameter_surface_result$parameter_name)) 2 else 3
  log_likelihood <- if (is.null(opt$lnL)) NULL else as.numeric(opt$lnL)
  aic <- if (!is.null(fit$aic)) {
    as.numeric(fit$aic)
  } else if (!is.null(opt$aic)) {
    as.numeric(opt$aic)
  } else if (!is.null(log_likelihood)) {
    as.numeric((2 * parameter_count) - (2 * log_likelihood))
  } else {
    NULL
  }
  aicc <- if (!is.null(opt$aicc)) {
    as.numeric(opt$aicc)
  } else if (!is.null(aic)) {
    fit_aicc(aic, length(trait_values), parameter_count)
  } else {
    NULL
  }
  summary <- list(
    taxon_count = length(trait_values),
    trait_name = case_payload$trait_name,
    model_name = case_payload$model_name,
    excluded_taxon_count = length(excluded_taxa),
    excluded_taxa = json_array(excluded_taxa),
    root_state = if (is.null(opt$z0)) NULL else as.numeric(opt$z0),
    rate = if (is.null(opt$sigsq)) NULL else as.numeric(opt$sigsq),
    log_likelihood = log_likelihood,
    aic = aic,
    aicc = aicc,
    likelihood_constant_policy = "full-gaussian-loglikelihood-includes-normalizing-constant",
    likelihood_comparison_policy = "raw-loglikelihood-and-derived-aic-are-directly-comparable-when-the-shared-gaussian-constant-policy-matches",
    missing_value_policy = "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values",
    standard_error_policy = "fitcontinuous-standard-error-explicitly-excluded-this-round",
    parameter_bound_policy = if (is.null(bounds)) {
      "reference-default-without-explicit-bounds"
    } else {
      "governed-bounded-grid-search"
    },
    hit_lower_parameter_boundary = boundary_hits$hit_lower,
    hit_upper_parameter_boundary = boundary_hits$hit_upper,
    parameter_name = parameter_surface_result$parameter_name,
    parameter_value = parameter_surface_result$parameter_value,
    optimizer_settings = case_payload$optimizer_settings,
    optimizer_result = normalize_optimizer_result(fit)
  )
  rows <- Filter(
    f = function(row) !is.null(row$value),
    x = list(
      list(parameter = "root_state", value = summary$root_state),
      list(parameter = "rate", value = summary$rate),
      list(parameter = "log_likelihood", value = summary$log_likelihood),
      list(parameter = "aic", value = summary$aic),
      list(parameter = "aicc", value = summary$aicc),
      if (!is.null(summary$parameter_name)) {
        list(parameter = summary$parameter_name, value = summary$parameter_value)
      }
    )
  )
  list(summary = summary, rows = rows)
}

build_fitcontinuous_model_comparison_payload <- function(tree, trait_values, excluded_taxa, case_payload) {
  candidate_model_names <- unname(unlist(case_payload$candidate_model_names))
  fitted_rows <- lapply(candidate_model_names, function(model_name) {
    model_case_payload <- case_payload
    model_case_payload$model_name <- model_name
    fit <- do.call(
      geiger::fitContinuous,
      list(
        phy = tree,
        dat = trait_values,
        model = model_name,
        bounds = fitcontinuous_bounds(model_case_payload),
        control = fitcontinuous_control(model_case_payload)
      )
    )
    opt <- fit$opt
    parameter_count <- fitcontinuous_parameter_count(model_name)
    log_likelihood <- if (is.null(opt$lnL)) NULL else as.numeric(opt$lnL)
    aic <- if (!is.null(fit$aic)) {
      as.numeric(fit$aic)
    } else if (!is.null(opt$aic)) {
      as.numeric(opt$aic)
    } else if (!is.null(log_likelihood)) {
      as.numeric((2 * parameter_count) - (2 * log_likelihood))
    } else {
      NULL
    }
    aicc <- if (!is.null(opt$aicc)) {
      as.numeric(opt$aicc)
    } else if (!is.null(aic)) {
      fit_aicc(aic, length(trait_values), parameter_count)
    } else {
      NULL
    }
    list(
      model = public_mode_name(model_name),
      rank = NULL,
      parameter_count = parameter_count,
      log_likelihood = log_likelihood,
      aic = aic,
      aicc = aicc,
      delta_aic = NULL,
      delta_aicc = NULL,
      selected = FALSE,
      comparable = !is.null(aic) && !is.null(aicc),
      likelihood_constant_policy = "full-gaussian-loglikelihood-includes-normalizing-constant",
      comparability_note = if (is.null(aicc)) {
        "sample size is too small to compute finite AICc for this parameter count"
      } else {
        ""
      }
    )
  })
  comparable_rows <- Filter(
    f = function(row) isTRUE(row$comparable),
    x = fitted_rows
  )
  if (length(comparable_rows) == 0) {
    stop("no comparable fitContinuous model retained a finite AICc surface")
  }
  best_aic <- min(vapply(comparable_rows, function(row) row$aic, numeric(1)))
  best_aicc <- min(vapply(comparable_rows, function(row) row$aicc, numeric(1)))
  ranked_rows <- comparable_rows[order(
    vapply(comparable_rows, function(row) row$aicc, numeric(1)),
    vapply(comparable_rows, function(row) row$aic, numeric(1)),
    vapply(comparable_rows, function(row) row$model, character(1))
  )]
  for (index in seq_along(ranked_rows)) {
    ranked_rows[[index]]$rank <- index
    ranked_rows[[index]]$delta_aic <- ranked_rows[[index]]$aic - best_aic
    ranked_rows[[index]]$delta_aicc <- ranked_rows[[index]]$aicc - best_aicc
    ranked_rows[[index]]$selected <- isTRUE(all.equal(
      ranked_rows[[index]]$aicc,
      best_aicc,
      tolerance = 1e-12
    ))
  }
  fitted_rows <- c(
    ranked_rows,
    Filter(
      f = function(row) !isTRUE(row$comparable),
      x = fitted_rows
    )
  )
  selected_row <- ranked_rows[[1]]
  runner_up_row <- if (length(ranked_rows) >= 2) ranked_rows[[2]] else NULL
  summary <- list(
    taxon_count = length(trait_values),
    trait_name = case_payload$trait_name,
    model_name = case_payload$model_name,
    selected_model = selected_row$model,
    model_ranking = json_array(vapply(fitted_rows, function(row) row$model, character(1))),
    comparable_model_count = length(ranked_rows),
    noncomparable_model_count = length(fitted_rows) - length(ranked_rows),
    runner_up_model = if (is.null(runner_up_row)) NULL else runner_up_row$model,
    runner_up_aicc_delta = if (is.null(runner_up_row)) NULL else runner_up_row$delta_aicc,
    warning_count = if (length(candidate_model_names) > 3) 1 else 0,
    likelihood_constant_policy = "full-gaussian-loglikelihood-includes-normalizing-constant",
    likelihood_comparison_policy = "relative-aic-and-aicc-ranking-is-permitted-only-when-all-candidate-modes-share-one-gaussian-likelihood-constant-policy",
    noncomparable_likelihood_models = json_array(character(0)),
    optimizer_settings = case_payload$optimizer_settings
  )
  list(summary = summary, rows = fitted_rows)
}

build_fitdiscrete_payload <- function(tree, trait_values, excluded_taxa, missing_value_taxa, missing_from_traits, extra_trait_taxa, case_payload) {
  bounds <- list()
  parameter_bounds <- NULL
  if (identical(case_payload$discrete_transform_name, "EB") && !is.null(case_payload$early_burst_bounds)) {
    parameter_bounds <- as.numeric(unlist(case_payload$early_burst_bounds))
    bounds <- list(a = parameter_bounds)
  }
  if (identical(case_payload$discrete_transform_name, "lambda") && !is.null(case_payload$lambda_bounds)) {
    parameter_bounds <- as.numeric(unlist(case_payload$lambda_bounds))
  }
  if (identical(case_payload$discrete_transform_name, "kappa") && !is.null(case_payload$kappa_bounds)) {
    parameter_bounds <- as.numeric(unlist(case_payload$kappa_bounds))
  }
  if (identical(case_payload$discrete_transform_name, "delta") && !is.null(case_payload$delta_bounds)) {
    parameter_bounds <- as.numeric(unlist(case_payload$delta_bounds))
  }
  fit <- geiger::fitDiscrete(
    phy = tree,
    dat = trait_values,
    model = case_payload$model_name,
    transform = if (is.null(case_payload$discrete_transform_name)) NULL else case_payload$discrete_transform_name,
    bounds = bounds
  )
  state_levels <- levels(trait_values)
  observed_state_count <- length(state_levels)
  state_counts <- vapply(
    state_levels,
    function(state) sum(as.character(trait_values) == state),
    integer(1)
  )
  sparse_states <- names(state_counts[state_counts < 2])
  transform_name <- if (is.null(case_payload$discrete_transform_name)) {
    NULL
  } else if (identical(case_payload$discrete_transform_name, "EB")) {
    "early-burst"
  } else {
    paste0("pagel-", case_payload$discrete_transform_name)
  }
  parameter_name <- NULL
  parameter_value <- NULL
  if (identical(case_payload$discrete_transform_name, "lambda")) {
    parameter_name <- "lambda"
    parameter_value <- as.numeric(fit$opt$lambda)
  } else if (identical(case_payload$discrete_transform_name, "kappa")) {
    parameter_name <- "kappa"
    parameter_value <- as.numeric(fit$opt$kappa)
  } else if (identical(case_payload$discrete_transform_name, "delta")) {
    parameter_name <- "delta"
    parameter_value <- as.numeric(fit$opt$delta)
  } else if (identical(case_payload$discrete_transform_name, "EB")) {
    parameter_name <- "a"
    parameter_value <- as.numeric(fit$opt$a)
  }
  boundary_hits <- hit_parameter_boundary(parameter_value, parameter_bounds)
  summary <- list(
    taxon_count = length(trait_values),
    trait_name = case_payload$trait_name,
    model_name = case_payload$model_name,
    transform_name = transform_name,
    observed_state_count = observed_state_count,
    state_order = json_array(state_levels),
    excluded_taxon_count = length(excluded_taxa),
    excluded_taxa = json_array(excluded_taxa),
    missing_value_taxa = json_array(missing_value_taxa),
    missing_from_traits = json_array(missing_from_traits),
    extra_trait_taxa = json_array(extra_trait_taxa),
    missing_value_policy = "prune-overlapping-missing-values",
    log_likelihood = as.numeric(fit$opt$lnL),
    parameter_count = as.integer(fit$opt$k),
    aic = as.numeric(fit$opt$aic),
    aicc = as.numeric(fit$opt$aicc),
    parameter_name = parameter_name,
    parameter_value = parameter_value,
    hit_lower_parameter_boundary = boundary_hits$hit_lower,
    hit_upper_parameter_boundary = boundary_hits$hit_upper,
    sparse_states = json_array(unname(sparse_states)),
    optimizer_settings = case_payload$optimizer_settings,
    optimizer_result = normalize_fitdiscrete_result(fit)
  )
  rows <- fitdiscrete_rate_rows(fit$opt, state_levels)
  list(summary = summary, rows = rows)
}

result <- tryCatch(
  {
    input_fixtures <- unname(unlist(case_payload$input_fixtures))
    tree <- ape::read.tree(input_fixtures[[1]])
    original_tree_tip_labels <- tree$tip.label
    traits_path <- input_fixtures[[2]]
    trait_table <- utils::read.table(
      traits_path,
      header = TRUE,
      sep = if (grepl("\\.csv$", traits_path, ignore.case = TRUE)) "," else "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    taxon_names <- trait_table[[case_payload$taxon_column]]
    raw_trait_column <- setNames(as.character(trait_table[[case_payload$trait_name]]), taxon_names)
    tree_only_taxa <- sort(setdiff(original_tree_tip_labels, taxon_names))
    extra_trait_taxa <- sort(setdiff(taxon_names, original_tree_tip_labels))
    overlapping_taxa <- intersect(original_tree_tip_labels, taxon_names)
    if (identical(case_payload$operation, "fit-discrete-mk")) {
      trimmed_trait_values <- trimws(raw_trait_column[overlapping_taxa])
      missing_value_taxa <- sort(
        overlapping_taxa[is.na(raw_trait_column[overlapping_taxa]) | !nzchar(trimmed_trait_values)]
      )
      kept_taxa <- overlapping_taxa[!(overlapping_taxa %in% missing_value_taxa)]
      trait_values <- stats::setNames(trimmed_trait_values[kept_taxa], kept_taxa)
      excluded_taxa <- sort(unique(c(tree_only_taxa, missing_value_taxa)))
      if (length(excluded_taxa) > 0) {
        tree <- ape::drop.tip(tree, excluded_taxa)
      }
      trait_values <- factor(trait_values)
      payload <- build_fitdiscrete_payload(
        tree,
        trait_values,
        excluded_taxa,
        missing_value_taxa,
        tree_only_taxa,
        extra_trait_taxa,
        case_payload
      )
    } else {
      raw_trait_values <- suppressWarnings(as.numeric(trait_table[[case_payload$trait_name]]))
      trait_values <- stats::setNames(raw_trait_values, taxon_names)
      empty_or_missing_values <- is.na(raw_trait_column[overlapping_taxa]) |
        !nzchar(trimws(raw_trait_column[overlapping_taxa]))
      missing_value_taxa <- sort(
        overlapping_taxa[is.na(trait_values[overlapping_taxa]) & empty_or_missing_values]
      )
      non_numeric_taxa <- sort(
        overlapping_taxa[is.na(trait_values[overlapping_taxa]) & !empty_or_missing_values]
      )
      kept_taxa <- names(trait_values)[!is.na(trait_values) & names(trait_values) %in% original_tree_tip_labels]
      trait_values <- trait_values[kept_taxa]
      excluded_taxa <- sort(unique(c(tree_only_taxa, missing_value_taxa, non_numeric_taxa)))
      if (length(excluded_taxa) > 0) {
        tree <- ape::drop.tip(tree, excluded_taxa)
      }
      payload <- if (identical(case_payload$operation, "compare-fitcontinuous-models")) {
        build_fitcontinuous_model_comparison_payload(tree, trait_values, excluded_taxa, case_payload)
      } else {
        build_fitcontinuous_payload(tree, trait_values, excluded_taxa, case_payload)
      }
      payload$summary$missing_from_traits <- json_array(tree_only_taxa)
      payload$summary$missing_value_taxa <- json_array(missing_value_taxa)
      payload$summary$non_numeric_taxa <- json_array(non_numeric_taxa)
      payload$summary$extra_trait_taxa <- json_array(extra_trait_taxa)
    }
    write_payload(summary_path, payload$summary)
    write_table(rows_path, payload$rows)
    write_payload(
      execution_path,
      list(
        status = "ok",
        r_version = r_version,
        geiger_version = as.character(packageVersion("geiger"))
      )
    )
    TRUE
  },
  error = function(error) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = class(error)[[1]],
        message = conditionMessage(error),
        r_version = r_version,
        geiger_version = as.character(packageVersion("geiger"))
      )
    )
    FALSE
  }
)

if (!isTRUE(result)) {
  quit(save = "no", status = 0)
}
