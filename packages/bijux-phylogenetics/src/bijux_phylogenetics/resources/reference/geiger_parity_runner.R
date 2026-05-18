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

normalize_optimizer_result <- function(fit) {
  result <- list()
  if (!is.null(fit$opt$method)) {
    result$best_method <- as.character(fit$opt$method)
  }
  if (is.matrix(fit$res) || is.data.frame(fit$res)) {
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

build_fitcontinuous_payload <- function(tree, trait_values, excluded_taxa, case_payload) {
  fit <- do.call(
    geiger::fitContinuous,
    list(
      phy = tree,
      dat = trait_values,
      model = case_payload$model_name,
      bounds = fitcontinuous_bounds(case_payload)
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
    missing_value_policy = "prune-tree-tip-overlap-with-missing-or-nonnumeric-trait-values",
    standard_error_policy = "tip-standard-errors-not-supported",
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
    raw_trait_values <- suppressWarnings(as.numeric(trait_table[[case_payload$trait_name]]))
    trait_values <- stats::setNames(raw_trait_values, taxon_names)
    tree_only_taxa <- sort(setdiff(original_tree_tip_labels, taxon_names))
    extra_trait_taxa <- sort(setdiff(taxon_names, original_tree_tip_labels))
    overlapping_taxa <- intersect(original_tree_tip_labels, taxon_names)
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
    payload <- build_fitcontinuous_payload(tree, trait_values, excluded_taxa, case_payload)
    payload$summary$missing_from_traits <- json_array(tree_only_taxa)
    payload$summary$missing_value_taxa <- json_array(missing_value_taxa)
    payload$summary$non_numeric_taxa <- json_array(non_numeric_taxa)
    payload$summary$extra_trait_taxa <- json_array(extra_trait_taxa)
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
