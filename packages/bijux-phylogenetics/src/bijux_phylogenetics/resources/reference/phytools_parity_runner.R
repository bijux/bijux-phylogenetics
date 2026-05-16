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
traits_path <- input_fixtures[[2]]
trait_name <- case_payload$trait_name
taxon_column <- case_payload$taxon_column
trait_table <- utils::read.table(
  traits_path,
  header = TRUE,
  sep = if (grepl("\\.csv$", traits_path, ignore.case = TRUE)) "," else "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
taxon_names <- trait_table[[taxon_column]]
is_discrete_operation <- case_payload$operation %in% c(
  "discrete-fit-mk",
  "discrete-stochastic-map",
  "discrete-stochastic-map-description",
  "discrete-ancestral-rerooting"
)
if (is_discrete_operation) {
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
} else {
  trait_values <- stats::setNames(as.numeric(trait_table[[trait_name]]), taxon_names)
  kept_taxa <- names(trait_values)[
    !is.na(trait_values) &
      names(trait_values) %in% tree$tip.label
  ]
  trait_values <- trait_values[kept_taxa]
}
excluded_taxa <- sort(setdiff(tree$tip.label, kept_taxa))
tree <- ape::drop.tip(tree, excluded_taxa)

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
  "discrete-stochastic-map-description" = build_describe_simmap_result(
    tree,
    trait_values,
    trait_name,
    excluded_taxa,
    case_payload$discrete_model
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
    "discrete-stochastic-map-description"
  )) {
    write_table(stochastic_map_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "discrete-ancestral-rerooting")) {
    write_table(rerooting_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "continuous-ancestral-fast-anc")) {
    write_table(fast_anc_rows_path, result_payload$rows)
  } else if (identical(case_payload$operation, "continuous-ancestral-anc-ml")) {
    write_table(anc_ml_rows_path, result_payload$rows)
  }
}
