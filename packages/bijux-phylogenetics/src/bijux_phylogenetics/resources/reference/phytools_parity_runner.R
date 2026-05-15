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
  if (length(value) > 1) {
    return(jsonlite::toJSON(unname(value), auto_unbox = TRUE))
  }
  value
}

execution_path <- file.path(output_root, "reference-execution.json")
summary_path <- file.path(output_root, "reference-summary.json")
summary_table_path <- file.path(output_root, "reference-summary.tsv")
fast_anc_rows_path <- file.path(output_root, "fast-anc-node-estimates.tsv")
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
trait_values <- trait_table[[trait_name]]
taxon_names <- trait_table[[taxon_column]]
names(trait_values) <- taxon_names
trait_values <- stats::setNames(as.numeric(trait_values), taxon_names)
kept_taxa <- names(trait_values)[
  !is.na(trait_values) &
    names(trait_values) %in% tree$tip.label
]
excluded_taxa <- sort(setdiff(tree$tip.label, kept_taxa))
trait_values <- trait_values[kept_taxa]
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
      excluded_taxa = excluded_taxa,
      tree_is_ultrametric = isTRUE(ape::is.ultrametric(tree))
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
  "continuous-ancestral-fast-anc" = build_fast_anc_result(
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
  write_table(fast_anc_rows_path, result_payload$rows)
}
