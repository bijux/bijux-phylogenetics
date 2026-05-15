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

write_table <- function(path, rows) {
  utils::write.table(
    as.data.frame(rows, stringsAsFactors = FALSE),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

execution_path <- file.path(output_root, "reference-execution.json")
summary_path <- file.path(output_root, "reference-summary.json")
summary_table_path <- file.path(output_root, "reference-summary.tsv")
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
  fit <- phytools::phylosig(
    tree,
    trait_values,
    method = "K",
    test = TRUE
  )
  simulated_values <- as.numeric(fit$sim.K)
  list(
    taxon_count = length(trait_values),
    trait_name = trait_name,
    k = as.numeric(fit$K),
    p_value = as.numeric(fit$P),
    permutation_count = length(simulated_values),
    simulated_k_minimum = as.numeric(min(simulated_values)),
    simulated_k_mean = as.numeric(mean(simulated_values)),
    simulated_k_maximum = as.numeric(max(simulated_values))
  )
}

summary_payload <- switch(
  case_payload$operation,
  "phylogenetic-signal-lambda" = build_lambda_summary(
    tree,
    trait_values,
    trait_name
  ),
  "phylogenetic-signal-k" = build_k_summary(
    tree,
    trait_values,
    trait_name
  ),
  stop(paste("unsupported phytools parity operation:", case_payload$operation))
)

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
    list(metric = key, value = summary_payload[[key]])
  })
)
