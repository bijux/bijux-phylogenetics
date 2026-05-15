args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("expected case-json path and output-root path")
}

case_path <- args[[1]]
output_root <- args[[2]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required for ape parity execution")
}

write_payload <- function(path, payload) {
  writeLines(
    jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE),
    con = path
  )
}

execution_path <- file.path(output_root, "reference-execution.json")
case_payload <- jsonlite::fromJSON(case_path)
r_version <- as.character(getRversion())

if (!requireNamespace("ape", quietly = TRUE)) {
  write_payload(
    execution_path,
    list(
      status = "unavailable",
      mismatch_reason = "ape_package_unavailable",
      message = "ape is not installed in the active R environment",
      case_id = case_payload$case_id,
      function_name = case_payload$function_name,
      input_fixture = case_payload$input_fixture,
      r_version = r_version,
      ape_version = NULL
    )
  )
  quit(save = "no", status = 0)
}

tree_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(
    ape::read.tree(case_payload$input_fixture),
    error = function(error) error
  )
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        message = conditionMessage(tree),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  summary_path <- file.path(output_root, "summary.json")
  tips_path <- file.path(output_root, "tips.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")

  summary_payload <- list(
    tip_count = length(tree$tip.label),
    internal_node_count = tree$Nnode,
    edge_count = nrow(tree$edge),
    rooted = ape::is.rooted(tree),
    tip_labels = unname(tree$tip.label),
    branch_length_count = if (is.null(tree$edge.length)) 0 else sum(!is.na(tree$edge.length))
  )
  write_payload(summary_path, summary_payload)

  tip_table <- data.frame(
    position = seq_along(tree$tip.label),
    label = unname(tree$tip.label),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    tip_table,
    file = tips_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  ape::write.tree(tree, file = newick_path)

  write_payload(
    execution_path,
    list(
      status = "ok",
      case_id = case_payload$case_id,
      function_name = case_payload$function_name,
      input_fixture = case_payload$input_fixture,
      r_version = r_version,
      ape_version = as.character(utils::packageVersion("ape")),
      outputs = list(
        summary_json = summary_path,
        tip_table = tips_path,
        normalized_tree = newick_path
      )
    )
  )
}

if (identical(case_payload$operation, "read-tree-summary")) {
  tree_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

write_payload(
  execution_path,
  list(
    status = "failed",
    mismatch_reason = "unsupported_operation",
    message = paste("unsupported ape parity operation:", case_payload$operation),
    case_id = case_payload$case_id,
    function_name = case_payload$function_name,
    input_fixture = case_payload$input_fixture,
    r_version = r_version,
    ape_version = as.character(utils::packageVersion("ape"))
  )
)
