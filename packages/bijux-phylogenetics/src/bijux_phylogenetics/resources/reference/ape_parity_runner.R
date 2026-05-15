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

write_table <- function(path, table) {
  utils::write.table(
    table,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

execution_path <- file.path(output_root, "reference-execution.json")
case_payload <- jsonlite::fromJSON(case_path)
r_version <- as.character(getRversion())

normalize_node_label <- function(label) {
  if (is.null(label) || is.na(label)) {
    return("")
  }
  as.character(label)
}

parse_support_label <- function(label) {
  text <- normalize_node_label(label)
  if (grepl("^'.*'$", text)) {
    text <- gsub("^'|'$", "", text)
    text <- gsub("''", "'", text, fixed = TRUE)
  }
  if (!nzchar(text)) {
    return("")
  }
  if (grepl("^[0-9.]+/[0-9.]+$", text)) {
    return(as.numeric(sub(".*/", "", text)))
  }
  numeric_value <- suppressWarnings(as.numeric(text))
  if (!is.na(numeric_value)) {
    return(numeric_value)
  }
  ""
}

descendant_taxa <- function(tree, node) {
  if (node <= length(tree$tip.label)) {
    return(normalize_node_label(tree$tip.label[[node]]))
  }
  children <- tree$edge[tree$edge[, 1] == node, 2]
  sort(unique(unlist(lapply(children, function(child) descendant_taxa(tree, child)))))
}

root_node <- function(tree) {
  internal_nodes <- sort(unique(tree$edge[, 1]))
  root_candidates <- setdiff(internal_nodes, tree$edge[, 2])
  if (length(root_candidates) == 0) {
    return(internal_nodes[[1]])
  }
  root_candidates[[1]]
}

node_kind <- function(tree, node, root_id) {
  if (node == root_id) {
    return("root")
  }
  if (node <= length(tree$tip.label)) {
    return("tip")
  }
  "internal"
}

node_label <- function(tree, node) {
  if (node <= length(tree$tip.label)) {
    return(normalize_node_label(tree$tip.label[[node]]))
  }
  node_labels <- tree$node.label
  index <- node - length(tree$tip.label)
  if (is.null(node_labels) || index < 1 || index > length(node_labels)) {
    return("")
  }
  normalize_node_label(node_labels[[index]])
}

node_branch_length <- function(tree, node) {
  match_index <- match(node, tree$edge[, 2])
  if (is.na(match_index) || is.null(tree$edge.length)) {
    return("")
  }
  value <- tree$edge.length[[match_index]]
  if (is.na(value)) {
    return("")
  }
  as.numeric(value)
}

tree_structure_rows <- function(tree, tree_index_value) {
  root_id <- root_node(tree)
  internal_nodes <- sort(unique(tree$edge[, 1]))
  nodes <- c(
    root_id,
    setdiff(internal_nodes, root_id),
    seq_len(length(tree$tip.label))
  )
  rows <- lapply(nodes, function(node) {
    taxa <- descendant_taxa(tree, node)
    label <- node_label(tree, node)
    list(
      tree_index = tree_index_value,
      node_kind = node_kind(tree, node, root_id),
      clade_id = paste(taxa, collapse = "|"),
      node_label = label,
      taxon_count = length(taxa),
      taxa = paste(taxa, collapse = "|"),
      support = parse_support_label(label),
      branch_length = node_branch_length(tree, node)
    )
  })
  order_key <- function(row) {
    node_rank <- switch(
      row$node_kind,
      root = 0L,
      internal = 1L,
      tip = 2L,
      9L
    )
    list(
      if (identical(tree_index_value, "")) 0L else as.integer(tree_index_value),
      node_rank,
      row$clade_id,
      row$node_label
    )
  }
  order_matrix <- do.call(
    rbind,
    lapply(rows, function(row) unlist(order_key(row), use.names = FALSE))
  )
  rows[do.call(order, as.data.frame(order_matrix, stringsAsFactors = FALSE))]
}

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
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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
  clades_path <- file.path(output_root, "clades.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  clade_rows <- tree_structure_rows(tree, "")

  summary_payload <- list(
    tree_count = 1,
    tip_count = length(tree$tip.label),
    internal_node_count = tree$Nnode,
    edge_count = nrow(tree$edge),
    rooted = ape::is.rooted(tree),
    tip_labels = unname(tree$tip.label),
    branch_length_count = if (is.null(tree$edge.length)) 0 else sum(!is.na(tree$edge.length))
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(clade_rows, stringsAsFactors = FALSE)))
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
        clades = clades_path,
        normalized_tree = newick_path
      )
    )
  )
}

root_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  rooted_tree <- tryCatch(
    ape::root(
      tree,
      outgroup = unlist(case_payload$outgroup_taxa),
      resolve.root = TRUE
    ),
    error = function(error) error
  )
  if (inherits(rooted_tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeRootingError",
        message = conditionMessage(rooted_tree),
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
  clades_path <- file.path(output_root, "clades.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  clade_rows <- tree_structure_rows(rooted_tree, "")

  summary_payload <- list(
    tree_count = 1,
    tip_count = length(rooted_tree$tip.label),
    internal_node_count = rooted_tree$Nnode,
    edge_count = nrow(rooted_tree$edge),
    rooted = ape::is.rooted(rooted_tree),
    tip_labels = unname(rooted_tree$tip.label),
    branch_length_count = if (is.null(rooted_tree$edge.length)) 0 else sum(!is.na(rooted_tree$edge.length))
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(clade_rows, stringsAsFactors = FALSE)))
  ape::write.tree(rooted_tree, file = newick_path)

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
        clades = clades_path,
        normalized_tree = newick_path
      )
    )
  )
}

unroot_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  unrooted_tree <- tryCatch(ape::unroot(tree), error = function(error) error)
  if (inherits(unrooted_tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
        message = conditionMessage(unrooted_tree),
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
  clades_path <- file.path(output_root, "clades.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  clade_rows <- tree_structure_rows(unrooted_tree, "")

  summary_payload <- list(
    tree_count = 1,
    tip_count = length(unrooted_tree$tip.label),
    internal_node_count = unrooted_tree$Nnode,
    edge_count = nrow(unrooted_tree$edge),
    rooted = ape::is.rooted(unrooted_tree),
    tip_labels = unname(unrooted_tree$tip.label),
    branch_length_count = if (is.null(unrooted_tree$edge.length)) 0 else sum(!is.na(unrooted_tree$edge.length))
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(clade_rows, stringsAsFactors = FALSE)))
  ape::write.tree(unrooted_tree, file = newick_path)

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
        clades = clades_path,
        normalized_tree = newick_path
      )
    )
  )
}

drop_tip_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  excluded_taxa <- unlist(case_payload$excluded_taxa)
  original_tip_labels <- unname(tree$tip.label)
  warnings <- character()
  pruned_tree <- withCallingHandlers(
    tryCatch(ape::drop.tip(tree, excluded_taxa), error = function(error) error),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(pruned_tree, "error") || is.null(pruned_tree)) {
    message <- if (inherits(pruned_tree, "error")) {
      conditionMessage(pruned_tree)
    } else if (length(warnings) > 0) {
      warnings[[1]]
    } else {
      "drop.tip did not return a valid tree"
    }
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreePruningError",
        message = message,
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
  clades_path <- file.path(output_root, "clades.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  clade_rows <- tree_structure_rows(pruned_tree, "")
  retained_tip_labels <- unname(pruned_tree$tip.label)

  summary_payload <- list(
    tree_count = 1,
    tip_count = length(retained_tip_labels),
    internal_node_count = pruned_tree$Nnode,
    edge_count = nrow(pruned_tree$edge),
    rooted = ape::is.rooted(pruned_tree),
    tip_labels = as.list(retained_tip_labels),
    branch_length_count = if (is.null(pruned_tree$edge.length)) 0 else sum(!is.na(pruned_tree$edge.length)),
    dropped_taxa = as.list(sort(setdiff(original_tip_labels, retained_tip_labels))),
    absent_requested_taxa = as.list(sort(setdiff(excluded_taxa, original_tip_labels)))
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(clade_rows, stringsAsFactors = FALSE)))
  ape::write.tree(pruned_tree, file = newick_path)

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
        clades = clades_path,
        normalized_tree = newick_path
      )
    )
  )
}

keep_tip_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  requested_taxa <- unlist(case_payload$requested_taxa)
  kept_tree <- withCallingHandlers(
    tryCatch(ape::keep.tip(tree, requested_taxa), error = function(error) error),
    warning = function(warning) {
      invokeRestart("muffleWarning")
    }
  )
  if (inherits(kept_tree, "error") || is.null(kept_tree)) {
    message <- if (inherits(kept_tree, "error")) {
      conditionMessage(kept_tree)
    } else {
      "keep.tip did not return a valid tree"
    }
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreePruningError",
        message = message,
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
  clades_path <- file.path(output_root, "clades.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  clade_rows <- tree_structure_rows(kept_tree, "")
  retained_tip_labels <- unname(kept_tree$tip.label)
  original_tip_labels <- unname(tree$tip.label)

  summary_payload <- list(
    tree_count = 1,
    tip_count = length(retained_tip_labels),
    internal_node_count = kept_tree$Nnode,
    edge_count = nrow(kept_tree$edge),
    rooted = ape::is.rooted(kept_tree),
    tip_labels = as.list(retained_tip_labels),
    branch_length_count = if (is.null(kept_tree$edge.length)) 0 else sum(!is.na(kept_tree$edge.length)),
    requested_taxa = as.list(sort(unique(requested_taxa))),
    dropped_taxa = as.list(sort(setdiff(original_tip_labels, retained_tip_labels)))
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(clade_rows, stringsAsFactors = FALSE)))
  ape::write.tree(kept_tree, file = newick_path)

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
        clades = clades_path,
        normalized_tree = newick_path
      )
    )
  )
}

extract_clade_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  node_id <- as.integer(case_payload$node_id)
  subtree <- tryCatch(ape::extract.clade(tree, node_id), error = function(error) error)
  if (inherits(subtree, "error") || is.null(subtree)) {
    message <- if (inherits(subtree, "error")) {
      conditionMessage(subtree)
    } else {
      "extract.clade did not return a valid tree"
    }
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeCladeExtractionError",
        message = message,
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
  clades_path <- file.path(output_root, "clades.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  clade_rows <- tree_structure_rows(subtree, "")
  retained_tip_labels <- unname(subtree$tip.label)
  source_node_label <- ""
  if (!is.null(tree$node.label)) {
    internal_node_offset <- node_id - ape::Ntip(tree)
    if (internal_node_offset >= 1 && internal_node_offset <= length(tree$node.label)) {
      source_node_label <- tree$node.label[[internal_node_offset]]
    }
  }

  summary_payload <- list(
    tree_count = 1,
    tip_count = length(retained_tip_labels),
    internal_node_count = subtree$Nnode,
    edge_count = nrow(subtree$edge),
    rooted = ape::is.rooted(subtree),
    tip_labels = as.list(retained_tip_labels),
    branch_length_count = if (is.null(subtree$edge.length)) 0 else sum(!is.na(subtree$edge.length)),
    requested_node_id = node_id,
    matched_node_id = node_id,
    matched_node_name = source_node_label
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(clade_rows, stringsAsFactors = FALSE)))
  ape::write.tree(subtree, file = newick_path)

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
        clades = clades_path,
        normalized_tree = newick_path
      )
    )
  )
}

mrca_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  requested_taxa <- as.character(unlist(case_payload$mrca_taxa))
  duplicate_requested_taxa <- sort(unique(requested_taxa[duplicated(requested_taxa)]))
  unique_requested_taxa <- sort(unique(requested_taxa))
  matched_node_id <- tryCatch(ape::getMRCA(tree, requested_taxa), error = function(error) error)
  if (inherits(matched_node_id, "error") || is.null(matched_node_id)) {
    message <- if (inherits(matched_node_id, "error")) {
      conditionMessage(matched_node_id)
    } else {
      "getMRCA did not return an internal node"
    }
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeMrcaError",
        message = message,
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  matched_node_id <- as.integer(matched_node_id)
  matched_taxa <- descendant_taxa(tree, matched_node_id)
  summary_path <- file.path(output_root, "summary.json")
  write_payload(
    summary_path,
    list(
      requested_taxa = as.list(sort(requested_taxa)),
      unique_requested_taxa = as.list(unique_requested_taxa),
      duplicate_requested_taxa = as.list(duplicate_requested_taxa),
      matched_node_id = matched_node_id,
      matched_node_name = node_label(tree, matched_node_id),
      matched_taxa = as.list(matched_taxa),
      matched_extra_taxa = as.list(sort(setdiff(matched_taxa, unique_requested_taxa))),
      matched_tip_count = length(matched_taxa),
      is_root = identical(matched_node_id, root_node(tree))
    )
  )
  write_payload(
    execution_path,
    list(
      status = "ok",
      case_id = case_payload$case_id,
      function_name = case_payload$function_name,
      input_fixture = case_payload$input_fixture,
      r_version = r_version,
      ape_version = as.character(utils::packageVersion("ape")),
      outputs = list(summary_json = summary_path)
    )
  )
}

monophyly_case <- function(case_payload, output_root, execution_path, r_version) {
  tree <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
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

  requested_taxa <- as.character(unlist(case_payload$requested_taxa))
  duplicate_requested_taxa <- sort(unique(requested_taxa[duplicated(requested_taxa)]))
  unique_requested_taxa <- sort(unique(requested_taxa))
  present_requested_taxa <- sort(intersect(unique_requested_taxa, tree$tip.label))
  missing_requested_taxa <- sort(setdiff(unique_requested_taxa, tree$tip.label))
  reroot <- isTRUE(case_payload$monophyly_reroot)
  monophyletic <- tryCatch(
    ape::is.monophyletic(tree, requested_taxa, reroot = reroot),
    error = function(error) error
  )
  if (inherits(monophyletic, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeMonophylyError",
        message = conditionMessage(monophyletic),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  matched_node_id <- NULL
  matched_taxa <- character()
  matched_extra_taxa <- character()
  matched_node_name <- ""
  if (length(present_requested_taxa) == 1) {
    matched_node_id <- match(present_requested_taxa[[1]], tree$tip.label)
    matched_taxa <- present_requested_taxa
  } else if (length(present_requested_taxa) >= 2) {
    matched_node_id <- ape::getMRCA(tree, present_requested_taxa)
    matched_taxa <- descendant_taxa(tree, matched_node_id)
    matched_extra_taxa <- sort(setdiff(matched_taxa, present_requested_taxa))
    matched_node_name <- node_label(tree, matched_node_id)
  }

  summary_path <- file.path(output_root, "summary.json")
  write_payload(
    summary_path,
    list(
      requested_taxa = as.list(sort(requested_taxa)),
      unique_requested_taxa = as.list(unique_requested_taxa),
      duplicate_requested_taxa = as.list(duplicate_requested_taxa),
      missing_requested_taxa = as.list(missing_requested_taxa),
      present_requested_taxa = as.list(present_requested_taxa),
      reroot = reroot,
      rooted = ape::is.rooted(tree),
      monophyletic = isTRUE(monophyletic),
      complementary_clade_used = isTRUE(monophyletic) && length(matched_extra_taxa) > 0,
      matched_node_id = matched_node_id,
      matched_node_name = matched_node_name,
      matched_taxa = as.list(matched_taxa),
      matched_extra_taxa = as.list(matched_extra_taxa),
      matched_tip_count = length(matched_taxa),
      is_root = if (is.null(matched_node_id)) NULL else identical(matched_node_id, root_node(tree))
    )
  )
  write_payload(
    execution_path,
    list(
      status = "ok",
      case_id = case_payload$case_id,
      function_name = case_payload$function_name,
      input_fixture = case_payload$input_fixture,
      r_version = r_version,
      ape_version = as.character(utils::packageVersion("ape")),
      outputs = list(summary_json = summary_path)
    )
  )
}

tree_set_case <- function(case_payload, output_root, execution_path, r_version) {
  tree_set <- tryCatch(ape::read.tree(case_payload$input_fixture), error = function(error) error)
  if (inherits(tree_set, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
        message = conditionMessage(tree_set),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }
  if (!inherits(tree_set, "multiPhylo")) {
    tree_set <- structure(list(tree_set), class = "multiPhylo")
  }

  summary_path <- file.path(output_root, "summary.json")
  clades_path <- file.path(output_root, "clades.tsv")
  tree_set_path <- file.path(output_root, "normalized-tree-set.nwk")
  tree_rows <- unlist(
    lapply(seq_along(tree_set), function(index) tree_structure_rows(tree_set[[index]], index)),
    recursive = FALSE
  )
  shared_tip_labels <- Reduce(
    intersect,
    lapply(tree_set, function(tree) sort(unname(tree$tip.label)))
  )
  summary_payload <- list(
    tree_count = length(tree_set),
    source_format = "newick",
    tree_indices = as.list(seq_along(tree_set)),
    shared_tip_labels = shared_tip_labels,
    unique_tip_label_count = length(shared_tip_labels)
  )
  write_payload(summary_path, summary_payload)
  write_table(clades_path, do.call(rbind.data.frame, c(tree_rows, stringsAsFactors = FALSE)))
  ape::write.tree(tree_set, file = tree_set_path)
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
        clades = clades_path,
        normalized_tree_set = tree_set_path
      )
    )
  )
}

dna_base_frequency_case <- function(case_payload, output_root, execution_path, r_version) {
  alignment <- ape::read.dna(case_payload$input_fixture, format = "fasta")
  alignment_matrix <- as.matrix(alignment)
  frequencies <- ape::base.freq(alignment, all = TRUE)
  summary_path <- file.path(output_root, "summary.json")
  table_path <- file.path(output_root, "base-frequency.tsv")

  write_payload(
    summary_path,
    list(
      sequence_count = nrow(alignment_matrix),
      alignment_length = ncol(alignment_matrix),
      state_count = length(frequencies)
    )
  )
  write_table(
    table_path,
    data.frame(
      state = names(frequencies),
      frequency = as.numeric(frequencies),
      stringsAsFactors = FALSE
    )
  )
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
        base_frequency = table_path
      )
    )
  )
}

dna_distance_case <- function(case_payload, output_root, execution_path, r_version) {
  alignment <- ape::read.dna(case_payload$input_fixture, format = "fasta")
  alignment_matrix <- as.matrix(alignment)
  distances <- as.matrix(
    ape::dist.dna(
      alignment,
      model = "RAW",
      pairwise.deletion = isTRUE(case_payload$pairwise_deletion),
      as.matrix = TRUE
    )
  )
  summary_path <- file.path(output_root, "summary.json")
  table_path <- file.path(output_root, "distance-matrix.tsv")

  rows <- data.frame(
    left_identifier = character(),
    right_identifier = character(),
    distance = numeric(),
    stringsAsFactors = FALSE
  )
  for (left in rownames(distances)) {
    for (right in colnames(distances)) {
      rows[nrow(rows) + 1, ] <- list(left, right, as.numeric(distances[left, right]))
    }
  }

  write_payload(
    summary_path,
    list(
      sequence_count = nrow(alignment_matrix),
      alignment_length = ncol(alignment_matrix),
      pairwise_deletion = isTRUE(case_payload$pairwise_deletion)
    )
  )
  write_table(table_path, rows)
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
        distance_matrix = table_path
      )
    )
  )
}

dna_translation_case <- function(case_payload, output_root, execution_path, r_version) {
  alignment <- ape::read.dna(case_payload$input_fixture, format = "fasta")
  genetic_code_id <- if (is.null(case_payload$genetic_code_id)) 1 else case_payload$genetic_code_id
  translated <- ape::trans(alignment, code = genetic_code_id)
  translated_matrix <- as.character(translated)
  identifiers <- rownames(translated_matrix)
  amino_acid_sequences <- apply(translated_matrix, 1, paste0, collapse = "")
  if (is.null(names(amino_acid_sequences)) && length(identifiers) == length(amino_acid_sequences)) {
    names(amino_acid_sequences) <- identifiers
  }
  summary_path <- file.path(output_root, "summary.json")
  table_path <- file.path(output_root, "translation.tsv")

  write_payload(
    summary_path,
    list(
      sequence_count = length(amino_acid_sequences),
      translated_length = if (length(amino_acid_sequences) == 0) 0 else ncol(translated_matrix),
      stop_codon_count = sum(grepl("\\*", amino_acid_sequences))
    )
  )
  write_table(
    table_path,
    data.frame(
      identifier = names(amino_acid_sequences),
      amino_acid_sequence = unname(amino_acid_sequences),
      stringsAsFactors = FALSE
    )
  )
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
        translation = table_path
      )
    )
  )
}

if (identical(case_payload$operation, "read-tree-structure") ||
    identical(case_payload$operation, "write-tree-structure")) {
  tree_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "root-tree-outgroup")) {
  root_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "unroot-tree")) {
  unroot_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "drop-tree-taxa")) {
  drop_tip_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "keep-tree-taxa")) {
  keep_tip_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "extract-tree-clade")) {
  extract_clade_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "get-tree-mrca")) {
  mrca_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "assess-tree-monophyly")) {
  monophyly_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "read-tree-set-structure") ||
    identical(case_payload$operation, "write-tree-set-structure")) {
  tree_set_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "dna-base-frequency")) {
  dna_base_frequency_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "dna-raw-distance")) {
  dna_distance_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "dna-translation")) {
  dna_translation_case(case_payload, output_root, execution_path, r_version)
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
