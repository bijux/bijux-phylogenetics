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

matrix_positive_definite <- function(matrix) {
  !inherits(tryCatch(chol(matrix), error = function(error) error), "error")
}

matrix_log_determinant <- function(matrix) {
  determinant_report <- determinant(matrix, logarithm = TRUE)
  if (!isTRUE(determinant_report$sign > 0)) {
    return(NULL)
  }
  as.numeric(determinant_report$modulus)
}

matrix_numeric_rank <- function(matrix, tolerance = 1e-12) {
  qr(matrix, tol = tolerance)$rank
}

branch_length_range <- function(tree) {
  if (is.null(tree$edge.length) || any(is.na(tree$edge.length))) {
    return(list(minimum = NULL, maximum = NULL))
  }
  list(
    minimum = as.numeric(min(tree$edge.length)),
    maximum = as.numeric(max(tree$edge.length))
  )
}

tree_has_polytomy <- function(tree) {
  any(as.integer(table(tree$edge[, 1])) > 2L)
}

signature_id <- function(taxa) {
  paste(sort(unique(as.character(taxa))), collapse = "|")
}

rooted_topology_signatures <- function(tree) {
  root_id <- root_node(tree)
  total_tip_count <- length(tree$tip.label)
  signatures <- list()
  for (node in sort(unique(tree$edge[, 1]))) {
    if (identical(node, root_id)) {
      next
    }
    taxa <- descendant_taxa(tree, node)
    if (length(taxa) <= 1L || length(taxa) >= total_tip_count) {
      next
    }
    signatures[[signature_id(taxa)]] <- sort(taxa)
  }
  signatures
}

canonical_unrooted_signature <- function(taxa, all_taxa) {
  selected <- sort(unique(as.character(taxa)))
  complement <- sort(setdiff(all_taxa, selected))
  if (length(selected) < length(complement)) {
    return(selected)
  }
  if (length(complement) < length(selected)) {
    return(complement)
  }
  if (length(selected) == 0L) {
    return(selected)
  }
  if (paste(selected, collapse = "|") <= paste(complement, collapse = "|")) {
    return(selected)
  }
  complement
}

unrooted_topology_signatures <- function(tree) {
  all_taxa <- sort(unname(tree$tip.label))
  signatures <- list()
  if (length(all_taxa) < 4L) {
    return(signatures)
  }
  partitions <- ape::prop.part(tree)
  for (partition in partitions) {
    taxa <- sort(unname(tree$tip.label[as.integer(partition)]))
    selected <- canonical_unrooted_signature(taxa, all_taxa)
    if (length(selected) <= 1L || length(selected) >= length(all_taxa)) {
      next
    }
    signatures[[signature_id(selected)]] <- selected
  }
  signatures
}

topology_distance_rows <- function(left_tree, right_tree, rf_mode) {
  left_signatures <- if (identical(rf_mode, "rooted")) {
    rooted_topology_signatures(left_tree)
  } else {
    unrooted_topology_signatures(left_tree)
  }
  right_signatures <- if (identical(rf_mode, "rooted")) {
    rooted_topology_signatures(right_tree)
  } else {
    unrooted_topology_signatures(right_tree)
  }
  left_ids <- names(left_signatures)
  right_ids <- names(right_signatures)
  shared_ids <- intersect(left_ids, right_ids)
  left_only_ids <- setdiff(left_ids, right_ids)
  right_only_ids <- setdiff(right_ids, left_ids)
  all_ids <- union(left_ids, right_ids)
  signature_lookup <- c(left_signatures, right_signatures[setdiff(right_ids, left_ids)])
  split_kind <- if (identical(rf_mode, "rooted")) "clade" else "split"
  rows <- lapply(all_ids, function(split_id) {
    taxa <- signature_lookup[[split_id]]
    list(
      split_id = split_id,
      split_kind = split_kind,
      comparison_status = if (split_id %in% shared_ids) {
        "shared"
      } else if (split_id %in% left_only_ids) {
        "left_only"
      } else {
        "right_only"
      },
      taxon_count = length(taxa),
      descendant_taxa = paste(taxa, collapse = "|"),
      left_present = split_id %in% left_ids,
      right_present = split_id %in% right_ids
    )
  })
  order_frame <- do.call(
    rbind.data.frame,
    c(rows, stringsAsFactors = FALSE)
  )
  order_frame <- order_frame[order(order_frame$taxon_count, order_frame$descendant_taxa), ]
  list(
    rows = order_frame,
    left_count = length(left_ids),
    right_count = length(right_ids),
    shared_count = length(shared_ids),
    left_only_count = length(left_only_ids),
    right_only_count = length(right_only_ids)
  )
}

require_identical_tree_set_taxa <- function(tree_set, message_text) {
  reference_taxa <- sort(unname(tree_set[[1]]$tip.label))
  for (tree in tree_set[-1]) {
    if (!identical(sort(unname(tree$tip.label)), reference_taxa)) {
      stop(message_text)
    }
  }
  reference_taxa
}

consensus_clade_frequency_rows <- function(tree_set) {
  counts <- list()
  tree_count <- length(tree_set)
  for (tree in tree_set) {
    root_id <- root_node(tree)
    total_tip_count <- length(tree$tip.label)
    for (node in sort(unique(tree$edge[, 1]))) {
      if (identical(node, root_id)) {
        next
      }
      taxa <- descendant_taxa(tree, node)
      if (length(taxa) <= 1L || length(taxa) >= total_tip_count) {
        next
      }
      clade_id <- paste(taxa, collapse = "|")
      if (is.null(counts[[clade_id]])) {
        counts[[clade_id]] <- list(clade = clade_id, tree_count = 0L)
      }
      counts[[clade_id]]$tree_count <- counts[[clade_id]]$tree_count + 1L
    }
  }
  if (length(counts) == 0L) {
    return(data.frame(
      clade = character(),
      tree_count = integer(),
      frequency = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(sort(names(counts)), function(clade_id) {
    count <- counts[[clade_id]]$tree_count
    list(
      clade = clade_id,
      tree_count = count,
      frequency = count / tree_count
    )
  })
  do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
}

clade_support_status <- function(supporting_tree_count, tree_count, node_kind, unscored_reason = NULL) {
  if (identical(node_kind, "root")) {
    return(list(
      support_status = "fixed",
      explanation = "the root spans the full compatible taxon set and is present in every comparison tree"
    ))
  }
  if (is.na(supporting_tree_count)) {
    if (identical(unscored_reason, "absent-root-split")) {
      return(list(
        support_status = "not-counted",
        explanation = "ape::prop.clades leaves this root-adjacent split unscored when the comparison tree set never realizes the matching bipartition"
      ))
    }
    return(list(
      support_status = "not-counted",
      explanation = "ape::prop.clades leaves this root-adjacent clade unscored because its complement is a singleton tip"
    ))
  }
  if (identical(supporting_tree_count, 0L)) {
    return(list(
      support_status = "absent",
      explanation = "the reference clade is absent from the comparison tree set"
    ))
  }
  if (identical(supporting_tree_count, tree_count)) {
    return(list(
      support_status = "fixed",
      explanation = "the reference clade is present in every comparison tree"
    ))
  }
  list(
    support_status = "partial-support",
    explanation = "the reference clade is present in only a subset of comparison trees"
  )
}

prop_clades_case <- function(case_payload, output_root, execution_path, r_version) {
  reference_tree <- tryCatch(
    ape::read.tree(case_payload$reference_tree_path),
    error = function(error) error
  )
  if (inherits(reference_tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
        message = conditionMessage(reference_tree),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  comparison_tree_set <- tryCatch(
    ape::read.tree(case_payload$input_fixture),
    error = function(error) error
  )
  if (inherits(comparison_tree_set, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeParseError",
        message = conditionMessage(comparison_tree_set),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }
  if (inherits(comparison_tree_set, "phylo")) {
    comparison_tree_set <- structure(list(comparison_tree_set), class = "multiPhylo")
  }

  compatible_taxa <- tryCatch(
    require_identical_tree_set_taxa(
      comparison_tree_set,
      "reference tree support mapping requires all comparison trees to share the exact same taxon set"
    ),
    error = function(error) error
  )
  if (inherits(compatible_taxa, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "PropCladesError",
        message = conditionMessage(compatible_taxa),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  reference_taxa <- sort(vapply(reference_tree$tip.label, normalize_node_label, character(1)))
  compatible_taxa <- sort(vapply(compatible_taxa, normalize_node_label, character(1)))
  if (!identical(reference_taxa, compatible_taxa)) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "PropCladesError",
        message = "reference tree and comparison tree set must share the exact same taxon set",
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  prop_clades <- tryCatch(
    ape::prop.clades(reference_tree, comparison_tree_set),
    error = function(error) error
  )
  if (inherits(prop_clades, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "PropCladesError",
        message = conditionMessage(prop_clades),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  tree_count <- length(comparison_tree_set)
  root_id <- root_node(reference_tree)
  root_children <- reference_tree$edge[reference_tree$edge[, 1] == root_id, 2]
  depths <- as.numeric(ape::node.depth.edgelength(reference_tree))
  internal_node_ids <- seq.int(
    length(reference_tree$tip.label) + 1,
    length(reference_tree$tip.label) + reference_tree$Nnode
  )
  rows <- data.frame(
    node_id = integer(),
    node_kind = character(),
    node_label = character(),
    descendant_taxa = character(),
    supporting_tree_count = character(),
    clade_frequency = character(),
    support_percent = character(),
    support_status = character(),
    explanation = character(),
    reference_branch_length = character(),
    reference_root_depth = character(),
    stringsAsFactors = FALSE
  )
  supported_clade_count <- 0L
  absent_clade_count <- 0L
  unscored_clade_count <- 0L

  for (index in seq_along(internal_node_ids)) {
    node_id <- internal_node_ids[[index]]
    taxa <- descendant_taxa(reference_tree, node_id)
    node_kind_value <- node_kind(reference_tree, node_id, root_id)
    supporting_tree_count <- prop_clades[[index]]
    unscored_reason <- NULL
    if (is.na(supporting_tree_count)) {
      supporting_tree_count_value <- NA_integer_
      clade_frequency <- NA_real_
      support_percent <- NA_real_
      if (length(taxa) == length(reference_taxa) - 1L) {
        unscored_reason <- "singleton-complement"
      } else if (node_id %in% root_children) {
        unscored_reason <- "absent-root-split"
      }
      unscored_clade_count <- unscored_clade_count + 1L
    } else {
      supporting_tree_count_value <- as.integer(supporting_tree_count)
      clade_frequency <- as.numeric(supporting_tree_count_value / tree_count)
      support_percent <- as.numeric(clade_frequency * 100.0)
      if (identical(node_kind_value, "root")) {
        NULL
      } else if (identical(supporting_tree_count_value, 0L)) {
        absent_clade_count <- absent_clade_count + 1L
      } else {
        supported_clade_count <- supported_clade_count + 1L
      }
    }
    status_payload <- clade_support_status(
      supporting_tree_count = supporting_tree_count_value,
      tree_count = tree_count,
      node_kind = node_kind_value,
      unscored_reason = unscored_reason
    )
    branch_length <- node_branch_length(reference_tree, node_id)
    rows[nrow(rows) + 1, ] <- list(
      as.integer(node_id),
      node_kind_value,
      node_label(reference_tree, node_id),
      paste(taxa, collapse = "|"),
      if (is.na(supporting_tree_count_value)) "" else as.character(supporting_tree_count_value),
      if (is.na(clade_frequency)) "" else as.character(clade_frequency),
      if (is.na(support_percent)) "" else as.character(support_percent),
      status_payload$support_status,
      status_payload$explanation,
      if (identical(branch_length, "")) "" else as.character(branch_length),
      if (is.na(depths[[node_id]])) "" else as.character(depths[[node_id]])
    )
  }

  summary_path <- file.path(output_root, "summary.json")
  rows_path <- file.path(output_root, "support-table.tsv")
  write_payload(
    summary_path,
    list(
      tree_count = tree_count,
      shared_taxa = unname(reference_taxa),
      shared_taxon_count = length(reference_taxa),
      internal_node_count = length(internal_node_ids),
      supported_clade_count = supported_clade_count,
      absent_clade_count = absent_clade_count,
      unscored_clade_count = unscored_clade_count
    )
  )
  write_table(rows_path, rows)
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
        support_table = rows_path
      )
    )
  )
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
    matched_node_name <- node_label(tree, matched_node_id)
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

consensus_case <- function(case_payload, output_root, execution_path, r_version) {
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

  shared_taxa <- tryCatch(
    require_identical_tree_set_taxa(
      tree_set,
      "consensus requires all trees to share the exact same taxon set"
    ),
    error = function(error) error
  )
  if (inherits(shared_taxa, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "ConsensusTreeError",
        message = conditionMessage(shared_taxa),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  consensus_method <- if (is.null(case_payload$consensus_method)) {
    "majority-rule"
  } else {
    as.character(case_payload$consensus_method)
  }
  consensus_threshold <- if (identical(consensus_method, "strict")) {
    1.0
  } else if (identical(consensus_method, "majority-rule")) {
    0.5
  } else {
    NA_real_
  }
  if (is.na(consensus_threshold)) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "ConsensusTreeError",
        message = paste("unsupported consensus method:", consensus_method),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  consensus_tree <- tryCatch(
    ape::consensus(tree_set, p = consensus_threshold),
    error = function(error) error
  )
  if (inherits(consensus_tree, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "ConsensusTreeError",
        message = conditionMessage(consensus_tree),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  frequency_rows <- consensus_clade_frequency_rows(tree_set)
  summary_path <- file.path(output_root, "summary.json")
  frequency_path <- file.path(output_root, "clade-frequencies.tsv")
  newick_path <- file.path(output_root, "normalized-tree.nwk")
  write_payload(
    summary_path,
    list(
      tree_count = length(tree_set),
      shared_taxa = as.list(shared_taxa),
      shared_taxon_count = length(shared_taxa),
      tip_count = length(consensus_tree$tip.label),
      rooted = ape::is.rooted(consensus_tree),
      consensus_method = consensus_method,
      consensus_threshold = consensus_threshold,
      included_clade_count = length(rooted_topology_signatures(consensus_tree)),
      clade_frequency_count = nrow(frequency_rows)
    )
  )
  write_table(frequency_path, frequency_rows)
  ape::write.tree(consensus_tree, file = newick_path)
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
        clade_frequencies = frequency_path,
        normalized_tree = newick_path
      )
    )
  )
}

topology_distance_case <- function(case_payload, output_root, execution_path, r_version) {
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
  if (length(tree_set) != 2L) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TopologyDistanceError",
        message = "ape topology-distance parity fixtures must contain exactly two trees",
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  rf_mode <- if (is.null(case_payload$rf_mode)) "rooted" else as.character(case_payload$rf_mode)
  left_tree <- tree_set[[1]]
  right_tree <- tree_set[[2]]
  split_report <- topology_distance_rows(left_tree, right_tree, rf_mode)
  shared_taxa <- sort(intersect(unname(left_tree$tip.label), unname(right_tree$tip.label)))
  distance <- as.numeric(ape::dist.topo(left_tree, right_tree))
  denominator <- split_report$left_count + split_report$right_count
  normalized_distance <- if (identical(denominator, 0L)) 0.0 else distance / denominator

  summary_path <- file.path(output_root, "summary.json")
  split_path <- file.path(output_root, "split-table.tsv")
  write_payload(
    summary_path,
    list(
      tip_count = length(shared_taxa),
      shared_taxa = as.list(shared_taxa),
      left_only_taxa = as.list(sort(setdiff(unname(left_tree$tip.label), unname(right_tree$tip.label)))),
      right_only_taxa = as.list(sort(setdiff(unname(right_tree$tip.label), unname(left_tree$tip.label)))),
      taxon_overlap_policy = "require-identical",
      rf_mode = rf_mode,
      rooted_left = ape::is.rooted(left_tree),
      rooted_right = ape::is.rooted(right_tree),
      polytomy_present_left = tree_has_polytomy(left_tree),
      polytomy_present_right = tree_has_polytomy(right_tree),
      left_split_count = split_report$left_count,
      right_split_count = split_report$right_count,
      shared_split_count = split_report$shared_count,
      left_only_split_count = split_report$left_only_count,
      right_only_split_count = split_report$right_only_count,
      robinson_foulds_distance = distance,
      normalized_robinson_foulds = normalized_distance,
      topology_equal = identical(distance, 0)
    )
  )
  write_table(split_path, split_report$rows)
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
        split_table = split_path
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
  ape_distance_model <- function(model_name) {
    normalized <- tolower(as.character(model_name))
    if (identical(normalized, "raw")) {
      return("RAW")
    }
    if (identical(normalized, "jc69")) {
      return("JC69")
    }
    if (identical(normalized, "k80")) {
      return("K80")
    }
    as.character(model_name)
  }
  canonical_distance_model <- function(model_name) {
    normalized <- tolower(as.character(model_name))
    if (identical(normalized, "raw")) {
      return("p-distance")
    }
    if (identical(normalized, "jc69")) {
      return("jukes-cantor")
    }
    if (identical(normalized, "k80")) {
      return("kimura-2-parameter")
    }
    normalized
  }
  alignment <- tryCatch(
    ape::read.dna(case_payload$input_fixture, format = "fasta"),
    error = function(error) error
  )
  if (inherits(alignment, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        error_type = "AlignmentParseError",
        message = conditionMessage(alignment),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }
  distance_result <- tryCatch(
    ape::dist.dna(
      alignment,
      model = ape_distance_model(case_payload$distance_model),
      pairwise.deletion = isTRUE(case_payload$pairwise_deletion),
      as.matrix = TRUE
    ),
    error = function(error) error
  )
  if (inherits(distance_result, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        error_type = "DnaDistanceError",
        message = conditionMessage(distance_result),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }
  alignment_matrix <- as.matrix(alignment)
  distances <- as.matrix(distance_result)
  summary_path <- file.path(output_root, "summary.json")
  table_path <- file.path(output_root, "distance-matrix.tsv")

  rows <- data.frame(
    left_identifier = character(),
    right_identifier = character(),
    distance = character(),
    distance_status = character(),
    stringsAsFactors = FALSE
  )
  finite_distance_count <- 0L
  undefined_distance_count <- 0L
  infinite_distance_count <- 0L
  for (left in rownames(distances)) {
    for (right in colnames(distances)) {
      distance_value <- as.numeric(distances[left, right])
      if (is.nan(distance_value)) {
        distance_text <- ""
        distance_status <- "undefined"
        undefined_distance_count <- undefined_distance_count + 1L
      } else if (is.infinite(distance_value)) {
        distance_text <- ""
        distance_status <- "infinite"
        infinite_distance_count <- infinite_distance_count + 1L
      } else {
        distance_text <- format(distance_value, digits = 16, scientific = FALSE, trim = TRUE)
        distance_status <- "finite"
        finite_distance_count <- finite_distance_count + 1L
      }
      rows[nrow(rows) + 1, ] <- list(left, right, distance_text, distance_status)
    }
  }

  write_payload(
    summary_path,
    list(
      sequence_count = nrow(alignment_matrix),
      alignment_length = ncol(alignment_matrix),
      pairwise_deletion = isTRUE(case_payload$pairwise_deletion),
      distance_model = canonical_distance_model(case_payload$distance_model),
      finite_distance_count = finite_distance_count,
      undefined_distance_count = undefined_distance_count,
      infinite_distance_count = infinite_distance_count
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

cophenetic_case <- function(case_payload, output_root, execution_path, r_version) {
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

  distances <- as.matrix(ape::cophenetic.phylo(tree))
  summary_path <- file.path(output_root, "summary.json")
  matrix_path <- file.path(output_root, "tip-distance-matrix.tsv")
  long_path <- file.path(output_root, "tip-distance-long.tsv")

  matrix_rows <- data.frame(
    taxon = rownames(distances),
    distances,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  long_rows <- data.frame(
    left_identifier = character(),
    right_identifier = character(),
    distance = numeric(),
    stringsAsFactors = FALSE
  )
  for (left in rownames(distances)) {
    for (right in colnames(distances)) {
      long_rows[nrow(long_rows) + 1, ] <- list(
        left,
        right,
        as.numeric(distances[left, right])
      )
    }
  }

  write_payload(
    summary_path,
    list(
      tip_count = length(tree$tip.label),
      rooted = ape::is.rooted(tree),
      tip_labels = as.list(unname(tree$tip.label)),
      pair_count = nrow(long_rows),
      diagonal_zero = all(diag(distances) == 0),
      symmetric = isTRUE(all.equal(distances, t(distances))),
      complete_branch_lengths = !is.null(tree$edge.length) && all(!is.na(tree$edge.length)),
      missing_branch_length_policy = "error"
    )
  )
  write_table(matrix_path, matrix_rows)
  write_table(long_path, long_rows)
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
        tip_distance_matrix = matrix_path,
        tip_distance_long = long_path
      )
    )
  )
}

vcv_case <- function(case_payload, output_root, execution_path, r_version) {
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

  covariance <- as.matrix(ape::vcv.phylo(tree))
  root_depths <- diag(covariance)
  branch_lengths <- branch_length_range(tree)
  covariance_rank <- matrix_numeric_rank(covariance)
  singular <- covariance_rank < nrow(covariance)
  positive_definite <- matrix_positive_definite(covariance)
  raw_log_determinant <- if (positive_definite) matrix_log_determinant(covariance) else NULL
  if (is.null(raw_log_determinant)) {
    positive_definite <- FALSE
  }
  condition_number <- if (singular) NULL else as.numeric(kappa(covariance, exact = TRUE))
  near_singular <- isTRUE(singular) || (!is.null(condition_number) && condition_number >= 1e12)

  summary_path <- file.path(output_root, "summary.json")
  matrix_path <- file.path(output_root, "covariance-matrix.tsv")
  long_path <- file.path(output_root, "covariance-long.tsv")

  matrix_rows <- data.frame(
    taxon = rownames(covariance),
    covariance,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  long_rows <- data.frame(
    left_taxon = character(),
    right_taxon = character(),
    shared_ancestry_covariance = numeric(),
    stringsAsFactors = FALSE
  )
  for (left in rownames(covariance)) {
    for (right in colnames(covariance)) {
      long_rows[nrow(long_rows) + 1, ] <- list(
        left,
        right,
        as.numeric(covariance[left, right])
      )
    }
  }

  write_payload(
    summary_path,
    list(
      tip_count = length(tree$tip.label),
      rooted = ape::is.rooted(tree),
      tip_labels = as.list(unname(tree$tip.label)),
      pair_count = nrow(long_rows),
      tree_is_ultrametric = diff(range(root_depths)) <= 1e-12,
      minimum_root_to_tip_depth = as.numeric(min(root_depths)),
      maximum_root_to_tip_depth = as.numeric(max(root_depths)),
      minimum_branch_length = branch_lengths$minimum,
      maximum_branch_length = branch_lengths$maximum,
      matrix_dimension = nrow(covariance),
      matrix_rank = covariance_rank,
      singular = singular,
      near_singular = near_singular,
      positive_definite = positive_definite,
      condition_number = condition_number,
      raw_log_determinant = raw_log_determinant
    )
  )
  write_table(matrix_path, matrix_rows)
  write_table(long_path, long_rows)
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
        covariance_matrix = matrix_path,
        covariance_long = long_path
      )
    )
  )
}

node_depth_case <- function(case_payload, output_root, execution_path, r_version) {
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

  tip_count <- length(tree$tip.label)
  internal_count <- tree$Nnode
  root_id <- root_node(tree)
  depths <- as.numeric(ape::node.depth.edgelength(tree))
  tip_depths <- depths[seq_len(tip_count)]
  internal_node_ids <- seq.int(tip_count + 1, tip_count + internal_count)
  internal_depths <- depths[internal_node_ids]
  zero_branch_length_count <- if (is.null(tree$edge.length)) 0L else sum(tree$edge.length == 0)

  rows <- data.frame(
    node_id = integer(),
    node_kind = character(),
    node_label = character(),
    descendant_taxa = character(),
    branch_length_depth = numeric(),
    branch_length = character(),
    stringsAsFactors = FALSE
  )
  for (node_id in c(seq_len(tip_count), internal_node_ids)) {
    branch_length <- node_branch_length(tree, node_id)
    rows[nrow(rows) + 1, ] <- list(
      as.integer(node_id),
      node_kind(tree, node_id, root_id),
      node_label(tree, node_id),
      paste(descendant_taxa(tree, node_id), collapse = "|"),
      as.numeric(depths[[node_id]]),
      if (identical(branch_length, "")) "" else as.character(branch_length)
    )
  }

  summary_path <- file.path(output_root, "summary.json")
  rows_path <- file.path(output_root, "node-depths.tsv")
  write_payload(
    summary_path,
    list(
      node_count = nrow(rows),
      tip_count = tip_count,
      internal_node_count = internal_count,
      rooted = ape::is.rooted(tree),
      tip_labels = as.list(unname(tree$tip.label)),
      tree_is_ultrametric = diff(range(tip_depths)) <= 1e-12,
      zero_branch_length_count = as.integer(zero_branch_length_count),
      minimum_tip_depth = as.numeric(min(tip_depths)),
      maximum_tip_depth = as.numeric(max(tip_depths)),
      minimum_internal_depth = as.numeric(min(internal_depths)),
      maximum_internal_depth = as.numeric(max(internal_depths))
    )
  )
  write_table(rows_path, rows)
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
        node_depths = rows_path
      )
    )
  )
}

branching_times_case <- function(case_payload, output_root, execution_path, r_version) {
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

  tip_count <- length(tree$tip.label)
  depths <- as.numeric(ape::node.depth.edgelength(tree))
  tip_depths <- depths[seq_len(tip_count)]
  internal_node_ids <- seq.int(tip_count + 1, tip_count + tree$Nnode)
  branching_times <- as.numeric(ape::branching.times(tree))
  zero_branch_length_count <- if (is.null(tree$edge.length)) 0L else sum(tree$edge.length == 0)
  root_age <- max(branching_times)

  rows <- data.frame(
    node_id = integer(),
    node_kind = character(),
    node_label = character(),
    descendant_taxa = character(),
    node_depth = numeric(),
    branching_time = numeric(),
    stringsAsFactors = FALSE
  )
  for (index in seq_along(internal_node_ids)) {
    node_id <- internal_node_ids[[index]]
    rows[nrow(rows) + 1, ] <- list(
      as.integer(node_id),
      node_kind(tree, node_id, root_node(tree)),
      node_label(tree, node_id),
      paste(descendant_taxa(tree, node_id), collapse = "|"),
      as.numeric(depths[[node_id]]),
      as.numeric(branching_times[[index]])
    )
  }

  summary_path <- file.path(output_root, "summary.json")
  rows_path <- file.path(output_root, "branching-times.tsv")
  write_payload(
    summary_path,
    list(
      internal_node_count = tree$Nnode,
      rooted = ape::is.rooted(tree),
      tip_labels = as.list(unname(tree$tip.label)),
      tree_is_ultrametric = diff(range(tip_depths)) <= 1e-12,
      root_age = as.numeric(root_age),
      zero_branch_length_count = as.integer(zero_branch_length_count),
      minimum_tip_depth = as.numeric(min(tip_depths)),
      maximum_tip_depth = as.numeric(max(tip_depths)),
      max_tip_depth_deviation = as.numeric(diff(range(tip_depths))),
      tolerance = 1e-12
    )
  )
  write_table(rows_path, rows)
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
        branching_times = rows_path
      )
    )
  )
}

ultrametric_case <- function(case_payload, output_root, execution_path, r_version) {
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

  tip_count <- length(tree$tip.label)
  tip_labels <- vapply(tree$tip.label, normalize_node_label, character(1))
  tip_depths <- as.numeric(ape::node.depth.edgelength(tree))[seq_len(tip_count)]
  tolerance <- if (is.null(case_payload$tolerance)) .Machine$double.eps^0.5 else as.numeric(case_payload$tolerance)
  option <- if (is.null(case_payload$ultrametric_option)) 1L else as.integer(case_payload$ultrametric_option)
  ultrametric <- tryCatch(
    ape::is.ultrametric(tree, tol = tolerance, option = option),
    error = function(error) error
  )
  if (inherits(ultrametric, "error")) {
    write_payload(
      execution_path,
      list(
        status = "failed",
        mismatch_reason = "reference_execution_failed",
        error_type = "TreeUltrametricError",
        message = conditionMessage(ultrametric),
        case_id = case_payload$case_id,
        function_name = case_payload$function_name,
        input_fixture = case_payload$input_fixture,
        r_version = r_version,
        ape_version = as.character(utils::packageVersion("ape"))
      )
    )
    quit(save = "no", status = 0)
  }

  minimum_tip_depth <- min(tip_depths)
  maximum_tip_depth <- max(tip_depths)
  mean_tip_depth <- mean(tip_depths)
  max_tip_depth_deviation <- maximum_tip_depth - minimum_tip_depth
  criterion_name <- if (identical(option, 1L)) "scaled-range" else "variance"
  criterion_value <- if (identical(option, 1L)) {
    if (isTRUE(all.equal(maximum_tip_depth, 0.0, tolerance = 1e-15))) {
      if (isTRUE(all.equal(max_tip_depth_deviation, 0.0, tolerance = 1e-15))) {
        0.0
      } else {
        Inf
      }
    } else {
      max_tip_depth_deviation / maximum_tip_depth
    }
  } else if (length(tip_depths) <= 1) {
    0.0
  } else {
    stats::var(tip_depths)
  }
  offending_taxa <- if (isTRUE(all.equal(max_tip_depth_deviation, 0.0, tolerance = 1e-12))) {
    character()
  } else {
    sort(
      unique(
        vapply(
          tree$tip.label[
            abs(tip_depths - minimum_tip_depth) <= 1e-12
            | abs(tip_depths - maximum_tip_depth) <= 1e-12
          ],
          normalize_node_label,
          character(1)
        )
      )
    )
  }

  summary_path <- file.path(output_root, "summary.json")
  rows_path <- file.path(output_root, "ultrametric-diagnostics.tsv")
  write_payload(
    summary_path,
    list(
      tip_count = tip_count,
      rooted = ape::is.rooted(tree),
      tip_labels = as.list(unname(tip_labels)),
      ultrametric = isTRUE(ultrametric),
      criterion_name = criterion_name,
      criterion_value = as.numeric(criterion_value),
      tolerance = tolerance,
      option = as.integer(option),
      minimum_tip_depth = minimum_tip_depth,
      maximum_tip_depth = maximum_tip_depth,
      mean_tip_depth = mean_tip_depth,
      max_tip_depth_deviation = max_tip_depth_deviation,
      root_age = maximum_tip_depth,
      offending_taxa = as.list(offending_taxa)
    )
  )
  write_table(
    rows_path,
    data.frame(
      node_id = seq_len(tip_count),
      tip_label = tip_labels,
      root_to_tip_depth = tip_depths,
      deviation_from_mean_depth = abs(tip_depths - mean_tip_depth),
      deviation_from_min_depth = tip_depths - minimum_tip_depth,
      deviation_from_max_depth = maximum_tip_depth - tip_depths,
      is_offending_taxon = tip_labels %in% offending_taxa,
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
        ultrametric_diagnostics = rows_path
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

if (identical(case_payload$operation, "tree-consensus")) {
  consensus_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-clade-support")) {
  prop_clades_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "dna-base-frequency")) {
  dna_base_frequency_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "dna-distance")) {
  dna_distance_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-tip-distance")) {
  cophenetic_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-topology-distance")) {
  topology_distance_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-brownian-covariance")) {
  vcv_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-node-depth")) {
  node_depth_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-branching-times")) {
  branching_times_case(case_payload, output_root, execution_path, r_version)
  quit(save = "no", status = 0)
}

if (identical(case_payload$operation, "tree-ultrametricity")) {
  ultrametric_case(case_payload, output_root, execution_path, r_version)
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
