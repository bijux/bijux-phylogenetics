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
  write_table(tips_path, tip_table)
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
  translated_matrix <- as.matrix(translated)
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

if (identical(case_payload$operation, "read-tree-summary")) {
  tree_case(case_payload, output_root, execution_path, r_version)
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
