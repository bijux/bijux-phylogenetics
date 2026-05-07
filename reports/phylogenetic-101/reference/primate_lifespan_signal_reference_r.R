#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(jsonlite)
  library(openxlsx)
  library(phytools)
  library(tidytree)
  library(treeio)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: primate_lifespan_signal_reference_r.R <r_repo_root> <out_dir>")
}

r_repo_root <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_dir <- file.path(r_repo_root, "PCM1_plots_signal", "Lecture", "R", "Data")
script_path <- file.path(
  r_repo_root,
  "PCM1_plots_signal",
  "Lecture",
  "R",
  "Rcripts",
  "PCM1_plots_signal.R"
)

checked_in_primate <- read.csv(file.path(data_dir, "primate.csv"), stringsAsFactors = FALSE)
if (names(checked_in_primate)[1] %in% c("", "X")) {
  checked_in_primate <- checked_in_primate[, -1, drop = FALSE]
}

primate_raw <- readWorkbook(file.path(data_dir, "primate_raw.xlsx"))
primate_raw$sex_dimorphism <- as.factor(primate_raw$sex_dimorphism)
primate_raw <- primate_raw %>% mutate(across(body_mass:social_group_size, ~as.numeric(.)))

primate <- primate_raw %>%
  group_by(family, species) %>%
  summarise(
    across(
      body_mass:social_group_size,
      .fns = ~mean(.x, na.rm = TRUE),
      .names = "{col}"
    ),
    sex_dimorphism = sex_dimorphism[!is.na(sex_dimorphism)],
    mating_system = mating_system[!is.na(mating_system)],
    .groups = "drop"
  )

reference_primate <- primate %>% arrange(species)
checked_in_primate_sorted <- checked_in_primate %>% arrange(species)
data_columns <- names(reference_primate)
checked_in_primate_sorted <- checked_in_primate_sorted[, data_columns]
numeric_columns <- c(
  "body_mass",
  "gestation",
  "home_range",
  "longevity",
  "social_group_size"
)
processed_match <- all(vapply(seq_len(nrow(reference_primate)), function(index) {
  all(vapply(data_columns, function(column_name) {
    left_value <- reference_primate[[column_name]][index]
    right_value <- checked_in_primate_sorted[[column_name]][index]
    if (column_name %in% numeric_columns) {
      return(abs(as.numeric(left_value) - as.numeric(right_value)) <= 1e-9)
    }
    identical(as.character(left_value), as.character(right_value))
  }, logical(1)))
}, logical(1)))

write.csv(reference_primate, file.path(out_dir, "reference_primate.csv"), row.names = FALSE)

original_tree <- read.tree(file.path(data_dir, "primatetree.nex"))
original_tree <- makeNodeLabel(original_tree)
missing_tips <- original_tree$tip.label[original_tree$tip.label %in% reference_primate$species == "FALSE"]
trimmed_tree <- ape::drop.tip(original_tree, missing_tips)
checked_in_trimmed_tree <- read.tree(file.path(data_dir, "trimmed_primatetree.nex"))

write.tree(trimmed_tree, file.path(out_dir, "reference_trimmed_primatetree.nwk"))

node_label_from_numeric <- function(tree, numeric_node) {
  tree$node.label[numeric_node - length(tree$tip.label)]
}

tip_id_lookup <- setNames(seq_along(trimmed_tree$tip.label), trimmed_tree$tip.label)
nodeid_examples <- list(
  Pan_paniscus = unname(treeio::nodeid(trimmed_tree, "Pan_paniscus")),
  Hylobates_lar = unname(treeio::nodeid(trimmed_tree, "Hylobates_lar")),
  Node32 = unname(treeio::nodeid(trimmed_tree, "Node32"))
)

tip_order <- data.frame(
  tip = trimmed_tree$tip.label,
  order = seq_along(trimmed_tree$tip.label),
  stringsAsFactors = FALSE
)
primate_aligned <- checked_in_primate %>%
  mutate(tip_order = tip_order$order[match(species, tip_order$tip)]) %>%
  arrange(tip_order)
primate_aligned$node <- treeio::nodeid(trimmed_tree, primate_aligned$species)
p_tree_data <- treeio::full_join(trimmed_tree, primate_aligned, by = "node")

tt77 <- extract.clade(trimmed_tree, 77)
rt130 <- rotateNodes(trimmed_tree, 130)
rt_all <- rotateNodes(trimmed_tree, "all")
unrooted_tree <- unroot(trimmed_tree)

set.seed(1)
random_tree <- rcoal(30)
write.tree(random_tree, file.path(out_dir, "random_tree_seed1.nwk"))

random_examples <- list(
  random_data = rTraitCont(random_tree, model = "BM", sigma = 0.5, root.value = 1),
  random_data2 = rTraitCont(random_tree, model = "BM", sigma = 0.5, root.value = 10),
  random_data3 = rTraitCont(random_tree, model = "BM", sigma = 5, root.value = 1),
  random_data4 = rTraitCont(random_tree, model = "OU", sigma = 5, root.value = 1, alpha = 0),
  random_data5 = rTraitCont(random_tree, model = "OU", sigma = 5, root.value = 1, alpha = 5)
)

random_results <- lapply(names(random_examples), function(name) {
  values <- random_examples[[name]]
  traits_path <- file.path(out_dir, paste0(name, ".csv"))
  write.csv(
    data.frame(species = names(values), value = as.numeric(values), stringsAsFactors = FALSE),
    traits_path,
    row.names = FALSE
  )
  fit <- if (name == "random_data5") {
    fitContinuous(random_tree, values, model = "lambda", control = list(niter = 1000))
  } else {
    fitContinuous(random_tree, values, model = "lambda")
  }
  list(
    name = name,
    traits_path = traits_path,
    lambda_value = unname(fit$opt$lambda),
    log_likelihood = unname(fit$opt$lnL),
    tip_count = length(values)
  )
})

primate_ll <- as.numeric(primate_aligned$longevity[match(trimmed_tree$tip.label, primate_aligned$species)])
names(primate_ll) <- trimmed_tree$tip.label
lambda_ll <- fitContinuous(trimmed_tree, primate_ll, model = "lambda")
trimmed_tree_lambda0 <- rescale(trimmed_tree, "lambda", 0)
lambda_ll0 <- fitContinuous(trimmed_tree_lambda0, primate_ll, model = "lambda")
ll_diff0 <- -2 * (lambda_ll0$opt$lnL - lambda_ll$opt$lnL)
p_value <- pchisq(ll_diff0, df = 1, lower.tail = FALSE)

long_ace <- ace(primate_aligned$longevity, trimmed_tree, type = "continuous", method = "pic")
primate_ace_long <- data.frame(
  species = c(primate_aligned$species, names(long_ace$ace)),
  longevity = c(primate_aligned$longevity, long_ace$ace),
  stringsAsFactors = FALSE
)
primate_ace_long$node <- treeio::nodeid(trimmed_tree, primate_ace_long$species)

get_desc <- function(tree, node) {
  kids <- tree$edge[tree$edge[, 1] == node, 2]
  out <- c()
  for (kid in kids) {
    if (kid <= length(tree$tip.label)) {
      out <- c(out, tree$tip.label[kid])
    } else {
      out <- c(out, get_desc(tree, kid))
    }
  }
  sort(out)
}

ancestral_nodes <- lapply(seq_along(long_ace$ace), function(index) {
  species <- names(long_ace$ace)[index]
  node <- treeio::nodeid(trimmed_tree, species)
  descendants <- get_desc(trimmed_tree, node)
  list(
    species = species,
    node = unname(node),
    descendant_taxa = descendants,
    signature = paste(descendants, collapse = "|"),
    estimate = unname(long_ace$ace[index]),
    lower_95 = unname(long_ace$CI95[index, 1]),
    upper_95 = unname(long_ace$CI95[index, 2])
  )
})

mrca_node <- getMRCA(trimmed_tree, tip = c("Pan_paniscus", "Hylobates_lar"))
tree_df <- tidytree::as_tibble(trimmed_tree) %>%
  mutate(
    ancestor_long = primate_ace_long$longevity[match(parent, primate_ace_long$node)],
    descendant_long = primate_ace_long$longevity[match(node, primate_ace_long$node)],
    diff_long = descendant_long - ancestor_long
  )

result <- list(
  script_path = script_path,
  versions = list(
    ape = as.character(packageVersion("ape")),
    geiger = as.character(packageVersion("geiger")),
    phytools = as.character(packageVersion("phytools")),
    treeio = as.character(packageVersion("treeio")),
    tidytree = as.character(packageVersion("tidytree")),
    openxlsx = as.character(packageVersion("openxlsx"))
  ),
  data_processing = list(
    raw_row_count = nrow(primate_raw),
    processed_row_count = nrow(reference_primate),
    processed_species_count = length(unique(reference_primate$species)),
    duplicate_species_after_grouping = max(reference_primate %>% count(species) %>% pull(n)),
    checked_in_processed_matches_reference = processed_match,
    checked_in_processed_path = file.path(data_dir, "primate.csv"),
    reference_processed_path = file.path(out_dir, "reference_primate.csv")
  ),
  tree_processing = list(
    original_tip_count = length(original_tree$tip.label),
    trimmed_tip_count = length(trimmed_tree$tip.label),
    missing_tips = missing_tips,
    rooted = is.rooted(trimmed_tree),
    binary = is.binary(trimmed_tree),
    ultrametric = is.ultrametric(trimmed_tree),
    checked_in_trimmed_tip_order_matches_reference = identical(
      trimmed_tree$tip.label,
      checked_in_trimmed_tree$tip.label
    ),
    checked_in_trimmed_tip_set_matches_reference = identical(
      sort(trimmed_tree$tip.label),
      sort(checked_in_trimmed_tree$tip.label)
    )
  ),
  tree_examples = list(
    extract_clade = list(
      source_node_numeric = 77,
      source_node_label = node_label_from_numeric(trimmed_tree, 77),
      tip_count = length(tt77$tip.label),
      taxa = sort(tt77$tip.label)
    ),
    unroot_tree = list(
      rooted = is.rooted(unrooted_tree),
      tip_count = length(unrooted_tree$tip.label)
    ),
    rotate_node = list(
      source_node_numeric = 130,
      source_node_label = node_label_from_numeric(trimmed_tree, 130),
      tip_count = length(rt130$tip.label),
      tip_order = rt130$tip.label,
      same_tip_set = identical(sort(rt130$tip.label), sort(trimmed_tree$tip.label)),
      all_equal_message = paste(all.equal(trimmed_tree, rt130), collapse = " | ")
    ),
    rotate_all = list(
      tip_count = length(rt_all$tip.label),
      tip_order = rt_all$tip.label,
      same_tip_set = identical(sort(rt_all$tip.label), sort(trimmed_tree$tip.label))
    )
  ),
  data_tree_alignment = list(
    aligned_species_equals_tip_order = identical(primate_aligned$species, trimmed_tree$tip.label),
    aligned_species_first_6 = primate_aligned$species[1:6],
    tip_order_first_6 = trimmed_tree$tip.label[1:6],
    nodeid_examples = nodeid_examples,
    joined_tip_count = length(p_tree_data@phylo$tip.label),
    joined_extra_rows = nrow(p_tree_data@extraInfo)
  ),
  random_signal = list(
    seed = 1,
    random_tree_path = file.path(out_dir, "random_tree_seed1.nwk"),
    examples = random_results
  ),
  primate_lambda = list(
    lambda_value = unname(lambda_ll$opt$lambda),
    log_likelihood = unname(lambda_ll$opt$lnL)
  ),
  primate_lambda_zero = list(
    lambda0_log_likelihood = unname(lambda_ll0$opt$lnL),
    likelihood_ratio = unname(ll_diff0),
    p_value = unname(p_value),
    lambda0_vcv_top3 = unname(vcv.phylo(trimmed_tree_lambda0)[1:3, 1:3]),
    real_vcv_top3 = unname(vcv.phylo(trimmed_tree)[1:3, 1:3])
  ),
  ancestral = list(
    nodewise = ancestral_nodes,
    mrca_node = unname(mrca_node),
    mrca_estimate = Filter(function(x) x$node == mrca_node, ancestral_nodes),
    increase_count = unname(sum(tree_df$diff_long > 0, na.rm = TRUE)),
    increase_gt12_count = unname(sum(tree_df$diff_long > 12, na.rm = TRUE))
  )
)

writeLines(
  toJSON(result, auto_unbox = TRUE, pretty = TRUE, digits = 16),
  file.path(out_dir, "r_reference_results.json")
)
