#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
bundle_root <- dirname(normalizePath(script_path, mustWork = TRUE))
results_root <- file.path(bundle_root, "results")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else results_root
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

payload <- list(
  study_id = "cichlid-bayesian-linear-and-phylogenetic-mixed-models",
  evidence_id = "evidence-005",
  evidence_title = "Cichlid phylogenetic fixed-effect posterior bundle",
  comparison_mode = "direct_parity",
  execution_mode = "study_reference_wrapper",
  source_intake_policy = "repository-owned-source",
  source_basis_locators = c(
    "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/provenance/lund-course-sources.json",
    "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/datasets/registry.json",
    "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/reference/README.md",
  ),
  reference_scripts = c(
    "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/reference/lund_pcm4_mixed_model_source.R",
  )
)

write_json(payload, file.path(out_dir, "reference-contract.json"), auto_unbox = TRUE, pretty = TRUE)
cat(toJSON(payload, auto_unbox = TRUE, pretty = TRUE))
cat("\n")

# To run the governed study reference, invoke the study-owned R script
# listed in `reference_scripts` with the source context required by that
# study. This bundle-local wrapper stays focused on one evidence unit and
# records the exact comparison contract without editing the untouched
# Lund materials in place.
