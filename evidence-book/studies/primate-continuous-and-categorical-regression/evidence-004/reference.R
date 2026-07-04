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
  study_id = "primate-continuous-and-categorical-regression",
  evidence_id = "evidence-004",
  evidence_title = "Primate intercept-free contrast bundle",
  comparison_mode = "direct_parity",
  execution_mode = "study_reference_wrapper",
  source_intake_policy = "repository-owned-source",
  source_basis_locators = c(
    "evidence-book/studies/primate-continuous-and-categorical-regression/provenance/lund-course-sources.json",
    "evidence-book/studies/primate-continuous-and-categorical-regression/datasets/registry.json",
    "evidence-book/studies/primate-continuous-and-categorical-regression/reference/README.md",
    "packages/bijux-phylogenetics/tests/test_pgls_interaction_coefficients.py",
    "packages/bijux-phylogenetics/tests/test_regression_model_selection.py",
  ),
  reference_scripts = c(
    "evidence-book/studies/primate-continuous-and-categorical-regression/reference/lund_pcm3_contcat_source.R",
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
