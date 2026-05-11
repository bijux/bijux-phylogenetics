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
  study_id = "primate-pgls-and-signal",
  evidence_id = "evidence-005",
  evidence_title = "Primate residual diagnostics parity bundle",
  comparison_mode = "direct_parity",
  execution_mode = "study_reference_wrapper",
  source_intake_policy = "read-only-external-source",
  source_basis_locators = c(
    "evidence-book/studies/primate-pgls-and-signal/provenance/lund-course-sources.json",
    "evidence-book/studies/primate-longevity-signal/datasets/reference_primate.csv",
    "evidence-book/studies/primate-longevity-signal/datasets/reference_trimmed_primatetree.nwk",
    "evidence-book/studies/primate-pgls-and-signal/reference/reference_results.json",
    "evidence-book/studies/primate-pgls-and-signal/evidence-005/residual-diagnostics-parity.json",
  ),
  reference_scripts = c(
    "evidence-book/studies/primate-pgls-and-signal/reference/primate_pgls_and_signal_reference_r.R",
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
