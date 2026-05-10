#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
bundle_root <- dirname(normalizePath(script_path, mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else tempdir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

payload <- list(
  study_id = "comparative-trust-boundaries",
  evidence_id = "evidence-003",
  evidence_title = "OU identifiability warning bundle",
  comparison_mode = "bijux_native_reinterpretation",
  execution_mode = "not_applicable",
  source_intake_policy = "repository-owned-source",
  source_basis_locators = c(
    "evidence-book/studies/comparative-trust-boundaries/provenance/runtime-sources.json",
    "packages/bijux-phylogenetics/tests/fixtures",
    "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative/models.py",
  ),
  reference_scripts = c(
  )
)

write_json(payload, file.path(out_dir, "reference-contract.json"), auto_unbox = TRUE, pretty = TRUE)
cat(toJSON(payload, auto_unbox = TRUE, pretty = TRUE))
cat("\n")

# This evidence unit is a Bijux-native or fixture-backed surface with no
# canonical external R comparison program. The wrapper therefore records
# the bundle-local contract instead of pretending an R parity run exists.
