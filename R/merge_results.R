argv <- commandArgs(trailingOnly = TRUE)
out_flag <- which(argv == "--out")
out_dir <- if (length(out_flag)) argv[out_flag + 1] else "results_combined"
dirs <- if (length(out_flag)) argv[-c(out_flag, out_flag + 1)] else argv
if (length(dirs) < 2) stop("usage: merge_results.R dir1 dir2 [dir3 ...] [--out combined_dir]")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

merge_csv <- function(filename) {
  parts <- lapply(dirs, function(d) {
    f <- file.path(d, filename)
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(invisible(NULL))
  combined <- do.call(rbind, parts)
  if ("model" %in% names(combined)) combined <- combined[!duplicated(combined[setdiff(names(combined), character(0))]), ]
  write.csv(combined, file.path(out_dir, filename), row.names = FALSE)
  cat(sprintf("%s: %d rows from %d source(s), models: %s\n", filename, nrow(combined), length(parts),
              paste(unique(combined$model), collapse = ", ")))
}

for (f in c("metrics_overall.csv", "metrics_by_horizon.csv", "runtime.csv")) merge_csv(f)

for (d in dirs) {
  cals <- list.files(d, pattern = "^calibration_.*\\.csv$", full.names = TRUE)
  file.copy(cals, out_dir, overwrite = TRUE)
}
cat(sprintf("\ncombined results written to %s/\n", out_dir))