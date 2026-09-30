source("R/utils/python_backend.R")
source("R/utils/data.R")
source("R/utils/calibration.R")
source("R/utils/metrics.R")
source("R/utils/forecast.R")

source("R/models/persistence.R")
source("R/models/arima.R")
source("R/models/rf.R") # backend_lagmodel/lag helpers live here
source("R/models/xgboost.R") # depends on rf.R being sourced first
source("R/models/chronos.R")
source("R/models/timesfm.R")

HORIZON <- 35L
LEVELS <- seq(0.1, 0.9, by = 0.1)
Z90 <- qnorm(0.9)
VARIABLE <- "Chla_ugL_mean"


get_backend <- function(name, ctx_len) {
  switch(name,
    persistence = backend_persistence(),
    arima       = backend_arima(),
    rf          = backend_rf(),
    xgboost     = backend_xgboost(),
    chronos     = backend_chronos(ctx_len),
    timesfm     = backend_timesfm(ctx_len)
  )
}

parse_args <- function(argv) {
  cmd <- argv[1]
  rest <- argv[-1]
  a <- list()
  i <- 1
  while (i <= length(rest)) {
    k <- gsub("-", "_", sub("^--", "", rest[i]))
    j <- i + 1
    vals <- character(0)
    while (j <= length(rest) && !grepl("^--", rest[j])) {
      vals <- c(vals, rest[j])
      j <- j + 1
    }
    a[[k]] <- if (length(vals)) vals else TRUE
    i <- j
  }
  list(cmd = cmd, a = a)
}

with_defaults <- function(a, defaults) {
  for (k in names(defaults)) if (is.null(a[[k]])) a[[k]] <- defaults[[k]]
  for (k in intersect(names(a), c("depth", "ctx_len", "members", "stride", "max_stale"))) a[[k]] <- as.numeric(a[[k]])
  a
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!length(argv) || !argv[1] %in% c("evaluate", "forecast")) stop("usage: chlorocast_fm.R evaluate|forecast [--options]")
  p <- parse_args(argv)
  common <- list(target = VARIABLE, site = "fcre", depth = 1.6, ctx_len = 512, members = 199)
  if (p$cmd == "evaluate") {
    a <- with_defaults(p$a, c(common, list(
      models = c("persistence", "arima", "rf", "xgboost", "chronos", "timesfm"), stride = 3,
      val_start = "2025-01-01", val_end = "2025-12-31", test_start = "2026-01-01",
      test_end = "2026-09-20", out = "results"
    )))
    cmd_evaluate(a)
  } else {
    a <- with_defaults(p$a, c(common, list(max_stale = 7, out = "forecasts")))
    if (is.null(a$model)) stop("--model is required")
    a$model <- match.arg(a$model, c("persistence", "arima", "rf", "xgboost", "chronos", "timesfm"))
    cmd_forecast(a)
  }
}

if (!interactive() && sys.nframe() == 0) main()
