source("R/utils/ensemble.R")
source("R/utils/python_backend.R")
source("R/utils/data.R")
source("R/utils/calibration.R")
source("R/utils/metrics.R")
source("R/utils/forecast.R")
source("R/models/backend_ensemble.R")
source("R/models/persistence.R")
source("R/models/arima.R")
source("R/models/rf.R") # backend_lagmodel/lag helpers live here
source("R/models/rf_direct.R")
source("R/models/xgboost.R") # depends on rf.R being sourced first
source("R/models/chronos.R")
source("R/models/timesfm.R")
source("R/models/rf_tfm_resid.R")

HORIZON <- 35L
LEVELS <- seq(0.1, 0.9, by = 0.1)
Z90 <- qnorm(0.9)
VARIABLE <- "Chla_ugL_mean"


get_backend <- function(name, ctx_len) {
  switch(name,
    persistence = backend_persistence(),
    arima = backend_arima(),
    rf = backend_rf(),
    rf_direct = backend_rf_direct(),
    xgboost = backend_xgboost(),
    chronos = backend_chronos(ctx_len),
    timesfm = backend_timesfm(ctx_len),
    rf_tfm_resid = backend_rf_tfm_resid(),
    ens = {
      w <- as.numeric(Sys.getenv("ENS_W", "0.7"))
      backend_ensemble(backend_timesfm(ctx_len), backend_rf(), w)
    },
    ens_direct = {
      w <- as.numeric(Sys.getenv("ENS_W", "0.7"))
      backend_ensemble(backend_timesfm(ctx_len), backend_rf_direct(), w)
    },
    ens_h = {
      w0 <- as.numeric(Sys.getenv("ENS_W0", "1"))
      w1 <- as.numeric(Sys.getenv("ENS_W1", "0.75"))
      backend_ensemble_h(backend_arima(), backend_rf(), w0, w1)
    },
    ens3 = {
      backend_ens3(backend_arima(), backend_rf(), backend_timesfm(ctx_len),
        v0 = as.numeric(Sys.getenv("ENS_V0", "0.75")),
        v1 = as.numeric(Sys.getenv("ENS_V1", "0.5"))
      )
    },
    stop(paste("unknown model:", name))
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
  set.seed(as.integer(Sys.getenv("SEED", "42")))
  if (!length(argv) || !argv[1] %in% c("evaluate", "forecast")) stop("usage: chlorocast_fm.R evaluate|forecast [--options]")
  p <- parse_args(argv)
  common <- list(target = VARIABLE, site = "fcre", depth = 1.6, ctx_len = 512, members = 199)
  if (p$cmd == "evaluate") {
    a <- with_defaults(p$a, c(common, list(
      models = c(
        "persistence", "arima", "rf",
        "rf_direct", "xgboost",
        "chronos", "timesfm"
      ),
      stride = 3,
      val_start = "2025-01-01",
      val_end = "2025-12-31",
      test_start = "2026-01-01",
      test_end = "2026-09-20",
      out = "results"
    )))

    cmd_evaluate(a)
  } else {
    a <- with_defaults(
      p$a,
      c(common, list(
        max_stale = 7,
        out = "forecasts"
      ))
    )

    if (is.null(a$model)) {
      stop("--model is required")
    }

    a$model <- match.arg(
      a$model,
      c(
        "persistence",
        "arima",
        "rf",
        "rf_direct",
        "xgboost",
        "chronos",
        "timesfm",
        "ens",
        "ens_direct",
        "ens_h",
        "ens3"
      )
    )
  }

    cmd_forecast(a)
  }

if (!interactive() && sys.nframe() == 0) main() 
