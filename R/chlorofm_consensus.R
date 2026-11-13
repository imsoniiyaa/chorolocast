# ChloroFM seasonal consensus, driven from R.
# R owns data preparation, local analogs, seasonal decomposition, scoring, and output.
# reticulate is used only for the official TimesFM neural-network inference backend.

`%||%` <- function(x, y) if (is.null(x)) y else x

.source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.source_dir <- if (!is.null(.source_file)) {
  dirname(normalizePath(.source_file))
} else if (length(.script_arg)) {
  dirname(normalizePath(sub("^--file=", "", .script_arg[1])))
} else {
  normalizePath("R")
}
source(file.path(.source_dir, "fm.R"))

ENV_COLUMNS <- c(
  "lake_Temp_C_mean", "lake_SpCond_uScm_mean", "lake_fDOM_QSU_mean",
  "lake_DOsat_percent_mean", "inflow_Flow_cms_mean", "inflow_Temp_C_mean",
  "met_PAR_umolm2s_mean", "met_AirTemp_C_mean", "met_RH_percent_mean",
  "met_WindSpeed_ms_mean", "met_Rain_mm_sum"
)

script_dir <- function() {
  .source_dir
}

locf <- function(x) {
  for (i in seq_along(x)) if (is.na(x[i]) && i > 1L) x[i] <- x[i - 1L]
  x
}

fill_both <- function(x) {
  x <- locf(x)
  if (any(!is.na(x))) x[seq_len(min(which(!is.na(x))) - 1L)] <- x[min(which(!is.na(x)))]
  x
}

roll_mean <- function(x, width) {
  vapply(seq_along(x), function(i) mean(x[max(1L, i - width + 1L):i], na.rm = TRUE), numeric(1))
}

load_daily_frame <- function(path, date_col = "date") {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!date_col %in% names(raw)) stop("missing date column: ", date_col)
  raw[[date_col]] <- as.Date(substr(as.character(raw[[date_col]]), 1, 10))
  need <- c(VARIABLE, ENV_COLUMNS)
  missing <- setdiff(need, names(raw))
  if (length(missing)) stop("missing required columns: ", paste(missing, collapse = ", "))
  days <- seq(min(raw[[date_col]]), max(raw[[date_col]]), by = "day")
  out <- data.frame(date = days)
  for (nm in need) {
    values <- tapply(as.numeric(raw[[nm]]), raw[[date_col]], function(v) {
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    })
    out[[nm]] <- as.numeric(values[as.character(days)])
  }
  out
}

make_features_r <- function(d) {
  ch <- locf(d[[VARIABLE]])
  doy <- as.integer(format(d$date, "%j"))
  f <- data.frame(sin = sin(2 * pi * doy / 365.25), cos = cos(2 * pi * doy / 365.25))
  for (w in c(1L, 7L, 14L, 30L, 60L)) f[[paste0("ch_mean", w)]] <- roll_mean(ch, w)
  f$ch_last <- ch
  f$ch_slope7 <- c(rep(NA_real_, 7L), ch[8:length(ch)] - ch[1:(length(ch) - 7L)])
  for (nm in ENV_COLUMNS) f[[nm]] <- roll_mean(locf(d[[nm]]), 7L)
  f
}

feature_columns_r <- function(group, all_names) {
  switch(group,
    season = c("sin", "cos"),
    season_chla = c("sin", "cos", "ch_last", "ch_mean7", "ch_mean30", "ch_slope7"),
    season_env = all_names,
    stop("unknown analog group: ", group)
  )
}

analog_samples_r <- function(d, f, origin, group, neighbors, recency,
                             members = 199L, horizon = HORIZON) {
  pos <- match(as.Date(origin), d$date)
  if (is.na(pos)) stop("origin outside data: ", origin)
  candidates <- seq.int(61L, pos - horizon)
  if (length(candidates) < 2L) stop("not enough historical analogs before ", origin)
  cols <- feature_columns_r(group, names(f))
  train <- as.matrix(f[candidates, cols, drop = FALSE])
  x0 <- as.numeric(f[pos, cols, drop = FALSE])
  scale <- apply(train, 2, sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  delta <- sweep(sweep(train, 2, x0, "-"), 2, scale, "/")
  distance <- sqrt(rowMeans(delta^2, na.rm = TRUE))
  distance[!is.finite(distance)] <- Inf
  if (recency != 0) {
    age_years <- as.numeric(d$date[pos] - d$date[candidates]) / 365.25
    distance <- distance * (1 + recency * pmax(age_years - 1, 0))
  }
  nearest <- candidates[head(order(distance), min(as.integer(neighbors), length(distance)))]
  target <- fill_both(d[[VARIABLE]])
  trajectories <- t(vapply(nearest, function(i) target[i + seq_len(horizon)], numeric(horizon)))
  pick <- floor((0:(members - 1L)) * nrow(trajectories) / members) + 1L
  trajectories[pick, , drop = FALSE]
}

foundation_samples_r <- function(backend, series, origins, members = 199L, context = 1024L) {
  contexts <- lapply(origins, function(o) {
    x <- make_context(series, as.Date(o), context, TRUE)
    if (is.null(x)) stop("unusable TimesFM context at ", o)
    x
  })
  q <- backend$predict(contexts, HORIZON)
  out <- array(NA_real_, c(length(origins), members, HORIZON))
  for (i in seq_along(origins)) {
    s <- quantiles_to_samples(matrix(q[i, , ], nrow = HORIZON), members)
    out[i, , ] <- inv(s, TRUE)
  }
  pmax(out, 0)
}

season_design_r <- function(dates) {
  day <- as.integer(format(as.Date(dates), "%j")) / 365.25
  cbind(1, sin(2*pi*day), cos(2*pi*day), sin(4*pi*day), cos(4*pi*day),
        sin(6*pi*day), cos(6*pi*day))
}

seasonal_foundation_r <- function(backend, series, origins, context = 1024L,
                                  strength = 0.5, members = 199L) {
  contexts <- vector("list", length(origins)); seasons <- vector("list", length(origins))
  for (i in seq_along(origins)) {
    pos <- match(as.Date(origins[i]), series$date)
    idx <- max(1L, pos - 1460L):pos
    valid <- !is.na(series$y[idx])
    X <- season_design_r(series$date[idx][valid])
    y <- log1p(pmax(series$y[idx][valid], 0))
    penalty <- diag(10, ncol(X)); penalty[1, 1] <- 1e-6
    coef <- solve(crossprod(X) + penalty, crossprod(X, y))
    ctx_idx <- tail(idx, context)
    ctx <- approx(seq_along(ctx_idx), series$y[ctx_idx], xout = seq_along(ctx_idx), rule = 2)$y
    contexts[[i]] <- as.numeric(log1p(pmax(ctx, 0)) - strength * season_design_r(series$date[ctx_idx]) %*% coef)
    future <- as.Date(origins[i]) + seq_len(HORIZON)
    seasons[[i]] <- as.numeric(strength * season_design_r(future) %*% coef)
  }
  q <- backend$predict(contexts, HORIZON)
  out <- array(NA_real_, c(length(origins), members, HORIZON))
  for (i in seq_along(origins)) {
    s <- quantiles_to_samples(matrix(q[i, , ], nrow = HORIZON), members)
    out[i, , ] <- pmax(expm1(sweep(s, 2, seasons[[i]], "+")), 0)
  }
  out
}

hybrid_samples_r <- function(d, origins, foundation, config, members = 199L) {
  f <- make_features_r(d)
  out <- array(NA_real_, dim(foundation)); cache <- new.env(parent = emptyenv())
  for (band in config$bands) {
    key <- paste(band$analog, band$neighbors, band$recency, sep = "_")
    if (!exists(key, cache, inherits = FALSE)) {
      a <- array(NA_real_, c(length(origins), members, HORIZON))
      for (i in seq_along(origins))
        a[i, , ] <- analog_samples_r(d, f, origins[i], band$analog,
                                     band$neighbors, band$recency, members)
      assign(key, a, cache)
    }
    analog <- get(key, cache, inherits = FALSE)
    horizons <- seq.int(band$start, band$end); w <- band$timesfm_weight
    for (i in seq_along(origins)) for (h in horizons)
      out[i, , h] <- w * sort(foundation[i, , h]) + (1 - w) * sort(analog[i, , h])
  }
  pmax(out, 0)
}

predict_consensus_r <- function(d, origins, config, members = 199L, include_baseline = FALSE) {
  series <- data.frame(date = d$date, y = d[[VARIABLE]])
  backend <- backend_python("timesfm", 1024L)
  foundation <- foundation_samples_r(backend, series, origins, members, 1024L)
  original <- hybrid_samples_r(d, origins, foundation, config, members)
  seasonal <- seasonal_foundation_r(backend, series, origins, 1024L, 0.5, members)
  result <- list(chlorofm_consensus = 0.5 * original + 0.5 * seasonal,
                 chlorofm_v1 = original, timesfm_log1024 = foundation)
  if (include_baseline)
    result$timesfm_log512 <- foundation_samples_r(backend, series, origins, members, 512L)
  result
}

outcomes_r <- function(d, origins) {
  y <- matrix(NA_real_, length(origins), HORIZON)
  for (i in seq_along(origins)) y[i, ] <- d[[VARIABLE]][match(as.Date(origins[i]), d$date) + seq_len(HORIZON)]
  y
}

scores_r <- function(samples, y, origins, model) {
  rows <- list(); k <- 1L
  for (i in seq_along(origins)) {
    cr <- crps_ens(samples[i, , ], ifelse(is.na(y[i, ]), 0, y[i, ]))
    ci <- apply(samples[i, , ], 2, quantile, probs = c(.025, .975), names = FALSE)
    for (h in which(!is.na(y[i, ]))) {
      rows[[k]] <- data.frame(model = model, origin = as.character(origins[i]), h = h,
        CRPS = cr[h], coverage95 = as.numeric(ci[1,h] <= y[i,h] && y[i,h] <= ci[2,h]),
        width95 = ci[2,h] - ci[1,h]); k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

read_consensus_config <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("install.packages('jsonlite')")
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

cmd_consensus_evaluate <- function(a) {
  d <- load_daily_frame(a$data); config <- read_consensus_config(a$config)
  origins <- seq(as.Date(a$start), as.Date(a$end), by = paste(a$stride, "days"))
  if (tail(origins, 1) + HORIZON > max(d$date)) stop("evaluation end must have a complete 35-day outcome")
  models <- predict_consensus_r(d, origins, config, a$members, TRUE); y <- outcomes_r(d, origins)
  detail <- do.call(rbind, lapply(names(models), function(nm) scores_r(models[[nm]], y, origins, nm)))
  overall <- do.call(rbind, lapply(split(detail, detail$model), function(x) data.frame(
    model = x$model[1], CRPS = mean(x$CRPS), n = nrow(x), coverage95 = mean(x$coverage95),
    width95 = mean(x$width95))))
  by_h <- do.call(rbind, lapply(split(detail, list(detail$model, detail$h), drop = TRUE), function(x)
    data.frame(model = x$model[1], h = x$h[1], CRPS = mean(x$CRPS), n = nrow(x))))
  dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
  write.csv(detail, file.path(a$out, "scores.csv"), row.names = FALSE)
  write.csv(overall, file.path(a$out, "metrics_overall.csv"), row.names = FALSE)
  write.csv(by_h, file.path(a$out, "metrics_by_horizon.csv"), row.names = FALSE)
  print(overall, row.names = FALSE)
}

cmd_consensus_forecast <- function(a) {
  d <- load_daily_frame(a$data); config <- read_consensus_config(a$config)
  ref <- if (!is.null(a$reference_date)) as.Date(a$reference_date) else max(d$date[!is.na(d[[VARIABLE]])])
  pos <- match(ref, d$date)
  if (is.na(pos) || is.na(d[[VARIABLE]][pos])) stop("reference date must have an observed target")
  models <- predict_consensus_r(d, ref, config, a$members, FALSE)
  df <- to_vera(models$chlorofm_consensus[1, , ], ref,
                "chlorocast_consensus_r_experimental", a$site, a$depth, VARIABLE)
  dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(a$out, sprintf("daily-%s-chlorocast_consensus_r_experimental.csv.gz", ref))
  con <- gzfile(path, "wt"); on.exit(close(con), add = TRUE)
  write.csv(df, con, row.names = FALSE, quote = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", path, nrow(df)))
}

# Friendly RStudio entry points. Run setwd() at the repository root first.
chlorofm_forecast <- function(data = "chla_daily.csv", reference_date = NULL,
                              out = "forecasts", members = 199L,
                              config = file.path(dirname(script_dir()), "config", "chlorofm.json"),
                              site = "fcre", depth = 1.6) {
  cmd_consensus_forecast(list(data = data, reference_date = reference_date,
    out = out, members = as.integer(members), config = config, site = site, depth = depth))
}

chlorofm_evaluate <- function(data = "chla_daily.csv", start = "2026-01-01",
                              end = "2026-08-16", stride = 7L,
                              out = "results_chlorofm_consensus_r", members = 199L,
                              config = file.path(dirname(script_dir()), "config", "chlorofm.json")) {
  cmd_consensus_evaluate(list(data = data, start = start, end = end,
    stride = as.integer(stride), out = out, members = as.integer(members), config = config))
}

main_consensus <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!length(argv) || !argv[1] %in% c("evaluate", "forecast"))
    stop("usage: Rscript R/chlorofm_consensus.R evaluate|forecast --data chla_daily.csv")
  p <- parse_args(argv); default_config <- file.path(dirname(script_dir()), "config", "chlorofm.json")
  common <- list(config = default_config, members = 199, site = "fcre", depth = 1.6)
  if (p$cmd == "evaluate") {
    a <- with_defaults(p$a, c(common, list(data = "chla_daily.csv", start = "2026-01-01",
      end = "2026-08-16", stride = 7, out = "results_chlorofm_consensus_r")))
    cmd_consensus_evaluate(a)
  } else {
    a <- with_defaults(p$a, c(common, list(data = "chla_daily.csv", out = "forecasts")))
    cmd_consensus_forecast(a)
  }
}

if (!interactive() && sys.nframe() == 0) main_consensus()
