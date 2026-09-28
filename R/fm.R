# set the value
HORIZON  <- 35L
LEVELS   <- seq(0.1, 0.9, by = 0.1)      
Z90      <- qnorm(0.9)
BLOOM    <- 20                          
VARIABLE <- "Chla_ugL_mean"

# find the path to the python backend
PY_DIR <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
  cands <- c(Sys.getenv("CHLOROCAST_PY_DIR", unset = NA), file.path(here, "..", "py"),
             file.path(getwd(), "py"), file.path(getwd(), "..", "py"))
  cands <- cands[!is.na(cands) & file.exists(file.path(cands, "fm_backend.py"))]
  if (length(cands)) normalizePath(cands[1]) else NA_character_
})
# finish a time series by filling in missing days and averaging duplicates
finish_series <- function(d, y) {
  m <- tapply(y, d, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
  days <- seq(min(as.Date(names(m))), max(as.Date(names(m))), by = "day")
  data.frame(date = days, y = as.numeric(m[as.character(days)]))
}
# load a time series from a CSV file or from VERA
load_series <- function(path, date_col = NULL, target = VARIABLE) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(date_col))
    date_col <- intersect(c("Date", "date", "datetime", "DateTime", "time", "timestamp"), names(df))[1]
  if (is.na(date_col)) stop("no date column found in: ", paste(names(df), collapse = ", "), " (use --date-col)")
  finish_series(as.Date(substr(as.character(df[[date_col]]), 1, 10)), as.numeric(df[[target]]))
}

# converting to and from log1p space, with a floor at 0
fwd <- function(x, log) { x <- pmax(x, 0); if (log) log1p(x) else x }
inv <- function(x, log) pmax(if (log) expm1(x) else x, 0)

# make a context vector for a given origin date, or return NULL if not enough data
make_context <- function(series, origin, ctx_len, log, max_stale = 7, min_obs = 30) {
  i <- match(origin, series$date)
  if (is.na(i)) return(NULL)
  s <- series$y[max(1, i - ctx_len + 1):i]
  if (sum(!is.na(s)) < min_obs) return(NULL)
  if (length(s) - max(which(!is.na(s))) > max_stale) return(NULL)
  s <- approx(seq_along(s), s, xout = seq_along(s), rule = 2)$y
  fwd(s, log)
}
# persistence model: predict that the future will be the same as the last observed value, 
# plus a quantile of the recent changes
backend_persistence <- function() {
  list(name = "persistence", predict = function(ctxs, horizon) {
    out <- array(0, c(length(ctxs), horizon, 9))
    for (i in seq_along(ctxs)) {
      x <- ctxs[[i]]; n <- length(x)
      for (h in seq_len(horizon)) {
        d <- if (n - h >= 5) x[(h + 1):n] - x[1:(n - h)] else NULL
        out[i, h, ] <- x[n] + if (is.null(d)) 0 else quantile(d, LEVELS, names = FALSE)
      }
    }
    out
  })
}
# Python backend: load the model and call its predict function
backend_python <- function(model, ctx_len) {
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("install.packages('reticulate')")
  if (is.na(PY_DIR)) stop("cannot find py/fm_backend.py; set CHLOROCAST_PY_DIR")
  m <- reticulate::import_from_path("fm_backend", path = PY_DIR)
  b <- m$load_backend(model, as.integer(ctx_len))
  list(name = model, predict = function(ctxs, horizon) {
    flat <- m$predict_flat(b, lapply(ctxs, as.numeric), as.integer(horizon))   
    aperm(array(as.numeric(unlist(flat)), dim = c(9, horizon, length(ctxs))), c(3, 2, 1))
  })
}

LAGS <- c(1, 2, 3, 5, 7, 10, 14, 21, 28, 35, 49, 63, 90)
ROLL <- c(7, 14, 30)

# build a lag table for a given time series, with specified lags and rolling windows
build_lag_table <- function(x, lags = LAGS, roll = ROLL) {
  n <- length(x); maxlag <- max(lags, roll)
  idx <- (maxlag + 1):(n - 1)
  if (length(idx) < 60) return(NULL)
  feats <- lapply(lags, function(L) x[idx - L])
  names(feats) <- paste0("lag", lags)
  for (w in roll) {
    feats[[paste0("ma", w)]] <- vapply(idx, function(t) mean(x[(t - w + 1):t]), numeric(1))
    feats[[paste0("sd", w)]] <- vapply(idx, function(t) sd(x[(t - w + 1):t]), numeric(1))
  }
  df <- as.data.frame(feats); df$y <- x[idx + 1]
  df
}
# compute the lagged features for a single row (the last row of a time series)
lag_state <- function(x, lags = LAGS, roll = ROLL) {
  n <- length(x)
  row <- as.list(setNames(x[n - lags + 1], paste0("lag", lags)))
  for (w in roll) {
    row[[paste0("ma", w)]] <- mean(x[(n - w + 1):n])
    row[[paste0("sd", w)]] <- sd(x[(n - w + 1):n])
  }
  as.data.frame(row)
}
# sample from a matrix of quantiles, using linear interpolation for the inner quantiles and normal approximation for the tails
sample_row_quantiles <- function(Q, u) {
  n <- nrow(Q); out <- numeric(n)
  inner <- u >= 0.1 & u <= 0.9; lo <- u < 0.1; hi <- u > 0.9
  s_lo <- pmax((Q[, 5] - Q[, 1]) / Z90, 1e-6); s_hi <- pmax((Q[, 9] - Q[, 5]) / Z90, 1e-6)
  if (any(inner)) out[inner] <- vapply(which(inner),
    function(i) approx(LEVELS, Q[i, ], xout = u[i], rule = 2)$y, numeric(1))
  out[lo] <- Q[lo, 1] + s_lo[lo] * (qnorm(u[lo]) - qnorm(0.1))
  out[hi] <- Q[hi, 9] + s_hi[hi] * (qnorm(u[hi]) - qnorm(0.9))
  out
}
# simulate multiple paths from a fitted model, using the lagged features and the step function
simulate_paths <- function(x0, step_fn, horizon, n) {
  mat <- matrix(rep(x0, n), nrow = n, byrow = TRUE)  
  out <- matrix(NA_real_, n, horizon)
  u_grid <- (seq_len(n) - 0.5) / n
  for (h in seq_len(horizon)) {
    state <- do.call(rbind, lapply(seq_len(n), function(i) lag_state(mat[i, ])))
    Q <- step_fn(state)                                 
    out[, h] <- sample_row_quantiles(Q, sample(u_grid))
    mat <- cbind(mat, out[, h])
  }
  out
}
# fit a lagged model and return a function that predicts the next value given the lagged features 
backend_lagmodel <- function(name, fit_and_predict) {
  list(name = name, predict = function(ctxs, horizon) {
    out <- array(NA_real_, c(length(ctxs), horizon, 9))
    for (i in seq_along(ctxs)) {
      x <- ctxs[[i]]
      tab <- build_lag_table(x)
      if (is.null(tab)) { out[i, , ] <- x[length(x)]; next }
      step_fn <- fit_and_predict(tab)
      paths <- simulate_paths(x, step_fn, horizon, 199)
      out[i, , ] <- t(apply(paths, 2, quantile, probs = LEVELS, names = FALSE))
    }
    out
  })
}
# random forest backend: fit a random forest to the lagged features and predict the next value
backend_rf <- function() {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("install.packages('ranger')")
  backend_lagmodel("rf", function(tab) {
    fit <- ranger::ranger(y ~ ., data = tab, quantreg = TRUE, num.trees = 200)
    function(newdata) predict(fit, data = newdata, type = "quantiles", quantiles = LEVELS)$predictions
  })
}
# xgboost backend: fit an xgboost model to the lagged features and predict the next value
backend_xgboost <- function() {
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("install.packages('xgboost')")
  backend_lagmodel("xgboost", function(tab) {
    x <- as.matrix(tab[, setdiff(names(tab), "y")]); y <- tab$y
    fit <- xgboost::xgboost(data = x, label = y, nrounds = 150, max_depth = 4, eta = 0.1,
                            objective = "reg:quantileerror", quantile_alpha = LEVELS, verbose = 0)
    function(newdata) matrix(predict(fit, as.matrix(newdata), reshape = TRUE), ncol = length(LEVELS))
  })
}
# ARIMA backend: fit an ARIMA model to the time series and predict the next values,
# using Fourier terms for seasonality if enough data is available
backend_arima <- function() {
  if (!requireNamespace("forecast", quietly = TRUE)) stop("install.packages('forecast')")
  order_cache <- NULL
  fourier <- function(t, K = 4, period = 365.25) {
    m <- do.call(cbind, lapply(1:K, function(k) cbind(sin(2 * pi * k * t / period), cos(2 * pi * k * t / period))))
    colnames(m) <- paste0(c("S", "C"), rep(1:K, each = 2)); m
  }
  list(name = "arima", predict = function(ctxs, horizon) {
    out <- array(NA_real_, c(length(ctxs), horizon, 9))
    for (i in seq_along(ctxs)) {
      x <- ctxs[[i]]; n <- length(x)
      xreg <- if (n >= 2 * 365) fourier(seq_len(n)) else NULL
      fit <- tryCatch({
        if (is.null(order_cache)) {
          f <- forecast::auto.arima(x, xreg = xreg, seasonal = FALSE, stepwise = TRUE, approximation = TRUE)
          order_cache <<- as.integer(forecast::arimaorder(f)[1:3])
          f
        } else {
          forecast::Arima(x, order = order_cache, xreg = xreg)
        }
      }, error = function(e) forecast::auto.arima(x, seasonal = FALSE))
      newxreg <- if (!is.null(xreg)) fourier(n + seq_len(horizon)) else NULL
      fc <- forecast::forecast(fit, h = horizon, xreg = newxreg, level = 95)
      se <- pmax((fc$upper[, 1] - fc$mean) / qnorm(0.975), 1e-6)
      for (h in seq_len(horizon)) out[i, h, ] <- as.numeric(fc$mean[h]) + se[h] * qnorm(LEVELS)
    }
    out
  })
}
# get the backend for a given model name and context length
get_backend <- function(name, ctx_len) {
  switch(name,
    persistence = backend_persistence(),
    arima       = backend_arima(),
    rf          = backend_rf(),
    xgboost     = backend_xgboost(),
    backend_python(name, ctx_len))
}
# convert a matrix of quantiles to samples, 
# using linear interpolation for the inner quantiles and normal approximation for the tails
quantiles_to_samples <- function(q, n) {
  q <- matrix(q, ncol = 9)
  q <- t(apply(q, 1, sort))
  H <- nrow(q); u <- (seq_len(n) - 0.5) / n
  inner <- u >= 0.1 & u <= 0.9; lo <- u < 0.1; hi <- u > 0.9
  s_lo <- pmax((q[, 5] - q[, 1]) / Z90, 1e-6)
  s_hi <- pmax((q[, 9] - q[, 5]) / Z90, 1e-6)
  out <- matrix(0, n, H)
  for (h in seq_len(H)) {
    out[inner, h] <- approx(LEVELS, q[h, ], xout = u[inner], rule = 2)$y
    out[lo, h]    <- q[h, 1] + s_lo[h] * (qnorm(u[lo]) - qnorm(0.1))
    out[hi, h]    <- q[h, 9] + s_hi[h] * (qnorm(u[hi]) - qnorm(0.9))
  }
  out
}
# compute the CRPS for an ensemble of samples and a vector of observations
apply_scale <- function(samples, scale) {
  med <- apply(samples, 2, median)
  scale <- scale[seq_len(ncol(samples))]
  sweep(sweep(sweep(samples, 2, med, "-"), 2, scale, "*"), 2, med, "+")
}

fit_scale <- function(Q, Y, n = 199, cover = 0.95, smooth = 5) {
  O <- dim(Q)[1]; H <- dim(Q)[2]
  scores <- matrix(NA_real_, O, H)
  for (o in seq_len(O)) {
    s <- quantiles_to_samples(matrix(Q[o, , ], nrow = H), n)
    med <- apply(s, 2, median)
    ci <- apply(s, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
    scores[o, ] <- abs(Y[o, ] - med) / pmax((ci[2, ] - ci[1, ]) / 2, 1e-6)
  }
  cs <- vapply(seq_len(H), function(h) {
    v <- scores[, h]; v <- v[!is.na(v)]
    if (length(v) >= 10) unname(quantile(v, cover)) else NA_real_
  }, numeric(1))
  ok <- !is.na(cs)
  if (!any(ok)) return(rep(1, H))
  if (sum(ok) == 1) cs <- rep(cs[ok], H) else cs <- approx(which(ok), cs[ok], xout = seq_len(H), rule = 2)$y
  half <- smooth %/% 2
  cs <- vapply(seq_len(H), function(h) median(cs[max(1, h - half):min(H, h + half)]), numeric(1))
  pmin(pmax(cs, 0.25), 6)
}

# compute the CRPS for an ensemble of samples and a vector of observations
crps_ens <- function(S, y) {
  n <- nrow(S); Ss <- matrix(apply(S, 2, sort), nrow = n); i <- seq_len(n)
  colMeans(abs(sweep(Ss, 2, y, "-"))) - colSums((2 * i - n - 1) * Ss) / n^2
}

# compute the forecast scores for a set of quantiles and observations, optionally applying a scale
score_forecasts <- function(Q, Y, log, scale = NULL, n = 199) {
  O <- dim(Q)[1]; H <- dim(Q)[2]
  rows <- vector("list", O)
  for (o in seq_len(O)) {
    s <- quantiles_to_samples(matrix(Q[o, , ], nrow = H), n)
    if (!is.null(scale)) s <- apply_scale(s, scale)
    S <- inv(s, log); y <- Y[o, ]
    ci <- apply(S, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
    med <- apply(S, 2, median)
    hh <- which(!is.na(y))
    cr <- crps_ens(S, ifelse(is.na(y), 0, y))
    rows[[o]] <- data.frame(origin = o, h = hh, y = y[hh], med = med[hh], crps = cr[hh],
                            cover95 = as.numeric(ci[1, hh] <= y[hh] & y[hh] <= ci[2, hh]),
                            width95 = ci[2, hh] - ci[1, hh])
  }
  d <- do.call(rbind, rows)
  d$abs_err <- abs(d$y - d$med); d$sq_err <- (d$y - d$med)^2; d$bloom <- d$y >= BLOOM
  d
}
# additional summary functions for the forecast scores, by horizon or overall
summarise_h <- function(d) {
  do.call(rbind, lapply(split(d, d$h), function(x)
    data.frame(h = x$h[1], n = nrow(x), MAE = mean(x$abs_err), RMSE = sqrt(mean(x$sq_err)),
               CRPS = mean(x$crps), coverage95 = mean(x$cover95), width95 = mean(x$width95))))
}
# summarize the forecast scores overall, including separate metrics for bloom events
summarise_all <- function(d) {
  b <- d[d$bloom, ]
  data.frame(n = nrow(d), MAE = mean(d$abs_err), RMSE = sqrt(mean(d$sq_err)), CRPS = mean(d$crps),
             coverage95 = mean(d$cover95), width95 = mean(d$width95), bloom_n = nrow(b),
             bloom_MAE = if (nrow(b)) mean(b$abs_err) else NA_real_,
             bloom_CRPS = if (nrow(b)) mean(b$crps) else NA_real_)
}

# collect the contexts and actual values for a set of forecast origins, using a given backend
collect <- function(series, backend, origins, ctx_len, log, horizon = HORIZON) {
  ctxs <- list(); keep <- integer(0)
  for (k in seq_along(origins)) {
    ctx <- make_context(series, origins[k], ctx_len, log)
    if (!is.null(ctx)) { ctxs[[length(ctxs) + 1]] <- ctx; keep <- c(keep, k) }
  }
  if (!length(ctxs)) stop("no usable forecast origins in this window")
  secs <- system.time(Q <- backend$predict(ctxs, horizon))[["elapsed"]]
  Y <- matrix(NA_real_, length(ctxs), horizon)
  for (j in seq_along(ctxs)) Y[j, ] <- series$y[match(origins[keep[j]], series$date) + seq_len(horizon)]
  good <- rowSums(!is.na(Y)) > 0
  list(Q = Q[good, , , drop = FALSE], Y = Y[good, , drop = FALSE],
       per_fc = secs / length(ctxs), n = sum(good))
}

#evaluate the models on the validation and test sets, saving the results to CSV files
cmd_evaluate <- function(a) {
  series <- if (!is.null(a$data)) load_series(a$data, a$date_col, a$target) else load_vera(a$site, a$depth)
  dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
  log <- !isTRUE(a$no_log)
  by <- paste(a$stride, "days")
  splits <- list(val  = seq(as.Date(a$val_start),  as.Date(a$val_end),  by = by),
                 test = seq(as.Date(a$test_start), as.Date(a$test_end), by = by))
  by_h <- list(); overall <- list(); runtime <- list()
  for (name in a$models) {
    t0 <- Sys.time(); backend <- get_backend(name, a$ctx_len)
    load_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    res <- lapply(splits, function(origins) collect(series, backend, origins, a$ctx_len, log))
    for (sp in names(res))
      cat(sprintf("[%s] %s: %d origins, %.3fs per 35-day forecast\n", name, sp, res[[sp]]$n, res[[sp]]$per_fc))
    runtime[[name]] <- data.frame(model = name, load_s = round(load_s, 2),
                                  per_35day_forecast_s = round(res$test$per_fc, 3))
    scale <- fit_scale(res$val$Q, fwd(res$val$Y, log), n = a$members)
    write.csv(data.frame(h = seq_along(scale), scale = scale, log1p = log, ctx_len = a$ctx_len),
              file.path(a$out, paste0("calibration_", name, ".csv")), row.names = FALSE)
    for (sp in names(res)) {
      variants <- if (sp == "val") list(raw = NULL) else list(raw = NULL, calibrated = scale)
      for (tag in names(variants)) {
        d <- score_forecasts(res[[sp]]$Q, res[[sp]]$Y, log, variants[[tag]], a$members)
        key <- paste(name, sp, tag)
        by_h[[key]]    <- cbind(model = name, split = sp, variant = tag, summarise_h(d))
        overall[[key]] <- cbind(model = name, split = sp, variant = tag, summarise_all(d))
      }
    }
  }
  by_h <- do.call(rbind, by_h); overall <- do.call(rbind, overall)
  write.csv(by_h,    file.path(a$out, "metrics_by_horizon.csv"), row.names = FALSE)
  write.csv(overall, file.path(a$out, "metrics_overall.csv"),    row.names = FALSE)
  write.csv(do.call(rbind, runtime), file.path(a$out, "runtime.csv"), row.names = FALSE)
  num <- function(x) { i <- vapply(x, is.numeric, logical(1)); x[i] <- round(x[i], 3); x }
  cat("\n=== overall (all horizons pooled) ===\n"); print(num(overall), row.names = FALSE)
  sel <- by_h[by_h$split == "test" & by_h$variant == "calibrated" & by_h$h %in% c(1, 7, 14, 21, 28, 35), ]
  cat("\n=== test, calibrated: selected horizons ===\n"); print(num(sel), row.names = FALSE)
  cat(sprintf("\nsaved to %s/  (target coverage95 = 0.95, runtime budget = 20 min)\n", a$out))
}
# for submission to VERA, convert a matrix of samples to a data frame with the required columns
to_vera <- function(S, ref, model_id, site, depth, variable) {
  n <- nrow(S); H <- ncol(S); fmt <- "%Y-%m-%d %H:%M:%S"
  data.frame(project_id = "vera4cast", model_id = model_id,
             datetime = rep(format(ref + seq_len(H), fmt), times = n),
             reference_datetime = format(ref, fmt), duration = "P1D", site_id = site,
             depth_m = depth, family = "ensemble", parameter = rep(seq_len(n), each = H),
             variable = variable, prediction = as.vector(t(S)), stringsAsFactors = FALSE)
}

# run a forecast for a given model and save the results to a CSV file
cmd_forecast <- function(a) {
  series <- if (!is.null(a$data)) load_series(a$data, a$date_col, a$target) else load_vera(a$site, a$depth)
  log <- !isTRUE(a$no_log)
  ref <- if (!is.null(a$reference_date)) as.Date(a$reference_date) else as.Date(Sys.time(), tz = "UTC")
  ok <- which(series$date <= ref & !is.na(series$y))
  last <- series$date[max(ok)]; lag <- as.integer(ref - last)
  ctx <- make_context(series, last, a$ctx_len, log, max_stale = a$max_stale)
  if (is.null(ctx)) stop(sprintf("no usable context: last observation %s is %d days before %s", last, lag, ref))
  t0 <- Sys.time(); backend <- get_backend(a$model, a$ctx_len)
  load_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  total <- HORIZON + lag
  s <- quantiles_to_samples(matrix(backend$predict(list(ctx), total)[1, , ], nrow = total), a$members)
  if (!is.null(a$calib)) {
    cal <- read.csv(a$calib)
    if (cal$log1p[1] != log) stop("calibration was fit with a different log1p setting")
    sc <- c(cal$scale, rep(tail(cal$scale, 1), max(0, total - nrow(cal))))
    s <- apply_scale(s, sc[seq_len(total)])
  }
  S <- inv(s, log)[, lag + seq_len(HORIZON), drop = FALSE]     # ref+1 ... ref+35
  model_id <- if (!is.null(a$model_id)) a$model_id else paste0("chlorocast_", a$model)
  df <- to_vera(S, ref, model_id, a$site, a$depth, a$target)
  stopifnot(!anyNA(df$prediction), all(df$prediction >= 0))
  dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(a$out, sprintf("daily-%s-%s.csv", ref, model_id))
  write.csv(df, f, row.names = FALSE, quote = FALSE)
  cat(sprintf("wrote %s (%d rows); context ends %s (lag %dd); load %.1fs, forecast %.2fs\n",
              f, nrow(df), last, lag, load_s, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# parse the command line arguments into a list of options
parse_args <- function(argv) {
  cmd <- argv[1]; rest <- argv[-1]; a <- list(); i <- 1
  while (i <= length(rest)) {
    k <- gsub("-", "_", sub("^--", "", rest[i])); j <- i + 1; vals <- character(0)
    while (j <= length(rest) && !grepl("^--", rest[j])) { vals <- c(vals, rest[j]); j <- j + 1 }
    a[[k]] <- if (length(vals)) vals else TRUE
    i <- j
  }
  list(cmd = cmd, a = a)
}

# default values for options, and convert certain options to numeric
with_defaults <- function(a, defaults) {
  for (k in names(defaults)) if (is.null(a[[k]])) a[[k]] <- defaults[[k]]
  for (k in intersect(names(a), c("depth", "ctx_len", "members", "stride", "max_stale"))) a[[k]] <- as.numeric(a[[k]])
  a
}
#main entry point: parse the command line arguments and call the appropriate function
main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (!length(argv) || !argv[1] %in% c("evaluate", "forecast")) stop("usage: chlorocast_fm.R evaluate|forecast [--options]")
  p <- parse_args(argv)
  common <- list(target = VARIABLE, site = "fcre", depth = 1.6, ctx_len = 512, members = 199)
  if (p$cmd == "evaluate") {
    a <- with_defaults(p$a, c(common, list(models = c("persistence", "arima", "rf", "xgboost", "chronos", "timesfm"), stride = 3,
      val_start = "2025-01-01", val_end = "2025-12-31", test_start = "2026-01-01",
      test_end = "2026-09-20", out = "results")))
    cmd_evaluate(a)
  } else {
    a <- with_defaults(p$a, c(common, list(max_stale = 7, out = "forecasts")))
    if (is.null(a$model)) stop("--model is required")
    a$model <- match.arg(a$model, c("persistence", "arima", "rf", "xgboost", "chronos", "timesfm"))
    cmd_forecast(a)
  }
}

if (!interactive() && sys.nframe() == 0) main()