# collect the contexts and actual values for a set of forecast origins, using a given backend
collect <- function(series, backend, origins, ctx_len, log, horizon = HORIZON) {
    ctxs <- list()
    keep <- integer(0)
    for (k in seq_along(origins)) {
        ctx <- make_context(series, origins[k], ctx_len, log)
        if (!is.null(ctx)) {
            ctxs[[length(ctxs) + 1]] <- ctx
            keep <- c(keep, k)
        }
    }
    if (!length(ctxs)) stop("no usable forecast origins in this window")
    secs <- system.time(Q <- backend$predict(ctxs, horizon))[["elapsed"]]
    Y <- matrix(NA_real_, length(ctxs), horizon)
    for (j in seq_along(ctxs)) Y[j, ] <- series$y[match(origins[keep[j]], series$date) + seq_len(horizon)]
    good <- rowSums(!is.na(Y)) > 0
    list(
        Q = Q[good, , , drop = FALSE], Y = Y[good, , drop = FALSE],
        per_fc = secs / length(ctxs), n = sum(good)
    )
}

# evaluate the models on the validation and test sets, saving the results to CSV files
cmd_evaluate <- function(a) {
    series <- if (!is.null(a$data)) load_series(a$data, a$date_col, a$target) else load_vera(a$site, a$depth)
    dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
    log <- !isTRUE(a$no_log)
    by <- paste(a$stride, "days")
    splits <- list(
        val = seq(as.Date(a$val_start), as.Date(a$val_end), by = by),
        test = seq(as.Date(a$test_start), as.Date(a$test_end), by = by)
    )
    by_h <- list()
    overall <- list()
    runtime <- list()
    for (name in a$models) {
        t0 <- Sys.time()
        backend <- get_backend(name, a$ctx_len)
        load_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        res <- lapply(splits, function(origins) collect(series, backend, origins, a$ctx_len, log))
        for (sp in names(res)) {
            cat(sprintf("[%s] %s: %d origins, %.3fs per 35-day forecast\n", name, sp, res[[sp]]$n, res[[sp]]$per_fc))
        }
        runtime[[name]] <- data.frame(
            model = name, load_s = round(load_s, 2),
            per_35day_forecast_s = round(res$test$per_fc, 3)
        )
        scale <- fit_scale(res$val$Q, fwd(res$val$Y, log), n = a$members)
        write.csv(data.frame(h = seq_along(scale), scale = scale, log1p = log, ctx_len = a$ctx_len),
            file.path(a$out, paste0("calibration_", name, ".csv")),
            row.names = FALSE
        )
        for (sp in names(res)) {
            variants <- if (sp == "val") list(raw = NULL) else list(raw = NULL, calibrated = scale)
            for (tag in names(variants)) {
                d <- score_forecasts(res[[sp]]$Q, res[[sp]]$Y, log, variants[[tag]], a$members)
                key <- paste(name, sp, tag)
                by_h[[key]] <- cbind(model = name, split = sp, variant = tag, summarise_h(d))
                overall[[key]] <- cbind(model = name, split = sp, variant = tag, summarise_all(d))
            }
        }
    }
    by_h <- do.call(rbind, by_h)
    overall <- do.call(rbind, overall)
    write.csv(by_h, file.path(a$out, "metrics_by_horizon.csv"), row.names = FALSE)
    write.csv(overall, file.path(a$out, "metrics_overall.csv"), row.names = FALSE)
    write.csv(do.call(rbind, runtime), file.path(a$out, "runtime.csv"), row.names = FALSE)
    num <- function(x) {
        i <- vapply(x, is.numeric, logical(1))
        x[i] <- round(x[i], 3)
        x
    }
    cat("\n=== overall (all horizons pooled) ===\n")
    print(num(overall), row.names = FALSE)
    sel <- by_h[by_h$split == "test" & by_h$variant == "calibrated" & by_h$h %in% c(1, 7, 14, 21, 28, 35), ]
    cat("\n=== test, calibrated: selected horizons ===\n")
    print(num(sel), row.names = FALSE)
    cat(sprintf("\nsaved to %s/  (target coverage95 = 0.95, runtime budget = 20 min)\n", a$out))
}
# for submission to VERA, convert a matrix of samples to a data frame with the required columns
to_vera <- function(S, ref, model_id, site, depth, variable) {
    n <- nrow(S)
    H <- ncol(S)
    fmt <- "%Y-%m-%d %H:%M:%S"
    data.frame(
        project_id = "vera4cast", model_id = model_id,
        datetime = rep(format(ref + seq_len(H), fmt), times = n),
        reference_datetime = format(ref, fmt), duration = "P1D", site_id = site,
        depth_m = depth, family = "ensemble", parameter = rep(seq_len(n), each = H),
        variable = variable, prediction = as.vector(t(S)), stringsAsFactors = FALSE
    )
}

# run a forecast for a given model and save the results to a CSV file
cmd_forecast <- function(a) {
    series <- if (!is.null(a$data)) load_series(a$data, a$date_col, a$target) else load_vera(a$site, a$depth)
    log <- !isTRUE(a$no_log)
    ref <- if (!is.null(a$reference_date)) as.Date(a$reference_date) else as.Date(Sys.time(), tz = "UTC")
    ok <- which(series$date <= ref & !is.na(series$y))
    last <- series$date[max(ok)]
    lag <- as.integer(ref - last)
    ctx <- make_context(series, last, a$ctx_len, log, max_stale = a$max_stale)
    if (is.null(ctx)) stop(sprintf("no usable context: last observation %s is %d days before %s", last, lag, ref))
    t0 <- Sys.time()
    backend <- get_backend(a$model, a$ctx_len)
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
    S <- inv(s, log)[, lag + seq_len(HORIZON), drop = FALSE]
    model_id <- if (!is.null(a$model_id)) a$model_id else paste0("chlorocast_", a$model)
    df <- to_vera(S, ref, model_id, a$site, a$depth, a$target)
    stopifnot(!anyNA(df$prediction), all(df$prediction >= 0))
    dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
    f <- file.path(a$out, sprintf("%s.csv", model_id))
    write.csv(df, f, row.names = FALSE, quote = FALSE)
    cat(sprintf(
        "wrote %s (%d rows); context ends %s (lag %dd); load %.1fs, forecast %.2fs\n",
        f, nrow(df), last, lag, load_s, as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ))
}